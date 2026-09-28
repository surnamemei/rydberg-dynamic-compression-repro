from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import sys
import time
import traceback
from pathlib import Path

import numpy as np
from scipy.signal import fftconvolve, resample_poly
from scipy.special import logsumexp

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "results" / "p1db_waveform_stage05"
sys.path.insert(0, str(ROOT))
import p1db_comm_stage0 as old

FS = 1e9
FS_SC = 20e6
UPSAMPLE = 50
E1 = old.E1
SIGMA = old.SIGMA
ND = 1501
MODS = ("CE", "QPSK", "16QAM", "OFDM")
MODELS = ("atomic", "static", "static_lti")
POINTS = (-15, -10, -6, -3, 0, 3, 6)
PILOT = (-6, 0, 3)
REPS = 4
N_SYMBOLS = 1020
SPS = 4
N_SC = N_SYMBOLS * SPS
NFFT = 256
N_ACTIVE = 68
CP = 16
WOLA = 16
N_FRAMES = N_SYMBOLS // N_ACTIVE
BLOCK = NFFT + CP
FRAME = BLOCK + WOLA
OFDM_BINS = np.r_[np.arange(1, N_ACTIVE // 2 + 1), np.arange(NFFT - N_ACTIVE // 2, NFFT)]
RRC = old.SingleCarrier._rrcosdesign(None, beta=.2, span=10, sps=SPS)
ALPHA_CE = .36
CONFIG_VERSION = "fair_v4_randomized_balanced_ofdm_training"

def config_for(mod):
    return CONFIG_VERSION if mod == "OFDM" else "fair_v2"


def _norm(x: np.ndarray) -> np.ndarray:
    return x / np.sqrt(np.mean(np.abs(x) ** 2))


def waveform(mod: str, seed: int):
    rng = np.random.default_rng(seed)
    if mod in ("QPSK", "16QAM"):
        alphabet = old.QPSK if mod == "QPSK" else old.QAM
        labels = rng.integers(0, len(alphabet), N_SYMBOLS)
        up = np.zeros(N_SC, complex)
        up[::SPS] = alphabet[labels]
        tx20 = fftconvolve(up, RRC, mode="same")
        return _norm(resample_poly(_norm(tx20), UPSAMPLE, 1)), labels, alphabet
    if mod == "CE":
        labels = rng.integers(0, 2, N_SYMBOLS)
        symbols = 2 * labels - 1
        impulses = np.zeros(N_SC)
        impulses[::SPS] = symbols
        t = np.arange(-8 * SPS, 8 * SPS + 1) / SPS
        alpha = np.sqrt(np.log(2)) / ALPHA_CE
        h = np.exp(-2 * (np.pi * t / alpha) ** 2)
        h /= h.sum()
        freq = fftconvolve(impulses, h, mode="same") * SPS
        freq_ns = np.repeat(freq, UPSAMPLE)
        phase = np.cumsum(freq_ns * (np.pi * 5e6 / (2 * FS)))
        tx = np.exp(1j * phase)
        return tx, labels, np.array([-1.0, 1.0])
    if mod == "OFDM":
        alphabet = old.QPSK
        labels = rng.integers(0, len(alphabet), (N_FRAMES, N_ACTIVE))
        # Four-frame known preamble uses an independent random permutation of all
        # QPSK points on each carrier; it identifies gain without correlating carriers.
        for k in range(N_ACTIVE):
            labels[:4, k] = rng.permutation(4)
        freq = np.zeros((N_FRAMES, NFFT), complex)
        freq[:, OFDM_BINS] = alphabet[labels]
        time = np.fft.ifft(freq, axis=1) * np.sqrt(NFFT)
        window = np.ones(FRAME)
        ramp = np.sin(np.pi / 2 * (np.arange(WOLA) + .5) / WOLA) ** 2
        window[:WOLA] = ramp
        window[-WOLA:] = ramp[::-1]
        tx20 = np.zeros(N_FRAMES * BLOCK + WOLA, complex)
        for j in range(N_FRAMES):
            block = np.r_[time[j, -CP:], time[j], time[j, :WOLA]]
            tx20[j * BLOCK:j * BLOCK + FRAME] += block * window
        tx20 = tx20[:N_FRAMES * BLOCK]
        return _norm(resample_poly(_norm(tx20), UPSAMPLE, 1)), labels, alphabet
    raise ValueError(mod)


def occupied99(x: np.ndarray) -> float:
    # One-sided 0.5% and 99.5% cumulative-power quantiles, with 8x zero padding.
    nfft = 1 << int(np.ceil(np.log2(len(x) * 8)))
    spec = np.fft.fftshift(np.fft.fft(x, nfft))
    p = np.abs(spec) ** 2
    c = np.cumsum(p) / np.sum(p)
    f = np.fft.fftshift(np.fft.fftfreq(nfft, 1 / FS))
    return float(f[np.searchsorted(c, .995)] - f[np.searchsorted(c, .005)])


def extract(intensity: np.ndarray, mod: str) -> np.ndarray:
    t = np.arange(len(intensity)) / FS
    bb = old.sosfilt(old.sos, (intensity - old.dc) * np.exp(-2j * np.pi * old.IF * t))[::UPSAMPLE]
    if mod == "OFDM":
        frames = bb[:N_FRAMES * BLOCK].reshape(N_FRAMES, BLOCK)[:, CP:]
        return (np.fft.fft(frames, axis=1) / np.sqrt(NFFT))[:, OFDM_BINS]
    if mod == "CE":
        z = bb
        dphi = np.angle(z[1:] * np.conj(z[:-1]))
        return dphi * (FS_SC / (2 * np.pi))
    return fftconvolve(bb, RRC, mode="same")[6::SPS]


def static_intensity(env: np.ndarray) -> np.ndarray:
    a = np.abs(env)
    c = np.interp(a, old.AMPS, old.HARM.real) + 1j * np.interp(a, old.AMPS, old.HARM.imag)
    out = np.zeros(len(env), dtype=complex)
    np.divide(c * env, a, out=out, where=a > 0)
    t = np.arange(len(env)) / FS
    return old.dc + np.real(out * np.exp(2j * np.pi * old.IF * t))


def static_lti_intensity(env: np.ndarray) -> np.ndarray:
    # CW fundamental divided by the small-signal LTI response gives an equivalent
    # complex envelope transfer, as in the validated Stage-0 control.
    if not hasattr(static_lti_intensity, "gain"):
        n = 40000
        t = np.arange(n) / FS
        probe = old.linear_intensity(np.full(n, .01, complex))
        hif = 2 * np.mean((probe[-4000:] - old.dc) * np.exp(-2j * np.pi * old.IF * t[-4000:])) / .01
        gain = np.zeros(len(old.AMPS), complex)
        gain[1:] = old.HARM[1:] / (old.AMPS[1:] * hif)
        gain[0] = gain[1]
        static_lti_intensity.gain = gain
    a = np.abs(env)
    gain = static_lti_intensity.gain
    envelope = (np.interp(a, old.AMPS, gain.real) + 1j * np.interp(a, old.AMPS, gain.imag)) * env
    return old.linear_intensity(envelope)


def atomic_intensity(env: np.ndarray, dt: float = 1e-9, nd: int = ND) -> np.ndarray:
    n = len(env)
    t = np.arange(n) * dt
    cfg = old.SimConfig(300, n, dt, n * dt, t)
    res = old.sim.run(old.raqr, cfg, old.raqr.A_LO + env * np.exp(2j * np.pi * old.IF * t), Nd=nd, device="cuda")
    return 10 ** (res.probeResponse / 10)


def fit_gain(x: np.ndarray, y: np.ndarray):
    xm, ym = np.mean(x), np.mean(y)
    den = np.vdot(x - xm, x - xm)
    if den.real < 1e-12:
        raise ValueError("training symbols do not identify a receiver gain")
    g = np.vdot(x - xm, y - ym) / den
    if not np.isfinite(g):
        raise ValueError("nonfinite receiver gain estimate")
    return g, ym - g * xm


def gaussian_air(z, labels, alphabet, var):
    z = np.asarray(z).ravel()
    labels = np.asarray(labels).ravel().astype(int)
    var = max(float(var), 1e-12)
    logq = -np.abs(z[:, None] - alphabet[None, :]) ** 2 / var
    vals = np.log2(len(alphabet)) + (logq[np.arange(len(z)), labels] - logsumexp(logq, axis=1)) / np.log(2)
    return float(np.mean(vals))


def split_eval(raw, labels, alphabet, mod, ce_offset=None):
    if mod == "OFDM":
        labels = np.asarray(labels)
        x = alphabet[labels]
        g = np.empty(N_ACTIVE, complex)
        b = np.empty(N_ACTIVE, complex)
        for k in range(N_ACTIVE):
            g[k], b[k] = fit_gain(x[:5, k], raw[:5, k])
        z = (raw - b) / g
        cal, cal_lab = z[5:10], labels[5:10]
        test, test_lab = z[10:], labels[10:]
        var = max(float(np.mean(np.abs(cal - alphabet[cal_lab]) ** 2)), 1e-12)
        air = gaussian_air(test, test_lab, alphabet, var)
        decisions = np.argmin(np.abs(test.reshape(-1, 1) - alphabet[None, :]) ** 2, axis=1)
        truth = test_lab.ravel()
        ber = float(np.mean(old.BITS_QPSK[decisions] != old.BITS_QPSK[truth]))
        return air, float(np.mean(decisions != truth)), ber, float(np.sqrt(np.mean(np.abs(test - alphabet[test_lab]) ** 2))), var, None

    labels = np.asarray(labels).ravel()
    raw = np.asarray(raw).ravel()
    if mod == "CE":
        n = len(labels)
        best = None
        offsets = [ce_offset] if ce_offset is not None else range(SPS)
        for off in offsets:
            ids = off + SPS * np.arange(n)
            ids = ids[ids + SPS <= len(raw)]
            count = len(ids)
            vals = np.array([np.mean(raw[i:i + SPS]) for i in ids])
            labs = labels[:count]
            g, b = fit_gain(alphabet[labs[:340]], vals[:340])
            resid = vals[:340] - (g * alphabet[labs[:340]] + b)
            score = float(np.mean(np.abs(resid) ** 2))
            if best is None or score < best[0]:
                best = score, int(off), vals, labs, g, b
        _, ce_offset, vals, labels, g, b = best
        z = (vals - b) / g
        cal, cal_lab = z[340:680], labels[340:680]
        test, test_lab = z[680:], labels[680:]
        var = max(float(np.mean((cal - alphabet[cal_lab]) ** 2)), 1e-12)
        air = gaussian_air(test.astype(complex), test_lab, alphabet.astype(complex), var)
        decisions = (test > 0).astype(int)
        return air, float(np.mean(decisions != test_lab)), float(np.mean(decisions != test_lab)), float(np.sqrt(np.mean((test - alphabet[test_lab]) ** 2))), var, ce_offset

    labels = labels[:len(raw)]
    x = alphabet[labels]
    g, b = fit_gain(x[:340], raw[:340])
    z = (raw - b) / g
    cal, cal_lab = z[340:680], labels[340:680]
    test, test_lab = z[680:], labels[680:]
    var = max(float(np.mean(np.abs(cal - alphabet[cal_lab]) ** 2)), 1e-12)
    air = gaussian_air(test, test_lab, alphabet, var)
    decisions = np.argmin(np.abs(test[:, None] - alphabet[None, :]) ** 2, axis=1)
    ser = float(np.mean(decisions != test_lab))
    if mod == "16QAM":
        ber = float(np.mean(old.BITS_QAM[decisions] != old.BITS_QAM[test_lab]))
    else:
        ber = float(np.mean(old.BITS_QPSK[decisions] != old.BITS_QPSK[test_lab]))
    return air, ser, ber, float(np.sqrt(np.mean(np.abs(test - alphabet[test_lab]) ** 2))), var, ce_offset


def rate_hz(mod: str) -> float:
    # Rate counts information-bearing QAM/CPFSK symbols and excludes CP samples,
    # while including the cyclic prefix in each OFDM block duration.
    return FS_SC / SPS if mod != "OFDM" else N_ACTIVE * FS_SC / BLOCK


def stats_row(mod, seed, tx, db):
    mag = np.abs(tx)
    papr = float(np.max(mag ** 2) / np.mean(mag ** 2))
    norm = mag / np.sqrt(np.mean(mag ** 2))
    scale = E1 * 10 ** (db / 20)
    actual = scale * mag
    rms = float(np.sqrt(np.mean(actual ** 2)))
    peak = float(np.max(actual))
    dt = 1 / FS
    row = {"config_version": config_for(mod), "modulation": mod, "seed": seed, "Pavg_over_P1dB_dB": db, "sample_rate_Hz": FS, "information_symbol_rate_Hz": rate_hz(mod), "n_information_symbols": N_SYMBOLS, "E_rms_Vpm": rms, "E_peak_Vpm": peak, "E_peak_over_E1dB": peak/E1, "PAPR_linear": papr, "PAPR_dB": 10*np.log10(papr), "B_occ_99_Hz": occupied99(tx), "amplitude_mean_norm": float(np.mean(norm)), "amplitude_std_norm": float(np.std(norm)), "amplitude_q01_norm": float(np.quantile(norm,.01)), "amplitude_q50_norm": float(np.quantile(norm,.5)), "amplitude_q90_norm": float(np.quantile(norm,.9)), "amplitude_q99_norm": float(np.quantile(norm,.99)), "amplitude_q999_norm": float(np.quantile(norm,.999))}
    for theta in (.5,.75,1.,1.25,1.5):
        mask = actual > theta * E1 * (1 + 1e-12)
        edge = np.diff(np.r_[False, mask, False].astype(np.int8))
        starts, ends = np.flatnonzero(edge == 1), np.flatnonzero(edge == -1)
        dwell = (ends - starts) * dt
        stem = f"theta_{theta:.2f}"
        row[f"exceedance_{stem}"] = float(np.mean(mask))
        row[f"occupancy_{stem}"] = float(np.mean(mask))
        row[f"dwell_mean_ns_{stem}"] = float(np.mean(dwell)*1e9) if len(dwell) else 0.
        row[f"dwell_median_ns_{stem}"] = float(np.median(dwell)*1e9) if len(dwell) else 0.
        row[f"dwell_p90_ns_{stem}"] = float(np.quantile(dwell,.9)*1e9) if len(dwell) else 0.
        row[f"dwell_max_ns_{stem}"] = float(np.max(dwell)*1e9) if len(dwell) else 0.
        row[f"crest_event_rate_Hz_{stem}"] = float(len(starts) / (len(tx)*dt))
    ac = mag - np.mean(mag)
    if float(np.mean(ac**2)) < 1e-24:
        row["envelope_corr_1e_ns"] = float("nan")
        row["envelope_corr_1e_s"] = float("nan")
    else:
        nfft = 1 << int(np.ceil(np.log2(2*len(ac)-1)))
        spectrum = np.fft.rfft(ac, nfft)
        corr = np.fft.irfft(spectrum*np.conj(spectrum), nfft)[:len(ac)]
        corr /= max(corr[0], 1e-30)
        cut = np.flatnonzero(corr <= np.exp(-1))
        row["envelope_corr_1e_ns"] = float(cut[0]*dt*1e9) if len(cut) else float(len(tx)*dt*1e9)
        row["envelope_corr_1e_s"] = row["envelope_corr_1e_ns"]*1e-9
    return row


def job_id(mod, seed, db, nd):
    return f"{config_for(mod)}_{mod}_s{seed}_p{db:+03d}_Nd{nd}"


def append_csv(path: Path, row: dict):
    path.parent.mkdir(parents=True, exist_ok=True)
    exists = path.exists() and path.stat().st_size > 0
    with path.open("a", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(row))
        if not exists:
            w.writeheader()
        w.writerow(row)
        f.flush()
        os.fsync(f.fileno())


def atomic_save_json(path: Path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(value, indent=2, sort_keys=True), encoding="utf-8")
    tmp.replace(path)


def file_hash(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def run_job(mod, seed, db, nd, tx, labels, alphabet, stats, ce_offset):
    jid = job_id(mod, seed, db, nd)
    jp = OUT / "jobs" / f"{jid}.json"
    if jp.exists():
        try:
            saved = json.loads(jp.read_text(encoding="utf-8"))
            if saved.get("status") == "complete" and saved.get("modulation") == mod and saved.get("seed") == seed and saved.get("power_db") == db and saved.get("Nd") == nd and np.isfinite(saved.get("AIR_native", np.nan)):
                csv_path = OUT / "02_all_results.csv"
                existing = set()
                if csv_path.exists():
                    with csv_path.open(newline="", encoding="utf-8") as f:
                        existing = {(r.get("job_id"), r.get("model")) for r in csv.DictReader(f)}
                for row in saved.get("rows", []):
                    if (jid, row.get("model")) not in existing:
                        append_csv(csv_path, row)
                return saved, saved.get("ce_timing_offset", ce_offset)
        except Exception:
            pass
    E = E1 * 10 ** (db / 20)
    env = E * tx
    noise = np.random.default_rng(seed + 9991).normal(0, SIGMA, len(tx))
    start = time.perf_counter()
    outputs = {"atomic": atomic_intensity(env, nd=nd), "static": static_intensity(env), "static_lti": static_lti_intensity(env)}
    if not np.isfinite(outputs["atomic"]).all():
        raise RuntimeError(f"Nonfinite atomic output for {jid}")
    out_rows = []
    for model in MODELS:
        raw = extract(outputs[model] + noise, mod)
        timing = ce_offset
        if mod == "CE" and timing is None and model == "atomic":
            # Fixed synchronization offset is selected from the designated low-power
            # training segment once, then reused at every level and control.
            _, _, _, _, _, timing = split_eval(raw, labels, alphabet, mod, None)
        air, ser, ber, evm, var, timing = split_eval(raw, labels, alphabet, mod, timing)
        if not np.isfinite([air, ser, ber, evm, var]).all():
            raise RuntimeError(f"Nonfinite receiver metric for {jid}, model={model}")
        rate = rate_hz(mod)
        b_occ = stats["B_occ_99_Hz"]
        row = {"job_id": jid, "config_version": config_for(mod), "modulation": mod, "seed": seed, "model": model, "Pavg_over_P1dB_dB": db, "Ppeak_over_P1dB_dB": db + stats["PAPR_dB"], "E_rms_Vpm": E, "E_peak_Vpm": stats["E_peak_Vpm"], "PAPR_dB": stats["PAPR_dB"], "B_occ_99_Hz": b_occ, "information_symbol_rate_Hz": rate, "AIR_native": air, "R_AIR_bit_s": air*rate, "eta_AIR_bit_s_Hz": air*rate/b_occ, "SER": ser, "BER": ber, "EVM": evm, "cal_var": var, "gmi_s": 1.0, "ce_timing_offset": timing, "n_test_symbols": int(np.asarray(labels).size - 681) if mod != "OFDM" else 5*N_ACTIVE, "Nd": nd, "dt_s": 1e-9, "noise_sigma_intensity": SIGMA, "runtime_s_all_models": time.perf_counter()-start}
        out_rows.append(row)
    result = {"job_id": jid, "status": "complete", "modulation": mod, "seed": seed, "power_db": db, "Nd": nd, "ce_timing_offset": timing, "AIR_native": out_rows[0]["AIR_native"], "rows": out_rows, "runtime_s": time.perf_counter()-start}
    atomic_save_json(jp, result)
    for row in out_rows:
        append_csv(OUT / "02_all_results.csv", row)
    log = OUT / "run.log"
    with log.open("a", encoding="utf-8") as f:
        f.write(f"DONE {jid} models={len(out_rows)} seconds={result['runtime_s']:.3f}\n")
        f.flush()
        os.fsync(f.fileno())
    return result, timing


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--grid", choices=("pilot", "high", "full"), default="pilot")
    ap.add_argument("--replicates", type=int, default=REPS)
    ap.add_argument("--nd", type=int, default=ND)
    ap.add_argument("--only-mod", choices=MODS)
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    points = PILOT if args.grid == "pilot" else ((-6, -3, 0, 3, 6) if args.grid == "high" else POINTS)
    mods = (args.only_mod,) if args.only_mod else MODS
    seeds = [20262001 + 100*r for r in range(args.replicates)]
    manifest = {"commit": "409572499ef54825b66c532e0aecdcffa89e7be1", "branch": "p1db-waveform-fairness", "stage": "0.5", "config_versions": {m: config_for(m) for m in mods}, "grid": args.grid, "replicates": args.replicates, "Nd": args.nd, "dt_s": 1e-9, "modulations": mods, "models": MODELS, "power_db": points, "seed_base": seeds, "unique_job_key": ["modulation", "B_occ_99_Hz", "information_symbol_rate_Hz", "power_db", "seed", "model", "Nd"], "run_started": time.strftime("%Y-%m-%dT%H:%M:%S%z"), "python": sys.version, "numpy": np.__version__, "cuda_dll_dir": os.environ.get("RYDBERG_CUDA_DLL_DIR", ".cuda-env/Library/bin"), "sigma_intensity": SIGMA}
    # Keep this invocation-level record separate from the consolidated manifest
    # that analyze_stage05.py writes after combining every run.
    atomic_save_json(OUT / "04_job_manifest_last_run.json", manifest)
    stat_path = OUT / "03_waveform_statistics_extended.csv"
    timing_by_seed = {}
    total = len(mods)*len(seeds)*len(points)
    done = 0
    for mod in mods:
        for rep, base in enumerate(seeds):
            seed = base + {"CE": 1, "QPSK": 2, "16QAM": 3, "OFDM": 4}[mod]
            tx, labels, alphabet = waveform(mod, seed)
            ce_timing = timing_by_seed.get((mod, seed))
            for db in points:
                done += 1
                try:
                    stats = stats_row(mod, seed, tx, db)
                    stats["seed_group"] = rep
                    existing_stats = set()
                    if stat_path.exists():
                        with stat_path.open(newline="", encoding="utf-8") as f:
                            existing_stats = {(r.get("config_version"), r.get("modulation"), int(r.get("seed", -1)), int(float(r.get("Pavg_over_P1dB_dB", -999)))) for r in csv.DictReader(f)}
                    if (config_for(mod), mod, seed, db) not in existing_stats:
                        append_csv(stat_path, stats)
                    result, ce_timing = run_job(mod, seed, db, args.nd, tx, labels, alphabet, stats, ce_timing)
                    timing_by_seed[(mod, seed)] = ce_timing
                    print(f"[{done}/{total}] {result['job_id']} status={result['status']} AIR={result['AIR_native']:.6f} CE_offset={ce_timing}", flush=True)
                except Exception as exc:
                    jid = job_id(mod, seed, db, args.nd)
                    atomic_save_json(OUT / "jobs" / f"{jid}.json", {"job_id": jid, "status": "failed", "error": repr(exc), "traceback": traceback.format_exc(), "modulation": mod, "seed": seed, "power_db": db, "Nd": args.nd})
                    with (OUT / "run.log").open("a", encoding="utf-8") as f:
                        f.write(f"FAILED {jid} {exc!r}\n{traceback.format_exc()}\n")
                        f.flush()
                    print(f"FAILED {jid}: {exc!r}", flush=True)


if __name__ == "__main__":
    main()


