"""Stage-0.6 controlled-waveform runner (dwell family and temporal-shuffle test).

Every job = (design, realization, variant, power, Nd, dt).  Within a realization all
variants share the amplitude multiset, the QPSK phase carrier, and the intensity-noise
draw, so A/B differences are paired.  One row per (job, model) goes to <tag>_rows.csv.
"""
from __future__ import annotations

import argparse
import json
import os
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.signal import butter, resample_poly, sosfilt
from scipy.special import logsumexp

from common import OUT, s
import models06 as M
import waveforms06 as W

CFG = json.loads((OUT / "stage06_config.json").read_text())
TAU_S = CFG["tau_atom_s"]
TAU = TAU_S * 1e9  # samples at 1 ns
DEC = 50  # 1 GHz -> 20 MHz
N = CFG["n_samples"]
CP = CFG["cyclic_prefix_samples"]
TS, TP = CFG["symbol_samples"], CFG["phase_ramp_samples"]
SKIP = CFG["analysis_skip_samples"]  # after the cyclic prefix
SOS = s.old.sos
SOS_DC = butter(4, 1e6, fs=1e9, output="sos")
QPSK = np.exp(1j * (np.pi / 4 + np.pi / 2 * np.arange(4)))
LAGS = np.arange(-10, 31)
CW_CARRIER = False
XIS = CFG["xis"]


_CTRL = {}


def controls(_nd=None):
    if not _CTRL:
        f1 = M.Controls.from_npz(OUT / "artifacts" / f"controls_Nd{CFG['controls_Nd']}.npz")
        _CTRL["f1"], _CTRL["mh"] = f1, M.ControlsMH(OUT / "artifacts" / f"controls_MH_Nd{CFG['controls_Nd']}.npz", f1)
    return _CTRL["f1"]


def lpf_delay():
    imp = np.zeros(20000)
    imp[0] = 1
    h = sosfilt(SOS, imp)
    return float(np.sum(np.arange(len(h)) * h) / np.sum(h))


DELAY = lpf_delay()


def realization(design, r):
    rng = np.random.default_rng(CFG["seed_base"] + 1000 * r + (0 if design == "dwell" else 500))
    if design == "dwell":
        amps, meta = W.dwell_family(rng, N, TAU, XIS, CFG["xi_ref"], CFG["edge_samples"], CFG["occupancy"], CFG["level_ratio"])
        amps = {f"xi_{k:g}": v for k, v in amps.items()}
        meta = {f"xi_{k:g}": {**v, "xi_nominal": k} for k, v in meta.items()}
    else:
        amps, info = W.shuffle_set(rng, N, TAU, **CFG["shuffle"])
        meta = {k: dict(info) for k in amps}
    car, labels = W.phase_carrier(np.random.default_rng(CFG["seed_base"] + 7 + 1000 * r), N, TS, TP)
    if CW_CARRIER:
        # Ablation: unmodulated carrier (physics metrics only; symbol metrics are degenerate).
        labels = np.zeros_like(labels)
        car = np.full(N, QPSK[0])
    noise = np.random.default_rng(CFG["seed_base"] + 9991 + 1000 * r).normal(0, s.SIGMA, N + CP)
    return amps, meta, car, labels, noise


def bb(y, dc):
    t = np.arange(len(y)) / 1e9
    return sosfilt(SOS, (y - dc) * np.exp(-2j * np.pi * 5e6 * t))[::DEC]


def ls_nmse(y, cols, powers=False):
    X = np.column_stack(cols)
    coef, *_ = np.linalg.lstsq(X, y, rcond=None)
    fit = X @ coef
    res_pow = float(np.mean(np.abs(y - fit) ** 2))
    fit_pow = float(np.mean(np.abs(fit - np.mean(fit)) ** 2))
    if powers:
        return res_pow / fit_pow, coef, res_pow, fit_pow
    return res_pow / fit_pow, coef


def fir_cols(x):
    n = len(x)
    cols = []
    for lag in LAGS:
        c = np.zeros(n, complex)
        if lag >= 0:
            c[lag:] = x[:n - lag]
        else:
            c[:lag] = x[-lag:]
        cols.append(c)
    return cols


def comm(z, xk, labels):
    k = np.arange(len(z))
    tr, ca, te = k % 3 == 0, k % 3 == 1, k % 3 == 2
    g, b = s.fit_gain(xk[tr], z[tr])
    zz = (z - b) / g
    var = max(float(np.mean(np.abs(zz[ca] - xk[ca]) ** 2)), 1e-18)
    unit = xk / QPSK[labels]
    cand = unit[:, None] * QPSK[None, :]
    logq = -np.abs(zz[:, None] - cand) ** 2 / var
    dens = 2 + (logq[k, labels] - logsumexp(logq, axis=1)) / np.log(2)
    dec = np.argmin(np.abs(zz[:, None] - cand) ** 2, axis=1)
    return {"AIR": float(np.mean(dens[te])), "SER": float(np.mean(dec[te] != labels[te])),
            "EVM": float(np.sqrt(np.mean(np.abs(zz[te] - xk[te]) ** 2) / np.mean(np.abs(xk[te]) ** 2))), "cal_var": var}


def analyze(outputs, env_cp, noise, labels, ctrl, a_unit):
    """outputs: {model: intensity incl. cyclic prefix}."""
    x = sosfilt(SOS, env_cp)[::DEC]
    w0 = (CP + SKIP) // DEC
    w1 = len(x) - 2000 // DEC
    xs = x[w0:w1]
    nsym = N // TS
    centers = CP + (np.arange(nsym) + 0.5) * TS
    idx = np.round(centers / DEC + DELAY / DEC).astype(int)
    keep = (idx >= w0) & (idx < len(x))
    idx, lab = idx[keep], labels[keep]
    res, ys = {}, {}
    for m, y in outputs.items():
        yb = bb(y, ctrl.dc)
        ys[m] = yb
        yn = bb(y + noise, ctrl.dc)
        d_gain, coef = ls_nmse(yb[w0:w1], [xs, np.ones_like(xs)])
        d_bla, _, r_pow, f_pow = ls_nmse(yb[w0:w1], fir_cols(xs) + [np.ones_like(xs)], powers=True)
        dcs = sosfilt(SOS_DC, y - ctrl.dc)[::DEC][w0:w1]
        res[m] = {"D_gain": d_gain, "D_BLA": d_bla, "BLA_residual_power": r_pow, "BLA_fit_ac_power": f_pow,
                  "gain_abs": float(abs(coef[0])), "gain_phase": float(np.angle(coef[0])),
                  "out_bb_rms": float(np.sqrt(np.mean(np.abs(yb[w0:w1]) ** 2))),
                  "dc_shift_mean": float(np.mean(dcs)), "dc_shift_std": float(np.std(dcs)),
                  **comm(yn[idx], x[idx], lab)}
    if "linear" in res:
        # Residual (non-BLA) power relative to a FIXED reference: the linear model's BLA-fit AC power
        # for the same waveform. Unlike D_BLA it does not depend on each model's own compressed gain.
        for m in res:
            res[m]["D_ref_linear"] = res[m]["BLA_residual_power"] / res["linear"]["BLA_fit_ac_power"]
    for m in ("static", "static_lti", "static_mh", "static_lti_mh"):
        if "atomic" in ys and m in ys:
            a, c = ys["atomic"][w0:w1], ys[m][w0:w1]
            g = np.vdot(c, a) / np.vdot(c, c)
            res["atomic"][f"NMSE_full_vs_{m}"] = float(np.mean(np.abs(a - g * c) ** 2) / np.mean(np.abs(a) ** 2))
    return res, ys, x


def run_job(design, r, variant, pdb, nd, dt, tag, cache, save_trace):
    jid = f"{design}_r{r}_{variant}_p{pdb:+d}_Nd{nd}_dt{dt*1e9:g}"
    jp = OUT / "jobs" / tag / f"{jid}.json"
    if jp.exists():
        return json.loads(jp.read_text())["rows"]
    if (design, r) not in cache:
        cache.clear()
        cache[(design, r)] = realization(design, r)
    amps, meta, car, labels, noise = cache[(design, r)]
    ctrl = controls(nd)
    a = amps[variant]
    env = ctrl.e1db * 10 ** (pdb / 20) * a * car
    env_cp = np.r_[env[-CP:], env]
    t0 = time.perf_counter()
    if dt == 1e-9:
        full = M.atomic(env_cp, nd)
    else:
        f = int(round(1e-9 / dt))
        full = M.atomic(resample_poly(env_cp, f, 1), nd, dt=dt)[::f]
    rt = time.perf_counter() - t0
    mh = _CTRL["mh"]
    outputs = {"atomic": full, "static": ctrl.static(env_cp), "static_lti": ctrl.static_lti(env_cp),
               "static_mh": mh.static_mh(env_cp), "static_lti_mh": mh.static_lti_mh(env_cp), "linear": ctrl.linear(env_cp)}
    if not np.isfinite(full).all():
        raise RuntimeError(f"nonfinite atomic output {jid}")
    res, ys, x = analyze(outputs, env_cp, noise, labels, ctrl, a)
    base = {"job_id": jid, "tag": tag, "design": design, "realization": r, "variant": variant, "Pavg_over_P1dB_dB": pdb, "Nd": nd, "dt_ns": dt * 1e9,
            "tau_atom_s": TAU_S, "E1dB_Vpm": ctrl.e1db, "runtime_atomic_s": rt, **{f"meta_{k}": v for k, v in meta[variant].items()}}
    rows = [{**base, "model": m, **v, "AIR_linear_ref": res["linear"]["AIR"], "loss_vs_linear": res["linear"]["AIR"] - v["AIR"]} for m, v in res.items()]
    if save_trace:
        (OUT / "traces" / tag).mkdir(parents=True, exist_ok=True)
        np.savez_compressed(OUT / "traces" / tag / f"{jid}.npz", x=x, a=a[::DEC], cp=CP // DEC, **{f"y_{m}": v for m, v in ys.items()},
                            **{f"dc_{m}": sosfilt(SOS_DC, outputs[m] - ctrl.dc)[::DEC] for m in outputs})
    jp.parent.mkdir(parents=True, exist_ok=True)
    tmp = jp.with_suffix(".tmp")
    tmp.write_text(json.dumps({"rows": rows}, default=float))
    tmp.replace(jp)
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--design", choices=("dwell", "shuffle"), required=True)
    ap.add_argument("--reals", default="0-7")
    ap.add_argument("--powers", default="-6,0,3,6")
    ap.add_argument("--variants", default="all")
    ap.add_argument("--nd", type=int, default=CFG["decisive_Nd"])
    ap.add_argument("--dt", type=float, default=1e-9)
    ap.add_argument("--tag", default="main")
    ap.add_argument("--traces", default="0")
    ap.add_argument("--cw-carrier", action="store_true")
    ap.add_argument("--n-samples", type=int, default=None, help="override waveform length (follow-up long-dwell family)")
    ap.add_argument("--xis", default=None, help="override dwell-family members, e.g. 1,4,8,16,32")
    a = ap.parse_args()
    global CW_CARRIER, N, XIS
    CW_CARRIER = a.cw_carrier
    if a.n_samples:
        N = a.n_samples
    if a.xis:
        XIS = [float(v) for v in a.xis.split(",")]
    lo, hi = (int(v) for v in a.reals.split("-")) if "-" in a.reals else (int(a.reals), int(a.reals))
    powers = [int(p) for p in a.powers.split(",")]
    trace_reals = {int(v) for v in a.traces.split(",") if v != ""}
    cache = {}
    out_csv = OUT / f"rows_{a.tag}_{a.design}.csv"
    for r in range(lo, hi + 1):
        if (a.design, r) not in cache:
            cache.clear()
            cache[(a.design, r)] = realization(a.design, r)
        variants = list(cache[(a.design, r)][0]) if a.variants == "all" else a.variants.split(",")
        for pdb in powers:
            for v in variants:
                t0 = time.perf_counter()
                rows = run_job(a.design, r, v, pdb, a.nd, a.dt, a.tag, cache, r in trace_reals)
                df = pd.DataFrame(rows)
                df.to_csv(out_csv, mode="a", header=not out_csv.exists(), index=False)
                at = df[df.model == "atomic"].iloc[0]
                print(f"{rows[0]['job_id']} AIRfull={at.AIR:.4f} AIRlti={df[df.model=='static_lti'].AIR.iloc[0]:.4f} "
                      f"AIRltiMH={df[df.model=='static_lti_mh'].AIR.iloc[0]:.4f} AIRstMH={df[df.model=='static_mh'].AIR.iloc[0]:.4f} "
                      f"Dfull={at.D_BLA:.3e} DltiMH={df[df.model=='static_lti_mh'].D_BLA.iloc[0]:.3e} {time.perf_counter()-t0:.1f}s", flush=True)


if __name__ == "__main__":
    main()
