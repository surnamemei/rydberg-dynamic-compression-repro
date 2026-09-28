"""P3: one standard memory model (GMP) vs the CW-derived surrogates (rules: results/tqe_viability/00_decision_rules.json, P3).

Data (reference configuration, OP1, 5 MHz IF): fair-waveform set (P4 re-run baseband), dwell/shuffle families (stored revision P1
traces), constant-amplitude phase cases (revision P2 runs; stored probe output at 50 MHz, rebuilt on the receiver chain).
All model outputs are receiver-chain basebands (20 MHz) referenced to the same zero-field level.
"""
from __future__ import annotations

import itertools
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from scipy.signal import butter, lfilter, resample_poly, sosfiltfilt, sosfilt, sosfreqz

from tv_common import TV, ROOT, R, RP, s

sys.path.insert(0, str(ROOT / "python" / "experiments" / "final_validation"))
import models06 as M  # noqa: E402
from common import OUT as OUT06  # noqa: E402

F1 = M.Controls.from_npz(OUT06 / "artifacts" / "controls_Nd4001.npz")
MH = M.ControlsMH(OUT06 / "artifacts" / "controls_MH_Nd4001.npz", F1)
E1 = F1.e1db
DC = F1.dc
IF = 5e6
FS, FS_BB = 1e9, 20e6
DEC = 50
SOS_IN = butter(8, 8e6, fs=FS, output="sos")
TRACES = ROOT / "results" / "revision" / "p1" / "traces" / "p1"
REVP2 = ROOT / "results" / "revision" / "p2"
RULES = json.loads((TV / "00_decision_rules.json").read_text())["P3_memory_model"]
COMPARATORS = ("static_mh", "static_lti_mh", "dsh")


# ------------------------------------------------------------------ receiver chain helpers
def bb(y, dc=DC):
    t = np.arange(len(y)) / FS
    return sosfilt(R.SOS, (y - dc) * np.exp(-2j * np.pi * IF * t))[::DEC]


def x_in(env):
    """GMP input: 8 MHz zero-phase low-pass of the 1 GHz envelope, 20 MHz, normalized by E1dB."""
    return sosfiltfilt(SOS_IN, env)[::DEC] / E1


def x_fit(env):
    return sosfilt(R.SOS, env)[::DEC]


def h20(ntap=128, nfft=1024):
    """20 MHz FIR with the frequency response of the 1 GHz receiver low-pass over |f| < 10 MHz."""
    f = np.fft.fftfreq(nfft, 1 / FS_BB)
    _, H = sosfreqz(R.SOS, worN=2 * np.pi * f / FS)
    return np.real(np.fft.ifft(H))[:ntap]


H20 = h20()


# ------------------------------------------------------------------ GMP basis (maximal set; smaller models are column subsets)
KA_MAX, MA_MAX, KB, LB_MAX = 5, 20, 2, 100


def basis_max(x):
    n = len(x)
    ax = np.abs(x)
    t = np.arange(n) / FS_BB
    down, up = np.exp(-2j * np.pi * IF * t), np.exp(2j * np.pi * IF * t)
    cols, names = [], []

    def lag(v, m):
        out = np.zeros_like(v)
        if m < n:
            out[m:] = v[:n - m]
        return out
    for k in range(KA_MAX):
        base = x * ax ** k
        for m in range(MA_MAX):
            cols.append(lag(base, m)); names.append(("aligned", k, m))
    for k in range(1, KB + 1):
        axk = ax ** k
        for l_ in range(1, LB_MAX + 1):
            cols.append(x * lag(axk, l_)); names.append(("lagging", k, l_))
    for k in range(1, 5):
        base = (ax ** k).astype(complex)
        for m in range(MA_MAX):
            cols.append(lag(base, m) * down); names.append(("zone0", k, m))
    for k in range(0, 3):
        base = x ** 2 * ax ** k
        for m in range(MA_MAX):
            cols.append(lag(base, m) * up); names.append(("zone2", k, m))
    cols.append(np.ones(n, complex)); names.append(("const", 0, 0))
    cols.append(down.astype(complex)); names.append(("zone0", 0, 0))
    Phi = np.column_stack(cols)
    Phi = lfilter(H20, [1.0], Phi, axis=0)
    return Phi, names


def subset(names, Ka, Ma, Lb):
    idx = []
    for i, (kind, k, m) in enumerate(names):
        if kind == "aligned" and k < Ka and m < Ma:
            idx.append(i)
        elif kind == "lagging" and m <= Lb:
            idx.append(i)
        elif kind in ("zone0", "zone2") and (m < Ma):
            idx.append(i)
        elif kind == "const":
            idx.append(i)
    return np.array(idx)


# ------------------------------------------------------------------ data loaders (lazy; one waveform at a time)
def split_of(family, key):
    sp = {"fair": {0: "train", 1: "train", 2: "val", 3: "val"}, "dwell": {0: "train", 1: "train", 2: "val", 3: "val"},
          "shuffle": {0: "train", 1: "train", 2: "val", 3: "val"}}
    if family == "phase":
        return {20260705: "train", 20260706: "train", 20260704: "val"}.get(key, "test")
    return sp[family].get(key, "test")


def list_waveforms():
    items = []
    for p in sorted(TRACES.glob("*.npz")):
        design, rest = p.stem.split("_r", 1)
        r = int(rest.split("_")[0])
        items.append({"family": "dwell" if design == "dwell" else "shuffle", "path": p, "r": r, "split": split_of("dwell" if design == "dwell" else "shuffle", r)})
    for p in sorted((TV / "p4" / "jobs").glob("*.npz")):
        z = np.load(p)
        items.append({"family": "fair", "path": p, "r": int(z["r"]), "split": split_of("fair", int(z["r"]))})
    cases = ["REF", "FAST", "SLOW", "JUMP", "LONGRAMP"]
    items.append({"family": "phase", "path": REVP2 / "OP1_CW_H.npz", "seed": 0, "case": "CW", "split": "train"})
    for case, seed in itertools.product(cases, RP.SEEDS):
        items.append({"family": "phase", "path": REVP2 / f"OP1_{case}_s{seed}_H.npz", "seed": seed, "case": case, "split": split_of("phase", seed)})
    return items


_REAL = {}


def realization(design, r):
    key = (design, r)
    if key not in _REAL:
        _REAL.clear()
        _REAL[key] = R.realization(design, r)
    return _REAL[key]


def load(item, comps=False):
    """Returns dict: x (GMP input), y (full), comps {name: y}, lin (linear model), xs (fit input), w (window), extras."""
    fam = item["family"]
    if fam in ("dwell", "shuffle"):
        z = np.load(item["path"])
        stem = item["path"].stem
        design = "dwell" if fam == "dwell" else "shuffle"
        variant = stem.split(f"_r{item['r']}_", 1)[1].rsplit("_p", 1)[0]
        pdb = int(stem.rsplit("_p", 1)[1].split("_")[0])
        amps, meta, car, labels, noise = realization(design, item["r"])
        env = E1 * 10 ** (pdb / 20) * amps[variant] * car
        env_cp = np.r_[env[-R.CP:], env]
        xs = x_fit(env_cp)
        assert np.allclose(xs, z["x"], rtol=1e-9, atol=1e-12), f"envelope mismatch {stem}"
        comps = {"static_mh": z["y_static_mh"], "static_lti_mh": z["y_static_lti_mh"], "dsh": bb(M.dsh(env_cp, R.TAU_S, MH))} if comps else {}
        w0, w1 = (R.CP + R.SKIP) // R.DEC, len(xs) - 2000 // R.DEC
        return {"x": x_in(env_cp), "y": z["y_atomic"], "comps": comps, "lin": z["y_linear"], "xs": xs, "w": (w0, w1),
                "id": stem.replace("_Nd4001_dt1", ""), "variant": variant, "p": pdb}
    if fam == "fair":
        z = np.load(item["path"])
        mod, seed, p = str(z["mod"]), int(z["seed"]), int(z["p"])
        tx, labels, alphabet = s.waveform(mod, seed)
        env = E1 * 10 ** (p / 20) * tx
        one = bb(np.ones(len(env)), dc=0.0)
        shift = (s.old.dc - DC) * one                       # archived receiver used old.dc; convert to the common reference
        comps = {m: z[f"bb_{m}"] + shift for m in COMPARATORS} if comps else {}
        n = len(z["bb_atomic"])
        return {"x": x_in(env), "y": z["bb_atomic"] + shift, "comps": comps, "lin": z["bb_linear"] + shift, "xs": x_fit(env),
                "w": (40, n - 40), "id": item["path"].stem, "mod": mod, "seed": seed, "p": p, "labels": labels, "alphabet": alphabet,
                "bb_noise": z["bb_noise"], "shift": shift}
    # phase
    z = np.load(item["path"])
    spec = json.loads(str(z["spec"]))
    n = int(round((RP.T0 + 60e-6) * FS))
    t = np.arange(n) / FS
    phi = RP.phase({**spec, "t_mod": spec.get("t_mod", 60e-6)}, t)
    env = float(z["amp"]) * np.exp(1j * phi)
    y50 = 10 ** (z["probe_db"] / 10)
    y1g = resample_poly(y50, 20, 1)[:n]
    comps = {"static_mh": bb(MH.static_mh(env)), "static_lti_mh": bb(MH.static_lti_mh(env)), "dsh": bb(M.dsh(env, R.TAU_S, MH))} if comps else {}
    tb = np.arange(n // DEC) / FS_BB
    w = (int(np.searchsorted(tb, 51e-6)), int(np.searchsorted(tb, 109.5e-6)))
    xwin = (int(np.searchsorted(tb, 89.5e-6)), int(np.searchsorted(tb, 109.5e-6)))
    return {"x": x_in(env), "y": bb(y1g), "comps": comps, "lin": bb(F1.linear(env)), "xs": x_fit(env), "w": w, "xwin": xwin,
            "id": item["path"].stem, "case": item["case"], "seed": item["seed"]}


# ------------------------------------------------------------------ metrics
def nmse(yh, y, w):
    a, b = w
    e = y[a:b] - yh[a:b]
    return float(np.mean(np.abs(e) ** 2) / np.mean(np.abs(y[a:b] - y[a:b].mean()) ** 2))


def gain_and_dref(yb, xs, w, lin_fit_pow=None):
    a, b = w
    xw = xs[a:b]
    X = np.column_stack([xw, np.ones_like(xw)])
    coef, *_ = np.linalg.lstsq(X, yb[a:b], rcond=None)
    cols = R.fir_cols(xw) + [np.ones_like(xw)]
    d_bla, _, r_pow, f_pow = R.ls_nmse(yb[a:b], cols, powers=True)
    return complex(coef[0]), r_pow, f_pow


def x_bb(yb, lin, xwin):
    a, b = xwin
    return float(np.median(np.abs(yb[a:b]) / np.abs(lin[a:b])))


def post_bb(bbv, mod):
    if mod == "OFDM":
        frames = bbv[:s.N_FRAMES * s.BLOCK].reshape(s.N_FRAMES, s.BLOCK)[:, s.CP:]
        return (np.fft.fft(frames, axis=1) / np.sqrt(s.NFFT))[:, s.OFDM_BINS]
    if mod == "CE":
        dphi = np.angle(bbv[1:] * np.conj(bbv[:-1]))
        return dphi * (s.FS_SC / (2 * np.pi))
    from scipy.signal import fftconvolve
    return fftconvolve(bbv, s.RRC, mode="same")[6::s.SPS]


def air(bb_model_common, d, tm, c=1.0):
    bbv = bb_model_common - d["shift"] + c * d["bb_noise"]
    return float(s.split_eval(post_bb(bbv, d["mod"]), d["labels"], d["alphabet"], d["mod"], tm)[0])


# ------------------------------------------------------------------ fitting
def combos():
    hp = RULES["hyperparameters"]
    return [dict(Ka=a, Ma=b, Lb=c, lam=d) for a, b, c, d in itertools.product(hp["K_a"], hp["M_a"], hp["L_b"], hp["ridge_lambda_relative"])]


def main():
    items = list_waveforms()
    fam_n = {}
    for it in items:
        if it["split"] == "train":
            fam_n[it["family"]] = fam_n.get(it["family"], 0) + 1
    # pass 1: weighted normal equations on the training set
    A = b = names = None
    cnt = {f: 0 for f in fam_n}
    for it in items:
        if it["split"] != "train":
            continue
        d = load(it)
        Phi, names = basis_max(d["x"])
        a0, a1 = d["w"]
        P_, y = Phi[a0:a1], d["y"][a0:a1]
        wgt = 1.0 / (fam_n[it["family"]] * (a1 - a0))
        A = (0 if A is None else A) + wgt * (P_.conj().T @ P_)
        b = (0 if b is None else b) + wgt * (P_.conj().T @ y)
        cnt[it["family"]] += 1
    print("training waveforms per family:", cnt, flush=True)
    sols = {}
    for c in combos():
        idx = subset(names, c["Ka"], c["Ma"], c["Lb"])
        Asub, bsub = A[np.ix_(idx, idx)], b[idx]
        scale = 1 / np.sqrt(np.real(np.diag(Asub)))
        As = Asub * scale[:, None] * scale[None, :]
        th = np.linalg.solve(As + c["lam"] * np.trace(As).real / len(idx) * np.eye(len(idx)), bsub * scale) * scale
        sols[json.dumps(c)] = (idx, th)
    # pass 2: validation NMSE per combo (family-averaged)
    val = {k: {} for k in sols}
    for it in items:
        if it["split"] != "val":
            continue
        d = load(it)
        Phi, _ = basis_max(d["x"])
        for k, (idx, th) in sols.items():
            val[k].setdefault(it["family"], []).append(nmse(Phi[:, idx] @ th, d["y"], d["w"]))
    vrows = [{**json.loads(k), **{f"val_nmse_{f}": float(np.median(v)) for f, v in fam.items()},
              "val_nmse_family_mean": float(np.mean([np.median(v) for v in fam.values()])), "n_coef": len(sols[k][0])} for k, fam in val.items()]
    vdf = pd.DataFrame(vrows).sort_values("val_nmse_family_mean")
    vdf.to_csv(TV / "03_memory_model_validation.csv", index=False)
    top = vdf.iloc[0]
    best_key = [k for k in sols if json.loads(k) == {"Ka": int(top.Ka), "Ma": int(top.Ma), "Lb": int(top.Lb), "lam": float(top.lam)}][0]
    idx, th = sols[best_key]
    print("selected:", best_key, "n_coef", len(idx), flush=True)
    np.savez_compressed(TV / "03_gmp_coefficients.npz", idx=idx, theta=th, names=np.array([f"{a}|{k}|{m}" for a, k, m in names], dtype=object),
                        selected=best_key)
    # pass 3: test metrics
    rows = []
    for it in items:
        is_ref = it["family"] == "phase" and it.get("case") == "CW"
        if it["split"] != "test" and not is_ref:
            continue
        d = load(it, comps=True)
        Phi, _ = basis_max(d["x"])
        yg = Phi[:, idx] @ th
        models = {"full": d["y"], "gmp": yg, **d["comps"]}
        g_lin, r_lin, f_lin = gain_and_dref(d["lin"], d["xs"], d["w"])
        base = {"family": it["family"], "id": d["id"], "split": "reference" if is_ref else "test"}
        for k in ("variant", "p", "mod", "seed", "case"):
            if k in d:
                base[k] = d[k]
        for m, yb in models.items():
            g, r_pow, _ = gain_and_dref(yb, d["xs"], d["w"])
            row = {**base, "model": m, "nmse_vs_full": np.nan if m == "full" else nmse(yb, d["y"], d["w"]),
                   "gain_rel": abs(g) / abs(g_lin), "D_ref": r_pow / f_lin}
            if it["family"] == "phase":
                row["g_bb"] = x_bb(yb, d["lin"], d["xwin"])
            if it["family"] == "fair":
                tm = json.loads((ROOT / "results" / "final_validation" / "regen_stage05" / "jobs" / f"{d['id']}.json").read_text())["ce_timing_offset"]
                row["AIR"] = air(yb, d, tm)
            rows.append(row)
    df = pd.DataFrame(rows)
    df.to_csv(TV / "03_memory_model_per_waveform.csv", index=False)
    print(df.groupby(["family", "model"]).nmse_vs_full.median().round(4).to_string())


if __name__ == "__main__":
    main()
