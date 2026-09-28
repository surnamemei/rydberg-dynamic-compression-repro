"""P4: AIR noise sensitivity of the Section III comparison (rules: results/tqe_viability/00_decision_rules.json, P4).

AIR at noise scale c is computed from each model's stored noiseless receiver baseband plus c x the baseband of the archived noise
draw (the receiver is linear up to the modulation-specific step). The loss, surrogate-error and threshold computations copy
final_validation/analyze_regen.py; at c = 1 they must reproduce the archived tables.
"""
from __future__ import annotations

import json
import sys

import numpy as np
import pandas as pd
from scipy import stats

from tv_common import TV, ROOT, s

sys.path.insert(0, str(ROOT / "python" / "experiments" / "final_validation"))
import analyze_regen as AR  # noqa: E402
from tv_p3_gmp import post_bb  # noqa: E402

REGEN = ROOT / "results" / "final_validation" / "regen_stage05"
SCALES = (0.25, 0.5, 1.0, 2.0)
MODELS = ("atomic", "static", "static_lti", "static_mh", "static_lti_mh", "linear", "dsh")
MODS = ("CE", "QPSK", "16QAM", "OFDM")


def ci(x):
    x = np.asarray(x, float)
    h = stats.t.ppf(.975, len(x) - 1) * x.std(ddof=1) / np.sqrt(len(x))
    return x.mean(), x.mean() - h, x.mean() + h


def air_rows():
    rows, val = [], []
    cache = {}
    for p in sorted((TV / "p4" / "jobs").glob("*.npz")):
        z = np.load(p)
        mod, seed, pdb = str(z["mod"]), int(z["seed"]), int(z["p"])
        if (mod, seed) not in cache:
            cache.clear()
            cache[(mod, seed)] = s.waveform(mod, seed)
        tx, labels, alphabet = cache[(mod, seed)]
        arch = json.loads((REGEN / "jobs" / f"{p.stem}.json").read_text())
        tm = arch["ce_timing_offset"]
        arch_air = {r["model"]: r["AIR_native"] for r in arch["rows"]}
        val.append({"job_id": p.stem, "orig_path_minus_archived": float(z["air_atomic_orig_c1"]) - arch_air["atomic"]})
        for c in SCALES:
            for m in MODELS:
                a = float(s.split_eval(post_bb(z[f"bb_{m}"] + c * z["bb_noise"], mod), labels, alphabet, mod, tm)[0])
                rows.append({"job_id": p.stem, "modulation": mod, "seed": seed, "Pavg_over_P1dB_dB": pdb, "model": m, "noise_scale": c, "AIR_native": a})
                if c == 1.0 and m in arch_air:
                    val[-1][f"decomposed_minus_archived_{m}"] = a - arch_air[m]
    return pd.DataFrame(rows), pd.DataFrame(val)


def tables(d):
    """analyze_regen.main() logic for one noise scale (losses, surrogate excess, thresholds)."""
    w = d.pivot_table(index=["modulation", "seed", "Pavg_over_P1dB_dB"], columns="model", values="AIR_native")
    loss_rows, exc_rows, thr_rows = [], [], []
    for mod in MODS:
        wm = w.loc[mod]
        for p in sorted(wm.index.get_level_values(1).unique()):
            seeds = [s_ for s_ in wm.index.get_level_values(0).unique() if (s_, p) in wm.index and (s_, -6) in wm.index]
            ref = np.array([wm.loc[(s_, -6), "atomic"] for s_ in seeds])
            full = np.array([wm.loc[(s_, p), "atomic"] for s_ in seeds])
            m, lo, hi, n = AR.ci(ref - full)
            rm, rlo, rhi, _ = AR.ci(100 * (ref - full) / ref)
            loss_rows.append({"modulation": mod, "Pavg_over_P1dB_dB": p, "AIR_full_mean": full.mean(), "loss_mean": m, "loss_ci_lo": lo,
                              "loss_ci_hi": hi, "rel_loss_pct": rm, "rel_loss_ci_lo": rlo, "rel_loss_ci_hi": rhi, "n": n})
            for ctrl in ("static", "static_lti", "static_mh", "static_lti_mh", "linear", "dsh"):
                e = np.array([wm.loc[(s_, p), ctrl] - wm.loc[(s_, p), "atomic"] for s_ in seeds])
                em, elo, ehi, _ = AR.ci(e)
                exc_rows.append({"modulation": mod, "Pavg_over_P1dB_dB": p, "control": ctrl, "excess_mean": em, "ci_lo": elo, "ci_hi": ehi, "n": len(e)})
        pp = [-6, -3, 0, 3, 6]
        seeds8 = [s_ for s_ in wm.index.get_level_values(0).unique() if all((s_, q) in wm.index for q in pp)]
        Rm = np.array([[100 * (wm.loc[(s_, -6), "atomic"] - wm.loc[(s_, q), "atomic"]) / wm.loc[(s_, -6), "atomic"] for q in pp] for s_ in seeds8])
        rng = np.random.default_rng(0)
        for lim in (5, 10):
            est = AR.threshold(pp, Rm.mean(axis=0), lim)
            boots = [AR.threshold(pp, Rm[rng.integers(0, len(Rm), len(Rm))].mean(axis=0), lim) for _ in range(2000)]
            thr_rows.append({"modulation": mod, "limit_pct": lim, "P_minus_P1dB_dB": est, "boot_lo": np.nanpercentile(boots, 2.5),
                             "boot_hi": np.nanpercentile(boots, 97.5), "n_seeds": len(Rm)})
    return pd.DataFrame(loss_rows), pd.DataFrame(exc_rows), pd.DataFrame(thr_rows)


def main():
    rows, val = air_rows()
    rows.to_csv(TV / "04_air_by_job_and_scale.csv", index=False)
    val.to_csv(TV / "04_validation_c1.csv", index=False)
    out, verdict = [], {"validation_max_abs_diff_c1": float(val.drop(columns="job_id").abs().max().max())}
    per_scale = {}
    for c in SCALES:
        L, E, T = tables(rows[rows.noise_scale == c])
        if c == 1.0:
            for mine, arch in ((L, "02_losses.csv"), (E, "03_excess_over_controls.csv"), (T, "04_thresholds.csv")):
                a = pd.read_csv(REGEN / arch)
                key = [k for k in ("modulation", "Pavg_over_P1dB_dB", "control", "limit_pct") if k in a.columns]
                mm = mine.merge(a, on=key, suffixes=("", "_arch"))
                num = [col for col in a.columns if col not in key and col + "_arch" in mm.columns and np.issubdtype(mm[col].dtype, np.number)]
                verdict[f"c1_reproduces_{arch}"] = float(max((mm[col] - mm[col + "_arch"]).abs().max() for col in num))
        for _, r in L.iterrows():
            out.append({"noise_scale": c, "table": "loss_full", **r.to_dict()})
        for _, r in E[E.control.isin(["static_mh", "static_lti_mh", "dsh"])].iterrows():
            out.append({"noise_scale": c, "table": "surrogate_minus_full", **r.to_dict()})
        for _, r in T.iterrows():
            out.append({"noise_scale": c, "table": "threshold", **r.to_dict()})
        # conclusion test
        cells = E[E.control.isin(["static_mh", "static_lti_mh"]) & E.Pavg_over_P1dB_dB.isin([0, 3, 6])]
        sig = cells[(cells.ci_lo > 0) | (cells.ci_hi < 0)]
        a_ok = all(((sig.control == ctl) & (sig.excess_mean.abs() >= .1)).any() for ctl in ("static_mh", "static_lti_mh"))
        b_ok = bool((sig.excess_mean > 0).any() and (sig.excess_mean < 0).any())
        rank = L[L.Pavg_over_P1dB_dB == 3].sort_values("loss_mean", ascending=False).modulation.tolist()
        per_scale[c] = {"a_each_surrogate_misses_by_0.1": bool(a_ok), "b_errors_of_both_signs": b_ok, "ranking_+3dB_most_to_least_loss": rank,
                        "P5_dB": dict(zip(T[T.limit_pct == 5].modulation, T[T.limit_pct == 5].P_minus_P1dB_dB)),
                        "P10_dB": dict(zip(T[T.limit_pct == 10].modulation, T[T.limit_pct == 10].P_minus_P1dB_dB))}
    ranks = {tuple(v["ranking_+3dB_most_to_least_loss"]) for v in per_scale.values()}
    holds = {c: v["a_each_surrogate_misses_by_0.1"] and v["b_errors_of_both_signs"] for c, v in per_scale.items()}
    if all(holds.values()) and len(ranks) == 1:
        cls = "ROBUST_TO_NOISE"
    elif holds[0.5] and holds[1.0] and holds[2.0]:
        cls = "PARTIALLY_ROBUST"
    else:
        cls = "NOISE_DEPENDENT"
    verdict.update({"classification": cls, "per_scale": {str(k): v for k, v in per_scale.items()}})
    pd.DataFrame(out).to_csv(TV / "04_noise_sensitivity.csv", index=False)
    (TV / "04_noise_verdict.json").write_text(json.dumps(verdict, indent=1, default=float))
    print(json.dumps(verdict, indent=1, default=float))
    tab = pd.DataFrame(out)
    lf = tab[(tab.table == "loss_full") & tab.Pavg_over_P1dB_dB.isin([0, 3, 6])].pivot_table(index=["modulation", "Pavg_over_P1dB_dB"], columns="noise_scale", values="loss_mean")
    print(lf.round(3).to_string())
    ex = tab[(tab.table == "surrogate_minus_full") & tab.Pavg_over_P1dB_dB.isin([0, 3, 6])].pivot_table(index=["control", "modulation", "Pavg_over_P1dB_dB"], columns="noise_scale", values="excess_mean")
    print(ex.round(2).to_string())


if __name__ == "__main__":
    main()
