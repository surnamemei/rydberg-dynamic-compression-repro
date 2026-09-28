"""Part IV: AIR noise sensitivity of the Section III conclusion (plan: results/tqe_viability_final/00_plan.md, section D).

AIR at noise scale c = AIR of (stored noiseless receiver baseband + c x baseband of the archived noise draw); the draw is paired
across models and scales and the symbol partition (split_eval, archived CE timing) is unchanged. Loss / excess / threshold logic
copies final_validation/analyze_regen.py; at c = 1 the archived tables must be reproduced.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from tvf_common import FIN, ROOT, TV, git_state, now, s

sys.path.insert(0, str(ROOT / "python" / "experiments" / "final_validation"))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "tqe_viability"))
import analyze_regen as AR  # noqa: E402
from tv_p3_gmp import post_bb  # noqa: E402

OUT = FIN / "04_noise"
REGEN = ROOT / "results" / "final_validation" / "regen_stage05"
SCALES = (0.25, 0.5, 1.0, 2.0, 4.0)
MODELS = ("atomic", "static", "static_lti", "static_mh", "static_lti_mh", "linear", "dsh")
SURR = ("static_mh", "static_lti_mh")
MODS = ("CE", "QPSK", "16QAM", "OFDM")
FMT_COL = {"CE": ("#e87ba4", "o"), "QPSK": ("#008300", "s"), "16QAM": ("#4a3aa7", "^"), "OFDM": ("#e34948", "D")}


def air_rows():
    rows, val, snr, bits = [], [], [], {}
    cache = {}
    for p in sorted((TV / "p4" / "jobs").glob("*.npz")):
        z = np.load(p)
        mod, seed, pdb = str(z["mod"]), int(z["seed"]), int(z["p"])
        if (mod, seed) not in cache:
            cache.clear()
            cache[(mod, seed)] = s.waveform(mod, seed)
        tx, labels, alphabet = cache[(mod, seed)]
        bits[mod] = float(np.log2(len(alphabet)))
        arch = json.loads((REGEN / "jobs" / f"{p.stem}.json").read_text())
        tm = arch["ce_timing_offset"]
        arch_air = {r["model"]: r["AIR_native"] for r in arch["rows"]}
        v = {"job_id": p.stem, "orig_path_minus_archived_atomic": float(z["air_atomic_orig_c1"]) - arch_air["atomic"]}
        n = len(z["bb_atomic"])
        w = slice(40, n - 40)
        p_lin, p_noise = np.mean(np.abs(z["bb_linear"][w]) ** 2), np.mean(np.abs(z["bb_noise"][w]) ** 2)
        for c in SCALES:
            snr.append({"job_id": p.stem, "modulation": mod, "seed": seed, "Pavg_over_P1dB_dB": pdb, "noise_scale": c,
                        "eff_SNR_dB": float(10 * np.log10(p_lin / (c * c * p_noise)))})
            for m in MODELS:
                a = float(s.split_eval(post_bb(z[f"bb_{m}"] + c * z["bb_noise"], mod), labels, alphabet, mod, tm)[0])
                rows.append({"job_id": p.stem, "modulation": mod, "seed": seed, "Pavg_over_P1dB_dB": pdb, "model": m, "noise_scale": c, "AIR_native": a})
                if c == 1.0 and m in arch_air:
                    v[f"decomposed_minus_archived_{m}"] = a - arch_air[m]
        val.append(v)
    return pd.DataFrame(rows), pd.DataFrame(val), pd.DataFrame(snr), bits


def tables(d):
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
            loss_rows.append({"modulation": mod, "Pavg_over_P1dB_dB": p, "AIR_full_mean": full.mean(), "AIR_ref_mean": ref.mean(), "loss_mean": m,
                              "loss_ci_lo": lo, "loss_ci_hi": hi, "rel_loss_pct": rm, "rel_loss_ci_lo": rlo, "rel_loss_ci_hi": rhi, "n": n})
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
                             "boot_hi": np.nanpercentile(boots, 97.5), "n_seeds": len(Rm), "censored_at_grid_top": bool(est == pp[-1]),
                             "below_grid": bool(np.isnan(est))})
    return pd.DataFrame(loss_rows), pd.DataFrame(exc_rows), pd.DataFrame(thr_rows)


def main():
    rows, val, snr, bits = air_rows()
    OUT.mkdir(parents=True, exist_ok=True)
    rows.to_csv(OUT / "air_by_job_and_scale.csv", index=False)
    val.to_csv(OUT / "validation_c1.csv", index=False)
    snr.to_csv(OUT / "effective_snr_by_job.csv", index=False)
    vv = val.drop(columns="job_id").abs().max()
    surr_cols = [c for c in vv.index if c.startswith("decomposed_minus_archived_") and not c.endswith("_atomic")]
    validation = {"max_abs_surrogate_and_linear": float(vv[surr_cols].max()), "max_abs_atomic": float(vv["decomposed_minus_archived_atomic"]),
                  "max_abs_atomic_original_path": float(vv["orig_path_minus_archived_atomic"])}
    validation["pass"] = bool(validation["max_abs_surrogate_and_linear"] <= 1e-9 and validation["max_abs_atomic"] <= 1e-3)
    out, per_scale = [], {}
    for c in SCALES:
        L, E, T = tables(rows[rows.noise_scale == c])
        if c == 1.0:
            for mine, arch in ((L, "02_losses.csv"), (E, "03_excess_over_controls.csv"), (T, "04_thresholds.csv")):
                a = pd.read_csv(REGEN / arch)
                key = [k for k in ("modulation", "Pavg_over_P1dB_dB", "control", "limit_pct") if k in a.columns]
                mm = mine.merge(a, on=key, suffixes=("", "_arch"))
                num = [col for col in a.columns if col not in key and col + "_arch" in mm.columns and np.issubdtype(mm[col].dtype, np.number)]
                validation[f"c1_reproduces_{arch}"] = float(np.nanmax([(mm[col] - mm[col + "_arch"]).abs().max() for col in num]))
        for _, r in L.iterrows():
            out.append({"noise_scale": c, "table": "loss_full", **r.to_dict()})
        for _, r in E[E.control.isin(["static_mh", "static_lti_mh", "dsh"])].iterrows():
            out.append({"noise_scale": c, "table": "surrogate_minus_full", **r.to_dict()})
        for _, r in T.iterrows():
            out.append({"noise_scale": c, "table": "threshold", **r.to_dict()})
        sn = snr[snr.noise_scale == c].groupby(["modulation", "Pavg_over_P1dB_dB"]).eff_SNR_dB.mean()
        for (mod, p), v in sn.items():
            out.append({"noise_scale": c, "table": "effective_SNR_dB", "modulation": mod, "Pavg_over_P1dB_dB": p, "eff_SNR_dB": v})
        # tests
        cells = E[E.control.isin(SURR) & E.Pavg_over_P1dB_dB.isin([0, 3, 6])]
        sig = cells[(cells.ci_lo > 0) | (cells.ci_hi < 0)]
        a_ok = {ctl: bool(((sig.control == ctl) & (sig.excess_mean.abs() >= .1)).any()) for ctl in SURR}
        b_ok = bool((sig.excess_mean > 0).any() and (sig.excess_mean < 0).any())
        l6, l3 = L[L.Pavg_over_P1dB_dB == -6].set_index("modulation"), L[L.Pavg_over_P1dB_dB == 3].set_index("modulation")
        informative = []
        for mod in MODS:
            top = bits[mod]
            a6, a3 = l6.loc[mod, "AIR_full_mean"], l3.loc[mod, "AIR_full_mean"]
            near_ceiling = a6 >= .95 * top and a3 >= .95 * top
            near_floor = a6 <= .05 * top and a3 <= .05 * top
            if not (near_ceiling or near_floor):
                informative.append(mod)
        loss3 = l3.loss_mean
        per_scale[c] = {"a_static_mh": a_ok["static_mh"], "a_static_lti_mh": a_ok["static_lti_mh"], "b_both_signs": b_ok,
                        "surrogate_failure_holds": bool(all(a_ok.values()) and b_ok), "informative_formats": informative,
                        "loss_+3dB_bit": {m: float(loss3[m]) for m in MODS}, "ranking_+3dB_all": loss3.sort_values(ascending=False).index.tolist(),
                        "baseline_AIR_-6dB": {m: float(l6.loc[m, "AIR_ref_mean"]) for m in MODS},
                        "P5_dB": {r.modulation: (None if np.isnan(r.P_minus_P1dB_dB) else float(r.P_minus_P1dB_dB)) for r in T[T.limit_pct == 5].itertuples()},
                        "P10_dB": {r.modulation: (None if np.isnan(r.P_minus_P1dB_dB) else float(r.P_minus_P1dB_dB)) for r in T[T.limit_pct == 10].itertuples()},
                        "P5_censored_top": {r.modulation: bool(r.censored_at_grid_top) for r in T[T.limit_pct == 5].itertuples()}}
    ref = per_scale[1.0]["loss_+3dB_bit"]
    for c, v in per_scale.items():
        inf = v["informative_formats"]
        if len(inf) < 2:
            v["broad_ordering_testable"], v["broad_ordering_unchanged"] = False, None
            continue
        cur = {m: v["loss_+3dB_bit"][m] for m in inf}
        r1 = {m: ref[m] for m in inf}
        v["broad_ordering_testable"] = True
        v["most_least_affected"] = [max(cur, key=cur.get), min(cur, key=cur.get)]
        v["broad_ordering_unchanged"] = bool(max(cur, key=cur.get) == max(r1, key=r1.get) and min(cur, key=cur.get) == min(r1, key=r1.get))
        v["detailed_ranking_unchanged"] = bool(sorted(inf, key=lambda m: -cur[m]) == sorted(inf, key=lambda m: -r1[m]))
    holds = {c: v["surrogate_failure_holds"] for c, v in per_scale.items()}
    order_ok = all(v["broad_ordering_unchanged"] for v in per_scale.values() if v["broad_ordering_testable"])
    if all(holds.values()) and order_ok:
        cls = "ROBUST"
    elif holds[0.5] and holds[1.0] and holds[2.0]:
        cls = "PARTIALLY ROBUST"
    else:
        cls = "NOISE DEPENDENT"
    # prior-pass P4 rule (scales 0.25-2; full ranking at +3 dB)
    pr = {c: per_scale[c] for c in (0.25, 0.5, 1.0, 2.0)}
    ranks = {tuple(v["ranking_+3dB_all"]) for v in pr.values()}
    ph = {c: v["surrogate_failure_holds"] for c, v in pr.items()}
    prior = ("ROBUST_TO_NOISE" if all(ph.values()) and len(ranks) == 1 else ("PARTIALLY_ROBUST" if ph[0.5] and ph[1.0] and ph[2.0] else "NOISE_DEPENDENT"))
    rev, dirty = git_state()
    verdict = {"generated": now(), "code_version": rev, "code_dirty": dirty, "classification": cls, "validation_c1": validation,
               "per_scale": {str(k): v for k, v in per_scale.items()}, "prior_pass_P4_rule": prior, "bits_per_symbol": bits}
    tab = pd.DataFrame(out)
    tab.to_csv(OUT / "noise_sensitivity.csv", index=False)
    (OUT / "noise_verdict.json").write_text(json.dumps(verdict, indent=1, default=float))
    figures(tab, per_scale)
    summary(tab, verdict)
    print(json.dumps(verdict, indent=1, default=float))


def figures(tab, per_scale):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axs = plt.subplots(1, 3, figsize=(11.5, 3.4))
    lf = tab[tab.table == "loss_full"]
    for mod in MODS:
        col, mk = FMT_COL[mod]
        g = lf[(lf.modulation == mod) & (lf.Pavg_over_P1dB_dB == 3)].sort_values("noise_scale")
        axs[0].errorbar(g.noise_scale, g.loss_mean, yerr=[g.loss_mean - g.loss_ci_lo, g.loss_ci_hi - g.loss_mean], color=col, marker=mk, ms=4, capsize=2, label=mod)
        axs[1].plot(g.noise_scale, g.AIR_ref_mean, color=col, marker=mk, ms=4, label=mod)
        sn = tab[(tab.table == "effective_SNR_dB") & (tab.modulation == mod) & (tab.Pavg_over_P1dB_dB == 3)].sort_values("noise_scale")
        axs[2].plot(sn.noise_scale, sn.eff_SNR_dB, color=col, marker=mk, ms=4, label=mod)
    for ax, lab in zip(axs, ("full-model AIR loss at +3 dB (bit)", "baseline AIR at −6 dB (bit)", "effective SNR at +3 dB (dB)")):
        ax.set_xscale("log", base=2)
        ax.set_xlabel("noise scale c (× archived σ)")
        ax.set_ylabel(lab)
    axs[0].legend(fontsize=7, frameon=False)
    fig.suptitle("Part IV: full-model AIR vs noise scale (8 realizations, paired noise draw)", fontsize=9)
    fig.tight_layout()
    fig.savefig(OUT / "fig_air_vs_noise.png", dpi=160)
    plt.close(fig)
    ex = tab[tab.table == "surrogate_minus_full"]
    fig, axs = plt.subplots(2, 3, figsize=(11.5, 5.6), sharey=True)
    for i, ctl in enumerate(SURR):
        for j, p in enumerate((0, 3, 6)):
            ax = axs[i, j]
            for mod in MODS:
                col, mk = FMT_COL[mod]
                g = ex[(ex.control == ctl) & (ex.modulation == mod) & (ex.Pavg_over_P1dB_dB == p)].sort_values("noise_scale")
                ax.errorbar(g.noise_scale, g.excess_mean, yerr=[g.excess_mean - g.ci_lo, g.ci_hi - g.excess_mean], color=col, marker=mk, ms=3.5,
                            capsize=2, label=mod)
            ax.axhline(0, color="#898781", lw=0.6)
            for y in (-.1, .1):
                ax.axhline(y, color="#898781", lw=0.5, ls=":")
            ax.set_xscale("log", base=2)
            ax.set_title(f"{'all-zone static' if ctl == 'static_mh' else 'all-zone LTI+static'}, {p:+d} dB", fontsize=8)
            if j == 0:
                ax.set_ylabel("AIR(surrogate) − AIR(full) (bit)")
            if i == 1:
                ax.set_xlabel("noise scale c")
    axs[0, 0].legend(fontsize=7, frameon=False)
    fig.tight_layout()
    fig.savefig(OUT / "fig_surrogate_error_vs_noise.png", dpi=160)
    plt.close(fig)


def summary(tab, v):
    L = ["# Part IV — AIR noise sensitivity", "", f"**Classification: {v['classification']}** (rule: `00_plan.md`, section D; generated {v['generated']}, "
         f"code {v['code_version'][:7]}).", "",
         f"Validation at c = 1: {json.dumps(v['validation_c1'], default=float)}.", "",
         "| c | (a) static | (a) LTI+static | (b) both signs | failure holds | informative formats | broad ordering (most, least) unchanged | detailed ranking unchanged | loss +3 dB (CE/QPSK/16QAM/OFDM) | baseline AIR −6 dB | P5 (dB) |",
         "|---|---|---|---|---|---|---|---|---|---|---|"]
    for c, p in v["per_scale"].items():
        L.append(f"| {c} | {p['a_static_mh']} | {p['a_static_lti_mh']} | {p['b_both_signs']} | {p['surrogate_failure_holds']} | {', '.join(p['informative_formats'])} | "
                 f"{p.get('broad_ordering_unchanged')} {p.get('most_least_affected', '')} | {p.get('detailed_ranking_unchanged')} | "
                 + "/".join(f"{p['loss_+3dB_bit'][m]:.2f}" for m in MODS) + " | " + "/".join(f"{p['baseline_AIR_-6dB'][m]:.2f}" for m in MODS) + " | "
                 + "/".join("n/a" if p['P5_dB'][m] is None else f"{p['P5_dB'][m]:+.2f}{'+' if p['P5_censored_top'][m] else ''}" for m in MODS) + " |")
    ex = tab[(tab.table == "surrogate_minus_full") & tab.control.isin(SURR) & tab.Pavg_over_P1dB_dB.isin([0, 3, 6])]
    L += ["", "Surrogate − full AIR (bit), paired mean [95% CI]:", "", "| control | format | power | " + " | ".join(f"c = {c:g}" for c in SCALES) + " |",
          "|---|---|---|" + "---|" * len(SCALES)]
    for (ctl, mod, p), g in ex.groupby(["control", "modulation", "Pavg_over_P1dB_dB"]):
        g = g.set_index("noise_scale")
        L.append(f"| {ctl} | {mod} | {int(p):+d} | " + " | ".join(f"{g.loc[c, 'excess_mean']:+.2f} [{g.loc[c, 'ci_lo']:+.2f}, {g.loc[c, 'ci_hi']:+.2f}]" for c in SCALES) + " |")
    L += ["", f"Prior-pass P4 rule (scales 0.25–2, full ranking): {v['prior_pass_P4_rule']}.",
          "P5 marked '+' is censored at the top of the power grid (+6 dB). Formats near ceiling/floor are excluded from ordering statements."]
    (OUT / "noise_summary.md").write_text("\n".join(L) + "\n")


if __name__ == "__main__":
    main()
