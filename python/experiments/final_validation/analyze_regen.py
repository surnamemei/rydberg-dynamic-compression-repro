"""Analysis of the converged Stage-0.5 regeneration (no hypothesis test; mechanical recomputation)."""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from fv_common import FV, ROOT

OUT = FV / "regen_stage05"
ARCH = ROOT / "results" / "p1db_waveform_stage05" / "02_all_results.csv"
SHIFT_DB = 20 * np.log10(0.0877075 / 0.0732189)  # archived axis -> converged axis (+1.57 dB)
MODS = ("CE", "QPSK", "16QAM", "OFDM")
C = {"CE": "#2a78d6", "QPSK": "#eb6834", "16QAM": "#1baf7a", "OFDM": "#eda100"}
plt.rcParams.update({"axes.grid": True, "grid.alpha": .25, "axes.spines.top": False, "axes.spines.right": False, "font.size": 9})


def ci(x):
    x = np.asarray(x, float)
    h = stats.t.ppf(.975, len(x) - 1) * x.std(ddof=1) / np.sqrt(len(x)) if len(x) > 1 else np.nan
    return x.mean(), x.mean() - h, x.mean() + h, len(x)


def threshold(curve_p, curve_rel, lim):
    """Highest power on the (monotone-interpolated) mean relative-loss curve with loss <= lim."""
    p, r = np.asarray(curve_p), np.asarray(curve_rel)
    ok = np.flatnonzero(r <= lim)
    if not len(ok):
        return np.nan
    i = ok[-1]
    if i == len(p) - 1:
        return float(p[-1])
    return float(p[i] + (lim - r[i]) * (p[i + 1] - p[i]) / (r[i + 1] - r[i]))


def main():
    d = pd.read_csv(OUT / "rows.csv").drop_duplicates(["job_id", "model"], keep="last")
    w = d.pivot_table(index=["modulation", "seed", "Pavg_over_P1dB_dB"], columns="model", values="AIR_native")
    loss_rows, exc_rows, thr_rows = [], [], []
    for mod in MODS:
        wm = w.loc[mod]
        for p in sorted(wm.index.get_level_values(1).unique()):
            seeds = [s_ for s_ in wm.index.get_level_values(0).unique() if (s_, p) in wm.index and (s_, -6) in wm.index]
            ref = np.array([wm.loc[(s_, -6), "atomic"] for s_ in seeds])
            full = np.array([wm.loc[(s_, p), "atomic"] for s_ in seeds])
            m, lo, hi, n = ci(ref - full)
            rm, rlo, rhi, _ = ci(100 * (ref - full) / ref)
            loss_rows.append({"modulation": mod, "Pavg_over_P1dB_dB": p, "AIR_full_mean": full.mean(), "loss_mean": m, "loss_ci_lo": lo, "loss_ci_hi": hi,
                              "rel_loss_pct": rm, "rel_loss_ci_lo": rlo, "rel_loss_ci_hi": rhi, "n": n})
            for ctrl in ("static", "static_lti", "static_mh", "static_lti_mh", "linear"):
                e = np.array([wm.loc[(s_, p), ctrl] - wm.loc[(s_, p), "atomic"] for s_ in seeds])
                em, elo, ehi, _ = ci(e)
                exc_rows.append({"modulation": mod, "Pavg_over_P1dB_dB": p, "control": ctrl, "excess_mean": em, "ci_lo": elo, "ci_hi": ehi, "n": len(e)})
        pp = [-6, -3, 0, 3, 6]
        seeds8 = [s_ for s_ in wm.index.get_level_values(0).unique() if all((s_, q) in wm.index for q in pp)]
        R = np.array([[100 * (wm.loc[(s_, -6), "atomic"] - wm.loc[(s_, q), "atomic"]) / wm.loc[(s_, -6), "atomic"] for q in pp] for s_ in seeds8])
        rng = np.random.default_rng(0)
        for lim in (5, 10):
            est = threshold(pp, R.mean(axis=0), lim)
            boots = [threshold(pp, R[rng.integers(0, len(R), len(R))].mean(axis=0), lim) for _ in range(2000)]
            thr_rows.append({"modulation": mod, "limit_pct": lim, "P_minus_P1dB_dB": est, "boot_lo": np.nanpercentile(boots, 2.5),
                             "boot_hi": np.nanpercentile(boots, 97.5), "n_seeds": len(R)})
    L, E, T = pd.DataFrame(loss_rows), pd.DataFrame(exc_rows), pd.DataFrame(thr_rows)
    L.to_csv(OUT / "02_losses.csv", index=False)
    E.to_csv(OUT / "03_excess_over_controls.csv", index=False)
    T.to_csv(OUT / "04_thresholds.csv", index=False)

    # archived Stage-0.5 (Nd=1501) on its own axis, shifted to the converged axis for comparison
    a = pd.read_csv(ARCH)
    a = a[a.config_version.isin(["fair_v2", "fair_v4_randomized_balanced_ofdm_training"]) & ~((a.modulation == "OFDM") & (a.config_version == "fair_v2"))]
    a = a[a.model == "atomic"].drop_duplicates(["job_id"])
    aw = a.pivot_table(index=["modulation", "seed"], columns="Pavg_over_P1dB_dB", values="AIR_native")
    arch_rows = []
    for mod in MODS:
        am = aw.loc[mod]
        for p in [q for q in am.columns if q >= -6]:
            ok = am[[-6]].join(am[[p]], rsuffix="_p").dropna() if p != -6 else am[[-6]].dropna()
            loss = 0.0 if p == -6 else float((ok.iloc[:, 0] - ok.iloc[:, 1]).mean())
            arch_rows.append({"modulation": mod, "p_archived_axis": p, "p_converged_axis": p + SHIFT_DB, "loss_mean": loss, "n": len(ok)})
    A = pd.DataFrame(arch_rows)
    A.to_csv(OUT / "05_archived_on_converged_axis.csv", index=False)

    fig, ax = plt.subplots(1, 2, figsize=(12, 4.2))
    for mod in MODS:
        g = L[L.modulation == mod].sort_values("Pavg_over_P1dB_dB")
        ax[0].errorbar(g.Pavg_over_P1dB_dB, g.loss_mean, yerr=[g.loss_mean - g.loss_ci_lo, g.loss_ci_hi - g.loss_mean], fmt="o-", capsize=2, color=C[mod], label=f"{mod} (Nd 4001)")
        h = A[A.modulation == mod].sort_values("p_converged_axis")
        ax[0].plot(h.p_converged_axis, h.loss_mean, "x--", color=C[mod], lw=.9, alpha=.8)
    ax[0].plot([], [], "kx--", lw=.9, label="archived Nd 1501, shifted +1.57 dB")
    ax[0].set_xlabel("Pavg / P1dB (dB, converged E1dB = 0.0732 V/m)")
    ax[0].set_ylabel("AIR loss vs −6 dB point (bit/symbol)")
    ax[0].legend(fontsize=7, frameon=False)
    ax[0].set_title("full-model AIR loss, 8 seeds (4 at −15/−10 dB), paired 95% CI", fontsize=9)
    ex = E[(E.Pavg_over_P1dB_dB == 3)]
    ctrls = [("static_mh", "all-zone static"), ("static_lti_mh", "all-zone LTI+static"), ("static", "F1 static"), ("static_lti", "F1 LTI+static")]
    for j, mod in enumerate(MODS):
        for k, (c, lab) in enumerate(ctrls):
            r = ex[(ex.modulation == mod) & (ex.control == c)].iloc[0]
            x = j + (k - 1.5) * .2
            ax[1].bar(x, r.excess_mean, width=.18, color=["#eb6834", "#1baf7a", "#f3b08f", "#8fd6bd"][k], label=lab if j == 0 else None)
            ax[1].errorbar(x, r.excess_mean, yerr=[[r.excess_mean - r.ci_lo], [r.ci_hi - r.excess_mean]], fmt="none", ecolor="k", capsize=2, lw=.8)
    ax[1].axhline(0, color="k", lw=.7)
    ax[1].set_xticks(range(4))
    ax[1].set_xticklabels(MODS)
    ax[1].set_ylabel("AIR(control) − AIR(full) at +3 dB (bit/symbol)")
    ax[1].legend(fontsize=7, frameon=False)
    ax[1].set_title("excess degradation of the full model over CW-derived controls", fontsize=9)
    fig.suptitle("Converged regeneration of the fair-waveform communication result (candidate main Fig. 2)")
    fig.tight_layout()
    fig.savefig(OUT / "fig_regen_air_loss.png", dpi=170)
    pd.set_option("display.width", 220)
    print(L.round(3).to_string())
    print(E[E.Pavg_over_P1dB_dB.isin([0, 3, 6])].pivot_table(index=["modulation", "Pavg_over_P1dB_dB"], columns="control", values="excess_mean").round(3).to_string())
    print(T.round(2).to_string())


if __name__ == "__main__":
    main()
