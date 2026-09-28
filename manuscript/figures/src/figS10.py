"""Supplementary Figure S10: noise sensitivity of the Section 3 comparison (prespecified; results/tqe_viability_final/00_plan.md,
section D; results/tqe_viability_final/04_noise/noise_sensitivity.csv).

(a) Full-model AIR loss at +3 dB vs noise scale c (x the archived sigma), by format (8 realizations, 95% CI).
(b) AIR(all-zone static) - AIR(full) at +3 dB vs c, by format (paired, 95% CI).
(c) AIR(all-zone LTI+static) - AIR(full) at +3 dB vs c, by format.
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from figstyle import (fs, COL2, ROOT, INK2, FMT, FMT_MARKER, FMT_LABEL, EB, panel, rec, save, relpath, zero_line)

TVF = ROOT / "results" / "tqe_viability_final"
T = pd.read_csv(TVF / "04_noise" / "noise_sensitivity.csv")
SRC = relpath(TVF / "04_noise" / "noise_sensitivity.csv")
KEY = {"CE": "CE", "QPSK": "QPSK", "16QAM": "16-QAM", "OFDM": "OFDM"}

fig, axs = plt.subplots(1, 3, figsize=(COL2, 2.5))
lf = T[(T.table == "loss_full") & (T.Pavg_over_P1dB_dB == 3)]
ax = axs[0]
panel(ax, "a", "full model, +3 dB")
zero_line(ax)
for mod, k in KEY.items():
    g = lf[lf.modulation == mod].sort_values("noise_scale")
    ax.errorbar(g.noise_scale, g.loss_mean, yerr=[g.loss_mean - g.loss_ci_lo, g.loss_ci_hi - g.loss_mean], color=FMT[k], marker=FMT_MARKER[k],
                ms=3.6, mec="white", mew=0.5, lw=1.0, label=FMT_LABEL[k], **EB)
    rec("a", f"loss +3 dB {mod}", g.noise_scale, g.loss_mean, g.loss_ci_lo, g.loss_ci_hi, source=SRC)
ax.set_ylabel("AIR loss re −6 dB (bit)")
ax.legend(loc="upper center", fontsize=fs(6.0), borderaxespad=0.2, ncol=2, columnspacing=0.8, handlelength=1.2)
ex = T[(T.table == "surrogate_minus_full") & (T.Pavg_over_P1dB_dB == 3)]
for ax, ctl, letter, title in ((axs[1], "static_mh", "b", "all-zone static, +3 dB"), (axs[2], "static_lti_mh", "c", "all-zone LTI+static, +3 dB")):
    panel(ax, letter, title)
    zero_line(ax)
    for y in (-0.1, 0.1):
        ax.axhline(y, color=INK2, lw=0.5, ls=(0, (1.2, 1.4)), zorder=0)
    for mod, k in KEY.items():
        g = ex[(ex.control == ctl) & (ex.modulation == mod)].sort_values("noise_scale")
        ax.errorbar(g.noise_scale, g.excess_mean, yerr=[g.excess_mean - g.ci_lo, g.ci_hi - g.excess_mean], color=FMT[k], marker=FMT_MARKER[k],
                    ms=3.6, mec="white", mew=0.5, lw=1.0, **EB)
        rec(letter, f"{ctl} minus full, +3 dB, {mod}", g.noise_scale, g.excess_mean, g.ci_lo, g.ci_hi, source=SRC)
    ax.set_ylim(-1.1, 1.1)
    ax.set_ylabel("AIR(surrogate) − AIR(full) (bit)")
for ax in axs:
    ax.set_xscale("log", base=2)
    ax.set_xticks([0.25, 0.5, 1, 2, 4])
    ax.set_xticklabels(["0.25", "0.5", "1", "2", "4"])
    ax.set_xlim(0.2, 5)
    ax.set_xlabel("noise scale c (× nominal σ)")
axs[0].set_ylim(-0.25, 1.3)
save(fig, "figS10", supplement=True)
