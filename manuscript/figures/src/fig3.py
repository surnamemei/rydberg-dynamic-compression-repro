"""Figure 3: fair-waveform comparison (OP1, Nd = 4001).

(a) full-model AIR loss relative to -6 dB vs Pavg/P1dB, paired 95% CIs (regen_stage05/02_losses.csv);
    dashed: screening results (Nd = 1501) on the converged axis (05_archived_on_converged_axis.csv, +1.57 dB shift);
(b) AIR(surrogate) - AIR(full) for the two all-zone surrogates at 0/+3/+6 dB (03_excess_over_controls.csv).
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.transforms import blended_transform_factory

from figstyle import (fs, COL2, INK, INK2, GRID, MODEL, MODEL_LABEL, MODEL_MARKER, FMT, FMT_MARKER, FMT_LABEL, PP, ND, FV, EB,
                      panel, rec, save, relpath, zero_line)

R = FV / "regen_stage05"
L = pd.read_csv(R / "02_losses.csv")
S = pd.read_csv(R / "05_archived_on_converged_axis.csv")
X = pd.read_csv(R / "03_excess_over_controls.csv")
for df in (L, S, X):
    df["fmt"] = df.modulation.replace({"16QAM": "16-QAM"})
FMTS = ["CE", "QPSK", "16-QAM", "OFDM"]

fig = plt.figure(figsize=(COL2, 2.75))
gs = fig.add_gridspec(1, 2, width_ratios=[1, 1.25])
ax = fig.add_subplot(gs[0, 0])
axb = fig.add_subplot(gs[0, 1])

# ------------------------------------------------------------------ (a)
panel(ax, "a", f"full model, OP1, {ND} = 4001")
zero_line(ax)
for f in FMTS:
    g = L[L.fmt == f].sort_values("Pavg_over_P1dB_dB")
    ax.errorbar(g.Pavg_over_P1dB_dB, g.loss_mean, yerr=[g.loss_mean - g.loss_ci_lo, g.loss_ci_hi - g.loss_mean], color=FMT[f],
                marker=FMT_MARKER[f], ms=3.4, mec="white", mew=0.5, lw=1.1, label=FMT_LABEL[f], zorder=3, **EB)
    rec("a", f"loss {f}", g.Pavg_over_P1dB_dB.values, g.loss_mean.values, g.loss_ci_lo.values, g.loss_ci_hi.values, source=relpath(R / "02_losses.csv"))
    s = S[S.fmt == f].sort_values("p_converged_axis")
    ax.plot(s.p_converged_axis, s.loss_mean, color=FMT[f], lw=0.8, ls=(0, (3, 2)), zorder=2)
    rec("a", f"screening {f} (shifted)", s.p_converged_axis.values, s.loss_mean.values, source=relpath(R / "05_archived_on_converged_axis.csv"),
        note="p_converged_axis = p_archived + 1.568 dB")
ax.set_xlim(-16, 8.3)
ax.set_xticks([-15, -10, -6, -3, 0, 3, 6])
ax.set_xticklabels(["−15", "−10", "−6", "−3", "0", "+3", "+6"])
ax.set_ylim(-0.25, 1.5)
ax.set_xlabel(f"{PP} (dB)")
ax.set_ylabel("AIR loss re −6 dB (bit/symbol)")
h, lab = ax.get_legend_handles_labels()
h.append(Line2D([], [], color=INK2, lw=0.8, ls=(0, (3, 2))))
lab.append(f"screening ({ND} = 1501), +1.57 dB")
ax.legend(h, lab, loc="upper left", fontsize=fs(6.3))

# ------------------------------------------------------------------ (b)
panel(axb, "b", "CW-matched surrogates vs full model, OP1")
zero_line(axb)
powers = [0, 3, 6]
ctl = [("static_mh", "static"), ("static_lti_mh", "lti")]
ticks, tlabels = [], []
for gi, f in enumerate(FMTS):
    base = gi * 3.8
    for pi, p in enumerate(powers):
        x = base + pi
        ticks.append(x)
        tlabels.append({0: "0", 3: "+3", 6: "+6"}[p])
        for ci, (c, key) in enumerate(ctl):
            r = X[(X.fmt == f) & (X.Pavg_over_P1dB_dB == p) & (X.control == c)].iloc[0]
            xx = x + (ci - 0.5) * 0.36
            axb.errorbar([xx], [r.excess_mean], yerr=[[r.excess_mean - r.ci_lo], [r.ci_hi - r.excess_mean]], color=MODEL[key],
                         marker=MODEL_MARKER[key], ms=3.6, mec="white", mew=0.5, lw=0, zorder=3, **EB,
                         label=MODEL_LABEL[key] if (gi == 0 and pi == 0) else None)
            rec("b", f"excess {c} {f}", p, r.excess_mean, r.ci_lo, r.ci_hi, source=relpath(R / "03_excess_over_controls.csv"))
    if gi < len(FMTS) - 1:
        axb.axvline(base + 2.9, color=GRID, lw=0.6, zorder=0)
axb.set_xticks(ticks)
axb.set_xticklabels(tlabels)
axb.tick_params(axis="x", length=0, pad=2)
tr = blended_transform_factory(axb.transData, axb.transAxes)
for gi, f in enumerate(FMTS):
    axb.text(gi * 3.8 + 1, -0.13, FMT_LABEL[f], transform=tr, ha="center", va="top", fontsize=fs(7.2))
axb.text(-0.9, -0.055, f"{PP} (dB):", transform=tr, ha="right", va="top", fontsize=fs(6.3), color=INK2)
axb.set_xlim(-0.8, 3 * 3.8 + 2.8)
axb.set_ylim(-1.6, 1.7)
axb.set_ylabel("AIR(surrogate) − AIR(full model) (bit/symbol)")
axb.legend(loc="upper left", fontsize=fs(6.3))

save(fig, "fig3")
