"""Supplementary Figure S6: phase-effect checks of the revision (prespecified P3, P4, P5).

(a) X of the reference QPSK case vs level relative to P1dB at OP1 and OP2 (3 zero-drift sequences, mean and 95% CI),
    with the memoryless frequency-dependent static surrogate g(a, f) (P5) and the 0.05 materiality threshold.
(b) Large-signal gain at a constant frequency offset relative to on-IF, R = g(df)/g_CW, at a = 2.4468 E1dB: sweep at
    OP1 and OP2; five offsets at OP3 (A_LO = 0.425 V/m).
(c) X of the reference case (3 sequences, mean and 95% CI) vs LO field at OP2, OP3 and OP1.
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from figstyle import (fs, COL2, ROOT, INK, INK2, EB, panel, rec, save, relpath, zero_line)

RV = ROOT / "results" / "revision"
S3 = pd.read_csv(RV / "21_p3_summary.csv")
T3 = pd.read_csv(RV / "43_p5_T3.csv")
SW = pd.read_csv(RV / "30_p4_offset_sweep.csv")
O3 = pd.read_csv(RV / "31_p4_op3_offsets.csv")
OPS = pd.read_csv(RV / "32_p4_operating_points.csv").set_index("op")
OPC = {"OP1": "#2a78d6", "OP2": "#4a3aa7", "OP3": "#008300"}
MK = {"OP1": "o", "OP2": "s", "OP3": "D"}

fig = plt.figure(figsize=(COL2, 2.5))
axs = fig.subplots(1, 3, gridspec_kw=dict(wspace=0.36, width_ratios=[1.05, 1.15, 0.8]))

ax = axs[0]
panel(ax, "a", "level sweep, 5 MHz IF")
zero_line(ax)
ax.axhline(0.05, color=INK, lw=0.7, ls=(0, (1.2, 1.4)), zorder=1)
ax.text(-6.8, 0.06, "0.05", ha="left", va="bottom", fontsize=fs(5.8), color=INK)
for op in ("OP1", "OP2"):
    g = S3[S3.op == op].sort_values("level_dB")
    ax.errorbar(g.level_dB, g.X_mean, yerr=[g.X_mean - g.ci_lo, g.ci_hi - g.X_mean], color=OPC[op], marker=MK[op], ms=3.6,
                mec="white", mew=0.5, lw=1.1, label=f"{op}, full model", **EB)
    rec("a", f"X REF {op}", g.level_dB.values, g.X_mean.values, g.ci_lo.values, g.ci_hi.values, source=relpath(RV / "21_p3_summary.csv"))
    t = T3[T3.op == op].sort_values("level_dB")
    ax.plot(t.level_dB, t.X_sur_mean, color=OPC[op], marker=MK[op], ms=3.4, mfc="white", mew=0.8, lw=0.8, ls=(0, (3, 2)),
            label=f"{op}, $g(a,f)$ surrogate")
    rec("a", f"X REF g(a,f) surrogate {op}", t.level_dB.values, t.X_sur_mean.values, source=relpath(RV / "43_p5_T3.csv"))
ax.set_xlim(-7, 9)
ax.set_xticks([-6, -3, 0, 3, 6, 8])
ax.set_ylim(-0.1, 0.8)
ax.set_xlabel("level re $P_\\mathrm{1dB}$ (dB)")
ax.set_ylabel("$X = 1 - g/g_\\mathrm{CW}$, reference sequence")
hh, ll = ax.get_legend_handles_labels()
order = [ll.index(t) for t in ("OP1, full model", "OP1, $g(a,f)$ surrogate", "OP2, full model", "OP2, $g(a,f)$ surrogate")]
ax.legend([hh[k] for k in order], [ll[k] for k in order], loc="upper left", fontsize=fs(5.8), borderaxespad=0.2)

ax = axs[1]
panel(ax, "b", "offset sweep, 5 MHz IF")
ax.axhline(1.0, color=INK2, lw=0.6, zorder=0)
for op in ("OP1", "OP2"):
    g = SW[SW.op == op].sort_values("offset_MHz")
    ax.plot(g.offset_MHz, g.R, color=OPC[op], marker=MK[op], ms=3.0, mec="white", mew=0.4, lw=1.0, label=f"{op}")
    rec("b", f"R {op}", g.offset_MHz.values, g.R.values, source=relpath(RV / "30_p4_offset_sweep.csv"))
g = O3.sort_values("offset_MHz")
ax.plot(g.offset_MHz, g.R, color=OPC["OP3"], marker=MK["OP3"], ms=3.6, mec="white", mew=0.5, ls="none", label="OP3 (five offsets)")
rec("b", "R OP3", g.offset_MHz.values, g.R.values, source=relpath(RV / "31_p4_op3_offsets.csv"))
ax.set_xlim(-3.3, 3.3)
ax.set_xticks([-3, -2, -1, 0, 1, 2, 3])
ax.set_ylim(0, 3.6)
ax.set_xlabel("constant frequency offset (MHz)")
ax.set_ylabel("$R = g(\\delta f)/g_\\mathrm{CW}$")
ax.legend(loc="upper left", fontsize=fs(5.8), borderaxespad=0.2)

ax = axs[2]
panel(ax, "c")
# descriptor right-aligned: centred, it collided with the narrow panel's "(c)" label (layout only)
ax.set_title("LO field, 5 MHz IF", loc="right", fontsize=8, pad=3.5, color=INK)
for op in ("OP2", "OP3", "OP1"):
    r = OPS.loc[op]
    m, lo, hi = r.X_REF_mean_3seeds, r.X_REF_ci_lo, r.X_REF_ci_hi
    ax.errorbar([r.A_LO_Vpm], [m], yerr=[[m - lo], [hi - m]], color=OPC[op], marker=MK[op], ms=4.0, mec="white", mew=0.5, ls="none", **EB)
    ax.text(r.A_LO_Vpm, hi + 0.03, op, ha="center", va="bottom", fontsize=fs(6.3), color=INK)
    rec("c", f"X REF 3 sequences {op}", r.A_LO_Vpm, m, lo, hi, source=relpath(RV / "32_p4_operating_points.csv"))
ax.set_xlim(0.32, 0.53)
ax.set_xticks([0.35, 0.425, 0.5])
ax.set_xticklabels(["0.35", "0.425", "0.5"])
ax.set_ylim(0, 0.7)
ax.set_xlabel("$A_\\mathrm{LO}$ (V/m)")
ax.set_ylabel("$X$, reference sequence")

save(fig, "figS6", supplement=True)
