"""Figure 2: thermal-quadrature convergence (OP1).

(a) baseband NMSE vs Nd relative to Nd = 8001 (02_nd_convergence.csv) with the prespecified tolerance 1e-3;
(b) E1dB at Nd = 1501 / 4001 / 8001 (controls_Nd*.npz).
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from figstyle import (fs, COL1, INK, INK2, MUTED, MODEL, FMT, ND, E1, RES06, panel, rec, save, relpath)

GRID_D = "#c3c2b7"

d = pd.read_csv(RES06 / "02_nd_convergence.csv")
SRC = relpath(RES06 / "02_nd_convergence.csv")
cases = [("A_lowpower_QPSK_m10", "QPSK −10 dB", FMT["QPSK"], "s", "white"),
         ("B_P1dB_OFDM_0", "OFDM 0 dB", FMT["OFDM"], "D", None),
         ("C_plus3_16QAM", "16-QAM +3 dB", FMT["16-QAM"], "^", None),
         ("D_strongest_QPSK_p6", "QPSK +6 dB", FMT["QPSK"], "s", None),
         ("E_stage05_dwell_pair", "dwell pair", MUTED, "o", None)]

fig = plt.figure(figsize=(COL1, 3.95))
gs = fig.add_gridspec(2, 1, height_ratios=[2.45, 0.95])
ax = fig.add_subplot(gs[0, 0])
axb = fig.add_subplot(gs[1, 0])

# ------------------------------------------------------------------ (a)
panel(ax, "a", "OP1")
ax.axvspan(1400, 1610, color="#f3f2ee", zorder=0, lw=0)
ax.axvline(4001, color=MODEL["full"], lw=0.6, zorder=0)
for key, lab, col, mk, mfc in cases:
    g = d[(d.case == key) & (d.Nd < 8001)].sort_values("Nd")
    ax.plot(g.Nd, g.baseband_NMSE_vs_top, color=col, marker=mk, ms=3.2, mfc=mfc or col, mec=col if mfc else "white", mew=0.7 if mfc else 0.4,
            lw=0.9, label=lab)
    rec("a", f"NMSE {lab}", g.Nd.values, g.baseband_NMSE_vs_top.values, source=SRC)
ax.axhline(1e-3, color=INK, lw=0.7, ls=(0, (1.2, 1.4)))
ax.text(2350, 1.3e-3, "prespecified tolerance", fontsize=fs(6.3), color=INK, va="bottom")
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlim(700, 8600)
ax.set_ylim(3e-7, 0.25)
ax.set_xticks([751, 1501, 4001, 8001])
ax.set_xticklabels(["751", "1501", "4001", "8001"])
ax.xaxis.set_minor_formatter(plt.NullFormatter())
ax.set_xlabel(ND)
ax.set_ylabel(f"baseband NMSE vs {ND} = 8001")
ax.text(1500, 0.14, "screening\n(fails)", ha="center", va="top", fontsize=fs(6.3), color=INK2)
ax.axvline(8001, color=GRID_D, lw=0.6, zorder=0)
ax.text(4001 * 0.96, 0.14, "decisive", ha="right", va="top", fontsize=fs(6.3), color=MODEL["full"])
ax.text(8001 * 0.96, 0.14, "reference", ha="right", va="top", fontsize=fs(6.3), color=INK2)
ax.legend(loc="lower left", fontsize=fs(6.0), handlelength=1.4, ncol=1, borderaxespad=0.2)

# ------------------------------------------------------------------ (b)
panel(axb, "b", "OP1")
nds = [1501, 4001, 8001]
e = []
for n in nds:
    z = np.load(RES06 / "artifacts" / f"controls_Nd{n}.npz")
    e.append(float(z["e1db"]))
    rec("b", f"E1dB Nd={n}", n, float(z["e1db"]), source=relpath(RES06 / "artifacts" / f"controls_Nd{n}.npz"))
for i, (n, v) in enumerate(zip(nds, e)):
    fc = "white" if n == 1501 else MODEL["full"]
    axb.plot([v], [i], marker="o", ms=4.6, mfc=fc, mec=MODEL["full"], mew=0.8, zorder=2)
    axb.annotate(f"{v:.4f} V/m", (v, i), xytext=(6, 0), textcoords="offset points", ha="left", va="center", fontsize=fs(6.3))
axb.set_yticks([0, 1, 2])
axb.set_yticklabels([f"{ND} = {n}" for n in nds])
axb.set_ylim(-0.6, 2.6)
axb.invert_yaxis()
axb.set_xlim(0.068, 0.096)
axb.set_xlabel(f"{E1} (V/m)")
axb.tick_params(axis="y", length=0)

save(fig, "fig2")
