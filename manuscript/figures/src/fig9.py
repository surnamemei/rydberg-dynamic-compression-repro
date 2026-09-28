"""Figure 9: fair cross-configuration comparison and LO-scaled IF dependence (prespecified; results/tqe_viability_final/00_plan.md).

(a) CW fundamental gain re the weak-tone slope vs signal/LO ratio (signal <= LO) at IF 5/10/15 MHz (standard probe) and with the
    weaker probe (CW tables of results/revision2/tables).
(b) X of the reference QPSK case at matched signal/LO 0.1-0.4, six sequences (mean and 95% CI) (01_matching/matching_results.csv).
(c) X at matched CW gain: 5 and 15 MHz (the only compression-matchable pair) and the reference IF at A_LO 0.425 and 0.35 V/m.
(d) CW gain vs r2 = f_IF / (Omega_LO/4pi) at three LO fields, signal/LO 0.358 (02_resonance/resonance_scan.csv).
(e) X_coh of one reference sequence vs r2 at the same three LO fields.
"""
import json

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from figstyle import (fs, COL2, ROOT, INK2, MUTED, EB, panel, rec, save, relpath, zero_line)

TVF = ROOT / "results" / "tqe_viability_final"
MR = pd.read_csv(TVF / "01_matching" / "matching_results.csv")
RUN = MR[MR.status == "run"]
SC = pd.read_csv(TVF / "02_resonance" / "resonance_scan.csv")
SC = SC[SC.status == "ok"]
R2 = ROOT / "results" / "revision2"
CONF = {"C1": ("5 MHz (reference)", "#2a78d6", "o"), "C2": ("10 MHz", "#4a3aa7", "s"), "C3": ("15 MHz", "#008300", "D"),
        "C4": ("5 MHz, weaker probe", "#e34948", "^"), "C5": ("5 MHz, $A_\\mathrm{LO}$ 0.425 V/m", "#eda100", "v"),
        "C6": ("5 MHz, $A_\\mathrm{LO}$ 0.35 V/m", "#e87ba4", "P")}
TAB = {"C1": R2 / "tables" / "OP1_IF5.npz", "C2": R2 / "tables" / "OP1_IF10_ext.npz", "C3": R2 / "tables" / "OP1_IF15_ext.npz",
       "C4": R2 / "tables" / "OP1_IF5_weakprobe_ext.npz"}
LOC = {0.5: CONF["C1"][1], 0.425: CONF["C5"][1], 0.35: CONF["C6"][1]}
LOM = {0.5: "o", 0.425: "v", 0.35: "P"}
SRC = relpath(TVF / "01_matching" / "matching_results.csv")
SRC2 = relpath(TVF / "02_resonance" / "resonance_scan.csv")

fig = plt.figure(figsize=(COL2, 5.1))
gs = fig.add_gridspec(2, 6)
axs = [fig.add_subplot(gs[0, 0:2]), fig.add_subplot(gs[0, 2:4]), fig.add_subplot(gs[0, 4:6]),
       fig.add_subplot(gs[1, 0:3]), fig.add_subplot(gs[1, 3:6])]

# ------------------------------------------------------------------ (a) CW response inside the superheterodyne regime
ax = axs[0]
panel(ax, "a", "CW response, signal ≤ LO")
zero_line(ax)
for k in ("C1", "C2", "C3", "C4"):
    lab, c, m = CONF[k]
    t = np.load(TAB[k])
    a, g = t["amps"][1:], t["gain_db"][1:]
    keep = a <= 0.5 + 1e-12
    ax.plot(a[keep] / 0.5, g[keep], color=c, lw=1.1, label=lab)
    rec("a", f"CW gain dB {k}", a[keep] / 0.5, g[keep], source=relpath(TAB[k]))
for r in (0.2, 0.358, 0.4):
    ax.axvline(r, color=MUTED, lw=0.6, ls=(0, (1.2, 1.4)), zorder=0)
ax.set_xscale("log")
ax.set_xlim(0.009, 1.05)
ax.set_xticks([0.01, 0.1, 1])
ax.set_xticklabels(["0.01", "0.1", "1"])
ax.set_ylim(-15, 7)
ax.set_xlabel("signal / LO field")
ax.set_ylabel("CW gain re weak-tone slope (dB)")
ax.legend(loc="lower left", fontsize=fs(6.0), borderaxespad=0.2, handlelength=1.4)

# ------------------------------------------------------------------ (b) matched signal/LO
ax = axs[1]
panel(ax, "b", "matched signal/LO")
zero_line(ax)
for k in ("C1", "C2", "C3", "C4"):
    lab, c, m = CONF[k]
    s = RUN[(RUN.configuration == k) & (RUN.matching_rule == "matched_signal_to_LO") & (RUN.target <= 0.4)].sort_values("target")
    ax.errorbar(s.target, s.X, yerr=[s.X - s.X_ci_lo, s.X_ci_hi - s.X], color=c, marker=m, ms=3.8, mec="white", mew=0.5, lw=1.0, **EB)
    rec("b", f"X {k}", s.target, s.X, s.X_ci_lo, s.X_ci_hi, source=SRC)
ax.set_xlim(0.07, 0.43)
ax.set_xticks([0.1, 0.2, 0.3, 0.4])
ax.set_ylim(-0.1, 0.6)
ax.set_xlabel("matched signal / LO")
ax.set_ylabel("X, reference phase case")

# ------------------------------------------------------------------ (c) matched CW gain
ax = axs[2]
panel(ax, "c", "matched CW gain")
zero_line(ax)
for k in ("C1", "C3", "C5", "C6"):
    lab, c, m = CONF[k]
    s = RUN[(RUN.configuration == k) & (RUN.matching_rule == "matched_CW_gain")].sort_values("target")
    ax.errorbar(s.target, s.X, yerr=[s.X - s.X_ci_lo, s.X_ci_hi - s.X], color=c, marker=m, ms=3.8, mec="white", mew=0.5, lw=1.0,
                label=lab if k in ("C5", "C6") else None, **EB)
    rec("c", f"X {k}", s.target, s.X, s.X_ci_lo, s.X_ci_hi, source=SRC)
ax.set_xlim(1.0, 0.35)
ax.set_ylim(-0.1, 0.65)
ax.set_xlabel("matched CW gain / small-signal")
ax.set_ylabel("X, reference phase case")
ax.legend(loc="lower right", fontsize=fs(6.0), borderaxespad=0.2, handlelength=1.2)

# ------------------------------------------------------------------ (d, e) LO-scaled IF scan
for ax, col, letter, title, ylab, ylim in ((axs[3], "g_CW", "d", "IF scan: CW gain", "CW gain / small-signal", (0.15, 1.85)),
                                          (axs[4], "X_coh", "e", "IF scan: phase effect", "$X_\\mathrm{coh}$, one sequence", (-0.12, 0.62))):
    panel(ax, letter, title)
    if col == "X_coh":
        zero_line(ax)
    ax.axvline(1.0, color=MUTED, lw=0.6, ls=(0, (3, 2)), zorder=0)
    ax.axvline(2.0, color=MUTED, lw=0.6, ls=(0, (3, 2)), zorder=0)
    for lo in (0.35, 0.425, 0.5):
        s = SC[SC.A_LO_Vpm == lo].sort_values("r2")
        ax.plot(s.r2, s[col], color=LOC[lo], marker=LOM[lo], ms=2.6, mec="white", mew=0.3, lw=0.9,
                label=f"$A_\\mathrm{{LO}}$ = {lo:g} V/m ($\\Omega_\\mathrm{{LO}}/2\\pi$ = {s.LO_Rabi_MHz.iloc[0]:.2f} MHz)")
        rec(letter, f"{col} vs r2, LO {lo:g}", s.r2, s[col], source=SRC2)
    ax.set_xlim(0.4, 4.0)
    ax.set_ylim(*ylim)
    ax.set_xlabel("$r_2 = f_\\mathrm{IF}/(\\Omega_\\mathrm{LO}/4\\pi)$")
    ax.set_ylabel(ylab)
axs[3].legend(loc="lower right", fontsize=fs(6.0), borderaxespad=0.2, handlelength=1.4)
save(fig, "fig9")
