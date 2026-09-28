"""Supplementary Figure S9: IF scan at matched signal/LO 0.358, three LO fields (prespecified; results/tqe_viability_final/00_plan.md,
section B; results/tqe_viability_final/02_resonance/resonance_scan.csv).

Rows: CW gain / small-signal; selectivity S (+-1.31 MHz, from the CW scan); X_coh of one reference sequence.
Columns: absolute IF; r2 = f_IF / (Omega_LO/4pi).
"""
import pandas as pd
import matplotlib.pyplot as plt

from figstyle import (fs, COL2, ROOT, MUTED, panel, rec, save, relpath, zero_line)

TVF = ROOT / "results" / "tqe_viability_final"
SC = pd.read_csv(TVF / "02_resonance" / "resonance_scan.csv")
SC = SC[SC.status == "ok"]
SRC = relpath(TVF / "02_resonance" / "resonance_scan.csv")
LOC = {0.5: "#2a78d6", 0.425: "#eda100", 0.35: "#e87ba4"}
LOM = {0.5: "o", 0.425: "v", 0.35: "P"}

fig, axs = plt.subplots(3, 2, figsize=(COL2, 5.6), sharex="col")
rows = [("g_CW", "CW gain / small-signal", (0.15, 1.85)), ("S", "selectivity S (±1.31 MHz)", (0, 1.3)), ("X_coh", "$X_\\mathrm{coh}$, one sequence", (-0.12, 0.62))]
letters = iter("abcdef")
for i, (col, ylab, ylim) in enumerate(rows):
    for j, (xc, xlab) in enumerate((("IF_MHz", "IF (MHz)"), ("r2", "$r_2 = f_\\mathrm{IF}/(\\Omega_\\mathrm{LO}/4\\pi)$"))):
        ax = axs[i, j]
        letter = next(letters)
        panel(ax, letter)
        if col == "X_coh":
            zero_line(ax)
        if xc == "r2":
            ax.axvline(1.0, color=MUTED, lw=0.6, ls=(0, (3, 2)), zorder=0)
            ax.axvline(2.0, color=MUTED, lw=0.6, ls=(0, (3, 2)), zorder=0)
        for lo in (0.35, 0.425, 0.5):
            s = SC[SC.A_LO_Vpm == lo].sort_values(xc)
            ax.plot(s[xc], s[col], color=LOC[lo], marker=LOM[lo], ms=2.6, mec="white", mew=0.3, lw=0.9,
                    label=f"$A_\\mathrm{{LO}}$ = {lo:g} V/m")
            rec(letter, f"{col} vs {xc}, LO {lo:g}", s[xc], s[col], source=SRC)
        ax.set_ylim(*ylim)
        if j == 0:
            ax.set_ylabel(ylab)
        if i == 2:
            ax.set_xlabel(xlab)
axs[0, 0].legend(loc="upper left", fontsize=fs(6.0), borderaxespad=0.2)
axs[0, 1].set_xlim(0.4, 4.0)
save(fig, "figS9", supplement=True)
