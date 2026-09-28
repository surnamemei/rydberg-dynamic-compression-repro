"""Supplementary Figure S4: Rydberg populations and probe-transmission baseline by case (OP1, end of the modulated segment).

(a) change of the thermally averaged Rydberg population (|3> + |4>) relative to CW (%);
(b) probe-transmission baseline shift relative to CW over the last 20 us (dB).
Source: phase_mechanism/03_phase_mechanism_results.csv (Nd = 4001, dt = 1 ns).
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from figstyle import (fs, COL2, INK2, GRID, MODEL, PM, panel, rec, save, relpath, tm)

R = pd.read_csv(PM / "03_phase_mechanism_results.csv")
main = R[(R.Nd == 4001) & (R.dt_ns == 1)].set_index("case")
cw = main.loc["CW"]
cases = [("QPSK_REF", "reference QPSK"), ("QPSK_REF_s2", "second sequence"), ("QPSK_SLOW", "slow (4 µs symbols)"),
         ("QPSK_FAST", "fast (0.5 µs symbols)"), ("QPSK_JUMP", "steps"), ("QPSK_LONGRAMP", "900 ns ramps"), ("TOGGLE", "deterministic sequence"),
         ("OFF_+0.125", "offset +0.125 MHz"), ("OFF_+1.31", "offset +1.31 MHz"), ("OFF_-1.31", "offset −1.31 MHz"),
         ("OFF_+2.62", "offset +2.62 MHz"), ("OFF_-2.62", "offset −2.62 MHz")]
ys = np.arange(len(cases))[::-1].astype(float)
ys[:7] += 0.6
ryd = {c: 100 * ((main.loc[c, "pop_rydberg_3"] + main.loc[c, "pop_rydberg_4"]) / (cw.pop_rydberg_3 + cw.pop_rydberg_4) - 1) for c, _ in cases}
pdb = {c: main.loc[c, "probe_db_last20us"] - cw.probe_db_last20us for c, _ in cases}

fig = plt.figure(figsize=(COL2, 2.5))
axa, axb = fig.subplots(1, 2, sharey=True)
for ax, vals, letter, xlabel in ((axa, ryd, "a", "Rydberg population change vs CW (%)"), (axb, pdb, "b", "probe-baseline shift vs CW (dB)")):
    panel(ax, letter, "OP1, end of the modulated segment")
    ax.axvline(0, color=INK2, lw=0.6, zorder=0)
    ax.axhline((ys[6] + ys[7]) / 2, color=GRID, lw=0.6)
    for (c, lab), y in zip(cases, ys):
        ax.plot([0, vals[c]], [y, y], color=MODEL["full"], lw=0.8, zorder=1)
        ax.plot([vals[c]], [y], marker="o", ms=3.8, color=MODEL["full"], mec="white", mew=0.4, zorder=2)
        rec(letter, f"{xlabel} {c}", c, float(vals[c]), source=relpath(PM / "03_phase_mechanism_results.csv"))
    ax.set_xlabel(xlabel)
    ax.tick_params(axis="y", length=0)
axa.set_yticks(ys)
axa.set_yticklabels([lab for _, lab in cases], fontsize=fs(6.4))
ph = [c for c, _ in cases if not c.startswith("OFF")]
axa.text(0.03, 0.72, f"phase cases: {tm(f'{min(ryd[c] for c in ph):+.1f}')}% to {tm(f'{max(ryd[c] for c in ph):+.1f}')}%",
         transform=axa.transAxes, ha="left", va="center", fontsize=fs(6.0), color=INK2)
axb.text(0.03, 0.72, f"phase cases: |shift| ≤ {max(abs(pdb[c]) for c in ph):.2f} dB", transform=axb.transAxes, ha="left", va="center", fontsize=fs(6.0), color=INK2)
rec("a", "range phase cases (%)", "min,max", [min(ryd[c] for c in ph), max(ryd[c] for c in ph)], note="ledger F_ryd")
rec("b", "max |shift| phase cases (dB)", "max", max(abs(pdb[c]) for c in ph), note="ledger F_probe")
rec("b", "max |shift| offsets (dB)", "max", max(abs(pdb[c]) for c, _ in cases if c.startswith("OFF")), note="ledger F_probe_off")

save(fig, "figS4", supplement=True)
