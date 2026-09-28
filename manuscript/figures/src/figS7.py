"""Supplementary Figure S7: internal state around an isolated +pi/2 phase step (reference configuration, OP1, 5 MHz IF;
prespecified P4 check, results/revision2/00_decision_rules.json).

(a) Output IF magnitude, IF component of the probe observable Im(rho21) and of the RF coherence rho43 (+IF), and |slow rho43|,
    each relative to its pre-step value; grey: output of the small-signal control (p4/*_traces.npz, 40_p4_internal_tracking.csv).
(b) Phase relative to the drive after the step: output and rho43 (+IF) at a_H; rho43 (+IF) in the small-signal control.
"""
import numpy as np
import matplotlib.pyplot as plt

from figstyle import (fs, COL2, ROOT, INK, MUTED, panel, rec, save, relpath, zero_line)

R2 = ROOT / "results" / "revision2"
fig = plt.figure(figsize=(COL2, 2.55))
axs = fig.subplots(1, 2)

# ------------------------------------------------------------------ (e), (f) internal traces around an isolated step
TR = np.load(R2 / "p4" / "step_aH_traces.npz")
TS = np.load(R2 / "p4" / "step_smallsignal_traces.npz")
SRC4 = relpath(R2 / "p4" / "step_aH_traces.npz")


def win(z):
    return (z["t"] >= 49.5e-6) & (z["t"] < 54.75e-6)


ax = axs[0]
panel(ax, "a", "π/2 step at $a_H$: magnitudes")
m, ms_ = win(TR), win(TS)
t, ts = TR["t"][m] * 1e6, TS["t"][ms_] * 1e6
ax.axvspan(49.85, 50.15, color="#e1e0d9", lw=0, zorder=0)
ax.plot(ts, TS["out_mag_norm"][ms_], color=MUTED, lw=1.0, label="output, small-signal step")
ax.plot(t, TR["out_mag_norm"][m], color=INK, lw=1.3, label="output")
ax.plot(t, TR["thermal|rho21 (probe)|Im part, IF"][m], color="#eda100", lw=0.9, ls=(0, (3, 1.5)), label="Im $\\rho_{21}$, IF")
ax.plot(t, TR["thermal|rho43 (RF)|+IF"][m], color="#2a78d6", lw=1.1, label="$\\rho_{43}$, IF")
ax.plot(t, TR["thermal|rho43 (RF)|slow"][m], color="#e34948", lw=1.0, ls=(0, (1.2, 1.2)), label="$|\\rho_{43}|$, slow")
for key, lab in (("out_mag_norm", "output"), ("thermal|rho21 (probe)|Im part, IF", "Im rho21 IF"), ("thermal|rho43 (RF)|+IF", "rho43 +IF"),
                 ("thermal|rho43 (RF)|slow", "rho43 slow")):
    rec("a", f"a_H {lab} (rel. pre-step)", t[::25], TR[key][m][::25], source=SRC4)
srch = (TR["t"] >= 49.8e-6) & (TR["t"] < 54.75e-6)
rec("a", "a_H output minimum (rel. pre-step), 49.8-54.75 us", "min", float(TR["out_mag_norm"][srch].min()), source=SRC4)
rec("a", "small-signal output (rel. pre-step)", ts[::25], TS["out_mag_norm"][ms_][::25], source=relpath(R2 / "p4" / "step_smallsignal_traces.npz"))
ax.set_xlim(49.5, 54.75)
ax.set_ylim(0, 3.25)
ax.set_xlabel("time (µs); step centred at 50 µs")
ax.set_ylabel("relative to pre-step value")
ax.legend(loc="upper center", fontsize=fs(5.8), borderaxespad=0.2, handlelength=1.5, ncol=2, columnspacing=0.8)

ax = axs[1]
panel(ax, "b", "phase relative to the drive")
zero_line(ax)
ax.axvspan(49.85, 50.15, color="#e1e0d9", lw=0, zorder=0)
ax.plot(ts, TS["thermal|rho43 (RF)|+IF|phase"][ms_], color=MUTED, lw=1.0, label="$\\rho_{43}$, small-signal step")
ax.plot(t, TR["out_phase_rel"][m], color=INK, lw=1.3, label="output")
ax.plot(t, TR["thermal|rho43 (RF)|+IF|phase"][m], color="#2a78d6", lw=1.1, label="$\\rho_{43}$, IF")
ok = (TR["t"] >= 49.2e-6) & (TR["t"] < 54.75e-6)
rec("b", "a_H rho43 +IF phase at the window end (rad)", "end", float(TR["thermal|rho43 (RF)|+IF|phase"][ok][-1]), source=SRC4)
rec("b", "a_H output phase (rad)", t[::25], TR["out_phase_rel"][m][::25], source=SRC4)
rec("b", "a_H rho43 +IF phase (rad)", t[::25], TR["thermal|rho43 (RF)|+IF|phase"][m][::25], source=SRC4)
rec("b", "small-signal rho43 +IF phase (rad)", ts[::25], TS["thermal|rho43 (RF)|+IF|phase"][ms_][::25], source=relpath(R2 / "p4" / "step_smallsignal_traces.npz"))
ax.set_xlim(49.5, 54.75)
ax.set_ylim(-15, 3)
ax.set_xlabel("time (µs)")
ax.set_ylabel("phase relative to drive, re pre-step (rad)")
ax.legend(loc="lower left", fontsize=fs(5.8), borderaxespad=0.2, handlelength=1.6)

save(fig, "figS7", supplement=True)
