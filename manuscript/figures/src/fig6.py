"""Figure 6: microsecond memory of the large-signal response (OP1, Nd = 4001).

(a) normalized step responses e(t) = (y - y_final)/(y_before - y_final) (tau_atom_traces.npz): probe transmission after a
    +1% LO step; IF envelope after a linear-regime signal step (0.05 -> 0.10 x 0.0877 V/m) and after a +5% signal step at
    0.0877 V/m; 1/e times from 03_tau_atom_results.csv;
(b) gain |IF envelope| / (amplitude x small-signal gain) during and after an isolated 60 us plateau, unmodulated carrier,
    +3 dB levels (plateau_gain_traces.npz), with CW-static levels and window medians from 12_plateau_steady_state.csv.
    Samples within the 200 ns demodulation boxcar around each level change (-100 to +200 ns) are omitted; the stored
    window metrics likewise exclude the first IF period after each change.
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from figstyle import (fs, COL2, INK, INK2, MUTED, MODEL, ND, TAU, RES06, panel, rec, save, relpath, zero_line)

Z = np.load(RES06 / "tau_atom_traces.npz")
T = pd.read_csv(RES06 / "03_tau_atom_results.csv")
T = T[T.Nd == 4001].set_index("experiment")
PL = np.load(RES06 / "plateau_gain_traces.npz")
PS = pd.read_csv(RES06 / "12_plateau_steady_state.csv").set_index("experiment")
FS_TXT = fs(6.5)

fig = plt.figure(figsize=(COL2, 2.6))
gs = fig.add_gridspec(1, 2, width_ratios=[1, 1.3])
ax = fig.add_subplot(gs[0, 0])
axb = fig.add_subplot(gs[0, 1])

# ------------------------------------------------------------------ (a)
panel(ax, "a", f"step responses, OP1, {ND} = 4001")
zero_line(ax)
for yv in (np.exp(-1), -np.exp(-1)):
    ax.axhline(yv, color=MUTED, lw=0.5, ls=(0, (1.2, 1.4)), zorder=0)
ax.text(14.9, np.exp(-1) + 0.02, "±1/e", ha="right", va="bottom", fontsize=fs(6.2), color=INK2)
series = [("LO_step_+1pct", INK, 1.0, "+1% LO step (probe transmission)"),
          ("IF_step_small_0.05to0.10E1", MUTED, 1.0, "linear-regime signal step (IF envelope)"),
          ("IF_step_P1dB_1.00to1.05E1", MODEL["full"], 1.1, "+5% signal step at 0.0877 V/m (IF envelope)")]
for key, col, lw, lab in series:
    t = Z[f"4001__{key}__t"] * 1e6
    e = Z[f"4001__{key}__e"]
    w = t <= 15.0
    ax.plot(t[w], e[w], color=col, lw=lw, label=lab)
    rec("a", f"e(t) {key}", "tau_1e_us", float(T.loc[key, "tau_1e_s"]) * 1e6, source=relpath(RES06 / "03_tau_atom_results.csv"))
    rec("a", f"e(t) {key} extremes", "min,max", [float(e.min()), float(e.max())], source=relpath(RES06 / "tau_atom_traces.npz"))
t_atom = float(T.loc["IF_step_P1dB_1.00to1.05E1", "tau_1e_s"]) * 1e6
ax.annotate(f"{TAU} = {t_atom:.3f} µs\n(1/e time of the +5% step)", xy=(t_atom, -0.37), xytext=(4.3, -0.62), fontsize=FS_TXT,
            arrowprops=dict(arrowstyle="-", lw=0.5, color=INK2, shrinkA=0, shrinkB=1))
ax.set_xlim(0, 15)
ax.set_ylim(-0.9, 1.05)
ax.set_xlabel("time after step (µs)")
ax.set_ylabel("normalized excess response e(t)")
ax.legend(loc="upper right", fontsize=fs(6.0), handlelength=1.5)

# ------------------------------------------------------------------ (b)
panel(axb, "b", "60 µs plateau, unmodulated carrier, +3 dB levels, OP1")
g = PL["plateau_p+3_cw"].astype(float).copy()
dt = float(PL["dt_s"]) * 1e6
warm, hold = float(PL["warm_s"]) * 1e6, float(PL["hold_s"]) * 1e6
t = np.arange(len(g)) * dt - warm
mask = ((t >= -0.1) & (t < 0.2)) | ((t >= hold - 0.1) & (t < hold + 0.2))
g[mask] = np.nan
sel = (t >= -12) & (t <= hold + 50)
axb.axvspan(0, hold, color="#f3f2ee", lw=0, zorder=0)
r = PS.loc["plateau_p+3_cw"]
axb.plot([-12, 0], [r.cw_static_gain_aL] * 2, color=INK2, lw=0.8, ls=(0, (3, 2)), zorder=1)
axb.plot([0, hold], [r.cw_static_gain_aH] * 2, color=INK2, lw=0.8, ls=(0, (3, 2)), zorder=1, label="CW-static gain at each level")
axb.plot([hold, hold + 50], [r.cw_static_gain_aL] * 2, color=INK2, lw=0.8, ls=(0, (3, 2)), zorder=1)
axb.plot(t[sel], g[sel], color=MODEL["full"], lw=0.9, label="full model", zorder=2)
rec("b", "trace", "min,max (masked)", [float(np.nanmin(g[sel])), float(np.nanmax(g[sel]))], source=relpath(RES06 / "plateau_gain_traces.npz"))
marks = [(0.2, 1.0, r.gain_plateau_first_us, "0.56 (first µs)", (3.0, 0.72)),
         (hold - 10, hold - 0.2, r.gain_plateau_last10us, "0.43 (last 10 µs)", (hold - 25, 0.24)),
         (hold + 0.2, hold + 1.0, r.gain_low_first_us_after, "0.36 (first µs after)", (hold + 3.5, 0.16)),
         (hold + 1.0, hold + 5.0, r.gain_low_1_to_5us_after, "0.60 (1–5 µs after)", (hold + 7.5, 0.47))]
for x0, x1, yv, lab, xyt in marks:
    axb.plot([x0, x1], [yv, yv], color=INK, lw=1.6, solid_capstyle="butt", zorder=4)
    axb.annotate(lab, xy=((x0 + x1) / 2, yv), xytext=xyt, fontsize=FS_TXT, arrowprops=dict(arrowstyle="-", lw=0.5, color=INK2, shrinkA=0, shrinkB=1))
    rec("b", f"window median {lab}", f"{x0:.1f}-{x1:.1f} us", float(yv), source=relpath(RES06 / "12_plateau_steady_state.csv"))
rec("b", "CW-static gain a_H / a_L", "aH,aL", [float(r.cw_static_gain_aH), float(r.cw_static_gain_aL)], source=relpath(RES06 / "12_plateau_steady_state.csv"))
axb.text(hold / 2, 1.2, "high level (plateau)", ha="center", va="center", fontsize=FS_TXT, color=INK2)
axb.set_xlim(-12, hold + 50)
axb.set_ylim(0, 1.3)
axb.set_xlabel("time from plateau start (µs)")
axb.set_ylabel("gain / small-signal gain")
axb.legend(loc="upper right", fontsize=fs(6.0))

save(fig, "fig6")
