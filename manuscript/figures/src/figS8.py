"""Supplementary Figure S8 (formerly main Fig. 9): archived configuration checks at levels matched to each configuration's
E1dB (prespecified; results/revision2/00_decision_rules.json). The matched level lies above the LO field at 10 MHz and with the
weaker probe (regime flag: results/tqe_viability_final/05_numerical/vc_regime_table.csv); kept as a record.

(a) CW fundamental gain re the weak-tone slope vs field amplitude at IF 5/10/15 MHz (standard probe) and with the
    weaker probe (5 MHz); markers: the matched phase level 2.4468 E1dB of each configuration (tables/*.npz).
(b) X of the reference QPSK case (3 sequences, mean and 95% CI) by configuration (11_p1_if_summary.csv, 31_p3_verdict.json).
(c) Paired fitted-gain change, xi = 4 minus xi = 0.05, +3 dB (4 realizations, 95% CI): full model and the CW-matched
    static surrogate of the same configuration (11_p1_if_summary.csv, 31_p3_verdict.json).
(d) R = g(df)/g_CW at the five prespecified offsets, 5/10/15 MHz (11_p1_if_summary.csv).
(e) Normalized IF-envelope excess after a +5% step at 1.0 E1dB (5 MHz, OP1), standard and weaker probe, identical code
    (p3/steps/*_plus5pct_at_1.0.{npz,json}); annotation: slow time constants of the prespecified two-exponential fit.
"""
import json

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from figstyle import (fs, COL2, ROOT, INK, INK2, MUTED, EB, panel, rec, save, relpath, zero_line)

R2 = ROOT / "results" / "revision2"
S1 = pd.read_csv(R2 / "11_p1_if_summary.csv").set_index("IF_MHz")
V3 = json.loads((R2 / "31_p3_verdict.json").read_text())
W = V3["per_probe"]["weak"]
CONF = [("5 MHz", "#2a78d6", "o"), ("10 MHz", "#4a3aa7", "s"), ("15 MHz", "#008300", "D"), ("weak probe", "#e34948", "^")]
LABEL = {"5 MHz": "5 MHz (reference)", "10 MHz": "10 MHz", "15 MHz": "15 MHz", "weak probe": "5 MHz, weaker probe"}
COL = {k: c for k, c, _ in CONF}
MK = {k: m for k, _, m in CONF}
TAB = {"5 MHz": R2 / "tables" / "OP1_IF5.npz", "10 MHz": R2 / "tables" / "OP1_IF10_ext.npz", "15 MHz": R2 / "tables" / "OP1_IF15.npz",
       "weak probe": R2 / "tables" / "OP1_IF5_weakprobe_ext.npz"}
SRC1, SRC3 = relpath(R2 / "11_p1_if_summary.csv"), relpath(R2 / "31_p3_verdict.json")

fig = plt.figure(figsize=(COL2, 4.9))
gs = fig.add_gridspec(2, 6)
axs = [fig.add_subplot(gs[0, 0:2]), fig.add_subplot(gs[0, 2:4]), fig.add_subplot(gs[0, 4:6]),
       fig.add_subplot(gs[1, 0:3]), fig.add_subplot(gs[1, 3:6])]

# ------------------------------------------------------------------ (a) CW compression curves
ax = axs[0]
panel(ax, "a", "CW calibration")
zero_line(ax)
ax.axhline(-1, color=INK2, lw=0.6, ls=(0, (1.2, 1.4)), zorder=0)
ax.axvline(0.5, color=MUTED, lw=0.7, ls=(0, (3, 2)), zorder=0)
ax.text(0.52, 5.2, "$A_\\mathrm{LO}$", fontsize=fs(6.3), color=INK2, ha="left", va="top")
for k, c, m in CONF:
    t = np.load(TAB[k])
    a, g = t["amps"][1:], t["gain_db"][1:]
    ax.plot(a, g, color=c, lw=1.1, label=LABEL[k])
    rec("a", f"CW gain dB {k}", a, g, source=relpath(TAB[k]))
    ah = 2.4468 * float(t["e1db"])
    gh = float(np.interp(ah, a, g))
    ax.plot([ah], [gh], marker=m, color=c, ms=4.2, mec="white", mew=0.6, ls="none")
    rec("a", f"matched level {k}", ah, gh, source=relpath(TAB[k]), note="2.4468 E1dB")
ax.set_xscale("log")
ax.set_xlim(0.004, 2.5)
ax.set_xticks([0.01, 0.1, 1])
ax.set_xticklabels(["0.01", "0.1", "1"])
ax.set_ylim(-50, 7)
ax.set_xlabel("CW field amplitude (V/m)")
ax.set_ylabel("fundamental gain re weak-tone slope (dB)")
ax.legend(loc="lower left", fontsize=fs(6.0), borderaxespad=0.2, handlelength=1.4)

# ------------------------------------------------------------------ (b) X of the reference case
ax = axs[1]
panel(ax, "b", "phase modulation (QPSK)")
zero_line(ax)
ax.axhline(0.05, color=INK, lw=0.7, ls=(0, (1.2, 1.4)), zorder=1)
ax.text(3.45, 0.065, "0.05", ha="right", va="bottom", fontsize=fs(6.0), color=INK)
vals = [(f"{f} MHz", S1.X_REF_mean[f], S1.X_REF_lo[f], S1.X_REF_hi[f], SRC1) for f in (5, 10, 15)] + [("weak probe", W["X_REF_mean"], W["X_REF_lo"], W["X_REF_hi"], SRC3)]
for i, (k, m, lo, hi, src) in enumerate(vals):
    ax.errorbar([i], [m], yerr=[[m - lo], [hi - m]], color=COL[k], marker=MK[k], ms=4.2, mec="white", mew=0.6, ls="none", **EB)
    rec("b", f"X REF {k}", i, m, lo, hi, source=src)
ax.set_xticks(range(4))
ax.set_xticklabels(["5", "10", "15", "5, weaker\nprobe"])
ax.set_xlim(-0.5, 3.5)
ax.set_ylim(-0.2, 0.6)
ax.set_xlabel("IF (MHz)")
ax.set_ylabel("$X = 1 - g/g_\\mathrm{CW}$ at the matched level")

# ------------------------------------------------------------------ (c) dwell gain change
ax = axs[2]
panel(ax, "c", "dwell: ξ = 4 minus 0.05")
zero_line(ax)
dvals = [(f"{f} MHz", S1.d_gain_full[f], S1.d_gain_full_lo[f], S1.d_gain_full_hi[f], S1.d_gain_static[f], SRC1) for f in (5, 10, 15)]
dvals += [("weak probe", W["d_gain_full"], W["d_gain_full_lo"], W["d_gain_full_hi"], W["d_gain_static"], SRC3)]
for i, (k, m, lo, hi, st, src) in enumerate(dvals):
    ax.errorbar([i - 0.12], [m], yerr=[[m - lo], [hi - m]], color=COL[k], marker=MK[k], ms=4.2, mec="white", mew=0.6, ls="none", **EB)
    ax.plot([i + 0.12], [st], color=COL[k], marker=MK[k], ms=4.0, mfc="white", mew=0.9, ls="none")
    rec("c", f"d gain full {k}", i, m, lo, hi, source=src)
    rec("c", f"d gain static {k}", i, st, source=src)
ax.set_xticks(range(4))
ax.set_xticklabels(["5", "10", "15", "5, weaker\nprobe"])
ax.set_xlim(-0.5, 3.5)
ax.set_ylim(-0.2, 0.02)
ax.set_xlabel("IF (MHz)")
ax.set_ylabel("paired change in fitted gain / small-signal")
ax.plot([], [], color=INK2, marker="o", ms=4.0, mec="white", ls="none", label="full model")
ax.plot([], [], color=INK2, marker="o", ms=4.0, mfc="white", mew=0.9, ls="none", label="static surrogate")
ax.legend(loc="lower right", fontsize=fs(6.0), borderaxespad=0.2, handletextpad=0.3)

# ------------------------------------------------------------------ (d) frequency selectivity
ax = axs[3]
panel(ax, "d", "constant offsets")
ax.axhline(1.0, color=INK2, lw=0.6, zorder=0)
offs = sorted([c for c in S1.columns if c.startswith("R_")], key=lambda c: float(c[2:]))
for f in (5, 10, 15):
    k = f"{f} MHz"
    x = [float(c[2:]) for c in offs]
    y = [S1.loc[f, c] for c in offs]
    ax.plot(x, y, color=COL[k], marker=MK[k], ms=3.6, mec="white", mew=0.5, lw=1.0, label=k)
    rec("d", f"R {k}", x, y, source=SRC1)
ax.set_xlim(-3, 3)
ax.set_xticks([-2.62, -1.31, 0, 1.31, 2.62])
ax.set_xticklabels(["−2.62", "−1.31", "0", "1.31", "2.62"])
ax.set_ylim(0, 3.3)
ax.set_xlabel("constant frequency offset (MHz)")
ax.set_ylabel("$R = g(\\delta f)/g_\\mathrm{CW}$")
ax.legend(loc="upper left", fontsize=fs(6.0), borderaxespad=0.2)

# ------------------------------------------------------------------ (e) memory: step responses, standard vs weaker probe
ax = axs[4]
panel(ax, "e", "+5% step at $E_\\mathrm{1dB}$: memory")
zero_line(ax)
for k, probe, lab in (("5 MHz", "standard", "standard probe"), ("weak probe", "weak", "weaker probe")):
    z = np.load(R2 / "p3" / "steps" / f"{probe}_plus5pct_at_1.0.npz")
    fit = json.loads((R2 / "p3" / "steps" / f"{probe}_plus5pct_at_1.0.json").read_text())
    m = z["t"] <= 12e-6
    ax.plot(z["t"][m] * 1e6, z["e"][m], color=COL[k], lw=1.2,
            label=f"{lab}: slow component {fit['fit2_tau_slow_s'] * 1e6:.2f} µs")
    rec("e", f"normalized excess {probe} (every 0.25 us)", z["t"][m][::25] * 1e6, z["e"][m][::25], source=relpath(R2 / "p3" / "steps" / f"{probe}_plus5pct_at_1.0.npz"))
    rec("e", f"slow time constant {probe} (us)", probe, fit["fit2_tau_slow_s"] * 1e6, source=relpath(R2 / "p3" / "steps" / f"{probe}_plus5pct_at_1.0.json"))
ax.set_xlim(0, 12)
ax.set_ylim(-0.5, 1.1)
ax.set_xlabel("time after the step (µs)")
ax.set_ylabel("normalized IF-envelope excess")
ax.legend(loc="upper right", fontsize=fs(6.0), borderaxespad=0.2)

save(fig, "figS8", supplement=True)
