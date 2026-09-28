"""Figure 8: additional compression under phase modulation at constant amplitude.

(a) X = 1 - g/g_CW by case at OP1 (a_H = 0.179 V/m): full model; LTI+static surrogate (identical to the single-slow-state
    surrogate at constant amplitude), X from its stored pipeline gains; numerical checks (Nd 8001, dt 0.5 ns)
    - phase_mechanism/03_phase_mechanism_results.csv. Full range shown (no clipping).
(b) X vs transition rate: original sequences (phase_mechanism) and six zero-drift sequences per rate (revision P2:
    mean and 95% CI, results/revision/11_p2_summary.csv; x = nominal rate of the reference sequence).
(c) X vs peak |f_inst - f_IF|: steps / 300 ns / 900 ns ramps vs constant offsets (02_phase_mechanism_configs.csv + X).
(d) POST HOC: min IF magnitude in the first 0.5 us after isolated +/-pi/2 steps (full vs small-signal), maximum
    response-phase lag, and transient oscillation periods after power and phase steps at OP1 and OP2
    (manuscript/posthoc/*.csv, replication/05_posthoc_C9a_wideband.json).
(e) OP1 vs OP2 on one X axis: six zero-drift sequences per phase case (revision P2, mean and 95% CI) and the single
    constant-offset runs (phase_mechanism, replication). Full range shown (no clipping).
"""
import json

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from figstyle import (fs, COL2, INK, INK2, MUTED, GRID, MODEL, MODEL_MARKER, PM, FV, POSTHOC, BOLD,
                      panel, rec, save, relpath, zero_line, tm)

R = pd.read_csv(PM / "03_phase_mechanism_results.csv")
C = pd.read_csv(PM / "02_phase_mechanism_configs.csv").set_index("case")
V = pd.read_csv(FV / "replication" / "03_replication_results.csv").set_index("case")
main = R[(R.Nd == 4001) & (R.dt_ns == 1)].set_index("case")
SRC_R = relpath(PM / "03_phase_mechanism_results.csv")
SRC_V = relpath(FV / "replication" / "03_replication_results.csv")
RV2 = pd.read_csv(FV.parent / "revision" / "11_p2_summary.csv").set_index(["op", "case"])
SRC_RV = relpath(FV.parent / "revision" / "11_p2_summary.csv")
BLUE, LIGHT = MODEL["full"], "#b7d3f6"
FS_TXT = fs(6.3)

fig = plt.figure(figsize=(COL2, 5.3))
sf = fig.subfigures(2, 1, height_ratios=[1.05, 1], hspace=0.03)
top = sf[0].subfigures(1, 2, width_ratios=[1.08, 1], wspace=0.02)
bot = sf[1].subfigures(1, 2, width_ratios=[1.05, 1], wspace=0.02)

# ------------------------------------------------------------------ (a) X by case
ax = top[0].subplots()
panel(ax, "a")
ax.set_title("OP1, $a_H$ = 0.179 V/m", loc="right", fontsize=fs(7.2), pad=3)
cases = [("QPSK_REF", "reference QPSK (1 µs, 300 ns ramps)"), ("QPSK_REF_s2", "second sequence"),
         ("QPSK_SLOW", "slow (4 µs symbols)"), ("QPSK_FAST", "fast (0.5 µs symbols)"),
         ("QPSK_JUMP", "steps"), ("QPSK_LONGRAMP", "900 ns ramps"), ("TOGGLE", "deterministic sequence"),
         ("OFF_+0.125", "offset +0.125 MHz"), ("OFF_+1.31", "offset +1.31 MHz"), ("OFF_-1.31", "offset −1.31 MHz"),
         ("OFF_+2.62", "offset +2.62 MHz"), ("OFF_-2.62", "offset −2.62 MHz"), ("CW", "unmodulated (CW)")]
ys = np.arange(len(cases))[::-1].astype(float)
ys[:7] += 0.6  # gap between phase-modulated and offset groups
ys[-1] -= 0.4
g_lti_cw = main.loc["CW", "g_final_static_lti_mh"]
for (case, lab), y in zip(cases, ys):
    x = main.loc[case, "X"]
    ax.plot([x], [y], marker="o", ms=4.2, color=BLUE, mec="white", mew=0.5, zorder=4, ls="none",
            label="full model" if case == "QPSK_REF" else None)
    xl = 1 - main.loc[case, "g_final_static_lti_mh"] / g_lti_cw
    ax.plot([xl], [y], marker=MODEL_MARKER["lti"], ms=3.4, color=MODEL["lti"], mec="white", mew=0.4, zorder=3, ls="none",
            label="LTI+static = single-slow-state" if case == "QPSK_REF" else None)
    rec("a", f"X full {case}", case, x, source=SRC_R)
    rec("a", f"X LTI+static(=DSH) {case}", case, xl, source=SRC_R, note="1 - g_final_static_lti_mh/g_CW; g_final_dsh identical")
    for nd, dt, mk, lab2 in ((8001, 1.0, "s", f"check: $N_\\mathrm{{d}}$ = 8001"), (4001, 0.5, "D", "check: dt = 0.5 ns")):
        rr = R[(R.case == case) & (R.Nd == nd) & (R.dt_ns == dt)]
        if len(rr):
            ax.plot([rr.X.iloc[0]], [y], marker=mk, ms=5.4, mfc="none", mec=INK, mew=0.6, zorder=5, ls="none",
                    label=lab2 if case == "QPSK_REF" else None)
            rec("a", f"X check Nd{nd} dt{dt} {case}", case, float(rr.X.iloc[0]), source=SRC_R)
ax.axvline(0, color=INK2, lw=0.6, zorder=0)
ax.axhline((ys[6] + ys[7]) / 2, color=GRID, lw=0.6)
ax.axhline((ys[11] + ys[12]) / 2, color=GRID, lw=0.6)
ax.set_yticks(ys)
ax.set_yticklabels([lab for _, lab in cases], fontsize=fs(6.4))
ax.tick_params(axis="y", length=0)
ax.set_xlim(-2.2, 0.85)
ax.set_ylim(ys.min() - 0.8, ys.max() + 1.45)
ax.set_xlabel(r"excess gain reduction $X = 1 - g/g_\mathrm{CW}$")
ax.text(-2.12, ys[0] + 0.95, "symbol phase\ntransitions", ha="left", va="center", fontsize=FS_TXT, color=INK2, linespacing=1.0)
ax.text(-2.12, ys[7] + 0.62, "constant frequency\noffsets", ha="left", va="center", fontsize=FS_TXT, color=INK2, linespacing=1.0)
ax.legend(loc="upper center", fontsize=fs(5.8), handletextpad=0.3, borderaxespad=0.2, bbox_to_anchor=(0.45, -0.2), ncol=2, columnspacing=1.0)

# ------------------------------------------------------------------ (b) rate and (c) peak excursion
axb, axc = top[1].subplots(1, 2, gridspec_kw=dict(width_ratios=[1, 1.3]))
panel(axb, "b", "rate, OP1")
rate = [("QPSK_SLOW", "LO1_QPSK_SLOW_zd"), ("QPSK_REF", "LO1_QPSK_REF_zd"), ("QPSK_FAST", "LO1_QPSK_FAST_zd")]
r_x = [C.loc[a, "transition_rate_per_us"] for a, _ in rate]
r_o = [main.loc[a, "X"] for a, _ in rate]
k6 = ["SLOW", "REF", "FAST"]
r_z = np.array([RV2.loc[("OP1", k), "X_mean"] for k in k6])
r_lo = np.array([RV2.loc[("OP1", k), "ci_lo"] for k in k6])
r_hi = np.array([RV2.loc[("OP1", k), "ci_hi"] for k in k6])
axb.plot(r_x, r_o, color=BLUE, marker="o", ms=4, mec="white", mew=0.5, lw=1.1, label="original sequence")
axb.errorbar(r_x, r_z, yerr=[r_z - r_lo, r_hi - r_z], color=BLUE, marker="s", ms=3.6, mfc="white", mec=BLUE, mew=0.8, lw=0.9,
             ls=(0, (3, 2)), elinewidth=0.8, capsize=1.6, capthick=0.8, label="six zero-drift sequences")
rec("b", "X original", r_x, r_o, source=SRC_R)
rec("b", "X zero-drift six sequences", r_x, r_z, r_lo, r_hi, source=SRC_RV)
zero_line(axb)
axb.set_xlim(0, 1.6)
axb.set_ylim(-0.05, 0.75)
axb.set_xticks([0.15, 0.73, 1.45])
axb.set_xticklabels(["0.15", "0.73", "1.45"])
axb.set_xlabel("nonzero transitions per µs")
axb.set_ylabel("X")
axb.legend(loc="upper left", fontsize=fs(5.8), borderaxespad=0.2)

panel(axc, "c")
# descriptor right-aligned: centred, it collided with the narrow panel's "(c)" label (layout only)
axc.set_title("peak excursion, OP1", loc="right", fontsize=8, pad=3.5, color=INK)
shape = [("QPSK_LONGRAMP", "900 ns"), ("QPSK_REF", "300 ns"), ("QPSK_JUMP", "steps")]
px = [C.loc[c, "peak_abs_freq_offset_Hz"] / 1e6 for c, _ in shape]
py = [main.loc[c, "X"] for c, _ in shape]
axc.plot(px, py, color=BLUE, marker="o", ms=4, mec="white", mew=0.5, lw=1.0, label="phase transitions")
for (c, lab), xx, yy in zip(shape, px, py):
    off, ha, va = {"QPSK_LONGRAMP": ((-2, -9), "right", "top"), "QPSK_REF": ((2, 6), "left", "bottom"), "QPSK_JUMP": ((0, 6), "center", "bottom")}[c]
    axc.annotate(lab, (xx, yy), xytext=off, textcoords="offset points", ha=ha, va=va, fontsize=fs(5.9))
rec("c", "X vs peak |f_inst - f_IF| (phase)", px, py, source=relpath(PM / "02_phase_mechanism_configs.csv"))
offs = ["OFF_+0.125", "OFF_+1.31", "OFF_-1.31", "OFF_+2.62", "OFF_-2.62"]
ox = [C.loc[c, "peak_abs_freq_offset_Hz"] / 1e6 for c in offs]
oy = [main.loc[c, "X"] for c in offs]
for c, xx, yy in zip(offs, ox, oy):
    axc.plot([xx], [yy], ls="none", marker="^" if "+" in c else "v", ms=4.2, mfc="white", mec=BLUE, mew=0.8,
             label=("offsets (▲ +, ▼ −)" if c == "OFF_+0.125" else None))
rec("c", "X vs |offset| (offsets)", ox, oy, source=relpath(PM / "02_phase_mechanism_configs.csv"))
zero_line(axc)
axc.set_xscale("log")
axc.set_xlim(0.07, 3e4)
axc.set_ylim(-2.2, 0.85)
axc.set_xlabel("peak |$f_\\mathrm{inst} - f_\\mathrm{IF}$| (MHz)")
axc.text(40, 0.26, "symbol phase\ntransitions", ha="center", va="center", fontsize=fs(5.9), color=INK)
axc.text(8, -1.05, "constant offsets\n(▲ +δ, ▼ −δ)", ha="left", va="center", fontsize=fs(5.9), color=INK)

# ------------------------------------------------------------------ (d) post hoc: isolated transitions and periods
d1, d2, d3 = bot[0].subplots(1, 3, gridspec_kw=dict(width_ratios=[1, 1, 1.15]))
panel(d1, "d")
bot[0].suptitle("post hoc: isolated ±π/2 steps (OP1) and oscillation periods", fontsize=fs(7.5), x=0.55)
it = pd.read_csv(POSTHOC / "posthoc_isolated_transitions.csv")
h = it[np.isclose(it.dphi_rad.abs(), np.pi / 2, atol=0.1)].reset_index(drop=True)
xk = np.arange(len(h))
d1.plot(xk, h.min_absfull_rel_first_0p5us, ls="none", marker="o", ms=4, color=BLUE, mec="white", mew=0.5, label="full model")
d1.plot(xk, h.min_abslin_rel_first_0p5us, ls="none", marker="s", ms=3.6, mfc="white", mec=MUTED, mew=0.8, label="small-signal")
rec("d", "min |z| rel. pre, first 0.5 us, full", h.t_boundary_us.values, h.min_absfull_rel_first_0p5us.values, source=relpath(POSTHOC / "posthoc_isolated_transitions.csv"))
rec("d", "min |z| rel. pre, first 0.5 us, small-signal", h.t_boundary_us.values, h.min_abslin_rel_first_0p5us.values, source=relpath(POSTHOC / "posthoc_isolated_transitions.csv"))
d1.set_ylim(0, 1.1)
d1.set_xlim(-0.6, len(h) - 0.4)
d1.set_xticks(xk)
d1.set_xticklabels([f"{t:.0f}" for t in h.t_boundary_us], fontsize=fs(6))
d1.set_xlabel("step time (µs)", fontsize=fs(6.8))
d1.set_ylabel("min IF magnitude, first 0.5 µs\n(relative to before the step)", fontsize=fs(6.8))
d1.legend(loc="lower left", fontsize=fs(5.6), borderaxespad=0.1, handletextpad=0.2)

well = h.min_absfull_rel_first_0p5us > 0.4
d2.plot(xk[well], h.max_abs_dtheta_rad[well], ls="none", marker="o", ms=4, color=BLUE, mec="white", mew=0.5)
d2.plot(xk[~well], h.max_abs_dtheta_rad[~well], ls="none", marker="o", ms=4, mfc="white", mec=BLUE, mew=0.8)
for k in xk[~well]:
    d2.annotate("magnitude 0.35:\nphase ill-defined", (k, h.max_abs_dtheta_rad[k]), xytext=(1.2, 2.45), textcoords="data",
                ha="center", va="center", fontsize=fs(5.5), color=INK2, arrowprops=dict(arrowstyle="-", lw=0.4, color=INK2, shrinkB=3))
rec("d", "max response-phase lag (rad)", h.t_boundary_us.values, h.max_abs_dtheta_rad.values, source=relpath(POSTHOC / "posthoc_isolated_transitions.csv"),
    note="open marker = min magnitude <= 0.4 (excluded from the quoted 0.9-1.7 rad range)")
d2.set_ylim(0, 3.3)
d2.set_xlim(-0.6, len(h) - 0.4)
d2.set_xticks(xk)
d2.set_xticklabels([f"{t:.0f}" for t in h.t_boundary_us], fontsize=fs(6))
d2.set_xlabel("step time (µs)", fontsize=fs(6.8))
d2.set_ylabel("max response-phase lag (rad)", fontsize=fs(6.8))

ts = pd.read_csv(POSTHOC / "posthoc_timescales.csv")
pw = ts[ts.transient.str.startswith("power")].osc_period_us.values
ph = ts[ts.transient.str.startswith("phase")].osc_period_us.values
c9 = json.loads((FV / "replication" / "05_posthoc_C9a_wideband.json").read_text())
p2w, p2h = c9["LO2_power_step_period_us"], np.asarray(c9["LO2_phase_step_periods_us"])
d3.axhspan(0.8, 2.2, color="#f3f2ee", lw=0, zorder=0)
d3.text(1.45, 2.12, "prespecified\nsearch band\n0.8–4 µs", ha="right", va="top", fontsize=fs(5.5), color=INK2)
jit = np.linspace(-0.12, 0.12, 9)
d3.plot(np.zeros(len(pw)) - 0.22, pw, ls="none", marker="s", ms=3.8, mfc=BLUE, mec="white", mew=0.4, label="power step")
d3.plot(jit[:len(ph)] + 0.12, ph, ls="none", marker="o", ms=3.2, mfc="white", mec=BLUE, mew=0.7, label="phase steps")
d3.plot([0.78], [p2w], ls="none", marker="s", ms=3.8, mfc=BLUE, mec="white", mew=0.4)
d3.plot(jit[:len(p2h)] + 1.12, p2h, ls="none", marker="o", ms=3.2, mfc="white", mec=BLUE, mew=0.7)
rec("d", "osc period OP1 power (us)", "OP1", pw, source=relpath(POSTHOC / "posthoc_timescales.csv"))
rec("d", "osc period OP1 phase (us)", "OP1", ph, source=relpath(POSTHOC / "posthoc_timescales.csv"))
rec("d", "osc period OP2 power (us)", "OP2", p2w, source=relpath(FV / "replication" / "05_posthoc_C9a_wideband.json"))
rec("d", "osc period OP2 phase (us)", "OP2", p2h, source=relpath(FV / "replication" / "05_posthoc_C9a_wideband.json"))
d3.set_xticks([0, 1])
d3.set_xticklabels(["OP1", "OP2"])
d3.set_xlim(-0.55, 1.55)
d3.set_ylim(0, 2.2)
d3.set_ylabel("oscillation period (µs)", fontsize=fs(6.8))
d3.legend(loc="center right", fontsize=fs(5.6), borderaxespad=0.1, handletextpad=0.2, bbox_to_anchor=(1.0, 0.52))

# ------------------------------------------------------------------ (e) OP2 vs OP1
ax = bot[1].subplots()
panel(ax, "e", "OP1 vs OP2 (both 5 MHz IF)")
pairs = [("reference", "REF"), ("slow", "SLOW"), ("fast", "FAST"), ("steps", "JUMP"), ("900 ns\nramps", "LONGRAMP"),
         ("offset\n+0.125", ("OFF_+0.125", "LO2_OFF_+0.125")), ("offset\n+1.31", ("OFF_+1.31", "LO2_OFF_+1.31")),
         ("offset\n−1.31", ("OFF_-1.31", "LO2_OFF_-1.31"))]
bw = 0.36
for i, (lab, key) in enumerate(pairs):
    ll = lab.replace(chr(10), " ")
    if isinstance(key, str):
        for op, dx, fc in (("OP1", -bw / 2 - 0.01, BLUE), ("OP2", bw / 2 + 0.01, LIGHT)):
            m, lo, hi = (RV2.loc[(op, key), c] for c in ("X_mean", "ci_lo", "ci_hi"))
            ax.bar(i + dx, m, bw, color=fc, edgecolor=BLUE, lw=0.6, zorder=2, label=(f"{op}, six zero-drift sequences" if i == 0 else None))
            ax.errorbar(i + dx, m, yerr=[[m - lo], [hi - m]], color=INK, lw=0, elinewidth=0.8, capsize=1.6, capthick=0.8, zorder=3)
            rec("e", f"X {op} {ll}", key, m, lo, hi, source=SRC_RV)
    else:
        x1, x2 = main.loc[key[0], "X"], V.loc[key[1], "X"]
        ax.bar(i - bw / 2 - 0.01, x1, bw, color="white", edgecolor=BLUE, lw=0.8, zorder=2, hatch="//////",
               label="OP1, single run" if i == 5 else None)
        ax.bar(i + bw / 2 + 0.01, x2, bw, color="white", edgecolor=BLUE, lw=0.6, zorder=2, hatch="....",
               label="OP2, single run" if i == 5 else None)
        rec("e", f"X OP1 {ll}", key[0], x1, source=SRC_R)
        rec("e", f"X OP2 {ll}", key[1], x2, source=SRC_V)
zero_line(ax)
ax.axhline(0.20, color=INK, lw=0.7, ls=(0, (1.2, 1.4)), zorder=1)
ax.text(len(pairs) - 0.45, 0.215, "prespecified OP2 magnitude\ncriterion $X_\\mathrm{REF}$ ≥ 0.20 (failed: 0.159)", ha="right", va="bottom", fontsize=fs(5.8), color=INK)
ax.set_xticks(range(len(pairs)))
ax.set_xticklabels([p[0] for p in pairs], fontsize=fs(6.1))
ax.tick_params(axis="x", length=0)
ax.set_xlim(-0.6, len(pairs) - 0.4)
ax.set_ylim(-0.92, 0.72)
ax.set_ylabel(r"$X = 1 - g/g_\mathrm{CW}$ (same operating point)")
hh, ll_ = ax.get_legend_handles_labels()
order = [ll_.index(t) for t in ("OP1, six zero-drift sequences", "OP2, six zero-drift sequences", "OP1, single run", "OP2, single run")]
ax.legend([hh[k] for k in order], [ll_[k] for k in order], loc="lower left", fontsize=fs(5.8), borderaxespad=0.2)
rec("e", "OP2 reference mean, two prespecified sequences", "V_Xref_mean", float(np.mean([V.loc["LO2_QPSK_REF_zd", "X"], V.loc["LO2_QPSK_REF_zd_s2", "X"]])), source=SRC_V)

save(fig, "fig8")
