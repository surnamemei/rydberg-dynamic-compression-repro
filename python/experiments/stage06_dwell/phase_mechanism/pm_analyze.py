"""Phase-mechanism analysis: metrics, pre-registered verdicts, figures."""
from __future__ import annotations

import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE))

import numpy as np
import pandas as pd
from scipy import stats

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import pm_run as P

PM = P.PM
PRE = json.loads((PM / "00_preregistration.json").read_text())
MH = np.load(P.OUT06 / "artifacts" / "controls_MH_Nd4001.npz")
T0_US, TEND_US = P.T0 * 1e6, P.T_END * 1e6
plt.rcParams.update({"axes.grid": True, "grid.alpha": .25, "axes.spines.top": False, "axes.spines.right": False, "font.size": 9, "lines.linewidth": 1.5})
BLUE, ORANGE, AQUA, YELLOW, GREY = "#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#8a8a85"
PHASE_CASES = ["QPSK_SLOW", "QPSK_REF", "QPSK_FAST", "QPSK_JUMP", "QPSK_LONGRAMP", "TOGGLE", "QPSK_REF_s2"]
OFFSETS = ["OFF_-2.62", "OFF_-1.31", "OFF_+0.125", "OFF_+1.31", "OFF_+2.62"]


def load(case, nd=4001, dt=1):
    return np.load(PM / "runs" / f"{case}_Nd{nd}_dt{dt:g}.npz")


def smooth(x, n=50):  # traces are stored at 20 ns; 50 samples = 1 us
    return np.convolve(x, np.ones(n) / n, "same")


def gain(z, model="full"):
    return np.abs(z[f"z_{model}"]) / np.abs(z["z_linear"])


def metrics(case, nd=4001, dt=1, cw_traj=None):
    z = load(case, nd, dt)
    t = z["t"] * 1e6
    g = gain(z)
    gs = smooth(g)
    # windows exclude the final 0.5 us: the moving-average demodulator has edge artifacts at the array end
    pre, last, mod = (t >= T0_US - 10) & (t < T0_US), (t >= TEND_US - 20.5) & (t < TEND_US - .5), (t >= T0_US) & (t < TEND_US - .5)
    g_pre, g_fin = float(np.median(g[pre])), float(np.median(g[last]))
    zl = z["z_linear"][last]
    g_coh = float(np.abs(np.mean(z["z_full"][last] * np.conj(zl) / np.abs(zl))) / np.mean(np.abs(zl)))
    # settling: first time a centered 2 us average completes 90% (63%) of the CW -> modulated change
    g2 = smooth(g, 100)[mod]
    change = g_fin - g_pre
    t63 = t90 = np.nan
    if abs(change) >= .02:
        frac = (g2 - g_pre) / change
        for lvl in (.63, .9):
            i = np.flatnonzero(frac >= lvl)
            val = float(t[mod][i[0]] - T0_US) if len(i) else np.nan
            t63, t90 = (val, t90) if lvl == .63 else (t63, val)
    sgn = gs[mod] - g_fin
    k = int(np.argmax(np.abs(sgn[t[mod] > T0_US + .5])))
    peak = float(sgn[t[mod] > T0_US + .5][k])
    integ = float(np.sum((cw_traj[mod] - gs[mod])) * (t[1] - t[0])) if cw_traj is not None else np.nan
    w = z["pv"]
    pops = (z["pops"] * w).sum(axis=1)
    row = {"case": case, "Nd": nd, "dt_ns": dt, "g_pre": g_pre, "g_final": g_fin, "g_coherent_final": g_coh,
           "t63_us": t63, "t90_us": t90, "peak_transient_dev_vs_final": peak, "integrated_excess_compression_vs_CW_us": integ,
           "probe_db_last20us": float(np.mean(z["probe_db"][last])), "pop_ground": pops[0], "pop_intermediate": pops[1],
           "pop_rydberg_3": pops[2], "pop_rydberg_4": pops[3], "pop_trace": float(pops.sum()),
           "abs_out_last20": float(np.median(np.abs(z["z_full"][last]))), "abs_lin_last20": float(np.median(np.abs(z["z_linear"][last])))}
    for m in ("static_mh", "static_lti_mh", "dsh"):
        if f"z_{m}" in z:
            row[f"g_final_{m}"] = float(np.median(gain(z, m)[last]))
    return row, gs


def organizer_test(df, var):
    """Pre-registered: equal-v cases differ by <=0.10 in X and per-level mean X is monotone (tol 0.03)."""
    d = df[["case", var, "X"]].copy()
    d["lvl"] = d[var].round(6)
    spread = d.groupby("lvl").X.agg(lambda s: s.max() - s.min())
    means = d.groupby("lvl").X.mean().sort_index().values
    inc = np.all(np.diff(means) >= -.03)
    dec = np.all(np.diff(means) <= .03)
    return bool((spread <= .10).all() and (inc or dec)), float(spread.max()), bool(inc or dec), float(stats.spearmanr(d[var], d.X).statistic)


def main():
    cfg = pd.read_csv(PM / "02_phase_mechanism_configs.csv")
    _, cw_traj = metrics("CW")
    rows = []
    for case in P.CASES:
        r, _ = metrics(case, cw_traj=cw_traj)
        rows.append(r)
    for case in P.CHECKS:
        for nd, dt in ((8001, 1), (4001, .5)):
            r, _ = metrics(case, nd, dt, cw_traj=cw_traj)
            rows.append(r)
    res = pd.DataFrame(rows)
    g_cw = float(res[(res.case == "CW") & (res.Nd == 4001) & (res.dt_ns == 1)].g_final.iloc[0])
    # CW-static compression at a_H (all-zone table, same small-signal normalization as FU-B)
    g_static = float(np.interp(P.A_H, MH["amps"], np.abs(MH["C"][:, 1])) / P.A_H / (abs(MH["C"][1, 1]) / MH["amps"][1]))
    res["X"] = 1 - res.g_final / g_cw
    res["delta_g_vs_CW"] = res.g_final - g_cw
    res["delta_g_vs_static_prediction"] = res.g_final - g_static
    res["static_prediction_g"] = g_static
    res = res.merge(cfg, on="case", how="left")
    res.to_csv(PM / "03_phase_mechanism_results.csv", index=False)
    main_ = res[(res.Nd == 4001) & (res.dt_ns == 1)].set_index("case")
    X = main_.X

    v = []
    def add(item, value, rule, outcome):
        v.append({"item": item, "value": value, "rule": rule, "outcome": outcome})
    effect = X["QPSK_REF"] >= .20 and abs(X["QPSK_REF_s2"] - X["QPSK_REF"]) <= .25 * X["QPSK_REF"]
    add("EFFECT", f"X_REF={X['QPSK_REF']:.3f}, X_REF_s2={X['QPSK_REF_s2']:.3f}", PRE["rules"]["EFFECT"], "PASS" if effect else "FAIL")
    num = []
    for c in P.CHECKS:
        base = X[c]
        for nd, dt in ((8001, 1), (4001, .5)):
            num.append(abs(float(res[(res.case == c) & (res.Nd == nd) & (res.dt_ns == dt)].X.iloc[0]) - base))
    numerics = max(num) <= .02
    add("numerics_ok", f"max |dX| = {max(num):.4f} over {len(num)} checks", PRE["rules"]["numerics_ok"], "PASS" if numerics else "FAIL")
    h3s = X["OFF_+0.125"] >= .8 * X["QPSK_REF"]
    add("H3_strong_KILL (mean-drift detuning)", f"X(OFF_+0.125)={X['OFF_+0.125']:.3f} vs 0.8 X_REF={.8 * X['QPSK_REF']:.3f}", PRE["rules"]["H3_strong_KILL"], "TRIGGERED" if h3s else "not triggered")
    h3p = any(X[c] >= .8 * X["QPSK_REF"] for c in ("OFF_+1.31", "OFF_-1.31", "OFF_+2.62", "OFF_-2.62"))
    add("H3_partial (instantaneous-frequency detuning)", "X offsets: " + ", ".join(f"{c} {X[c]:+.3f}" for c in OFFSETS), PRE["rules"]["H3_partial"], "SUPPORTED" if h3p else "REJECTED")
    h1 = (X["QPSK_SLOW"] < X["QPSK_REF"] < X["QPSK_FAST"]) and (X["QPSK_FAST"] - X["QPSK_SLOW"] >= .1) and all(X[c] < .5 * X["QPSK_REF"] for c in ("OFF_+0.125", "OFF_+1.31", "OFF_-1.31"))
    add("H1_transition_rate", f"SLOW {X['QPSK_SLOW']:.3f} / REF {X['QPSK_REF']:.3f} / FAST {X['QPSK_FAST']:.3f}", PRE["rules"]["H1_transition_rate"], "SUPPORTED" if h1 else "NOT SUPPORTED (non-monotone)" if X["QPSK_SLOW"] < X["QPSK_REF"] else "NOT SUPPORTED")
    h2 = (X["QPSK_JUMP"] > X["QPSK_REF"] > X["QPSK_LONGRAMP"]) and (X["QPSK_JUMP"] - X["QPSK_LONGRAMP"] >= .1)
    add("H2_frequency_swing", f"JUMP {X['QPSK_JUMP']:.3f} / REF {X['QPSK_REF']:.3f} / LONGRAMP {X['QPSK_LONGRAMP']:.3f}", PRE["rules"]["H2_frequency_swing"], "SUPPORTED" if h2 else "REJECTED")
    h4 = X["QPSK_JUMP"] >= .2 and X["QPSK_REF"] < .5 * X["QPSK_JUMP"] and X["QPSK_LONGRAMP"] < .5 * X["QPSK_JUMP"]
    add("H4_discrete_jump", f"JUMP {X['QPSK_JUMP']:.3f} vs REF {X['QPSK_REF']:.3f}, LONGRAMP {X['QPSK_LONGRAMP']:.3f}", PRE["rules"]["H4_discrete_jump"], "SUPPORTED" if h4 else "REJECTED")
    st = abs(X["TOGGLE"] - X["QPSK_REF"]) <= .25 * X["QPSK_REF"]
    add("stochasticity_irrelevant", f"TOGGLE {X['TOGGLE']:.3f} vs REF {X['QPSK_REF']:.3f} (|d|={abs(X['TOGGLE'] - X['QPSK_REF']):.3f}, limit {.25 * X['QPSK_REF']:.3f})", PRE["rules"]["stochasticity_irrelevant"], "PASS" if st else "FAIL (pattern matters)")
    org_df = main_.drop(index="QPSK_REF_s2").reset_index()
    passing = []
    for var, label in (("transition_rate_per_us", "transition rate"), ("mean_abs_freq_offset_Hz", "mean |f_inst - f_IF|"), ("rms_freq_offset_Hz", "RMS(f_inst - f_IF)"),
                       ("peak_abs_freq_offset_Hz", "peak |f_inst - f_IF|"), ("mean_signed_freq_offset_Hz", "mean signed offset")):
        ok, spread, mono, rho = organizer_test(org_df, var)
        add(f"organizer: {label}", f"max within-level spread {spread:.3f}, monotone={mono}, Spearman {rho:+.2f}", PRE["rules"]["organizer"], "PASS" if ok else "FAIL")
        if ok:
            passing.append(label)
        ok2, spread2, mono2, rho2 = organizer_test(org_df[org_df.kind != "offset"], var)
        add(f"organizer (phase-modulated cases + CW only): {label}", f"spread {spread2:.3f}, monotone={mono2}, Spearman {rho2:+.2f}", "descriptive (not pre-registered)", "PASS" if ok2 else "FAIL")
    if not effect:
        dec = "HOLD"
    elif h3s or X["QPSK_REF"] < .05:
        dec = "KILL"
    elif numerics and len(passing) == 1:
        dec = "GO"
    else:
        dec = "CONDITIONAL GO"
    add("DECISION", f"organizers passing: {passing or 'none'}", PRE["rules"]["decision"], dec)
    vd = pd.DataFrame(v)
    vd.to_csv(PM / "04_phase_mechanism_verdicts.csv", index=False)
    print(vd[["item", "value", "outcome"]].to_string(), flush=True)
    print(main_[["g_pre", "g_final", "X", "g_coherent_final", "t63_us", "t90_us", "peak_transient_dev_vs_final", "integrated_excess_compression_vs_CW_us", "g_final_static_lti_mh", "g_final_dsh", "probe_db_last20us"]].round(4).to_string())
    figures(main_, res, g_cw, g_static)


def isolated_transition_response():
    """SLOW case: gain around each nonzero transition (4 us apart), aligned on the boundary."""
    z = load("QPSK_SLOW")
    t = z["t"] * 1e6
    g = smooth(gain(z), 10)
    d = P.transitions(P.CASES["QPSK_SLOW"])
    out = []
    for k, dk in enumerate(d, start=1):
        if dk == 0:
            continue
        tb = T0_US + k * 4.0
        m = (t >= tb - 1) & (t < tb + 3)
        if m.sum() < 190:
            continue
        out.append((dk, t[m] - tb, g[m]))
    return out


def figures(main_, res, g_cw, g_static):
    # pm01: RF signal (real part of E e^{jwt}) around the first nonzero transition for REF / JUMP / LONGRAMP
    fig, ax = plt.subplots(3, 1, figsize=(11, 6.4), sharex=True)
    d = P.transitions(P.CASES["QPSK_REF"])
    k1 = int(np.flatnonzero(d != 0)[0]) + 1
    tb = P.T0 + k1 * 1e-6
    t = np.arange(int((tb - 1.2e-6) * 1e9), int((tb + 1.2e-6) * 1e9)) * 1e-9
    for a, c in zip(ax, ("QPSK_JUMP", "QPSK_REF", "QPSK_LONGRAMP")):
        phi = P.phase(P.CASES[c], t)
        a.plot((t - tb) * 1e6, P.A_H * np.cos(2 * np.pi * P.IF * t + phi), color=BLUE, lw=.8, label="Re{E(t) e^{jωt}}")
        a.plot((t - tb) * 1e6, np.full_like(t, P.A_H), color="k", lw=.8, ls="--", label="|E| (constant)")
        a.set_ylabel("field (V/m)")
        a.set_title(f"{c}: phase step {d[k1 - 1] / np.pi:+.2f}π, ramp {P.CASES[c]['Tp'] * 1e9:.0f} ns", fontsize=9)
    ax[0].legend(fontsize=7, frameon=False, loc="upper right")
    ax[-1].set_xlabel("time from symbol boundary (µs)")
    fig.suptitle("fig_pm01 — constant-amplitude test signals: identical net phase step, different transition shapes")
    fig.tight_layout()
    fig.savefig(PM / "fig_pm01_case_waveforms.png", dpi=160)
    plt.close(fig)

    # pm02: phase and instantaneous frequency, first 12 us of modulation
    cases = ["CW", "OFF_+1.31", "OFF_-2.62", "QPSK_SLOW", "QPSK_REF", "QPSK_FAST", "QPSK_JUMP", "QPSK_LONGRAMP", "TOGGLE"]
    fig, ax = plt.subplots(len(cases), 2, figsize=(12, 13), sharex=True)
    t = np.arange(int((P.T0 - 1e-6) * 1e9), int((P.T0 + 12e-6) * 1e9)) * 1e-9
    for i, c in enumerate(cases):
        phi = P.phase(P.CASES[c], t)
        f = np.r_[0, np.diff(phi)] / (2 * np.pi * 1e-9)
        ax[i, 0].plot((t - P.T0) * 1e6, phi / np.pi, color=BLUE, lw=1)
        ax[i, 1].plot((t - P.T0) * 1e6, np.clip(f / 1e6, -6, 6), color=ORANGE, lw=1)
        ax[i, 0].set_ylabel(c, fontsize=7)
    ax[0, 0].set_title("phase φ(t) / π", fontsize=9)
    ax[0, 1].set_title("instantaneous frequency offset (MHz; jumps clipped at ±6)", fontsize=9)
    ax[-1, 0].set_xlabel("time from modulation start (µs)")
    ax[-1, 1].set_xlabel("time from modulation start (µs)")
    fig.suptitle("fig_pm02 — phase and instantaneous-frequency trajectories")
    fig.tight_layout()
    fig.savefig(PM / "fig_pm02_phase_and_instfreq.png", dpi=150)
    plt.close(fig)

    # pm03: steady compression by case
    order = ["CW"] + OFFSETS + ["QPSK_SLOW", "QPSK_REF", "QPSK_REF_s2", "QPSK_FAST", "QPSK_JUMP", "QPSK_LONGRAMP", "TOGGLE"]
    fig, ax = plt.subplots(figsize=(12, 4.2))
    xs = np.arange(len(order))
    col = [GREY if c == "CW" else AQUA if c.startswith("OFF") else BLUE for c in order]
    ax.bar(xs, [main_.loc[c, "g_final"] for c in order], color=col, width=.7)
    for c in P.CHECKS:
        i = order.index(c)
        for nd, dt, mk in ((8001, 1, "o"), (4001, .5, "s")):
            g = float(res[(res.case == c) & (res.Nd == nd) & (res.dt_ns == dt)].g_final.iloc[0])
            ax.plot(i, g, mk, mfc="none", mec="k", ms=6, label=f"check Nd={nd}, dt={dt:g} ns" if c == "CW" else None)
    ax.axhline(g_static, color=ORANGE, lw=1.5, label=f"static / Hammerstein / single-slow-state prediction (all cases) = {g_static:.3f}")
    ax.axhline(g_cw, color="k", lw=.8, ls=":", label=f"full model, CW = {g_cw:.3f}")
    ax.set_xticks(xs)
    ax.set_xticklabels(order, rotation=35, ha="right", fontsize=8)
    ax.set_ylabel("final gain g (large / small-signal)")
    ax.legend(fontsize=7, frameon=False, loc="upper left")
    ax.set_title("fig_pm03 — steady compression at constant |E| = %.3f V/m (+3 dB level); lower = more compressed" % P.A_H, fontsize=9)
    fig.tight_layout()
    fig.savefig(PM / "fig_pm03_steady_compression_by_case.png", dpi=160)
    plt.close(fig)

    # pm04: transition-rate dependence
    fig, ax = plt.subplots(figsize=(7, 4.2))
    for c in ["CW", "QPSK_SLOW", "QPSK_REF", "QPSK_REF_s2", "QPSK_FAST"]:
        ax.plot(main_.loc[c, "transition_rate_per_us"], main_.loc[c, "X"], "o", color=BLUE, ms=7)
        if c != "QPSK_REF_s2":
            ax.annotate(c.replace("QPSK_", ""), (main_.loc[c, "transition_rate_per_us"], main_.loc[c, "X"]), textcoords="offset points", xytext=(6, -12), fontsize=7)
    for c, mk in (("QPSK_JUMP", "^"), ("QPSK_LONGRAMP", "v"), ("TOGGLE", "D")):
        ax.plot(main_.loc[c, "transition_rate_per_us"], main_.loc[c, "X"], mk, color=YELLOW, ms=7, label=c)
    for c in OFFSETS:
        ax.plot(0, max(main_.loc[c, "X"], -.33), "v" if main_.loc[c, "X"] < -.33 else "x", color=AQUA, ms=7)
    ax.plot([], [], "x", color=AQUA, label="constant offsets, 0 transitions (X = %.2f … %.2f; ▼ off-scale, see fig_pm06)" % (main_.loc[OFFSETS, "X"].max(), main_.loc[OFFSETS, "X"].min()))
    ax.set_ylim(-.36, .8)
    xr = main_.loc[["CW", "QPSK_SLOW", "QPSK_REF", "QPSK_FAST"], ["transition_rate_per_us", "X"]].sort_values("transition_rate_per_us")
    ax.plot(xr.transition_rate_per_us, xr.X, color=BLUE, lw=1)
    ax.axhline(0, color="k", lw=.7)
    top = ax.secondary_xaxis("top", functions=(lambda r: r * 2.535, lambda x: x / 2.535))
    top.set_xlabel("transitions per τ_atom (2.535 µs)")
    ax.set_xlabel("nonzero phase transitions per µs")
    ax.set_ylabel("extra compression X = 1 − g/g_CW")
    ax.legend(fontsize=7, frameon=False, loc="lower right")
    ax.set_title("fig_pm04 — transition-rate dependence (300 ns ramps unless marked)", fontsize=9)
    fig.tight_layout()
    fig.savefig(PM / "fig_pm04_transition_rate_dependence.png", dpi=160)
    plt.close(fig)

    # pm05: ramp dependence at fixed rate and fixed net phase changes
    fig, ax = plt.subplots(1, 2, figsize=(11, 3.8))
    cs = ["QPSK_JUMP", "QPSK_REF", "QPSK_LONGRAMP"]
    ramps = [0, 300, 900]
    ax[0].plot(ramps, [main_.loc[c, "X"] for c in cs], "o-", color=BLUE)
    for r_, c in zip(ramps, cs):
        ax[0].annotate(c, (r_, main_.loc[c, "X"]), textcoords="offset points", xytext=(4, 6), fontsize=7)
    ax[0].set_xlabel("phase-ramp duration (ns)")
    ax[0].set_ylabel("extra compression X")
    ax[0].set_ylim(0, .8)
    ax[0].set_title("same symbols, rate, net phase steps; only ramp differs", fontsize=9)
    ax[1].semilogx([main_.loc[c, "peak_abs_freq_offset_Hz"] / 1e6 for c in cs], [main_.loc[c, "X"] for c in cs], "o-", color=BLUE, label="phase ramps (rate 0.73/µs)")
    ax[1].semilogx([main_.loc[c, "peak_abs_freq_offset_Hz"] / 1e6 for c in OFFSETS], [main_.loc[c, "X"] for c in OFFSETS], "x", color=AQUA, ms=8, label="constant offsets (|δ|)")
    ax[1].axhline(0, color="k", lw=.7)
    ax[1].set_xlabel("peak |instantaneous frequency offset| (MHz)")
    ax[1].set_ylabel("extra compression X")
    ax[1].legend(fontsize=7, frameon=False)
    ax[1].set_title("frequency swing does not organize X", fontsize=9)
    fig.suptitle("fig_pm05 — phase-ramp dependence (H2/H4)")
    fig.tight_layout()
    fig.savefig(PM / "fig_pm05_phase_ramp_dependence.png", dpi=160)
    plt.close(fig)

    # pm06: frequency-offset controls
    fig, ax = plt.subplots(1, 2, figsize=(11, 3.8))
    dl = [-2.618, -1.309, 0, .125, 1.309, 2.618]
    cs = ["OFF_-2.62", "OFF_-1.31", "CW", "OFF_+0.125", "OFF_+1.31", "OFF_+2.62"]
    ax[0].plot(dl, [main_.loc[c, "g_final"] for c in cs], "o-", color=AQUA, label="full model, constant offset")
    ax[0].axhline(main_.loc["QPSK_REF", "g_final"], color=BLUE, ls="--", label="full model, QPSK_REF")
    ax[0].axhline(g_static, color=ORANGE, label="static / Hammerstein / slow-state prediction")
    ax[0].set_xlabel("constant frequency offset δ (MHz)")
    ax[0].set_ylabel("final gain g (vs small-signal at the same frequency)")
    ax[0].legend(fontsize=7, frameon=False)
    ax[1].plot(dl, [main_.loc[c, "abs_out_last20"] / main_.loc["CW", "abs_out_last20"] for c in cs], "o-", color=BLUE, label="|large-signal output| / CW")
    ax[1].plot(dl, [main_.loc[c, "abs_lin_last20"] / main_.loc["CW", "abs_lin_last20"] for c in cs], "s--", color=GREY, label="|small-signal output| / CW")
    ax[1].set_xlabel("constant frequency offset δ (MHz)")
    ax[1].set_ylabel("relative to exact 5 MHz")
    ax[1].legend(fontsize=7, frameon=False)
    fig.suptitle("fig_pm06 — frequency-offset controls: detuning REDUCES compression; the nonlinear response is sharply frequency-selective at 5 MHz")
    fig.tight_layout()
    fig.savefig(PM / "fig_pm06_frequency_offset_controls.png", dpi=160)
    plt.close(fig)

    # pm07: transient trajectories + isolated-transition response
    fig, ax = plt.subplots(1, 3, figsize=(15, 4))
    for c, colr in (("CW", "k"), ("OFF_+0.125", AQUA), ("OFF_+1.31", AQUA), ("OFF_-2.62", AQUA)):
        z = load(c)
        ax[0].plot(z["t"] * 1e6 - T0_US, smooth(gain(z)), color=colr, lw=1.2, ls="-" if c != "OFF_-2.62" else "--", label=c)
    for c, colr, ls in (("QPSK_SLOW", BLUE, ":"), ("QPSK_REF", BLUE, "-"), ("QPSK_FAST", BLUE, "--"), ("QPSK_JUMP", YELLOW, "-"), ("QPSK_LONGRAMP", YELLOW, "--"), ("TOGGLE", ORANGE, "-")):
        z = load(c)
        ax[1].plot(z["t"] * 1e6 - T0_US, smooth(gain(z)), color=colr, lw=1.2, ls=ls, label=c)
    for a in ax[:2]:
        a.axhline(g_static, color=ORANGE, lw=.8, ls=":")
        a.set_xlim(-5, 59)
        a.set_xlabel("time from modulation start (µs)")
        a.legend(fontsize=7, frameon=False)
    ax[0].set_ylabel("gain g(t), 1 µs average")
    ax[0].set_title("CW and constant offsets", fontsize=9)
    ax[1].set_title("phase-modulated cases", fontsize=9)
    for dk, tt, gg in isolated_transition_response():
        ax[2].plot(tt, gg, color=BLUE if abs(dk) < 3 else ORANGE, lw=.9, alpha=.8)
    ax[2].plot([], [], color=BLUE, label="±π/2 steps")
    ax[2].plot([], [], color=ORANGE, label="π steps")
    ax[2].axvline(0, color="k", lw=.7)
    ax[2].set_xlabel("time from symbol boundary (µs)")
    ax[2].set_ylabel("gain g(t), 200 ns average")
    ax[2].set_title("isolated transitions (QPSK_SLOW, 4 µs apart)", fontsize=9)
    ax[2].legend(fontsize=7, frameon=False)
    fig.suptitle("fig_pm07 — each phase transition knocks the phase-locked response out of lock: brief spike (π) or small bump (±π/2), then a delayed dip that recovers over ~µs")
    fig.tight_layout()
    fig.savefig(PM / "fig_pm07_transient_response_comparison.png", dpi=160)
    plt.close(fig)

    # pm08: atomic-state evidence
    fig, ax = plt.subplots(1, 3, figsize=(15, 4))
    zc = load("CW")
    base_cw = smooth(zc["probe_db"])
    for c, colr, ls in (("QPSK_REF", BLUE, "-"), ("TOGGLE", ORANGE, "-"), ("QPSK_SLOW", BLUE, ":"), ("OFF_+1.31", AQUA, "-"), ("OFF_-2.62", AQUA, "--")):
        z = load(c)
        ax[0].plot(z["t"] * 1e6 - T0_US, smooth(z["probe_db"]) - base_cw, color=colr, ls=ls, lw=1.1, label=c)
    ax[0].set_xlim(-5, 59)
    ax[0].set_xlabel("time from modulation start (µs)")
    ax[0].set_ylabel("probe transmission baseline − CW (dB, 1 µs avg)")
    ax[0].legend(fontsize=7, frameon=False)
    order = ["CW"] + OFFSETS + ["QPSK_SLOW", "QPSK_REF", "QPSK_FAST", "QPSK_JUMP", "QPSK_LONGRAMP", "TOGGLE"]
    ryd = [(main_.loc[c, "pop_rydberg_3"] + main_.loc[c, "pop_rydberg_4"]) / (main_.loc["CW", "pop_rydberg_3"] + main_.loc["CW", "pop_rydberg_4"]) - 1 for c in order]
    ax[1].bar(range(len(order)), np.array(ryd) * 100, color=[GREY if c == "CW" else AQUA if c.startswith("OFF") else BLUE for c in order])
    ax[1].set_xticks(range(len(order)))
    ax[1].set_xticklabels(order, rotation=40, ha="right", fontsize=7)
    ax[1].set_ylabel("thermal-avg Rydberg population vs CW (%)")
    ax[1].set_title("populations at plateau end (110 µs)", fontsize=9)
    vx = zc["vx"]
    r_cw = zc["pops"][2] + zc["pops"][3]
    for c, colr in (("QPSK_REF", BLUE), ("TOGGLE", ORANGE), ("OFF_+1.31", AQUA)):
        z = load(c)
        ax[2].plot(vx, (z["pops"][2] + z["pops"][3]) - r_cw, color=colr, lw=.9, label=c)
    ax[2].set_xlim(-40, 40)
    ax[2].set_xlabel("velocity class v (m/s)")
    ax[2].set_ylabel("Rydberg population − CW (per class)")
    ax[2].legend(fontsize=7, frameon=False)
    ax[2].set_title("velocity-resolved (end-of-plateau snapshot)", fontsize=9)
    fig.suptitle("fig_pm08 — atomic-state proxies: populations change little; the extra compression lives in the coherent (phase-locked) response")
    fig.tight_layout()
    fig.savefig(PM / "fig_pm08_atomic_state_comparison.png", dpi=160)
    plt.close(fig)


if __name__ == "__main__":
    main()
