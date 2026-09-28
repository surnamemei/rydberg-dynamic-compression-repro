"""Figure 4: temporal order at identical amplitude samples (OP1, realization 0 for (a)-(c)).

(a) envelopes: stored decimated inputs of traces/main (dwell xi 0.05 / 4; shuffle original / cycle block-shuffle);
(b) amplitude histograms (identical) and time-weighted dwell survival: full-rate envelopes regenerated from the stored
    seeds (stimuli.py; bit-identical to the stored decimated inputs), thresholds as in 04_controlled_pair_configs.csv;
(c) |baseband output| at +3 dB, xi 0.05 and 4: stored traces (full model, all-zone static, all-zone LTI+static, linear);
(d) cycle block-shuffle minus original, paired change in D, 8 realizations, 95% CIs (06_temporal_shuffle_results.csv).
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from figstyle import (fs, COL2, INK, INK2, MUTED, GRID, MODEL, MODEL_LABEL, MODEL_MARKER, PP, TAU, RES06, EB,
                      panel, rec, save, relpath, zero_line)
import stimuli as S

TR = RES06 / "traces" / "main"
FS_TXT = fs(6.5)
fig = plt.figure(figsize=(COL2, 5.0))
sf = fig.subfigures(2, 2, width_ratios=[1.12, 1], height_ratios=[1, 1.05], wspace=0.02, hspace=0.02)

# ------------------------------------------------------------------ (a) envelopes (stored decimated inputs)
axs_a = list(sf[0, 0].subplots(4, 1, gridspec_kw=dict(hspace=0.0)))
traces_a = [("dwell", "xi_0.05", "dwell family: ξ = 0.05 (reference), OP1 stimuli"), ("dwell", "xi_4", "dwell family: ξ = 4"),
            ("shuffle", "ORIGINAL", "clustered envelope (original)"), ("shuffle", "BLOCK_SHUFFLED_CYCLES", "cycle block-shuffle")]
T_WIN = 200.0
for i, (ax, (des, v, lab)) in enumerate(zip(axs_a, traces_a)):
    f = TR / f"{des}_r0_{v}_p+3_Nd4001_dt1.npz"
    a = np.load(f)["a"]
    t = np.arange(len(a)) * 0.05
    sel = t < T_WIN
    ax.plot(t[sel], a[sel], color=INK if des == "dwell" else INK2, lw=0.45)
    ax.text(1.0, 1.0, lab, transform=ax.transAxes, ha="right", va="bottom", fontsize=FS_TXT)
    top = 2.2 if des == "dwell" else 1.08 * float(a[sel].max())
    ax.set_ylim(0, top)
    ax.set_xlim(0, T_WIN)
    ax.set_yticks([0, 1] if des == "dwell" else [0, 2, 4])
    ax.spines["bottom"].set_visible(i == 3)
    ax.tick_params(axis="x", bottom=(i == 3), labelbottom=(i == 3))
    rec("a", f"envelope {des} {v} max in window", "max", float(a[sel].max()), source=relpath(f), note="stored decimated input a[::50]")
axs_a[-1].set_xlabel("time (µs)")
sf[0, 0].supylabel("|E| / RMS", fontsize=fs(7.5))
panel(axs_a[0], "a")

# ------------------------------------------------------------------ (b) histograms and dwell survival (regenerated full-rate envelopes)
axh, axd = sf[0, 1].subplots(1, 2)
panel(axh, "b")
env = {"dwell": S.envelopes("dwell", 0), "shuffle": S.envelopes("shuffle", 0)}
ver = {d: S.verify(d, env[d], 0) for d in env}
assert all(v == 0.0 for d in ver for v in ver[d].values()), ver
cfg = pd.read_csv(RES06 / "04_controlled_pair_configs.csv")
cfg0 = cfg[cfg.realization == 0].set_index(["design", "variant"])
bins = np.linspace(0, 9.0, 181)
pairs = [("dwell", "xi_0.05", "xi_4", INK, "ξ = 0.05", "ξ = 4"), ("shuffle", "ORIGINAL", "BLOCK_SHUFFLED_CYCLES", MUTED, "original", "cycle block-shuffle")]
for des, v0, v1, col, l0, l1 in pairs:
    a0, a1 = env[des][v0], env[des][v1]
    assert np.array_equal(np.sort(a0), np.sort(a1))
    h0, _ = np.histogram(a0, bins=bins, density=True)
    h1, _ = np.histogram(a1, bins=bins, density=True)
    lw_h = 0.7 if des == "dwell" else 1.0
    axh.stairs(h0, bins, color=col, lw=lw_h, label=l0)
    axh.stairs(h1, bins, color=col, lw=lw_h, ls=(0, (2.2, 1.6)), label=l1, baseline=None)
    rec("b", f"histogram identical {des}", "max|h0-h1|", float(np.max(np.abs(h0 - h1))), note="regenerated full-rate envelopes; sorted samples identical")
    thr = float(np.median(a0)) if des == "shuffle" else float((a0.min() + a0.max()) / 2)
    for v, ls, lab in ((v0, "-", l0), (v1, (0, (2.2, 1.6)), l1)):
        runs = np.sort(S.dwell_runs(env[des][v], thr)).astype(float)
        wts = np.cumsum(runs[::-1])[::-1] / runs.sum()
        axd.step(runs / S.TAU_SAMPLES, wts, where="post", color=col, lw=1.0, ls=ls)
        tw = float(np.sum(runs ** 2) / np.sum(runs)) / S.TAU_SAMPLES
        rec("b", f"time-weighted dwell/tau {des} {v}", v, tw, source=relpath(RES06 / "04_controlled_pair_configs.csv"),
            note=f"stored T_dwell_over_tau = {cfg0.loc[(des, v), 'T_dwell_over_tau']:.6f}")
axh.set_yscale("log")
axh.set_xlim(0, 9)
axh.set_ylim(1e-4, 30)
axh.set_xlabel("|E| / RMS")
axh.set_ylabel("probability density")
axh.legend(loc="upper right", fontsize=fs(6.0), handlelength=2.0)
axd.set_xscale("log")
axd.set_xlim(2.5e-3, 30)
axd.set_ylim(0, 1.03)
axd.set_xlabel(f"high-field dwell / {TAU}")
axd.set_ylabel("fraction of high-field time\nin runs ≥ dwell")
axd.text(0.03, 0.03, "ξ = 0.05", transform=axd.transAxes, fontsize=fs(6.0), color=INK)
axd.annotate("ξ = 4", xy=(1.5, 0.93), xytext=(0.5, 0.80), fontsize=fs(6.0), color=INK, ha="center",
             arrowprops=dict(arrowstyle="-", lw=0.5, color=INK2))
axd.annotate("original ≈\nblock-shuffle", xy=(0.9, 0.55), xytext=(0.11, 0.42), fontsize=fs(6.0), color=INK2, ha="center",
             arrowprops=dict(arrowstyle="-", lw=0.5, color=INK2))

# ------------------------------------------------------------------ (c) baseband outputs at +3 dB (stored traces)
axc = list(sf[1, 0].subplots(2, 1))
panel(axc[0], "c")
W0, W1 = 8.0, 38.0
for i, v in enumerate(("xi_0.05", "xi_4")):
    f = TR / f"dwell_r0_{v}_p+3_Nd4001_dt1.npz"
    z = np.load(f)
    cp = int(z["cp"])
    t = np.arange(len(z["x"]) - cp) * 0.05
    w0, w1 = (S.CFG["cyclic_prefix_samples"] + S.CFG["analysis_skip_samples"]) // S.DEC - cp, len(t) - 2000 // S.DEC
    norm = float(np.sqrt(np.mean(np.abs(z["y_linear"][cp:][w0:w1]) ** 2)))
    sel = (t >= W0) & (t <= W1)
    ax = axc[i]
    ax.plot(t[sel], np.abs(z["y_linear"][cp:][sel]) / norm, color=MUTED, lw=0.8, label="small-signal (linear) model")
    for m, key in (("static_mh", "static"), ("static_lti_mh", "lti"), ("atomic", "full")):
        y = np.abs(z["y_" + m][cp:]) / norm
        ax.plot(t[sel], y[sel], color=MODEL[key], lw=0.9 if key != "full" else 1.1, label=MODEL_LABEL[key], zorder=3 if key == "full" else 2)
        rec("c", f"|y| {m} {v} window mean", "mean", float(np.mean(y[sel])), source=relpath(f), note="normalized by RMS |y_linear| over the analysis window")
    ax.set_xlim(W0, W1)
    ax.set_ylim(0, 2.8)
    ax.text(1.0, 1.0, {"xi_0.05": "ξ = 0.05", "xi_4": "ξ = 4"}[v] + ", +3 dB, OP1", transform=ax.transAxes, ha="right", va="bottom", fontsize=FS_TXT)
    if i == 0:
        ax.tick_params(axis="x", labelbottom=False)
axc[1].set_xlabel("time (µs)")
sf[1, 0].supylabel("|output| / RMS small-signal output", fontsize=fs(7.5))
axc[1].legend(loc="upper right", fontsize=fs(6.0), ncol=1, handlelength=1.4, borderaxespad=0.1)

# ------------------------------------------------------------------ (d) cycle block-shuffle change in D
axd2 = sf[1, 1].subplots()
panel(axd2, "d", "cycle block-shuffle minus original, OP1")
Sh = pd.read_csv(RES06 / "06_temporal_shuffle_results.csv")
g = Sh[(Sh.metric == "D_BLA") & (Sh.variant == "BLOCK_SHUFFLED_CYCLES")]
zero_line(axd2)
for k, (m, key) in enumerate((("atomic", "full"), ("static_mh", "static"), ("static_lti_mh", "lti"))):
    r = g[g.model == m].sort_values("Pavg_over_P1dB_dB")
    x = r.Pavg_over_P1dB_dB.values + (k - 1) * 0.45
    axd2.errorbar(x, r.Delta_mean, yerr=[r.Delta_mean - r.Delta_ci_lo, r.Delta_ci_hi - r.Delta_mean], color=MODEL[key], marker=MODEL_MARKER[key],
                  ms=3.8, mec="white", mew=0.5, lw=0, label=MODEL_LABEL[key], zorder=3, **EB)
    rec("d", f"dD cycle-shuffle {m}", r.Pavg_over_P1dB_dB.values, r.Delta_mean.values, r.Delta_ci_lo.values, r.Delta_ci_hi.values,
        source=relpath(RES06 / "06_temporal_shuffle_results.csv"))
axd2.set_xticks([-6, 0, 3, 6])
axd2.set_xticklabels(["−6", "0", "+3", "+6"])
axd2.set_xlim(-7.5, 7.5)
axd2.set_ylim(-0.2, 0.15)
axd2.set_xlabel(f"{PP} (dB)")
axd2.set_ylabel("paired change in D")
axd2.legend(loc="lower left", fontsize=fs(6.3))

save(fig, "fig4")
