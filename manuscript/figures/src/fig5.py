"""Figure 5: what temporal order changes at fixed amplitude statistics (OP1; gain decomposition post hoc).

(a) Fitted gain relative to the small-signal model vs dwell, QPSK phase carrier, +3 dB (n = 8; revision P1 re-run).
(b) The same for the unmodulated carrier (archived carrier-ablation rows, n = 4).
(c) Fixed-reference residual distortion D_ref vs dwell, QPSK carrier, +3 dB (2-us FIR, n = 8).
(d) Declustering: paired change in D_ref (cycle block-shuffle minus original) vs FIR span, +3 and +6 dB, full model and
    the two all-zone surrogates (n = 8; results/revision/p1/02_p1_contrasts.csv, lag columns taken from the periodic input).
The prespecified FIR-span/normalization robustness table is Supplementary Table S7 (results/revision/p1/04_*).
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy import stats

from figstyle import (fs, COL2, XI, ROOT, INK, INK2, MODEL, MODEL_MARKER, EB, panel, rec, save, relpath, zero_line, mline)

P1 = ROOT / "results" / "revision" / "p1"
G = pd.read_csv(P1 / "07_p1_posthoc_gain_vs_dwell.csv")
M = pd.read_csv(P1 / "01_p1_metric_rows.csv")
C = pd.read_csv(P1 / "02_p1_contrasts.csv")
KEY = {"atomic": "full", "static_mh": "static", "static_lti_mh": "lti", "dsh_tau2.535us": "dsh"}

fig = plt.figure(figsize=(COL2, 2.75))
axs = fig.subplots(1, 4, gridspec_kw=dict(width_ratios=[1, 1, 1, 1.08])).ravel()


def xi_axis(ax):
    ax.set_xscale("log")
    ax.set_xticks([0.05, 0.25, 1, 4])
    ax.set_xticklabels(["0.05", "0.25", "1", "4"])
    ax.xaxis.set_minor_formatter(plt.NullFormatter())
    ax.set_xlabel(f"{XI} (ξ)")


for ax, car, letter, title in ((axs[0], "qpsk", "a", "QPSK carrier"),
                               (axs[1], "cw", "b", "unmodulated")):
    panel(ax, letter, title)
    for m in ("atomic", "dsh_tau2.535us", "static_lti_mh", "static_mh"):
        g = G[(G.carrier == car) & (G.Pavg_over_P1dB_dB == 3) & (G.model == m)].sort_values("xi")
        mline(ax, g.xi, g.gain_rel_mean, KEY[m], yerr=[g.gain_rel_mean - g.gain_rel_lo, g.gain_rel_hi - g.gain_rel_mean],
              label=(letter == "a"), ms=3.0, lw=1.0)
        rec(letter, f"gain_rel {m} {car} +3 dB", g.xi.values, g.gain_rel_mean.values, g.gain_rel_lo.values, g.gain_rel_hi.values,
            source=relpath(P1 / "07_p1_posthoc_gain_vs_dwell.csv"))
    xi_axis(ax)
    ax.set_ylim(0, 1.0)
    ax.set_ylabel("fitted gain / small-signal gain")

ax = axs[2]
panel(ax, "c", "residual")
d = M[(M.design == "dwell") & (M.Pavg_over_P1dB_dB == 3) & (M.method == "context") & (M.span_us == 2)].copy()
d["xi"] = d.variant.str[3:].astype(float)
for m in ("atomic", "static_lti_mh", "static_mh"):
    s = d[d.model == m].groupby("xi").D_ref
    mu, sd, n = s.mean(), s.std(ddof=1), s.count()
    h = stats.t.ppf(.975, n - 1) * sd / np.sqrt(n)
    mline(ax, mu.index, mu.values, KEY[m], yerr=[h.values, h.values], label=False, ms=3.0, lw=1.0)
    rec("c", f"D_ref {m} +3 dB", mu.index.values, mu.values, (mu - h).values, (mu + h).values, source=relpath(P1 / "01_p1_metric_rows.csv"))
xi_axis(ax)
ax.set_yscale("log")
ax.set_yticks([0.02, 0.05, 0.1, 0.2, 0.5])
ax.set_yticklabels(["0.02", "0.05", "0.1", "0.2", "0.5"])
ax.yaxis.set_minor_formatter(plt.NullFormatter())
ax.set_ylim(0.015, 0.7)
ax.set_ylabel("$D_\\mathrm{ref}$ = residual / linear-model power")

ax = axs[3]
panel(ax, "d", "declustering")
zero_line(ax)
sh = C[(C.design == "shuffle") & (C.method == "context") & (C.metric == "D_ref")]
off = {"atomic": 0.0, "static_lti_mh": -0.09, "static_mh": 0.09}
col = {"atomic": "d_atomic", "static_lti_mh": "d_static_lti_mh", "static_mh": "d_static_mh"}
for pdb, filled in ((3, True), (6, False)):
    for m in ("atomic", "static_lti_mh", "static_mh"):
        g = sh[sh.Pavg_over_P1dB_dB == pdb].sort_values("span_us")
        x = np.log2(g.span_us.values) + off[m] + (0.0 if pdb == 3 else 0.03)
        y, lo, hi = g[col[m]].values, g[col[m] + "_lo"].values, g[col[m] + "_hi"].values
        ax.errorbar(x, y, yerr=[y - lo, hi - y], color=MODEL[KEY[m]], marker=MODEL_MARKER[KEY[m]], ms=3.2, ls="none",
                    mfc=MODEL[KEY[m]] if filled else "white", mec=MODEL[KEY[m]] if not filled else "white", mew=0.8 if not filled else 0.5, **EB)
        rec("d", f"dD_ref shuffle {m} {pdb:+d} dB", g.span_us.values, y, lo, hi, source=relpath(P1 / "02_p1_contrasts.csv"))
ax.set_xticks(np.log2([2, 4, 8, 16]))
ax.set_xticklabels(["2", "4", "8", "16"])
ax.set_xlim(0.6, 4.4)
ax.set_ylim(-0.0115, 0.0068)
ax.set_xlabel("FIR span (µs)")
ax.set_ylabel("paired change in $D_\\mathrm{ref}$")
hh, ll = axs[0].get_legend_handles_labels()
hh += [Line2D([], [], color=INK2, marker="o", ls="none", ms=3.2, mec="white", mew=0.5),
       Line2D([], [], color=INK2, marker="o", ls="none", ms=3.2, mfc="white", mec=INK2, mew=0.8)]
ll += ["(d) +3 dB", "(d) +6 dB"]
fig.legend(hh, ll, loc="outside lower center", ncol=6, fontsize=fs(6.2), handlelength=1.4, columnspacing=1.0, handletextpad=0.4)

save(fig, "fig5")
