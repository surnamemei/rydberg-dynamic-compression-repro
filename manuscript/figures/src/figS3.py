"""Supplementary Figure S3: dwell dependence at -6 and +6 dB (OP1, QPSK phase carrier).

Paired change in D (member minus xi = 0.05 reference), full model and the two all-zone surrogates, 8 realizations,
95% CIs (05_controlled_pair_results.csv).
"""
import pandas as pd
import matplotlib.pyplot as plt

from figstyle import (fs, COL1, XI, RES06, panel, rec, save, relpath, zero_line, mline)

P = pd.read_csv(RES06 / "05_controlled_pair_results.csv")
P = P[(P.metric == "D_BLA") & P.variant.str.startswith("xi_")].copy()
P["xi"] = P.variant.str[3:].astype(float)
fig = plt.figure(figsize=(COL1, 2.5))
axs = fig.subplots(1, 2)
for ax, p, letter in ((axs[0], -6, "a"), (axs[1], 6, "b")):
    panel(ax, letter, f"{'−6' if p < 0 else '+6'} dB, OP1")
    zero_line(ax)
    for m, key in (("atomic", "full"), ("static_mh", "static"), ("static_lti_mh", "lti")):
        g = P[(P.model == m) & (P.Pavg_over_P1dB_dB == p)].sort_values("xi")
        mline(ax, g.xi, g.Delta_mean, key, yerr=[g.Delta_mean - g.Delta_ci_lo, g.Delta_ci_hi - g.Delta_mean], label=(p > 0), ms=3.0, lw=1.0)
        rec(letter, f"dD {m} {p:+d} dB", g.xi.values, g.Delta_mean.values, g.Delta_ci_lo.values, g.Delta_ci_hi.values,
            source=relpath(RES06 / "05_controlled_pair_results.csv"))
    ax.set_xscale("log")
    ax.set_xlim(0.07, 5.6)
    ax.set_xticks([0.1, 1, 4])
    ax.set_xticklabels(["0.1", "1", "4"])
    ax.xaxis.set_minor_formatter(plt.NullFormatter())
    ax.set_xlabel(f"{XI} (ξ)")
axs[0].set_ylabel("paired change in D")
axs[0].set_ylim(-0.19, 0.03)
axs[1].set_ylim(-0.95, 0.62)
axs[1].legend(loc="lower right", fontsize=fs(5.6), borderaxespad=0.2, bbox_to_anchor=(1.0, 0.1))

save(fig, "figS3", supplement=True)
