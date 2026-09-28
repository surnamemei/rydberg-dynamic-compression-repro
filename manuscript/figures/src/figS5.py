"""Supplementary Figure S5: self-normalized distortion change vs dwell (OP1, QPSK phase carrier unless stated).

(a) paired change in D (member minus xi = 0.05 reference) vs nominal T_dwell/tau_atom at 0 and +3 dB, full model and the two
    all-zone surrogates, 8 realizations, 95% CIs; inset -6 dB; top axis x 2.25 (05_controlled_pair_results.csv);
(b) single-slow-state surrogate vs full model for the 12 prespecified dwell contrasts and 2 clustering contrasts, with the
    +/-25% relative-error tolerance (16_followup_single_slow_state.csv).
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from figstyle import (fs, COL2, INK, INK2, MUTED, MODEL, MODEL_LABEL, MODEL_MARKER, XI, TAU, E1, RES06, EB,
                      panel, rec, save, relpath, zero_line, mline)

P = pd.read_csv(RES06 / "05_controlled_pair_results.csv")
P = P[(P.metric == "D_BLA") & P.variant.str.startswith("xi_")].copy()
P["xi"] = P.variant.str[3:].astype(float)
SRC = relpath(RES06 / "05_controlled_pair_results.csv")
MODELS = [("atomic", "full"), ("static_mh", "static"), ("static_lti_mh", "lti")]
F225 = 2.25  # 2.535 us / 1.127 us (1/e time at the converged E1dB; 03b_tau_atom_sensitivity.csv)

fig = plt.figure(figsize=(COL2, 3.25))
gs = fig.add_gridspec(1, 3, width_ratios=[1, 1, 0.95])
axs = [fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1])]
axb = fig.add_subplot(gs[0, 2])


def draw(ax, p, ms=3.4, lw=1.1, label=True, record=True):
    for m, key in MODELS:
        g = P[(P.model == m) & (P.Pavg_over_P1dB_dB == p)].sort_values("xi")
        mline(ax, g.xi, g.Delta_mean, key, yerr=[g.Delta_mean - g.Delta_ci_lo, g.Delta_ci_hi - g.Delta_mean], label=label, ms=ms, lw=lw)
        if record:
            rec("a", f"dD {m} {p:+d} dB", g.xi.values, g.Delta_mean.values, g.Delta_ci_lo.values, g.Delta_ci_hi.values, source=SRC)


for ax, p, lab in ((axs[0], 0, "0 dB"), (axs[1], 3, "+3 dB")):
    zero_line(ax)
    draw(ax, p, label=(p == 3))
    ax.set_xscale("log")
    ax.set_xlim(0.07, 5.6)
    ax.set_xticks([0.1, 0.25, 0.5, 1, 2, 4])
    ax.set_xticklabels(["0.1", "0.25", "0.5", "1", "2", "4"])
    ax.xaxis.set_minor_formatter(plt.NullFormatter())
    ax.set_xlabel(f"{XI} (ξ)")
    sec = ax.secondary_xaxis("top", functions=(lambda x: x * F225, lambda x: x / F225))
    sec.set_xticks([0.25, 0.5, 1, 2, 4, 8])
    sec.set_xticklabels(["0.25", "0.5", "1", "2", "4", "8"])
    sec.xaxis.set_minor_formatter(plt.NullFormatter())
    sec.tick_params(labelsize=fs(6.3))
    sec.set_xlabel(r"$T_\mathrm{dwell}/\tau_{1/e}$ at converged " + E1 + " (= 2.25 ξ)", fontsize=fs(6.5), labelpad=2)
    ax.text(0.03, 0.97, f"{lab}, OP1", transform=ax.transAxes, ha="left", va="top", fontsize=fs(7))
panel(axs[0], "a")
axs[0].set_ylabel("paired change in D (member − reference)")
axs[0].set_ylim(-0.12, 0.23)
axs[1].set_ylim(-1.55, 0.52)
axs[1].legend(loc="center right", fontsize=fs(6.3), bbox_to_anchor=(1.0, 0.5))

# inset: -6 dB (in the 0 dB facet)
ins = axs[0].inset_axes([0.14, 0.50, 0.40, 0.30])
zero_line(ins, lw=0.5)
draw(ins, -6, ms=2.2, lw=0.8, label=False)
ins.set_xscale("log")
ins.set_xlim(0.07, 5.6)
ins.set_ylim(-0.19, 0.03)
ins.set_xticks([0.1, 1])
ins.set_xticklabels(["0.1", "1"])
ins.xaxis.set_minor_formatter(plt.NullFormatter())
ins.tick_params(labelsize=fs(5.8), length=1.8, pad=1)
ins.set_title("−6 dB", fontsize=fs(6.2), pad=1.5)
for s in ("top", "right"):
    ins.spines[s].set_visible(True)
    ins.spines[s].set_linewidth(0.5)

# ------------------------------------------------------------------ (b) single-slow-state surrogate
panel(axb, "b", "single-slow-state, OP1")
D = pd.read_csv(RES06 / "16_followup_single_slow_state.csv")
lim = (-0.26, 0.42)
xx = np.linspace(0, lim[1], 2)
axb.fill_between(xx, 0.75 * xx, 1.25 * xx, color="#f3f2ee", lw=0, zorder=0, gid="decor")
xn = np.linspace(lim[0], 0, 2)
axb.fill_between(xn, 1.25 * xn, 0.75 * xn, color="#f3f2ee", lw=0, zorder=0, gid="decor")
axb.plot(lim, lim, color=INK2, lw=0.6, zorder=1)
axb.axhline(0, color=MUTED, lw=0.4, zorder=0)
axb.axvline(0, color=MUTED, lw=0.4, zorder=0)
groups = [("QPSK carrier dwell", "o", MODEL["dsh"], "dwell, QPSK carrier (8 realizations)"),
          ("CW carrier dwell", "o", "white", "dwell, unmodulated carrier (4 realizations)"),
          ("shuffle cycles", "D", MODEL["dsh"], "clustering (cycle block-shuffle)")]
for key, mk, fc, lab in groups:
    g = D[D.set.str.startswith(key)]
    axb.plot(g.Delta_full, g.Delta_DSH, ls="none", marker=mk, ms=4.2, mfc=fc, mec=INK2, mew=0.6, label=lab, zorder=3)
    rec("b", f"DSH vs full: {key}", g.Delta_full.values, g.Delta_DSH.values, source=relpath(RES06 / "16_followup_single_slow_state.csv"),
        note="within_25pct=" + ",".join(str(bool(v)) for v in g.within_25pct))
n_dw = int(D[~D.set.str.startswith("shuffle")].within_25pct.sum())
n_cl = int(D[D.set.str.startswith("shuffle")].within_25pct.sum())
rec("b", "dwell contrasts within 25%", "count", n_dw, note="of 12")
rec("b", "clustering contrasts within 25%", "count", n_cl, note="of 2")
axb.text(0.40, 0.28, "±25% tolerance", fontsize=fs(6.2), color=INK2, ha="right", rotation=0)
axb.text(0.40, -0.02, f"dwell: {n_dw} of 12 within", fontsize=fs(6.3), ha="right", va="top", color=INK)
axb.text(0.03, -0.245, f"clustering (diamonds):\n{n_cl} of 2 within", fontsize=fs(6.3), ha="left", va="bottom", color=INK)
axb.set_xlim(*lim)
axb.set_ylim(*lim)
axb.set_aspect("equal")
axb.set_xlabel("paired change in D, full model")
axb.set_ylabel("paired change in D, single-slow-state")
axb.legend(loc="upper center", bbox_to_anchor=(0.5, -0.24), fontsize=fs(5.9), handletextpad=0.3, borderaxespad=0.0, ncol=1)

save(fig, "figS5", supplement=True)
