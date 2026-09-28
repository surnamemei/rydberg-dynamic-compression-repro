"""Figure 7: long-dwell behaviour (OP1, +3 dB).

(a) late-plateau gain (last 10 us of an isolated 60 us plateau), both carriers, 0 and +3 dB levels, vs CW-static gain
    (prespecified 2% test) - 12_plateau_steady_state.csv;
(b) D_ref (fixed-reference distortion; POST-HOC metric) vs T_dwell/tau_atom, both carriers, full model and surrogates,
    mean and 95% CI over 4 realizations - 15_followup_long_dwell_per_waveform.csv (metric D_ref_linear);
(c) fitted gain relative to the small-signal model vs dwell, per-realization ratio gain_abs(model)/gain_abs(linear),
    mean and t-based 95% CI (n = 4) - rows_fu_long_{cw,qpsk}_dwell.csv;
(d) PRESPECIFIED self-normalized D vs dwell (unmodulated carrier) and the failed criterion - 15_...csv (metric D_BLA).
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats

from figstyle import (fs, COL2, INK, INK2, MUTED, GRID, MODEL, MODEL_LABEL, MODEL_MARKER, XI, RES06, EB, BOLD,
                      panel, rec, save, relpath, zero_line, mline, tm)

Q = pd.read_csv(RES06 / "15_followup_long_dwell_per_waveform.csv")
PS = pd.read_csv(RES06 / "12_plateau_steady_state.csv").set_index("experiment")
SRC15 = relpath(RES06 / "15_followup_long_dwell_per_waveform.csv")
MODELS = [("atomic", "full"), ("static_mh", "static"), ("static_lti_mh", "lti")]
CAR = [("cw", "unmodulated carrier"), ("qpsk", "QPSK phase carrier")]
XT = [0.05, 1, 4, 8, 16, 32]
FS_TXT = fs(6.4)


def xaxis(ax, label=True, narrow=False):
    ax.set_xscale("log")
    ax.set_xlim(0.035, 48)
    xt = [0.05, 1, 4, 32] if narrow else XT
    ax.set_xticks(xt)
    ax.set_xticklabels([f"{v:g}" for v in xt])
    ax.xaxis.set_minor_formatter(plt.NullFormatter())
    if label:
        ax.set_xlabel(f"{XI} (ξ)")


def tci(x):
    x = np.asarray(x, float)
    n = len(x)
    m = x.mean()
    h = stats.t.ppf(0.975, n - 1) * x.std(ddof=1) / np.sqrt(n)
    return m, m - h, m + h


fig = plt.figure(figsize=(COL2, 4.7))
sf = fig.subfigures(2, 2, width_ratios=[0.78, 1.22], hspace=0.03, wspace=0.02)

# ------------------------------------------------------------------ (a) late-plateau gain (prespecified)
ax = sf[0, 0].subplots()
panel(ax, "a", "isolated plateau, OP1 (prespecified)")
for i, p in enumerate((0, 3)):
    ref = float(PS.loc[f"plateau_p{p:+d}_cw", "cw_static_gain_aH"])
    ax.plot([i - 0.32, i + 0.32], [ref, ref], color=INK2, lw=0.9, ls=(0, (3, 2)), zorder=1, label="CW-static gain" if i == 0 else None)
    for k, (car, mk, fc, lab) in enumerate((("cw", "o", MODEL["full"], "unmodulated carrier"), ("qpsk", "D", "white", "QPSK phase carrier"))):
        r = PS.loc[f"plateau_p{p:+d}_{car}"]
        x = i + (k - 0.5) * 0.3
        ax.plot([x], [r.gain_plateau_last10us], ls="none", marker=mk, ms=5, mfc=fc, mec=MODEL["full"], mew=0.9, zorder=3,
                label=lab if i == 0 else None)
        rel = 100 * r.plateau_vs_cw_static_rel
        ax.annotate(tm(f"{rel:+.2f}%" if abs(rel) < 1 else f"{rel:+.1f}%"), (x, r.gain_plateau_last10us),
                    xytext=(0, 6 if car == "cw" else -7), textcoords="offset points", ha="center", va="bottom" if car == "cw" else "top", fontsize=FS_TXT)
        rec("a", f"late-plateau gain {car} {p:+d} dB", p, float(r.gain_plateau_last10us), source=relpath(RES06 / "12_plateau_steady_state.csv"),
            note=f"rel to CW static {rel:+.3f}%")
    rec("a", f"CW-static gain {p:+d} dB", p, ref, source=relpath(RES06 / "12_plateau_steady_state.csv"))
ax.set_xticks([0, 1])
ax.set_xticklabels(["0 dB levels", "+3 dB levels"])
ax.set_xlim(-0.6, 1.6)
ax.set_ylim(0, 0.9)
ax.set_ylabel("late-plateau gain / small-signal gain")
ax.text(0.5, 0.97, "tolerance ±2% (unmodulated): passed", transform=ax.transAxes, ha="center", va="top", fontsize=FS_TXT, color=INK2)
ax.legend(loc="lower left", fontsize=fs(6.0), borderaxespad=0.2)

# ------------------------------------------------------------------ (b) D_ref (post hoc)
axs_b = sf[0, 1].subplots(1, 2, sharey=True)
panel(axs_b[0], "b")
sf[0, 1].suptitle("fixed-reference distortion $D_\\mathrm{ref}$ — post-hoc metric (OP1, +3 dB)", fontsize=fs(7.5), x=0.56, y=1.0)
q = Q[Q.metric == "D_ref_linear"]
for ax, (car, clab) in zip(axs_b, CAR):
    for m, key in MODELS:
        g = q[(q.carrier == car) & (q.model == m)].sort_values("xi")
        mline(ax, g.xi, g["mean"], key, yerr=[g["mean"] - g.ci_lo, g.ci_hi - g["mean"]], label=(car == "cw"), ms=3.3)
        rec("b", f"D_ref {m} {car}", g.xi.values, g["mean"].values, g.ci_lo.values, g.ci_hi.values, source=SRC15)
    ax.set_yscale("log")
    ax.set_ylim(2.5e-3, 3)
    xaxis(ax, narrow=True)
    ax.text(0.97, 0.97, clab, transform=ax.transAxes, ha="right", va="top", fontsize=fs(6.8))
axs_b[0].set_ylabel("$D_\\mathrm{ref}$ (post hoc)")
axs_b[0].annotate("+45%", xy=(4, 0.075), xytext=(1.6, 0.17), fontsize=FS_TXT, arrowprops=dict(arrowstyle="-", lw=0.5, color=INK2))
axs_b[0].annotate("+4.9% [+2.7%, +7.1%]", xy=(32, 0.035), xytext=(1.2, 0.012), fontsize=FS_TXT, arrowprops=dict(arrowstyle="-", lw=0.5, color=INK2))
h, lab = axs_b[0].get_legend_handles_labels()
axs_b[1].legend(h, lab, loc="lower right", fontsize=fs(5.9), borderaxespad=0.2)
rec("b", "gap annotations (ledger Q_refgap_cw_xi_4 / xi_32)", "text", 0.0, note="+45%; +4.9% [+2.7%, +7.1%]")

# ------------------------------------------------------------------ (c) fitted gain
axs_c = sf[1, 0].subplots(1, 2, sharey=True)
panel(axs_c[0], "c")
sf[1, 0].suptitle("fitted gain (OP1, +3 dB)", fontsize=fs(7.5), x=0.6, y=1.0)
for ax, (car, clab) in zip(axs_c, CAR):
    f = RES06 / f"rows_fu_long_{car}_dwell.csv"
    rows = pd.read_csv(f).drop_duplicates(["job_id", "model"], keep="last")
    lin = rows[rows.model == "linear"].set_index(["variant", "realization"]).gain_abs
    for m, key in MODELS:
        xs, ms_, lo, hi = [], [], [], []
        for v in ["xi_0.05", "xi_1", "xi_4", "xi_8", "xi_16", "xi_32"]:
            r = rows[(rows.model == m) & (rows.variant == v)].set_index("realization").gain_abs / lin.loc[v]
            mm, l, h = tci(r.values)
            xs.append(float(v[3:])); ms_.append(mm); lo.append(l); hi.append(h)
        xs, ms_, lo, hi = map(np.asarray, (xs, ms_, lo, hi))
        mline(ax, xs, ms_, key, yerr=[ms_ - lo, hi - ms_], label=False, ms=3.3)
        rec("c", f"gain {m} {car}", xs, ms_, lo, hi, source=relpath(f), note="per-realization gain_abs/linear gain_abs; t 95% CI, n=4")
    ax.set_ylim(0, 0.85)
    xaxis(ax, narrow=True)
    ax.text(0.97, 0.97, clab, transform=ax.transAxes, ha="right", va="top", fontsize=fs(6.8))
axs_c[0].set_ylabel("fitted gain / small-signal gain")

# ------------------------------------------------------------------ (d) prespecified D and the failed criterion
ax = sf[1, 1].subplots()
panel(ax, "d", "prespecified self-normalized D, unmodulated carrier (OP1, +3 dB)")
b = Q[(Q.metric == "D_BLA") & (Q.carrier == "cw")]
for m, key in MODELS:
    g = b[b.model == m].sort_values("xi")
    mline(ax, g.xi, g["mean"], key, yerr=[g["mean"] - g.ci_lo, g.ci_hi - g["mean"]], label=True, ms=3.3)
    rec("d", f"D (prespecified) {m} cw", g.xi.values, g["mean"].values, g.ci_lo.values, g.ci_hi.values, source=SRC15)
gap = b[b.model == "full_minus_static_mh"].set_index("xi")["mean"]
rec("d", "full-static gap", [4.0, 32.0], [gap[4.0], gap[32.0]], source=SRC15, note=f"fraction {gap[32.0] / gap[4.0]:.3f}")
ax.set_yscale("log")
ax.set_ylim(0.05, 150)
xaxis(ax)
ax.set_ylabel("D (prespecified, self-normalized)")
st = b[(b.model == "static_mh") & (b.xi == 0.05)]["mean"].iloc[0]
ax.annotate(f"static D ≈ {st:.0f}: denominator collapses\n(ill-conditioned for a compressed\ntwo-level envelope)", xy=(0.05, st), xytext=(0.12, 30),
            fontsize=FS_TXT, arrowprops=dict(arrowstyle="-", lw=0.5, color=INK2))
ax.text(0.99, 0.60, f"criterion: full − static gap at ξ = 32 < 25% of its ξ = 4 value\n" +
        tm(f"observed {gap[4.0]:+.3f}") + " → " + tm(f"{gap[32.0]:+.3f}") + f" ({100 * gap[32.0] / gap[4.0]:.0f}% remains): not met",
        transform=ax.transAxes, ha="right", va="center", fontsize=FS_TXT, color=INK,
        bbox=dict(boxstyle="round,pad=0.25", fc="white", ec=INK2, lw=0.5))
ax.legend(loc="upper right", fontsize=fs(6.0))

save(fig, "fig7")
