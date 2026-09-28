"""Figure 10: CW-derived surrogates vs one standard identified memory model (GMP), held-out test data (prespecified;
results/tqe_viability_final/00_plan.md, section C; results/tqe_viability_final/03_memory_model/gmp_verdict.json).

(a) Median test NMSE vs the full model by waveform family.
(b) Mean |AIR - AIR_full| at +3 dB by format (test realizations; archived noise).
(c) Held-out contrasts relative to the full model (model / full): dwell gain contrast, declustering D_ref contrast, and X, X_coh,
    X_mag of the reference phase case (test sequences).
(d) X_coh of the reference phase case at held-out drive levels (revision P3 runs, three sequences).
"""
import json

import numpy as np
import matplotlib.pyplot as plt

from figstyle import (fs, COL2, ROOT, INK2, MUTED, MODEL, MODEL_LABEL, MODEL_MARKER, panel, rec, save, relpath, zero_line)

TVF = ROOT / "results" / "tqe_viability_final"
V = json.loads((TVF / "03_memory_model" / "gmp_verdict.json").read_text())
SRC = relpath(TVF / "03_memory_model" / "gmp_verdict.json")
GMP = "#4a3aa7"
MODS = [("gmp", "GMP", GMP, "*"), ("static_mh", MODEL_LABEL["static"], MODEL["static"], MODEL_MARKER["static"]),
        ("static_lti_mh", MODEL_LABEL["lti"], MODEL["lti"], MODEL_MARKER["lti"]), ("dsh", MODEL_LABEL["dsh"], MODEL["dsh"], MODEL_MARKER["dsh"])]

fig, axs = plt.subplots(1, 4, figsize=(COL2, 2.55), gridspec_kw={"width_ratios": [1.0, 1.0, 1.25, 1.0]})

# ------------------------------------------------------------------ (a) NMSE by family
ax = axs[0]
panel(ax, "a", "waveform error")
fams = [("dwell", "dwell"), ("shuffle", "shuffle"), ("phase", "phase"), ("fair", "fair")]
w = 0.19
for j, (m, lab, c, _) in enumerate(MODS):
    x = np.arange(len(fams)) + (j - 1.5) * w
    y = [V["nmse"][f][m] for f, _ in fams]
    ax.bar(x, y, w * 0.92, color=c, label=lab, zorder=2)
    rec("a", f"median NMSE {m}", [f for f, _ in fams], y, source=SRC)
ax.set_yscale("log")
ax.set_ylim(0.02, 8.0)
ax.set_xticks(range(len(fams)))
ax.set_xticklabels([l_ for _, l_ in fams], fontsize=fs(6.5), rotation=30, ha="right", rotation_mode="anchor")
ax.set_ylabel("median test NMSE vs full model")
ax.legend(loc="upper center", fontsize=fs(5.8), borderaxespad=0.2, handlelength=0.9, ncol=2, columnspacing=0.6)

# ------------------------------------------------------------------ (b) AIR error by format
ax = axs[1]
panel(ax, "b", "AIR error")
fm = [("QPSK", "QPSK"), ("16QAM", "16-QAM"), ("OFDM", "OFDM"), ("CE", "CE")]
for j, (m, lab, c, _) in enumerate(MODS):
    x = np.arange(len(fm)) + (j - 1.5) * w
    y = [V["air_err"][f][m] for f, _ in fm]
    ax.bar(x, y, w * 0.92, color=c, zorder=2)
    rec("b", f"mean |AIR error| {m}", [f for f, _ in fm], y, source=SRC)
ax.axhline(0.05, color=INK2, lw=0.6, ls=(0, (1.2, 1.4)), zorder=0)
ax.set_xticks(range(len(fm)))
ax.set_xticklabels([l_ for _, l_ in fm], fontsize=fs(6.5), rotation=30, ha="right", rotation_mode="anchor")
ax.set_ylim(0, 1.0)
ax.set_ylabel("mean |AIR − AIR$_\\mathrm{full}$| (bit)")

# ------------------------------------------------------------------ (c) contrasts relative to the full model
ax = axs[2]
panel(ax, "c", "contrasts / full")
keys = [("dwell_gain", "dwell\nΔgain"), ("decluster_Dref", "declust.\nΔD$_\\mathrm{ref}$"), ("X_REF", "X"), ("X_coh_REF", "X$_\\mathrm{coh}$"),
        ("X_mag_REF", "X$_\\mathrm{mag}$")]
ax.axhspan(0.75, 1.25, color="#f3f2ee", lw=0, zorder=0)
ax.axhline(1.0, color=MODEL["full"], lw=0.8, zorder=1)
zero_line(ax)
for j, (m, lab, c, mk) in enumerate(MODS):
    x = np.arange(len(keys)) + (j - 1.5) * 0.14
    y = [V["contrasts"][k][m] / V["contrasts"][k]["full"] for k, _ in keys]
    ax.plot(x, y, marker=mk, color=c, ls="none", ms=4.2 if mk == "*" else 3.6, mec="white", mew=0.4, zorder=3)
    rec("c", f"model/full {m}", [k for k, _ in keys], y, source=SRC)
ax.set_xticks(range(len(keys)))
ax.set_xticklabels([l_ for _, l_ in keys], fontsize=fs(6.3))
ax.set_ylim(-1.5, 3.0)
ax.set_ylabel("model contrast / full-model contrast")

# ------------------------------------------------------------------ (d) held-out drive levels
ax = axs[3]
panel(ax, "d", "held-out levels")
zero_line(ax)
levs = sorted(int(k) for k in V["levels"])
for m, lab, c, mk in [("full", MODEL_LABEL["full"], MODEL["full"], MODEL_MARKER["full"])] + MODS:
    y = [V["levels"][str(l_)][m]["X_coh"] for l_ in levs]
    ax.plot(levs, y, marker=mk, color=c, ms=3.6 if mk != "*" else 4.4, mec="white", mew=0.4, lw=1.0, label=lab if m == "full" else None)
    rec("d", f"X_coh {m}", levs, y, source=SRC)
ax.set_xlim(-7, 9)
ax.set_xticks([-6, -3, 0, 3, 6])
ax.set_ylim(-0.12, 0.65)
ax.set_xlabel("$P_\\mathrm{avg}/P_\\mathrm{1dB}$ (dB)")
ax.set_ylabel("$X_\\mathrm{coh}$, reference phase case")
ax.legend(loc="upper left", fontsize=fs(5.8), borderaxespad=0.2, handlelength=1.2)
save(fig, "fig10")
