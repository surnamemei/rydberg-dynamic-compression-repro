"""Shared style, palette, paths and export helpers for the final manuscript figures.

Every figure script plots stored results only (CSV/JSON/NPZ under results/ and manuscript/posthoc/).
No simulation code is imported. Each script records the plotted values with `rec()`; `save()` writes
PDF (TrueType-embedded text), SVG (text as paths) and a 600-dpi PNG preview, plus
figures/audit/<name>_plotted_values.csv for the numerical cross-check in audit_figures.py.
"""
from __future__ import annotations

import csv
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib import font_manager as fm  # noqa: E402

ROOT = Path(__file__).resolve().parents[3]
MS = ROOT / "manuscript"
RES06 = ROOT / "results" / "stage06_dwell_physics"
PM = RES06 / "phase_mechanism"
FV = ROOT / "results" / "final_validation"
ST05 = ROOT / "results" / "p1db_waveform_stage05"
POSTHOC = MS / "posthoc"
OUT_MAIN = MS / "figures" / "final"
OUT_SUPP = MS / "figures" / "supplement"
AUDIT = MS / "figures" / "audit"
for d in (OUT_MAIN, OUT_SUPP, AUDIT):
    d.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------- page geometry (inches)
COL1 = 242.67355 / 72.27   # IEEE TQE template (ieeeaccess.cls) \columnwidth = 3.358 in (85.3 mm)
COL2 = 505.12177 / 72.27   # IEEE TQE template \textwidth = 6.989 in (177.5 mm)


def fs(x):
    """Map the original (pre-TQE) annotation sizes to the TQE sizes: nothing below 6.5 pt; labels ~8 pt."""
    x = float(x)
    if x < 6.0:
        return 6.5
    if x < 6.5:
        return 6.8
    if x < 7.0:
        return 7.0
    if x < 7.5:
        return 7.5
    return 8.0

# ---------------------------------------------------------------- font: FreeSans (Helvetica metric clone), TrueType files only
fm.fontManager.ttflist = [f for f in fm.fontManager.ttflist if not (f.name == "FreeSans" and f.fname.endswith(".otf"))]
BOLD = 600  # FreeSans Bold is registered at weight 600

# ---------------------------------------------------------------- ink and palette (dataviz reference palette, light mode, white print surface)
INK = "#0b0b0b"
INK2 = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
REF_GRAY = "#898781"
# Models: categorical slots 1-4 (validated adjacent pairs on #ffffff; single-slow-state never shares a panel with all-zone static)
MODEL = {"full": "#2a78d6", "static": "#eb6834", "lti": "#1baf7a", "dsh": "#eda100"}
MODEL_LABEL = {"full": "full model", "static": "all-zone static", "lti": "all-zone LTI+static", "dsh": "single-slow-state"}
MODEL_MARKER = {"full": "o", "static": "s", "lti": "^", "dsh": "D"}
# Formats: contiguous slots 5-8 (validated adjacent; secondary encoding by marker and direct label)
FMT = {"CE": "#e87ba4", "QPSK": "#008300", "16-QAM": "#4a3aa7", "OFDM": "#e34948"}
FMT_MARKER = {"CE": "o", "QPSK": "s", "16-QAM": "^", "OFDM": "D"}
FMT_LABEL = {"CE": "CE (GMSK)", "QPSK": "QPSK", "16-QAM": "16-QAM", "OFDM": "OFDM"}

plt.rcParams.update({
    "font.family": "FreeSans",
    "mathtext.fontset": "custom",
    "mathtext.rm": "FreeSans",
    "mathtext.it": "FreeSans:italic",
    "mathtext.bf": "FreeSans:semibold",
    "mathtext.sf": "FreeSans",
    "font.size": 8,
    "axes.titlesize": 8,
    "axes.labelsize": 8,
    "xtick.labelsize": 7.5,
    "ytick.labelsize": 7.5,
    "legend.fontsize": 7,
    "legend.frameon": False,
    "legend.handlelength": 1.8,
    "legend.borderaxespad": 0.3,
    "legend.labelspacing": 0.3,
    "text.color": INK,
    "axes.labelcolor": INK,
    "axes.edgecolor": INK2,
    "axes.linewidth": 0.6,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.titlepad": 3.0,
    "axes.grid": False,
    "grid.color": GRID,
    "grid.linewidth": 0.4,
    "grid.linestyle": "-",
    "xtick.color": INK2,
    "ytick.color": INK2,
    "xtick.labelcolor": INK,
    "ytick.labelcolor": INK,
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "xtick.minor.width": 0.4,
    "ytick.minor.width": 0.4,
    "xtick.major.size": 2.5,
    "ytick.major.size": 2.5,
    "xtick.minor.size": 1.5,
    "ytick.minor.size": 1.5,
    "lines.linewidth": 1.1,
    "lines.markersize": 3.6,
    "lines.markeredgewidth": 0.5,
    "lines.solid_capstyle": "round",
    "lines.solid_joinstyle": "round",
    "errorbar.capsize": 1.6,
    "patch.linewidth": 0.5,
    "hatch.linewidth": 0.5,
    "figure.dpi": 150,
    "savefig.dpi": 600,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "svg.fonttype": "path",
    "axes.unicode_minus": True,
    "figure.constrained_layout.use": True,
    "figure.constrained_layout.h_pad": 0.02,
    "figure.constrained_layout.w_pad": 0.02,
    "figure.constrained_layout.hspace": 0.04,
    "figure.constrained_layout.wspace": 0.04,
})

EB = dict(elinewidth=0.8, capsize=1.6, capthick=0.8)  # CI error-bar style

# ---------------------------------------------------------------- notation (identical in every figure)
ND = r"$N_\mathrm{d}$"
PP = r"$P_\mathrm{avg}/P_\mathrm{1dB}$"
E1 = r"$E_\mathrm{1dB}$"
TAU = r"$\tau_\mathrm{ref}$"
XI = r"$T_\mathrm{dwell}/\tau_\mathrm{ref}$"
ALO = r"$A_\mathrm{LO}$"
OP1 = r"OP1 ($A_\mathrm{LO}$ = 0.5 V/m)"
OP2 = r"OP2 ($A_\mathrm{LO}$ = 0.35 V/m)"


def panel(ax, letter, title=None):
    """Panel label '(a)' in bold at top-left; optional short descriptor in regular weight."""
    ax.set_title(f"({letter})", loc="left", fontweight=BOLD, fontsize=9.5, pad=3)
    if title:
        ax.set_title(title, loc="center", fontsize=8, pad=3.5, color=INK)


def mline(ax, x, y, key, yerr=None, label=True, ms=None, **kw):
    """Model series: colour + marker by model; optional CI error bars (yerr = [lo_err, hi_err])."""
    lab = MODEL_LABEL[key] if label is True else (label or None)
    kw.setdefault("lw", 1.1)
    if yerr is not None:
        return ax.errorbar(x, y, yerr=yerr, color=MODEL[key], marker=MODEL_MARKER[key], ms=ms or 3.6,
                           mec="white", mew=0.5, label=lab, **EB, **kw)
    return ax.plot(x, y, color=MODEL[key], marker=MODEL_MARKER[key], ms=ms or 3.6, mec="white", mew=0.5, label=lab, **kw)


def zero_line(ax, y=0.0, **kw):
    kw.setdefault("color", INK2)
    kw.setdefault("lw", 0.6)
    kw.setdefault("zorder", 0)
    ax.axhline(y, **kw)


# ---------------------------------------------------------------- plotted-value recorder
_REC: list[dict] = []


def rec(panel_id, series, x, y, lo=None, hi=None, source="", note=""):
    """Record plotted values (scalars or equal-length sequences) for the audit cross-check."""
    import numpy as np
    xs = np.atleast_1d(np.asarray(x, dtype=object))
    ys = np.atleast_1d(np.asarray(y, dtype=float))
    los = np.atleast_1d(np.asarray(lo, dtype=float)) if lo is not None else [None] * len(ys)
    his = np.atleast_1d(np.asarray(hi, dtype=float)) if hi is not None else [None] * len(ys)
    if len(xs) == 1 and len(ys) > 1:
        xs = [xs[0]] * len(ys)
    for xi, yi, l, h in zip(xs, ys, los, his):
        _REC.append({"panel": panel_id, "series": series, "x": xi, "y": float(yi),
                     "ci_lo": "" if l is None else float(l), "ci_hi": "" if h is None else float(h),
                     "source": source, "note": note})


def clip_report(fig):
    """Count data points (lines, error bars, bars, filled regions, step patches) that fall outside their axes' limits."""
    import numpy as np
    from matplotlib.collections import LineCollection, PolyCollection
    from matplotlib.patches import Rectangle, StepPatch
    out = []
    for ax in fig.get_axes():
        (x0, x1), (y0, y1) = sorted(ax.get_xlim()), sorted(ax.get_ylim())
        tx, ty = 1e-9 * max(abs(x0), abs(x1), 1e-12), 1e-9 * max(abs(y0), abs(y1), 1e-12)

        def n_out(xy):
            xy = np.asarray(xy, float).reshape(-1, 2)
            xy = xy[np.all(np.isfinite(xy), axis=1)]
            bad = (xy[:, 0] < x0 - tx) | (xy[:, 0] > x1 + tx) | (xy[:, 1] < y0 - ty) | (xy[:, 1] > y1 + ty)
            return int(bad.sum())
        n = 0
        for ln in ax.get_lines():
            if ln.get_transform() == ax.transData and ln.get_gid() != "decor":
                n += n_out(ln.get_xydata())
        for c in ax.collections:
            if c.get_transform() == ax.transData and isinstance(c, (LineCollection, PolyCollection)) and c.get_gid() != "decor":
                for p in c.get_paths():
                    n += n_out(p.vertices)
        for p in ax.patches:
            if p.get_transform() == ax.transData and isinstance(p, Rectangle) and p.get_height() != 0:
                n += n_out([[p.get_x(), p.get_y()], [p.get_x(), p.get_y() + p.get_height()]])
            elif isinstance(p, StepPatch):
                vals = np.asarray(p.get_data().values, float)
                vals = vals[np.isfinite(vals) & (vals > 0)] if ax.get_yscale() == "log" else vals[np.isfinite(vals)]
                n += int(((vals < y0 - ty) | (vals > y1 + ty)).sum())
        if n:
            out.append(f"{ax.get_title(loc='left') or ax.get_title() or 'axes'}: {n} points outside limits")
    return out


def save(fig, name, supplement=False):
    """Write <name>.pdf/.svg/.png into figures/final or figures/supplement and the plotted-value CSV."""
    fig.canvas.draw()
    clips = clip_report(fig)
    _REC.append({"panel": "_audit", "series": "points outside axis limits", "x": "count", "y": float(sum(int(c.split(": ")[-1].split()[0]) for c in clips)),
                 "ci_lo": "", "ci_hi": "", "source": "", "note": "; ".join(clips)})
    out = OUT_SUPP if supplement else OUT_MAIN
    fig.savefig(out / f"{name}.pdf", metadata={"Creator": None, "Producer": None, "CreationDate": None})
    fig.savefig(out / f"{name}.svg", metadata={"Creator": None, "Date": None, "Format": None, "Type": None})
    fig.savefig(out / f"{name}.png", dpi=600, metadata={"Software": None})
    with open(AUDIT / f"{name}_plotted_values.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["panel", "series", "x", "y", "ci_lo", "ci_hi", "source", "note"])
        w.writeheader()
        w.writerows(_REC)
    _REC.clear()
    plt.close(fig)
    print(f"wrote {out.relative_to(ROOT)}/{name}.{{pdf,svg,png}}", file=sys.stderr)


def relpath(p: Path) -> str:
    return str(Path(p).resolve().relative_to(ROOT))


def tm(s) -> str:
    """Typographic minus for formatted numbers in annotation text (tick labels already use U+2212)."""
    return str(s).replace("-", "−")
