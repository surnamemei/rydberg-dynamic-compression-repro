"""Figure 1: receiver model and CW-matched surrogate controls (OP1).

(a) schematic; (b) CW fundamental gain vs amplitude, Nd = 1501 / 4001 / 8001 (controls_Nd*.npz);
(c) |c0|, |c1|, |c2| vs amplitude, Nd = 4001, 0-3.2 E1dB (controls_MH_Nd4001.npz); (d) surrogate structure.
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

from figstyle import (fs, COL2, INK, INK2, MUTED, MODEL, ND, E1, TAU, RES06, BOLD, panel, rec, save, relpath, zero_line)

A = RES06 / "artifacts"
CW = {n: np.load(A / f"controls_Nd{n}.npz") for n in (1501, 4001, 8001)}
MH = np.load(A / "controls_MH_Nd4001.npz")
FS_TXT = fs(6.8)


def box(ax, x, y, w, h, text, fc="white", ec=INK2, fs=FS_TXT, lw=0.6, style="round,pad=0.004,rounding_size=0.012", **kw):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle=style, fc=fc, ec=ec, lw=lw, transform=ax.transData, clip_on=False))
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=fs, color=INK, **kw)


def arrow(ax, p0, p1, style="-|>", lw=0.7, color=INK2, ms=6, **kw):
    ax.add_patch(FancyArrowPatch(p0, p1, arrowstyle=style, mutation_scale=ms, lw=lw, color=color, shrinkA=0, shrinkB=0, clip_on=False, **kw))


fig = plt.figure(figsize=(COL2, 5.85))
gs = fig.add_gridspec(3, 2, height_ratios=[1.75, 2.0, 1.55])
axa = fig.add_subplot(gs[0, :])
axb = fig.add_subplot(gs[1, 0])
axc = fig.add_subplot(gs[1, 1])
axd = fig.add_subplot(gs[2, :])

# ------------------------------------------------------------------ (a) schematic
ax = axa
ax.set_xlim(0, 1)
ax.set_ylim(0, 1)
ax.axis("off")
panel(ax, "a")
lv = {1: 0.47, 2: 0.66, 3: 0.86, 4: 0.97}
names = {1: "|1⟩ ground", 2: "|2⟩ excited", 3: "|3⟩ Rydberg", 4: "|4⟩ Rydberg"}
for k, y in lv.items():
    ax.plot([0.088, 0.136], [y, y], color=INK, lw=1.0, solid_capstyle="butt")
    ax.text(0.083, y, names[k], ha="right", va="center", fontsize=FS_TXT)
arrow(ax, (0.100, lv[1]), (0.100, lv[2]), lw=0.9, color=INK)
arrow(ax, (0.112, lv[2]), (0.112, lv[3]), lw=0.9, color=INK)
arrow(ax, (0.125, lv[3] + 0.006), (0.125, lv[4] - 0.006), style="<|-|>", lw=0.9, color=INK, ms=4.5)
ax.text(0.145, (lv[1] + lv[2]) / 2, "probe, 852 nm, " + r"$\Omega_p/2\pi$ = 8.08 MHz", ha="left", va="center", fontsize=FS_TXT)
ax.text(0.145, (lv[2] + lv[3]) / 2, "coupling, 510 nm, " + r"$\Omega_c/2\pi$ = 2.05 MHz", ha="left", va="center", fontsize=FS_TXT)
ax.text(0.145, (lv[3] + lv[4]) / 2 + 0.005, "RF (LO + signal): " + r"$A_\mathrm{LO} + E(t)\,e^{j2\pi f_\mathrm{IF}t}$, " + r"$f_\mathrm{IF}$ = 5 MHz" + "\n"
        + r"$A_\mathrm{LO}$ = 0.5 V/m (OP1) or 0.35 V/m (OP2)", ha="left", va="center", fontsize=FS_TXT)

# thermal velocity classes
x0, x1, yb, yt = 0.47, 0.63, 0.50, 0.80
v = np.linspace(-3.5, 3.5, 29)
w = np.exp(-v ** 2 / 2)
for vi, wi in zip(v, w):
    xx = x0 + (vi + 3.5) / 7 * (x1 - x0)
    kept = abs(vi) <= 3.0 + 1e-9
    ax.plot([xx, xx], [yb, yb + wi * (yt - yb)], color=MODEL["full"] if kept else MUTED, lw=1.0 if kept else 0.6, solid_capstyle="butt")
ax.plot([x0 - 0.005, x1 + 0.005], [yb, yb], color=INK2, lw=0.6)
ax.text((x0 + x1) / 2, yb - 0.03, ND + r" velocity classes $v$ (±3σ kept of a ±3.5σ grid)", ha="center", va="top", fontsize=FS_TXT)
ax.text((x0 + x1) / 2, yt + 0.04, "one density matrix per class (RK4)", ha="center", va="bottom", fontsize=FS_TXT)
arrow(ax, (0.425, 0.62), (0.462, 0.62))

# thermal sum
bx, by, bw_, bh_ = 0.70, 0.54, 0.17, 0.28
arrow(ax, (0.645, 0.66), (bx - 0.004, 0.66))
box(ax, bx, by, bw_, bh_, "thermal sum of probe coherences\n→ probe transmission\n(1 GHz intensity trace)", fc="#f4f3f0")

# receiver chain (second row)
yc, ch = 0.15, 0.22
chain_boxes = [(0.155, 0.12, "DC removal"), (0.300, 0.15, "mix with " + r"$e^{-j2\pi f_\mathrm{IF}t}$"),
               (0.475, 0.19, "low-pass, 3 MHz\n(8th-order Butterworth)"), (0.690, 0.14, "decimate to 20 MHz")]
for i, (xx, ww, lab) in enumerate(chain_boxes):
    box(ax, xx, yc - ch / 2, ww, ch, lab)
    if i < len(chain_boxes) - 1:
        arrow(ax, (xx + ww + 0.002, yc), (chain_boxes[i + 1][0] - 0.004, yc))
arrow(ax, (0.832, yc), (0.862, yc))
ax.text(0.866, yc, "baseband output", ha="left", va="center", fontsize=FS_TXT)
ax.text(0.0, yc, "receiver chain\n(common to the\nfull model and\nall surrogates)", ha="left", va="center", fontsize=FS_TXT, color=INK2)
# elbow connector: thermal-sum box -> first chain box
ax.plot([bx + bw_ / 2, bx + bw_ / 2, 0.215], [by, 0.40, 0.40], color=INK2, lw=0.7, solid_capstyle="butt", solid_joinstyle="miter")
arrow(ax, (0.215, 0.40), (0.215, yc + ch / 2 + 0.004))

# ------------------------------------------------------------------ (b) CW fundamental gain vs amplitude
ax = axb
panel(ax, "b", "CW fundamental gain (OP1)")
sty = {4001: dict(color=MODEL["full"], lw=1.2, marker="o", ms=3.2, mec="white", mew=0.4, zorder=3),
       8001: dict(color="#104281", lw=0.9, ls=(0, (1.2, 1.2)), zorder=4),
       1501: dict(color=MUTED, lw=1.0, marker="o", ms=3.0, mfc="white", mec=MUTED, mew=0.7, zorder=2)}
for n in (4001, 8001, 1501):
    d = CW[n]
    a, g, e = d["amps"][1:], d["gain_db"][1:], float(d["e1db"])
    sel = a <= 0.16
    ax.plot(a[sel], g[sel], label=f"{ND} = {n}  ({E1} = {e:.4f} V/m)", **sty[n])
    ax.plot([e, e], [-5.2, -1], color=sty[n]["color"], lw=0.6, ls=(0, (2, 1.5)), zorder=1)
    rec("b", f"gain_dB Nd={n}", a[sel], g[sel], source=relpath(A / f"controls_Nd{n}.npz"))
    rec("b", f"E1dB Nd={n}", "E1dB", e, source=relpath(A / f"controls_Nd{n}.npz"), note="stored e1db")
ax.axhline(-1, color=INK2, lw=0.6, zorder=0)
ax.text(0.003, -0.92, "−1 dB", fontsize=FS_TXT, color=INK2, va="bottom")
ax.set_xlim(0, 0.16)
ax.set_ylim(-5.2, 0.4)
ax.set_xlabel("CW field amplitude $a$ (V/m)")
ax.set_ylabel("gain re weak-tone slope (dB)")
ax.legend(loc="lower left", fontsize=fs(6.3), frameon=True, facecolor="white", edgecolor="none", framealpha=1)

# ------------------------------------------------------------------ (c) harmonic zones
ax = axc
panel(ax, "c", f"CW harmonic amplitudes (OP1, {ND} = 4001)")
a = MH["amps"]
C = np.abs(MH["C"])
e1 = float(CW[4001]["e1db"])
sel = a <= 3.2 * e1 + 1e-9
zc = {0: ("#104281", "s", "|$c_0$| rectified baseline"), 1: (MODEL["full"], "o", "|$c_1$| fundamental"), 2: ("#86b6ef", "^", "|$c_2$| second harmonic")}
for k in (1, 0, 2):
    col, mk, lab = zc[k]
    ax.plot(a[sel], C[sel, k] * 1e3, color=col, marker=mk, ms=3.0, mec="white", mew=0.4, lw=1.0, label=lab)
    rec("c", f"|c{k}|x1e3", a[sel], C[sel, k] * 1e3, source=relpath(A / "controls_MH_Nd4001.npz"))
ax.axvline(e1, color=INK2, lw=0.6, ls=(0, (2, 1.5)), zorder=0)
ax.text(e1 + 0.003, 0.50, E1, fontsize=FS_TXT, color=INK2, va="top")
i_near = int(np.argmin(np.abs(a - e1)))
rng = (a >= 0.75 * e1) & (a <= 3.2 * e1)
ratio = C[rng, 0] / C[rng, 1]
i_max = np.flatnonzero(rng)[int(np.argmax(ratio))]
ax.annotate(f"|$c_0/c_1$| = {C[i_near, 0] / C[i_near, 1]:.2f}", xy=(a[i_near], C[i_near, 0] * 1e3), xytext=(0.086, 0.175),
            fontsize=FS_TXT, arrowprops=dict(arrowstyle="-", lw=0.5, color=INK2, shrinkA=0, shrinkB=2))
ax.annotate(f"|$c_0/c_1$| = {C[i_max, 0] / C[i_max, 1]:.1f}", xy=(a[i_max], C[i_max, 0] * 1e3), xytext=(0.150, 0.70),
            fontsize=FS_TXT, arrowprops=dict(arrowstyle="-", lw=0.5, color=INK2, shrinkA=0, shrinkB=2))
rec("c", "|c0/c1| nearest E1dB", a[i_near], C[i_near, 0] / C[i_near, 1], source=relpath(A / "controls_MH_Nd4001.npz"))
rec("c", "max |c0/c1| 0.75-3.2 E1dB", a[i_max], C[i_max, 0] / C[i_max, 1], source=relpath(A / "controls_MH_Nd4001.npz"))
ax.set_xlim(0, 3.2 * e1 + 0.004)
ax.set_ylim(0, 0.80)
ax.set_xlabel("CW field amplitude $a$ (V/m)")
ax.set_ylabel(r"|$c_k$| (probe transmission, $10^{-3}$)")
ax.legend(loc="upper left", fontsize=fs(6.3))

# ------------------------------------------------------------------ (d) surrogate structure
ax = axd
ax.set_xlim(0, 1)
ax.set_ylim(0, 1)
ax.axis("off")
panel(ax, "d")
col_x = {"lab": 0.0, 0: 0.175, 1: 0.385, 2: 0.865}
col_w = {0: 0.19, 1: 0.46, 2: 0.135}
hdr_y = 0.90
ax.text(col_x[0] + col_w[0] / 2, hdr_y, "zone $k$ = 0 (baseline)", ha="center", va="center", fontsize=FS_TXT, fontweight=BOLD)
ax.text(col_x[1] + col_w[1] / 2, hdr_y, "zone $k$ = 1 (fundamental)", ha="center", va="center", fontsize=FS_TXT, fontweight=BOLD)
ax.text(col_x[2] + col_w[2] / 2, hdr_y, "zones $k$ ≥ 2", ha="center", va="center", fontsize=FS_TXT, fontweight=BOLD)
rows = [("static", "all-zone static", 0.70), ("lti", "all-zone LTI+static\n(Hammerstein)", 0.46), ("dsh", "single-slow-state", 0.22)]
MEM, LTI_FC, SLOW_FC = "white", "#e8f1fb", "#fdf1d6"
bh = 0.17


def chain(y, x, items, total_w):
    """items: list of (text, facecolor, weight); drawn left to right inside total_w with arrows."""
    n = len(items)
    gap = 0.018
    ws = [wt for _, _, wt in items]
    unit = (total_w - gap * (n - 1)) / sum(ws)
    xx = x
    for i, (t, fc, wt) in enumerate(items):
        bw_ = wt * unit
        box(ax, xx, y - bh / 2, bw_, bh, t, fc=fc)
        if i < n - 1:
            arrow(ax, (xx + bw_ + 0.001, y), (xx + bw_ + gap - 0.001, y), ms=4.5, lw=0.6)
        xx += bw_ + gap


for key, lab, y in rows:
    ax.add_patch(FancyBboxPatch((col_x["lab"], y - 0.035), 0.012, 0.07, boxstyle="square,pad=0", fc=MODEL[key], ec="none"))
    ax.text(col_x["lab"] + 0.02, y, lab, ha="left", va="center", fontsize=FS_TXT)
    if key == "static":
        chain(y, col_x[0], [(r"$c_0(|E|)$", MEM, 1)], col_w[0])
        chain(y, col_x[1], [(r"$c_1(|E|)\,u$  (CW table)", MEM, 1)], col_w[1])
    else:
        chain(y, col_x[0], [(r"$c_0(|E|)$", MEM, 1), ("baseline\nkernel", LTI_FC, 1)], col_w[0])
        if key == "lti":
            chain(y, col_x[1], [(r"CW gain at $|E|$", MEM, 1.2), ("small-signal\nfield kernel", LTI_FC, 1)], col_w[1])
        else:
            chain(y, col_x[1], [(r"$|E|^2$ → low-pass," + "\n" + r"$\tau$ = " + TAU + r" → $P_s$", SLOW_FC, 1.55), ("CW gain at √" + r"$P_s$", MEM, 1.0),
                                ("small-signal\nfield kernel", LTI_FC, 1)], col_w[1])
    chain(y, col_x[2], [(r"$c_k(|E|)\,u^k$", MEM, 1)], col_w[2])
ax.text(0.0, -0.02, "White: memoryless CW table.  Blue: small-signal LTI kernel.  Yellow: slow power state.\n"
        "Every surrogate reproduces the full model's CW tones in every zone; its output passes through the receiver chain in (a).",
        ha="left", va="bottom", fontsize=fs(6.3), color=INK2)

save(fig, "fig1")
