"""P2 analysis: LO / dressed-state resonance hypothesis (rules: results/tqe_viability/00_decision_rules.json, P2)."""
from __future__ import annotations

import json

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from tv_common import TV, P_GRID_NS, SEEDS3, lo_rabi_hz
import rev2_p2 as P2

LOS = (0.35, 0.425, 0.5)
COLORS = {0.35: "#4a3aa7", 0.425: "#008300", 0.5: "#2a78d6"}


def selectivity(f, g, delta_max=1.31e6, step=0.01e6):
    lf = np.log(g)
    out = np.full(len(f), np.nan)
    for i, fi in enumerate(f):
        d = np.arange(-delta_max, delta_max + step / 2, step)
        fd = fi + d
        m = (fd >= f.min()) & (fd <= f.max())
        out[i] = np.max(np.abs(np.interp(fd[m], f, lf) - lf[i]))
    return out


def collapse_rms(curves):
    """curves: {lo: (abscissa, values)}; RMS over a common abscissa grid (overlap) of the across-LO standard deviation."""
    lo_ = max(a.min() for a, _ in curves.values())
    hi_ = min(a.max() for a, _ in curves.values())
    if hi_ <= lo_:
        return np.nan
    grid = np.linspace(lo_, hi_, 200)
    stack = np.array([np.interp(grid, a, v) for a, v in curves.values()])
    return float(np.sqrt(np.mean(np.var(stack, axis=0))))


def main():
    rows = []
    for lo in LOS:
        f_lo = lo_rabi_hz(lo)
        recs = []
        for P in P_GRID_NS:
            f = 1e9 / P
            d = TV / "p2" / f"LO{lo:g}"
            c = P2.metrics(d / f"P{P}_CW.npz")
            q = P2.metrics(d / f"P{P}_REF_s{SEEDS3[0]}.npz")
            recs.append({"A_LO_Vpm": lo, "LO_Rabi_MHz": f_lo / 1e6, "IF_MHz": f / 1e6, "r1": f / f_lo, "r2": f / (f_lo / 2),
                         "g_CW": c["g_med"], **{f"X_{k}": 1 - q[f"g_{k}"] / c[f"g_{k}"] for k in ("med", "coh", "mag", "rms")},
                         "phase_circ_sd_rad": q["phase_circ_sd_rad"]})
        recs.sort(key=lambda r: r["IF_MHz"])
        f = np.array([r["IF_MHz"] for r in recs]) * 1e6
        S = selectivity(f, np.array([r["g_CW"] for r in recs]))
        for r, sv in zip(recs, S):
            r["selectivity_S131"] = sv
        rows += recs
    df = pd.DataFrame(rows).rename(columns={"X_med": "X"})
    df.to_csv(TV / "02_resonance_scan.csv", index=False)
    # features and collapse
    feats = {}
    for lo, g in df.groupby("A_LO_Vpm"):
        feats[lo] = {"g_min": g.loc[g.g_CW.idxmin()], "X_max": g.loc[g.X.idxmax()], "S_max": g.loc[g.selectivity_S131.idxmax()]}
    spread = {}
    for ft in ("g_min", "X_max", "S_max"):
        for ax in ("IF_MHz", "r1", "r2"):
            v = [feats[lo][ft][ax] for lo in LOS]
            spread[(ft, ax)] = float(max(v) / min(v))
    rms = {}
    for q, col in (("ln_g", "g_CW"), ("X", "X")):
        for ax in ("IF_MHz", "r1", "r2"):
            cur = {}
            for lo, g in df.groupby("A_LO_Vpm"):
                g = g.sort_values(ax)
                vals = np.log(g[col].values) if q == "ln_g" else g[col].values
                cur[lo] = (g[ax].values, vals)
            rms[(q, ax)] = collapse_rms(cur)
    verdict = {}
    for ax in ("r1", "r2"):
        ok = [ft for ft in ("g_min", "X_max", "S_max") if spread[(ft, ax)] <= 1.08 and spread[(ft, "IF_MHz")] >= 1.25]
        red = {q: rms[(q, ax)] / rms[(q, "IF_MHz")] for q in ("ln_g", "X")}
        verdict[ax] = {"features_aligned": ok, "rms_ratio": red}
    clear = any(len(v["features_aligned"]) >= 2 and all(r <= .5 for r in v["rms_ratio"].values()) for v in verdict.values())
    sugg = any(len(v["features_aligned"]) >= 1 or any(r <= .75 for r in v["rms_ratio"].values()) for v in verdict.values())
    cls = "CLEAR_RESONANCE_SCALING" if clear else ("SUGGESTIVE_ONLY" if sugg else "NO_SIMPLE_SCALING")
    # transient oscillation frequencies (archived) vs candidate scales
    osc = {0.5: {"power_step_MHz": [1 / 1.62, 1 / 1.60], "phase_step_MHz": [1 / 1.74, 1 / 1.63]}, 0.35: {"power_step_MHz": [1 / 0.560], "phase_step_MHz": [1 / 0.558]}}
    osc_rows = []
    for lo, o in osc.items():
        fl = lo_rabi_hz(lo) / 1e6
        osc_rows.append({"A_LO_Vpm": lo, "observed_MHz": o, "Omega_LO_2pi_MHz": fl, "Omega_LO_4pi_MHz": fl / 2,
                         "IF_minus_OmegaLO_4pi_MHz": abs(5 - fl / 2), "2IF_minus_OmegaLO_2pi_MHz": abs(10 - fl)})
    out = {"classification": cls, "feature_locations": {str(lo): {ft: {ax: float(feats[lo][ft][ax]) for ax in ("IF_MHz", "r1", "r2")} for ft in feats[lo]} for lo in LOS},
           "spread": {f"{a}|{b}": v for (a, b), v in spread.items()}, "collapse_rms": {f"{a}|{b}": v for (a, b), v in rms.items()},
           "per_normalization": verdict, "transient_oscillation": osc_rows}
    (TV / "02_resonance_verdict.json").write_text(json.dumps(out, indent=1, default=float))
    # figure
    fig, axs = plt.subplots(3, 3, figsize=(11, 8.5), sharex="col")
    for j, (ax_name, lab) in enumerate((("IF_MHz", "IF (MHz)"), ("r1", r"$r_1 = f_\mathrm{IF}/(\Omega_\mathrm{LO}/2\pi)$"), ("r2", r"$r_2 = f_\mathrm{IF}/(\Omega_\mathrm{LO}/4\pi)$"))):
        for i, (col, ylab) in enumerate((("g_CW", "CW gain / small-signal"), ("X", "X, reference sequence"), ("selectivity_S131", "selectivity S (±1.31 MHz)"))):
            ax = axs[i, j]
            for lo, g in df.groupby("A_LO_Vpm"):
                g = g.sort_values(ax_name)
                ax.plot(g[ax_name], g[col], "o-", ms=3, color=COLORS[lo], label=f"$A_\\mathrm{{LO}}$ = {lo:g} V/m")
            if i == 1:
                ax.axhline(0, color="0.5", lw=0.6)
            if j > 0:
                ax.axvline(1.0, color="0.6", lw=0.7, ls=":")
            ax.set_ylabel(ylab if j == 0 else "")
            if i == 2:
                ax.set_xlabel(lab)
    axs[0, 0].legend(fontsize=8)
    fig.suptitle(f"Resonance scan at signal/LO = 0.358 (standard probe): {cls}", fontsize=10)
    fig.tight_layout()
    fig.savefig(TV / "fig_resonance_scaling.png", dpi=150)
    pd.set_option("display.width", 250)
    print(df.round(3).to_string())
    print(json.dumps({k: v for k, v in out.items() if k not in ("transient_oscillation",)}, indent=1, default=float))
    print(json.dumps(osc_rows, indent=1, default=float))


if __name__ == "__main__":
    main()
