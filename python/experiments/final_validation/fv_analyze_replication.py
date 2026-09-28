"""Analysis of the second-operating-point replication (preregistered rules, evaluated mechanically)."""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from fv_common import FV, OUT06, ROOT

sys.path.insert(0, str(ROOT / "manuscript" / "posthoc"))
from posthoc_timescales import dominant_period  # noqa: E402

REP = FV / "replication"
PM_RUNS = OUT06 / "phase_mechanism" / "runs"
PRE = json.loads((REP / "00_preregistration.json").read_text())
BLUE, ORANGE, AQUA, YELLOW, GREY = "#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#8a8a85"
plt.rcParams.update({"axes.grid": True, "grid.alpha": .25, "axes.spines.top": False, "axes.spines.right": False, "font.size": 9})


def metrics(path):
    """Identical windows/definitions to phase_mechanism/pm_analyze.metrics."""
    z = np.load(path)
    t = z["t"] * 1e6
    g = np.abs(z["z_full"]) / np.abs(z["z_linear"])
    pre, last = (t >= 40) & (t < 50), (t >= 89.5) & (t < 109.5)
    mod = (t >= 50) & (t < 109.5)
    g_pre, g_fin = float(np.median(g[pre])), float(np.median(g[last]))
    g2 = np.convolve(g, np.ones(100) / 100, "same")[mod]
    t63 = np.nan
    if abs(g_fin - g_pre) >= .02:
        i = np.flatnonzero((g2 - g_pre) / (g_fin - g_pre) >= .63)
        t63 = float(t[mod][i[0]] - 50) if len(i) else np.nan
    pops = (z["pops"] * z["pv"]).sum(axis=1)
    return {"g_pre": g_pre, "g_final": g_fin, "t63_us": t63, "probe_db_last20": float(np.mean(z["probe_db"][last])),
            "pop_rydberg": float(pops[2] + pops[3]), "abs_full_last20": float(np.median(np.abs(z["z_full"][last]))),
            "abs_lin_last20": float(np.median(np.abs(z["z_linear"][last])))}


def phase_step_periods(path, ts_us=4.0):
    z = np.load(path)
    t, phi = z["t"] * 1e6, z["phi"]
    u = z["z_full"] * np.conj(z["z_linear"])
    out = []
    for k in range(1, int(round(60 / ts_us))):
        tb = 50 + k * ts_us
        i0 = np.searchsorted(t, tb)
        if abs(phi[min(i0 + 50, len(phi) - 1)] - phi[max(i0 - 50, 0)]) < 1e-6:
            continue
        post, pre = (t >= tb + .3) & (t < tb + 3.7), (t >= tb - .6) & (t < tb - .2)
        if post.sum() < 150:
            continue
        out.append(dominant_period(np.angle(u[post] * np.conj(np.mean(u[pre]))), t[1] - t[0]))
    return out


def power_step_period(path):
    z = np.load(path)
    t = z["t"] * 1e6
    g = np.abs(z["z_full"]) / np.abs(z["z_linear"])
    fin = np.median(g[(t >= 90) & (t < 109.5)])
    m = (t >= 50.3) & (t <= 58)
    return dominant_period(g[m] - fin, t[1] - t[0])


def main():
    cfg = json.loads((REP / "02_run_configs.json").read_text())
    ctl2 = np.load(FV / "artifacts" / "controls_LO0.35_Nd4001.npz")
    a2 = cfg["a_prime_Vpm"]
    g_static2 = float(np.interp(a2, ctl2["amps"], np.abs(ctl2["harm"])) / a2 / ctl2["slope"])
    rows = []
    for f in sorted((REP / "runs").glob("*.npz")):
        name, nd = f.stem.rsplit("_Nd", 1)
        rows.append({"case": name, "Nd": int(nd), **metrics(f), **{k: v for k, v in cfg["cases"].get(name, {}).items() if isinstance(v, (int, float, str))}})
    pm_cw = metrics(PM_RUNS / "CW_Nd4001_dt1.npz")
    pm_ref = metrics(PM_RUNS / "QPSK_REF_Nd4001_dt1.npz")
    pm_off = metrics(PM_RUNS / "OFF_+0.125_Nd4001_dt1.npz")
    res = pd.DataFrame(rows)
    g_cw2 = float(res[res.case == "LO2_CW"].g_final.iloc[0])
    res["g_CW_same_point"] = np.where(res.case.str.startswith("LO2"), g_cw2, pm_cw["g_final"])
    res["X"] = 1 - res.g_final / res.g_CW_same_point
    res["g_static_prediction"] = np.where(res.case.str.startswith("LO2"), g_static2, np.nan)
    res.to_csv(REP / "03_replication_results.csv", index=False)
    X = res.set_index("case").X
    x_ref2 = float(np.mean([X["LO2_QPSK_REF_zd"], X["LO2_QPSK_REF_zd_s2"]]))
    x_pm_ref = 1 - pm_ref["g_final"] / pm_cw["g_final"]

    v = []
    add = lambda item, value, rule, outcome: v.append({"item": item, "value": value, "rule": rule, "outcome": outcome})  # noqa: E731
    h_gen = x_ref2 >= .20
    add("H_gen", f"X_REF' = {x_ref2:.3f} (seeds {X['LO2_QPSK_REF_zd']:.3f} / {X['LO2_QPSK_REF_zd_s2']:.3f}); original point X_REF = {x_pm_ref:.3f} (+pi), {X['LO1_QPSK_REF_zd']:.3f} (zero-drift)", PRE["hypotheses"]["H_gen"], "PASS" if h_gen else "FAIL")
    h_shape = abs(X["LO2_QPSK_JUMP_zd"] - X["LO2_QPSK_LONGRAMP_zd"]) <= .10
    add("H_shape", f"JUMP' {X['LO2_QPSK_JUMP_zd']:.3f}, REF' {X['LO2_QPSK_REF_zd']:.3f}, LONGRAMP' {X['LO2_QPSK_LONGRAMP_zd']:.3f}", PRE["hypotheses"]["H_shape"], "PASS" if h_shape else "FAIL")
    offs = ["LO2_OFF_+0.125", "LO2_OFF_+1.31", "LO2_OFF_-1.31"]
    h_det = all(X[c] < .5 * x_ref2 for c in offs)
    add("H_detune", ", ".join(f"{c} {X[c]:+.3f}" for c in offs) + f" vs 0.5 X_REF' = {.5 * x_ref2:.3f}", PRE["hypotheses"]["H_detune"], "PASS" if h_det else "FAIL")
    h_drift = X["LO1_QPSK_FAST_zd"] >= X["LO1_QPSK_REF_zd"] - .03
    add("H_drift (original point)", f"zero-drift SLOW/REF/FAST = {X['LO1_QPSK_SLOW_zd']:.3f} / {X['LO1_QPSK_REF_zd']:.3f} / {X['LO1_QPSK_FAST_zd']:.3f} (+pi convention: 0.159 / {x_pm_ref:.3f} / 0.438)", PRE["hypotheses"]["H_drift"], "PASS (turnover was drift)" if h_drift else "FAIL (turnover persists without drift)")
    add("rate triplet at LO' (descriptive)", f"SLOW'/REF'/FAST' = {X['LO2_QPSK_SLOW_zd']:.3f} / {X['LO2_QPSK_REF_zd']:.3f} / {X['LO2_QPSK_FAST_zd']:.3f}", "descriptive", "-")
    p_pow = power_step_period(REP / "runs" / "LO2_POWERSTEP_Nd4001.npz")
    p_ph = phase_step_periods(REP / "runs" / "LO2_QPSK_SLOW_zd_Nd4001.npz")
    c9 = abs(p_pow - np.median(p_ph)) <= .15
    add("C9a replication at LO'", f"power-step period {p_pow:.2f} us; phase-step periods {np.round(p_ph, 2).tolist()} (median {np.median(p_ph):.2f})", PRE["hypotheses"]["C9a_replication"], "PASS" if c9 else "FAIL")
    x8 = 1 - float(res[(res.case == "LO1_OFF_+0.125")].g_final.iloc[0]) / pm_cw["g_final"]
    x4 = 1 - pm_off["g_final"] / pm_cw["g_final"]
    lr = json.loads((REP / "long_dwell_xi32_Nd8001_rows.json").read_text())
    d8 = {r["model"]: r for r in lr}
    L = pd.read_csv(OUT06 / "rows_fu_long_cw_dwell.csv").drop_duplicates(["job_id", "model"], keep="last")
    d4 = L[(L.realization == 0) & (L.variant == "xi_32")].set_index("model")
    rel = abs(d8["atomic"]["D_ref_linear"] / d4.loc["atomic", "D_ref_linear"] - 1)
    num_ok = abs(x8 - x4) <= .02 and rel <= .02
    add("numerics (Nd=8001 closure)", f"OFF_+0.125: X {x4:+.4f} (4001) vs {x8:+.4f} (8001); xi=32 long dwell D_ref {d4.loc['atomic', 'D_ref_linear']:.5f} vs {d8['atomic']['D_ref_linear']:.5f} (rel {rel:.2%}); gain_abs rel change {abs(d8['atomic']['gain_abs'] / d4.loc['atomic', 'gain_abs'] - 1):.2%}", PRE["hypotheses"]["numerics"], "PASS" if num_ok else "FAIL")
    add("static prediction at LO'", f"g_static' = {g_static2:.3f}; full model CW g = {g_cw2:.3f}; E1dB' = {cfg['E1dB_prime_Vpm']:.4f} V/m, a' = {a2:.4f} V/m", "pipeline sanity (CW should match static)", "OK" if abs(g_cw2 - g_static2) < .02 else "CHECK")
    cls = "A_replicates" if (h_gen and h_shape and h_det) else ("B_weakens_but_qualitative" if x_ref2 >= .05 else "C_fails")
    add("CLASSIFICATION", f"X_REF' = {x_ref2:.3f}", json.dumps(PRE["classification"]), cls)
    vd = pd.DataFrame(v)
    vd.to_csv(REP / "04_replication_verdicts.csv", index=False)
    pd.set_option("display.width", 250)
    pd.set_option("display.max_colwidth", 200)
    print(vd[["item", "value", "outcome"]].to_string())
    print(res[["case", "Nd", "g_pre", "g_final", "X", "t63_us", "probe_db_last20", "pop_rydberg", "mean_signed_freq_offset_Hz"]].round(4).to_string())
    figures(res, pm_cw, x_pm_ref, g_static2, p_pow, p_ph)


def figures(res, pm_cw, x_pm_ref, g_static2, p_pow, p_ph):
    R = res.set_index("case")
    fig, ax = plt.subplots(1, 2, figsize=(12, 4), sharey=True)
    pm = pd.read_csv(OUT06 / "phase_mechanism" / "03_phase_mechanism_results.csv")
    pm = pm[(pm.Nd == 4001) & (pm.dt_ns == 1)].set_index("case")
    left = [("CW", "CW"), ("OFF_+0.125", "+0.125 MHz"), ("OFF_+1.31", "+1.31 MHz"), ("OFF_-1.31", "−1.31 MHz"), ("QPSK_REF", "REF (+π)"), ("QPSK_JUMP", "JUMP (+π)"),
            ("QPSK_LONGRAMP", "LONGRAMP (+π)"), ("QPSK_SLOW", "SLOW (+π)"), ("QPSK_FAST", "FAST (+π)")]
    zd1 = [("LO1_QPSK_SLOW_zd", "SLOW zd"), ("LO1_QPSK_REF_zd", "REF zd"), ("LO1_QPSK_FAST_zd", "FAST zd")]
    xs = [pm.loc[c, "X"] for c, _ in left] + [R.loc[c, "X"] for c, _ in zd1]
    labs = [l for _, l in left] + [l for _, l in zd1]
    cols = [GREY] + [AQUA] * 3 + [BLUE] * 5 + [YELLOW] * 3
    ax[0].bar(range(len(xs)), np.clip(xs, -1, None), color=cols)
    ax[0].set_xticks(range(len(xs)))
    ax[0].set_xticklabels(labs, rotation=40, ha="right", fontsize=7)
    ax[0].set_title("original point: LO 0.5 V/m, a_H = 0.179 V/m (offsets clipped at −1)", fontsize=9)
    right = [("LO2_CW", "CW"), ("LO2_OFF_+0.125", "+0.125 MHz"), ("LO2_OFF_+1.31", "+1.31 MHz"), ("LO2_OFF_-1.31", "−1.31 MHz"), ("LO2_QPSK_REF_zd", "REF zd"),
             ("LO2_QPSK_REF_zd_s2", "REF zd s2"), ("LO2_QPSK_JUMP_zd", "JUMP zd"), ("LO2_QPSK_LONGRAMP_zd", "LONGRAMP zd"), ("LO2_QPSK_SLOW_zd", "SLOW zd"), ("LO2_QPSK_FAST_zd", "FAST zd")]
    xs2 = [R.loc[c, "X"] for c, _ in right]
    ax[1].bar(range(len(xs2)), np.clip(xs2, -1, None), color=[GREY] + [AQUA] * 3 + [YELLOW] * 6)
    ax[1].set_xticks(range(len(xs2)))
    ax[1].set_xticklabels([l for _, l in right], rotation=40, ha="right", fontsize=7)
    ax[1].set_title("second point: LO 0.35 V/m, a' = 2.45 E1dB' (zero-drift QPSK)", fontsize=9)
    for a in ax:
        a.axhline(0, color="k", lw=.7)
    ax[0].set_ylabel("extra compression X = 1 − g/g_CW (same operating point)")
    fig.suptitle("fig_fv01 — phase-driven compression at the original and the second operating point (static / Hammerstein / slow-state predict X = 0)")
    fig.tight_layout()
    fig.savefig(REP / "fig_fv01_compression_by_case.png", dpi=160)
    plt.close(fig)

    fig, ax = plt.subplots(1, 2, figsize=(12, 3.8))
    z = np.load(REP / "runs" / "LO2_POWERSTEP_Nd4001.npz")
    t = z["t"] * 1e6
    g = np.abs(z["z_full"]) / np.abs(z["z_linear"])
    ax[0].plot(t - 50, np.convolve(g, np.ones(10) / 10, "same"), color=BLUE, lw=1)
    ax[0].axhline(g_static2, color=ORANGE, lw=1, label=f"CW-static gain at a' = {g_static2:.3f}")
    ax[0].set_xlim(-3, 25)
    ax[0].set_xlabel("time from power step a'/3 → a' (µs)")
    ax[0].set_ylabel("gain g(t)")
    ax[0].legend(fontsize=7, frameon=False)
    ax[0].set_title(f"power step at LO' (dominant period {p_pow:.2f} µs)", fontsize=9)
    zz = np.load(REP / "runs" / "LO2_QPSK_SLOW_zd_Nd4001.npz")
    t2 = zz["t"] * 1e6
    u = zz["z_full"] * np.conj(zz["z_linear"])
    for k in range(1, 15):
        tb = 50 + 4 * k
        m, pre = (t2 >= tb) & (t2 < tb + 3.7), (t2 >= tb - .6) & (t2 < tb - .2)
        th = np.angle(u[m] * np.conj(np.mean(u[pre])))
        if np.max(np.abs(th)) > .05:
            ax[1].plot(t2[m] - tb, th / np.pi, color=BLUE, lw=.8)
    ax[1].axhline(0, color="k", lw=.6)
    ax[1].set_xlabel("time after phase step (µs)")
    ax[1].set_ylabel("response phase change (units of π)")
    ax[1].set_title(f"isolated phase steps at LO' (median period {np.median(p_ph):.2f} µs)", fontsize=9)
    fig.suptitle("fig_fv02 — power-step vs phase-step transients at the second operating point")
    fig.tight_layout()
    fig.savefig(REP / "fig_fv02_transients_second_point.png", dpi=160)
    plt.close(fig)


if __name__ == "__main__":
    main()
