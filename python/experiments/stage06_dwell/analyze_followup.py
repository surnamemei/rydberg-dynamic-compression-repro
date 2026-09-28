"""Stage-0.6 follow-up analysis (FU-A..FU-E, pre-registered in 11_followup_preregistration.json)."""
from __future__ import annotations

import json

import numpy as np
import pandas as pd
from scipy import stats

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from common import OUT

C = {"atomic": "#2a78d6", "static_mh": "#eb6834", "static_lti_mh": "#1baf7a", "dsh_tau2.535us": "#eda100", "linear": "#8a8a85"}
LBL = {"atomic": "M_FULL", "static_mh": "M_STATIC (all-zone)", "static_lti_mh": "M_LTI_STATIC (all-zone)",
       "dsh_tau2.535us": "single slow state, τ=2.535 µs (not fitted)", "linear": "linear"}
MODELS = ["atomic", "static_mh", "static_lti_mh", "dsh_tau2.535us"]
TAU = json.loads((OUT / "stage06_config.json").read_text())["tau_atom_s"]
TAU_1E_CONV = float(pd.read_csv(OUT / "03b_tau_atom_sensitivity.csv").set_index("experiment").loc["tau_1.0xE1", "tau_1e_s"])
plt.rcParams.update({"axes.grid": True, "grid.alpha": .25, "axes.spines.top": False, "axes.spines.right": False, "font.size": 9, "lines.linewidth": 1.6})


def load(tag, design="dwell"):
    rows = pd.read_csv(OUT / f"rows_{tag}_{design}.csv").drop_duplicates(["job_id", "model"], keep="last")
    p = OUT / f"rows_{tag}_{design}_dsh.csv"
    if p.exists():
        rows = pd.concat([rows, pd.read_csv(p)], ignore_index=True)
    return rows


def ci(x):
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    if len(x) < 2:
        return np.nan, np.nan, np.nan, len(x)
    h = stats.t.ppf(.975, len(x) - 1) * x.std(ddof=1) / np.sqrt(len(x))
    return x.mean(), x.mean() - h, x.mean() + h, len(x)


def paired(df, ref, metric="D_BLA", reals=None):
    if reals is not None:
        df = df[df.realization.isin(reals)]
    w = df.pivot_table(index=["Pavg_over_P1dB_dB", "variant", "realization"], columns="model", values=metric)
    out = []
    for (p, v), g in w.groupby(level=[0, 1]):
        if v == ref:
            continue
        a = w.loc[(p, ref)]
        b = g.droplevel([0, 1])
        for m in w.columns:
            mu, lo, hi, n = ci((b[m] - a[m]).values)
            out.append({"Pavg_over_P1dB_dB": p, "variant": v, "xi": float(v[3:]) if v.startswith("xi_") else np.nan, "model": m, "metric": metric,
                        "Delta_mean": mu, "Delta_ci_lo": lo, "Delta_ci_hi": hi, "n": n})
        for m in ("static_mh", "static_lti_mh", "dsh_tau2.535us"):
            if m in w.columns:
                mu, lo, hi, n = ci(((b.atomic - a.atomic) - (b[m] - a[m])).values)
                out.append({"Pavg_over_P1dB_dB": p, "variant": v, "xi": float(v[3:]) if v.startswith("xi_") else np.nan, "model": f"full_minus_{m}",
                            "metric": metric, "Delta_mean": mu, "Delta_ci_lo": lo, "Delta_ci_hi": hi, "n": n})
    return pd.DataFrame(out)


def per_waveform(df, metric="D_BLA"):
    w = df.pivot_table(index=["Pavg_over_P1dB_dB", "variant", "realization"], columns="model", values=metric)
    out = []
    for (p, v), g in w.groupby(level=[0, 1]):
        for m in w.columns:
            mu, lo, hi, n = ci(g[m].values)
            out.append({"Pavg_over_P1dB_dB": p, "variant": v, "xi": float(v[3:]) if v.startswith("xi_") else np.nan, "model": m, "mean": mu, "ci_lo": lo, "ci_hi": hi, "n": n})
        for m in ("static_mh", "static_lti_mh", "dsh_tau2.535us"):
            if m in w.columns:
                mu, lo, hi, n = ci((g.atomic - g[m]).values)
                out.append({"Pavg_over_P1dB_dB": p, "variant": v, "xi": float(v[3:]) if v.startswith("xi_") else np.nan, "model": f"full_minus_{m}", "mean": mu, "ci_lo": lo, "ci_hi": hi, "n": n})
    return pd.DataFrame(out)


def get(t, p, xi, m, col="Delta_mean"):
    r = t[(t.Pavg_over_P1dB_dB == p) & (np.isclose(t.xi, xi)) & (t.model == m)]
    return float(r[col].iloc[0]) if len(r) else np.nan


def main():
    verdict = []
    # FU-A -------------------------------------------------------------------------------
    tau = pd.read_csv(OUT / "03b_tau_atom_sensitivity.csv")
    t10 = float(tau.set_index("experiment").loc["tau_1.0xE1", "tau_1e_s"])
    verdict.append({"item": "FU-A tau_1e at converged E1dB", "value": f"{t10 * 1e6:.3f} us (vs 2.535)", "rule": ">20% difference -> report both xi axes",
                    "outcome": "TRIGGERED" if abs(t10 / TAU - 1) > .2 else "not triggered"})
    # FU-B -------------------------------------------------------------------------------
    pl = pd.read_csv(OUT / "12_plateau_steady_state.csv")
    for r in pl.itertuples():
        if r.carrier == "cw":
            verdict.append({"item": f"FU-B H_QS_cw {r.Pavg_over_P1dB_dB:+d} dB", "value": f"{100 * r.plateau_vs_cw_static_rel:+.2f}% vs CW static",
                            "rule": "|dev| < 2%", "outcome": "SUPPORTED" if abs(r.plateau_vs_cw_static_rel) < .02 else "not supported"})
        else:
            verdict.append({"item": f"FU-B H_carrier {r.Pavg_over_P1dB_dB:+d} dB", "value": f"{100 * r.plateau_vs_cw_static_rel:+.1f}% vs CW static",
                            "rule": "< -10%", "outcome": "SUPPORTED" if r.plateau_vs_cw_static_rel < -.10 else "not supported"})
    # FU-C: CW-carrier dwell family vs the same realizations with the QPSK carrier --------
    cw = load("fu_cwcarrier")
    qp = load("main")
    pc = paired(cw, "xi_0.05")
    pq = paired(qp, "xi_0.05", reals=range(4))
    pc["carrier"], pq["carrier"] = "cw", "qpsk"
    both = pd.concat([pc, pq], ignore_index=True)
    both.to_csv(OUT / "14_followup_carrier_ablation_pairs.csv", index=False)
    lo, hi = [get(pc, 3, 4, "full_minus_static_lti_mh", c) for c in ("Delta_ci_lo", "Delta_ci_hi")]
    mu = get(pc, 3, 4, "full_minus_static_lti_mh")
    mu_q = get(pq, 3, 4, "full_minus_static_lti_mh")
    verdict.append({"item": "FU-C amplitude-dwell dynamics without phase modulation (+3 dB, xi=4)", "value": f"full-LTI_MH = {mu:+.3f} [{lo:+.3f},{hi:+.3f}] (QPSK carrier, same reals: {mu_q:+.3f})",
                    "rule": "CI excludes 0, same sign as QPSK", "outcome": "SUPPORTED" if (lo > 0 or hi < 0) and np.sign(mu) == np.sign(mu_q) else "not supported"})
    # FU-D: long dwell ------------------------------------------------------------------------
    longs = {}
    for car in ("qpsk", "cw"):
        tag = f"fu_long_{car}"
        if not (OUT / f"rows_{tag}_dwell.csv").exists():
            continue
        d = load(tag)
        pw = per_waveform(d)
        pw["carrier"] = car
        pr = per_waveform(d, "D_ref_linear")
        pr["carrier"], pr["metric"] = car, "D_ref_linear"
        pw["metric"] = "D_BLA"
        longs[car] = pd.concat([pw, pr], ignore_index=True)
        e4, e32 = get(pw, 3, 4, "full_minus_static_mh", "mean"), get(pw, 3, 32, "full_minus_static_mh", "mean")
        lo32, hi32 = get(pw, 3, 32, "full_minus_static_mh", "ci_lo"), get(pw, 3, 32, "full_minus_static_mh", "ci_hi")
        ok = (abs(e32) < .25 * abs(e4)) or (lo32 <= 0 <= hi32)
        verdict.append({"item": f"FU-D quasi-static convergence, {car} carrier (+3 dB)", "value": f"D_full-D_static: xi=4 {e4:+.3f}, xi=32 {e32:+.3f} [{lo32:+.3f},{hi32:+.3f}]",
                        "rule": "|xi=32| < 25% of |xi=4| or CI includes 0", "outcome": "SUPPORTED" if ok else "not supported"})
    if longs:
        pd.concat(longs.values(), ignore_index=True).to_csv(OUT / "15_followup_long_dwell_per_waveform.csv", index=False)
    # FU-E: single slow state -------------------------------------------------------------------
    rows_e = []
    for name, t in (("QPSK carrier dwell (main, 8 reals)", paired(qp, "xi_0.05")), ("CW carrier dwell (FU-C, 4 reals)", pc)):
        for p in (0, 3):
            for xi in (1, 2, 4):
                f, dsh = get(t, p, xi, "atomic"), get(t, p, xi, "dsh_tau2.535us")
                rows_e.append({"set": name, "Pavg_over_P1dB_dB": p, "contrast": f"xi={xi:g} vs 0.05", "Delta_full": f, "Delta_DSH": dsh,
                               "rel_error": abs(f - dsh) / abs(f) if f else np.nan})
    sh = load("main", "shuffle")
    ps = paired(sh, "ORIGINAL")
    for p in (3, 6):
        r = ps[(ps.Pavg_over_P1dB_dB == p) & (ps.variant == "BLOCK_SHUFFLED_CYCLES")].set_index("model")
        f, dsh = r.loc["atomic", "Delta_mean"], r.loc["dsh_tau2.535us", "Delta_mean"]
        rows_e.append({"set": "shuffle cycles (main, 8 reals)", "Pavg_over_P1dB_dB": p, "contrast": "cycles vs original", "Delta_full": f, "Delta_DSH": dsh,
                       "rel_error": abs(f - dsh) / abs(f)})
    fe = pd.DataFrame(rows_e)
    fe["within_25pct"] = fe.rel_error <= .25
    fe.to_csv(OUT / "16_followup_single_slow_state.csv", index=False)
    for name, g in fe.groupby("set"):
        verdict.append({"item": f"FU-E DSH explains ordering effect: {name}", "value": f"{int(g.within_25pct.sum())}/{len(g)} contrasts within 25%; median rel. error {g.rel_error.median():.2f}",
                        "rule": "all contrasts within 25%", "outcome": "SUPPORTED" if g.within_25pct.all() else "not supported"})
    nm = []
    for tag, design in (("main", "dwell"), ("main", "shuffle"), ("fu_cwcarrier", "dwell"), ("fu_long_qpsk", "dwell"), ("fu_long_cw", "dwell")):
        p = OUT / f"rows_{tag}_{design}_dsh.csv"
        if not p.exists():
            continue
        d = pd.read_csv(p)
        d = d[(d.model == "dsh_tau2.535us") & d.NMSE_vs_full_trace.notna()]
        if "NMSE_vs_full_trace_static_lti_mh" not in d:
            continue
        # the reference-model NMSEs are stored on the first model's row of each job
        ref = pd.read_csv(p)
        ref = ref[ref.NMSE_vs_full_trace_static_lti_mh.notna()][["job_id", "NMSE_vs_full_trace_static_lti_mh"]]
        m = d[["job_id", "Pavg_over_P1dB_dB", "variant", "NMSE_vs_full_trace"]].merge(ref, on="job_id")
        m["DSH_better"] = m.NMSE_vs_full_trace < m.NMSE_vs_full_trace_static_lti_mh
        m["set"] = f"{tag}/{design}"
        nm.append(m)
    if nm:
        nm = pd.concat(nm, ignore_index=True)
        nm.to_csv(OUT / "16b_followup_dsh_trace_nmse.csv", index=False)
        for s_, g in nm.groupby("set"):
            verdict.append({"item": f"FU-E trace NMSE DSH < LTI_MH ({s_}, realization 0)", "value": f"{g.DSH_better.mean():.0%} of {len(g)} waveforms; median NMSE DSH {g.NMSE_vs_full_trace.median():.3f} vs LTI_MH {g.NMSE_vs_full_trace_static_lti_mh.median():.3f}",
                            "rule": ">= 75% of waveforms", "outcome": "SUPPORTED" if g.DSH_better.mean() >= .75 else "not supported"})
    v = pd.DataFrame(verdict)
    v.to_csv(OUT / "14_followup_verdicts.csv", index=False)
    print(v.to_string(), flush=True)
    figures(pl, pc, pq, longs, fe)


def figures(pl, pc, pq, longs, fe):
    z = np.load(OUT / "plateau_gain_traces.npz")
    dt, warm, hold = float(z["dt_s"]), float(z["warm_s"]), float(z["hold_s"])
    fig, ax = plt.subplots(1, 2, figsize=(12, 3.8), sharey=True)
    for a, p in zip(ax, (0, 3)):
        for car, ls in (("cw", "-"), ("qpsk", "--")):
            g = z[f"plateau_p{p:+d}_{car}"]
            t = np.arange(len(g)) * dt * 1e6 - warm * 1e6
            k = np.ones(10) / 10
            a.plot(t, np.convolve(g, k, "same"), ls, color="#2a78d6", lw=1.3, label=f"M_FULL, {'unmodulated' if car == 'cw' else 'QPSK phase'} carrier")
        r = pl[(pl.Pavg_over_P1dB_dB == p) & (pl.carrier == "cw")].iloc[0]
        a.hlines([r.cw_static_gain_aH], 0, hold * 1e6, color="#eb6834", lw=1.4, label="CW-static gain at a_H")
        a.hlines([r.cw_static_gain_aL], -10, 0, color="#eb6834", lw=1.4, ls=":")
        a.hlines([r.cw_static_gain_aL], hold * 1e6, hold * 1e6 + 60, color="#eb6834", lw=1.4, ls=":", label="CW-static gain at a_L")
        a.set_xlim(-10, hold * 1e6 + 40)
        a.set_title(f"single 60 µs plateau at the dwell-family levels, {p:+d} dB", fontsize=9)
        a.set_xlabel("time from plateau start (µs)")
    ax[0].set_ylabel("|IF envelope| / (a · small-signal gain)")
    ax[1].legend(fontsize=7, frameon=False, loc="lower right")
    fig.suptitle("FU-B — plateau steady state: unmodulated carrier converges to CW-static; the QPSK phase carrier compresses 55–60% deeper")
    fig.tight_layout()
    fig.savefig(OUT / "fig15_plateau_carrier_effect.png", dpi=170)
    plt.close(fig)

    fig, ax = plt.subplots(1, 2, figsize=(12, 4), sharex=True)
    for a, p in zip(ax, (0, 3)):
        for t, ls, car in ((pq, "-", "QPSK carrier"), (pc, "--", "unmodulated carrier")):
            for m in MODELS:
                g = t[(t.Pavg_over_P1dB_dB == p) & (t.model == m)].sort_values("xi")
                if g.empty:
                    continue
                a.errorbar(g.xi, g.Delta_mean, yerr=[g.Delta_mean - g.Delta_ci_lo, g.Delta_ci_hi - g.Delta_mean], fmt="o" + ls, ms=3.5, capsize=2,
                           color=C[m], label=f"{LBL[m]}, {car}")
        a.axhline(0, color="k", lw=.7)
        a.set_xscale("log")
        a.set_title(f"Pavg/P1dB = {p:+d} dB (realizations 0–3, paired 95% CI)", fontsize=9)
        a.set_xlabel(f"T_dwell / τ_atom of B   (τ_atom = {TAU * 1e6:.3f} µs; ×{TAU / TAU_1E_CONV:.2f} for τ_1e at E1dB)")
    ax[0].set_ylabel("Δ D_BLA = D(B) − D(A = ξ 0.05)")
    ax[1].legend(fontsize=6.5, frameon=False, ncol=2)
    fig.suptitle("FU-C/E — dwell effect with and without phase modulation, and the single-slow-state control")
    fig.tight_layout()
    fig.savefig(OUT / "fig16_carrier_ablation_dwell.png", dpi=170)
    plt.close(fig)

    if longs:
        fig, ax = plt.subplots(1, 2, figsize=(12, 4), sharey=False)
        for a, car in zip(ax, ("qpsk", "cw")):
            if car not in longs:
                continue
            t = longs[car]
            t = t[t.metric == "D_ref_linear"]
            for m in MODELS:
                g = t[(t.model == m)].sort_values("xi")
                if g["mean"].isna().all():
                    continue
                a.errorbar(g.xi, g["mean"], yerr=[g["mean"] - g.ci_lo, g.ci_hi - g["mean"]], fmt="o-", ms=3.5, capsize=2, color=C[m], label=LBL[m])
            a.set_xscale("log")
            a.set_yscale("log")
            a.set_xlabel("T_dwell / τ_atom")
            a.set_title(f"long-dwell family (2.4 ms), +3 dB, {'QPSK phase' if car == 'qpsk' else 'unmodulated'} carrier", fontsize=9)
        ax[0].set_ylabel("residual distortion / linear-model signal power\n(mean, 95% CI over 4 realizations)")
        ax[0].legend(fontsize=7, frameon=False)
        fig.suptitle("FU-D — does the full model approach the quasi-static prediction at long dwell?")
        fig.tight_layout()
        fig.savefig(OUT / "fig17_long_dwell_quasistatic.png", dpi=170)
        plt.close(fig)

    fig, ax = plt.subplots(figsize=(7, 4.4))
    sets = list(fe.set.unique())
    for i, s_ in enumerate(sets):
        g = fe[fe.set == s_]
        ax.scatter(g.Delta_full, g.Delta_DSH, s=36, color=("#2a78d6", "#eb6834", "#1baf7a")[i % 3], label=s_, zorder=3)
    lim = np.nanmax(np.abs(fe[["Delta_full", "Delta_DSH"]].values)) * 1.1
    xx = np.linspace(-lim, lim, 10)
    ax.plot(xx, xx, color="k", lw=.8)
    ax.fill_between(xx, .75 * xx, 1.25 * xx, color="0.9", zorder=0, label="±25% band (pre-registered)")
    ax.axhline(0, color="k", lw=.5)
    ax.axvline(0, color="k", lw=.5)
    ax.set_xlabel("ordering effect in M_FULL (Δ D_BLA)")
    ax.set_ylabel("same contrast in single-slow-state control")
    ax.legend(fontsize=7, frameon=False)
    ax.set_title("FU-E — how much of the ordering effect does one |E|²-driven slow state (τ from Step 2) explain?", fontsize=9)
    fig.tight_layout()
    fig.savefig(OUT / "fig18_single_slow_state_control.png", dpi=170)
    plt.close(fig)


if __name__ == "__main__":
    main()
