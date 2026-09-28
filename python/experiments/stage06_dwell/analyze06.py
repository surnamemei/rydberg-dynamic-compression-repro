"""Stage-0.6 analysis: matching tables, paired A/B statistics, numerical margins, figures."""
from __future__ import annotations

import json

import numpy as np
import pandas as pd
from scipy import stats

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from common import OUT, ROOT
import run06 as R
import waveforms06 as W

C = {"atomic": "#2a78d6", "static_mh": "#eb6834", "static_lti_mh": "#1baf7a", "static": "#eb6834", "static_lti": "#1baf7a", "linear": "#8a8a85"}
LBL = {"atomic": "M_FULL (thermal transient)", "static_mh": "M_STATIC (all-zone, CW-matched)", "static_lti_mh": "M_LTI_STATIC (all-zone)",
       "static": "M_STATIC F1 (Stage-0.5, fundamental only)", "static_lti": "M_LTI_STATIC F1 (Stage-0.5)", "linear": "small-signal linear"}
CTRL = ("static_mh", "static_lti_mh", "static", "static_lti")
TAU = R.TAU_S
plt.rcParams.update({"axes.grid": True, "grid.alpha": .25, "axes.spines.top": False, "axes.spines.right": False, "font.size": 9, "lines.linewidth": 1.6})


def load(tag, design):
    p = OUT / f"rows_{tag}_{design}.csv"
    if not p.exists():
        return pd.DataFrame()
    return pd.read_csv(p).drop_duplicates(["job_id", "model"], keep="last")


def ci(x):
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    n = len(x)
    m = float(np.mean(x)) if n else np.nan
    sd = float(np.std(x, ddof=1)) if n > 1 else np.nan
    h = float(stats.t.ppf(.975, n - 1) * sd / np.sqrt(n)) if n > 1 else np.nan
    return m, sd, m - h, m + h, n


# ----------------------------------------------------------------------------- waveform matching
def waveform_tables():
    rows = []
    for design, ref in (("dwell", "xi_0.05"), ("shuffle", "ORIGINAL")):
        for r in range(8):
            amps, meta, car, _, _ = R.realization(design, r)
            a0 = amps[ref]
            thr = np.median(a0) if design == "shuffle" else float((a0.min() + a0.max()) / 2)
            base = None
            for v, a in amps.items():
                env = a * car
                runs = W.dwell_runs(a, thr)
                f = np.fft.fftfreq(len(env), 1e-9)
                P = np.abs(np.fft.fft(env)) ** 2
                row = {"design": design, "realization": r, "variant": v, "reference_variant": ref,
                       "RMS_unit": float(np.sqrt(np.mean(a ** 2))), "peak_unit": float(a.max()), "energy_unit": float(np.sum(a ** 2)),
                       "threshold_unit": thr, "occupancy_above_threshold": float(np.mean(a > thr)),
                       "KS_vs_reference": W.ks_distance(a, a0), "sorted_samples_identical": bool(np.array_equal(np.sort(a), np.sort(a0))),
                       "B99_complex_env_Hz": W.occupied99(env), "B99_amplitude_fluct_Hz": W.occupied99(a - a.mean()),
                       "inband_pm3MHz_power_fraction": float(P[np.abs(f) <= 3e6].sum() / P.sum()),
                       "dwell_mean_event_ns": float(np.mean(runs)), "dwell_time_weighted_ns": float(np.sum(runs.astype(float) ** 2) / np.sum(runs)),
                       "dwell_p90_ns": float(np.quantile(runs, .9)), "n_high_events": int(len(runs)),
                       "env_corr_1e_ns": W.env_corr_1e(a), **{f"meta_{k}": val for k, val in meta[v].items()}}
                row["T_dwell_over_tau"] = row["dwell_time_weighted_ns"] * 1e-9 / TAU
                row["C_env_tau_env_over_tau"] = row["env_corr_1e_ns"] * 1e-9 / TAU
                if v == ref:
                    base = row
                for k in ("RMS_unit", "peak_unit", "energy_unit", "B99_complex_env_Hz", "inband_pm3MHz_power_fraction"):
                    row[f"{k}_mismatch_pct"] = 100 * (row[k] / base[k] - 1)
                row["occupancy_mismatch_pct_points"] = 100 * (row["occupancy_above_threshold"] - base["occupancy_above_threshold"])
                rows.append(row)
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "04_controlled_pair_configs.csv", index=False)
    return df


# ----------------------------------------------------------------------------- paired statistics
def paired(df, ref, key_cols=("Pavg_over_P1dB_dB",)):
    out = []
    metrics = ("D_BLA", "AIR", "EVM", "SER", "loss_vs_linear")
    wide = df.pivot_table(index=["realization", *key_cols, "variant"], columns="model", values=list(metrics))
    for keys, g in df.groupby([*key_cols, "variant"]):
        if keys[-1] == ref:
            continue
        pk = keys[:-1]
        for met in metrics:
            for m in ("atomic", *CTRL, "linear"):
                vals, dyn = [], {c: [] for c in CTRL}
                for r in sorted(g.realization.unique()):
                    try:
                        b = wide.loc[(r, *pk, keys[-1]), (met, m)]
                        a = wide.loc[(r, *pk, ref), (met, m)]
                    except KeyError:
                        continue
                    vals.append(b - a)
                    if m == "atomic":
                        for c in CTRL:
                            dyn[c].append((b - a) - (wide.loc[(r, *pk, keys[-1]), (met, c)] - wide.loc[(r, *pk, ref), (met, c)]))
                mu, sd, lo, hi, n = ci(vals)
                row = {**dict(zip(key_cols, pk)), "variant": keys[-1], "reference": ref, "metric": met, "model": m,
                       "Delta_mean": mu, "Delta_sd": sd, "Delta_ci_lo": lo, "Delta_ci_hi": hi, "n": n}
                if m == "atomic":
                    for c in CTRL:
                        mu2, sd2, lo2, hi2, _ = ci(dyn[c])
                        row.update({f"full_minus_{c}_mean": mu2, f"full_minus_{c}_ci_lo": lo2, f"full_minus_{c}_ci_hi": hi2})
                out.append(row)
    return pd.DataFrame(out)


def per_waveform_excess(df):
    out = []
    for keys, g in df.groupby(["Pavg_over_P1dB_dB", "variant"]):
        w = g.pivot_table(index="realization", columns="model", values=["D_BLA", "AIR", "loss_vs_linear"])
        row = {"Pavg_over_P1dB_dB": keys[0], "variant": keys[1]}
        for m in ("atomic", *CTRL, "linear"):
            for met in ("D_BLA", "AIR", "loss_vs_linear"):
                mu, sd, lo, hi, n = ci(w[(met, m)])
                row.update({f"{met}_{m}_mean": mu, f"{met}_{m}_ci_lo": lo, f"{met}_{m}_ci_hi": hi})
        for c in CTRL:
            mu, _, lo, hi, n = ci(w[("D_BLA", "atomic")] - w[("D_BLA", c)])
            row.update({f"dynexcess_D_vs_{c}_mean": mu, f"dynexcess_D_vs_{c}_ci_lo": lo, f"dynexcess_D_vs_{c}_ci_hi": hi})
            mu, _, lo, hi, n = ci(w[("AIR", c)] - w[("AIR", "atomic")])
            row.update({f"dynexcess_loss_vs_{c}_mean": mu, f"dynexcess_loss_vs_{c}_ci_lo": lo, f"dynexcess_loss_vs_{c}_ci_hi": hi})
        row["n"] = n
        out.append(row)
    return pd.DataFrame(out)


# ----------------------------------------------------------------------------- figures
def fig_nd():
    d = pd.read_csv(OUT / "02_nd_convergence.csv")
    fig, ax = plt.subplots(1, 3, figsize=(13, 3.8))
    cols = ["#2a78d6", "#eb6834", "#1baf7a", "#4a3aa7"]
    for (case, g), c in zip(d[d.case != "E_stage05_dwell_pair"].groupby("case"), cols):
        g = g.sort_values("Nd")
        ax[0].plot(g.Nd, g.AIR, "o-", color=c, label=case)
        ax[1].plot(g.Nd, g.AIR_diff_vs_top, "o-", color=c, label=case)
        ax[2].loglog(g.Nd[g.Nd < g.Nd.max()], g.baseband_NMSE_vs_top[g.Nd < g.Nd.max()], "o-", color=c, label=case)
    e = d[d.case == "E_stage05_dwell_pair"].sort_values("Nd")
    ax[2].loglog(e.Nd[e.Nd < e.Nd.max()], e.baseband_NMSE_vs_top[e.Nd < e.Nd.max()], "s--", color="#8a8a85", label="Stage-0.5 Stage-K pair")
    ax[1].axhspan(-.005, .005, color="0.85", zorder=0, label="±0.005 tolerance")
    ax[2].axhline(1e-3, color="k", ls=":", lw=1, label="NMSE tolerance 1e-3")
    for a in ax:
        a.axvline(1501, color="0.5", lw=.8, ls="--")
        a.axvline(4001, color="k", lw=.8)
        a.set_xlabel("Nd (velocity classes before ±3σ truncation)")
    ax[0].set_ylabel("AIR (bit/native symbol)")
    ax[1].set_ylabel("AIR − AIR(Nd=8001)")
    ax[2].set_ylabel("baseband NMSE vs Nd=8001")
    ax[1].legend(fontsize=7, frameon=False)
    ax[2].legend(fontsize=7, frameon=False)
    fig.suptitle("Step 1 — thermal quadrature convergence (dashed: Nd=1501, solid: decisive Nd=4001)")
    fig.tight_layout()
    fig.savefig(OUT / "fig01_nd_convergence.png", dpi=170)
    plt.close(fig)


def fig_tau():
    z = np.load(OUT / "tau_atom_traces.npz")
    d = pd.read_csv(OUT / "03_tau_atom_results.csv")
    exps = [("LO_step_+1pct", "small LO step (+1%)"), ("IF_step_small_0.05to0.10E1", "IF envelope step, linear point"),
            ("IF_step_P1dB_1.00to1.05E1", "IF envelope step at P1dB (+5%) — defines τ_atom"), ("IF_plateau_0.5to1.5E1_off_recovery", "recovery after 30 µs plateau at 1.5·E1dB")]
    fig, ax = plt.subplots(1, 4, figsize=(15, 3.6))
    for a, (e, title) in zip(ax, exps):
        for nd, c, ls in ((1501, "#eb6834", "--"), (4001, "#2a78d6", "-"), (8001, "#1baf7a", ":")):
            k = f"{nd}__{e}"
            if k + "__t" in z:
                a.semilogy(z[k + "__t"] * 1e6, np.abs(z[k + "__e"]) + 1e-6, ls, color=c, lw=1.2, label=f"Nd={nd}")
        row = d[(d.experiment == e) & (d.Nd == 4001)].iloc[0]
        a.axhline(np.exp(-1), color="k", lw=.7, ls=":")
        a.axvline(row.tau_1e_s * 1e6, color="k", lw=.8)
        a.set_title(f"{title}\nτ_1e={row.tau_1e_s*1e6:.3f} µs, τ_90={row.tau_90_s*1e6:.2f} µs", fontsize=8)
        a.set_xlim(0, 20)
        a.set_ylim(1e-3, 3)
        a.set_xlabel("time after step (µs)")
    ax[0].set_ylabel("|normalized excess response|")
    ax[0].legend(fontsize=7, frameon=False)
    fig.suptitle("Step 2 — atomic relaxation (full thermal averaging); τ_atom ≡ τ_1e of the P1dB IF-envelope step = %.2f µs" % (TAU * 1e6))
    fig.tight_layout()
    fig.savefig(OUT / "fig02_atomic_relaxation_tau.png", dpi=170)
    plt.close(fig)


def fig_waveforms(cfgdf):
    ampd, _, _, _, _ = R.realization("dwell", 0)
    amps, _, _, _, _ = R.realization("shuffle", 0)
    fig, ax = plt.subplots(8, 1, figsize=(12, 12), sharex=True)
    t = np.arange(R.N) * 1e-3
    sel = [("xi_0.05", ampd), ("xi_0.25", ampd), ("xi_1", ampd), ("xi_4", ampd), ("ORIGINAL", amps), ("BLOCK_SHUFFLED_CYCLES", amps), ("BLOCK_SHUFFLED_FIXED", amps), ("SHUFFLED", amps)]
    for a, (v, src) in zip(ax, sel):
        a.plot(t[:200_000:5], src[v][:200_000:5], color="#2a78d6" if v.startswith("xi") else "#4a3aa7", lw=.6)
        a.set_ylabel(v.replace("_", " "), fontsize=7)
    ax[-1].set_xlabel("time (µs)")
    fig.suptitle("Step 3 — controlled amplitude envelopes (unit RMS; realization 0; first 200 µs). Dwell family: exact common multiset; shuffle set: exact common multiset")
    fig.tight_layout()
    fig.savefig(OUT / "fig03_controlled_envelopes.png", dpi=150)
    plt.close(fig)

    fig, ax = plt.subplots(1, 2, figsize=(11, 3.6))
    for a, (src, vs) in zip(ax, ((ampd, ["xi_0.05", "xi_1", "xi_4"]), (amps, ["ORIGINAL", "BLOCK_SHUFFLED_CYCLES", "SHUFFLED"]))):
        bins = np.linspace(0, max(src[v].max() for v in vs) * 1.02, 120)
        for v, c, ls in zip(vs, ("#2a78d6", "#eb6834", "#1baf7a"), ("-", "--", ":")):
            a.hist(src[v], bins=bins, histtype="step", color=c, ls=ls, lw=1.6, label=v, density=True)
        a.set_yscale("log")
        a.set_xlabel("|E| / RMS")
        a.legend(fontsize=7, frameon=False)
    ks = cfgdf.KS_vs_reference.max()
    fig.suptitle(f"Step 4 — amplitude PDFs overlap exactly (sorted samples identical in every realization; max KS = {ks:.1e})")
    ax[0].set_ylabel("density")
    fig.tight_layout()
    fig.savefig(OUT / "fig04_amplitude_pdf_match.png", dpi=170)
    plt.close(fig)

    fig, ax = plt.subplots(1, 2, figsize=(11, 3.6))
    cm = plt.get_cmap("viridis")
    for i, v in enumerate(["xi_0.05", "xi_0.1", "xi_0.25", "xi_0.5", "xi_1", "xi_2", "xi_4"]):
        a0 = ampd[v]
        runs = np.sort(W.dwell_runs(a0, (a0.min() + a0.max()) / 2)).astype(float)
        w = np.cumsum(runs[::-1])[::-1] / runs.sum()
        ax[0].loglog(runs / (TAU * 1e9), w, color=cm(i / 6), label=v)
    for v, c in zip(["ORIGINAL", "BLOCK_SHUFFLED_CYCLES", "BLOCK_SHUFFLED_FIXED", "SHUFFLED"], ("#2a78d6", "#eb6834", "#1baf7a", "#8a8a85")):
        a0 = amps[v]
        runs = np.sort(W.dwell_runs(a0, np.quantile(a0, .75))).astype(float)
        w = np.cumsum(runs[::-1])[::-1] / runs.sum()
        ax[1].loglog(runs / (TAU * 1e9), w, color=c, label=v)
    for a in ax:
        a.set_xlabel("high-field dwell / τ_atom")
        a.legend(fontsize=7, frameon=False)
    ax[0].set_ylabel("fraction of high-field time in events ≥ dwell")
    ax[0].set_title("dwell family (threshold = level midpoint)", fontsize=9)
    ax[1].set_title("shuffle set (threshold = 75th percentile)", fontsize=9)
    fig.suptitle("Step 4/5 — time-weighted dwell distributions (same occupancy, different temporal organization)")
    fig.tight_layout()
    fig.savefig(OUT / "fig05_dwell_distribution_match.png", dpi=170)
    plt.close(fig)

    fig, ax = plt.subplots(1, 2, figsize=(11, 3.6))
    _, _, car, _, _ = R.realization("dwell", 0)
    for a, (src, vs) in zip(ax, ((ampd, ["xi_0.05", "xi_0.25", "xi_1", "xi_4"]), (amps, ["ORIGINAL", "BLOCK_SHUFFLED_CYCLES", "BLOCK_SHUFFLED_FIXED", "SHUFFLED"]))):
        for v, c in zip(vs, ("#2a78d6", "#eb6834", "#1baf7a", "#8a8a85")):
            x = src[v] * car
            seg = 1 << 16
            k = len(x) // seg
            P = np.mean(np.abs(np.fft.fft(x[:k * seg].reshape(k, seg), axis=1)) ** 2, axis=0)
            f = np.fft.fftfreq(seg, 1e-9)
            o = np.argsort(f)
            sm = np.convolve(P[o], np.ones(9) / 9, "same")
            a.semilogy(f[o] / 1e6, sm / sm.max(), color=c, lw=1, label=v)
        a.axvspan(-3, 3, color="0.9", zorder=0, label="receiver band ±3 MHz")
        a.set_xlim(-30, 30)
        a.set_ylim(1e-7, 2)
        a.set_xlabel("offset from IF (MHz)")
        a.legend(fontsize=7, frameon=False)
    ax[0].set_ylabel("normalized PSD of E(t)")
    fig.suptitle("Step 4 — spectra of the complex envelope incl. QPSK phase carrier (realization 0)")
    fig.tight_layout()
    fig.savefig(OUT / "fig06_spectra_match.png", dpi=170)
    plt.close(fig)


def _xi(v):
    return float(v.split("_")[1])


def fig_pairs(pairs, pw):
    xs = sorted({_xi(v) for v in pairs[pairs.variant.str.startswith("xi")].variant})
    powers = sorted(pairs.Pavg_over_P1dB_dB.unique())
    for fname, models, title in (("fig07_full_atomic_pair_difference.png", ["atomic"], "M_FULL"),
                                 ("fig08_static_control_pair_difference.png", ["atomic", "static_mh", "static"], "M_STATIC (all-zone and Stage-0.5 F1) vs M_FULL"),
                                 ("fig09_lti_static_control_pair_difference.png", ["atomic", "static_lti_mh", "static_lti"], "M_LTI_STATIC (all-zone and Stage-0.5 F1) vs M_FULL")):
        fig, ax = plt.subplots(2, len(powers), figsize=(3.4 * len(powers), 6.2), sharex=True)
        for j, p in enumerate(powers):
            for i, met in enumerate(("D_BLA", "AIR")):
                a = ax[i, j]
                for m in models:
                    g = pairs[(pairs.metric == met) & (pairs.model == m) & (pairs.Pavg_over_P1dB_dB == p)].copy()
                    g["xi"] = g.variant.map(_xi)
                    g = g.sort_values("xi")
                    ls = "--" if m in ("static", "static_lti") else "-"
                    a.errorbar(g.xi, g.Delta_mean, yerr=[g.Delta_mean - g.Delta_ci_lo, g.Delta_ci_hi - g.Delta_mean], fmt="o" + ls, ms=4, capsize=2,
                               color=C[m], label=LBL[m])
                a.axhline(0, color="k", lw=.7)
                a.set_xscale("log")
                if i == 0:
                    a.set_title(f"Pavg/P1dB = {p:+d} dB", fontsize=9)
                if i == 1:
                    a.set_xlabel("nominal T_dwell / τ_atom of B (A: 0.05)")
            ax[0, 0].set_ylabel("Δ D_BLA = D(B) − D(A)")
            ax[1, 0].set_ylabel("Δ AIR = AIR(B) − AIR(A)")
        ax[0, -1].legend(fontsize=7, frameon=False)
        fig.suptitle(f"A/B pair difference vs dwell — {title} (mean and paired 95% t-CI over 8 realizations)")
        fig.tight_layout()
        fig.savefig(OUT / fname, dpi=170)
        plt.close(fig)

    fig, ax = plt.subplots(1, 2, figsize=(11, 4))
    cm = {p: c for p, c in zip(powers, ("#8a8a85", "#2a78d6", "#eb6834", "#e34948"))}
    for p in powers:
        g = pw[(pw.Pavg_over_P1dB_dB == p) & pw.variant.str.startswith("xi")].copy()
        g["xi"] = g.variant.map(_xi)
        g = g.sort_values("xi")
        for a, met in zip(ax, ("dynexcess_D_vs_static_lti_mh", "dynexcess_loss_vs_static_lti_mh")):
            a.errorbar(g.xi, g[met + "_mean"], yerr=[g[met + "_mean"] - g[met + "_ci_lo"], g[met + "_ci_hi"] - g[met + "_mean"]], fmt="o-", ms=4, capsize=2, color=cm[p], label=f"{p:+d} dB")
    for a in ax:
        a.set_xscale("log")
        a.axhline(0, color="k", lw=.7)
        a.set_xlabel("nominal T_dwell / τ_atom")
        a.legend(fontsize=7, frameon=False, title="Pavg/P1dB")
    ax[0].set_ylabel("D_BLA(full) − D_BLA(LTI+static, all-zone)")
    ax[1].set_ylabel("AIR(LTI+static) − AIR(full)  [bit/symbol]")
    fig.suptitle("Step 11 — dynamic excess (full minus the best tested LTI+static control) vs dwell")
    fig.tight_layout()
    fig.savefig(OUT / "fig10_dynamic_excess_loss.png", dpi=170)
    plt.close(fig)


def fig_loss_vs_dwell(dw, sh, cfgdf):
    tw = cfgdf.set_index(["design", "realization", "variant"])
    fig, ax = plt.subplots(2, 2, figsize=(11, 7))
    for col, (df, xkey, xl) in enumerate(((dw, "T_dwell_over_tau", "time-weighted T_dwell / τ_atom"), (sh, "C_env_tau_env_over_tau", "C_env = τ_env / τ_atom"))):
        if df.empty:
            continue
        d = df[df.Pavg_over_P1dB_dB == 3].copy()
        d["x"] = [tw.loc[(r.design, r.realization, r.variant), xkey] for r in d.itertuples()]
        for i, met in enumerate(("loss_vs_linear", "D_BLA")):
            a = ax[i, col]
            for m in ("atomic", "static_mh", "static_lti_mh"):
                g = d[d.model == m].groupby("variant").agg(x=("x", "mean"), y=(met, "mean"), s=(met, "std")).sort_values("x")
                a.errorbar(g.x, g.y, yerr=g.s, fmt="o-", ms=4, capsize=2, color=C[m], label=LBL[m])
            a.set_xscale("log")
            a.set_xlabel(xl)
            a.set_ylabel("AIR loss vs linear model (bit/symbol)" if met == "loss_vs_linear" else "D_BLA")
        ax[0, col].set_title(("dwell family" if col == 0 else "shuffle set") + " at +3 dB (mean ± SD over realizations)", fontsize=9)
    ax[0, 0].legend(fontsize=7, frameon=False)
    fig.suptitle("Step 11 — dimensionless organization of loss")
    fig.tight_layout()
    fig.savefig(OUT / "fig11_loss_vs_dwell_over_tau.png", dpi=170)
    plt.close(fig)


def fig_traces():
    fig, ax = plt.subplots(3, 2, figsize=(13, 8.5), sharex="col")
    for col, v in enumerate(("xi_0.05", "xi_4")):
        d = np.load(OUT / "traces" / "main" / f"dwell_r0_{v}_p+3_Nd4001_dt1.npz")
        cp = int(d["cp"])
        x = d["x"][cp:]
        t = np.arange(len(x)) * .05
        g = np.vdot(x, d["y_linear"][cp:]) / np.vdot(x, x)
        sl = (t > 100) & (t < 190)
        ax[0, col].plot(t[sl], d["a"][sl], color="k", lw=.7)
        ax[0, col].set_ylabel("input |E| / RMS")
        ax[0, col].set_title(f"{v.replace('_', ' = ')} (T_dwell/τ nominal), +3 dB, realization 0", fontsize=9)
        ax[1, col].plot(t[sl], np.abs(g * x)[sl], color=C["linear"], lw=.8, label="linear small-signal prediction")
        for m in ("atomic", "static_mh", "static_lti_mh"):
            ax[1, col].plot(t[sl], np.abs(d["y_" + m][cp:])[sl], color=C[m], lw=.8, label=LBL[m])
        ax[1, col].set_ylabel("|baseband output|")
        for m in ("atomic", "static_mh", "static_lti_mh"):
            gi = np.abs(d["y_" + m][cp:]) / np.maximum(np.abs(g * x), 1e-12)
            k = np.ones(40) / 40
            ax[2, col].plot(t[sl], np.convolve(gi, k, "same")[sl], color=C[m], lw=1)
        ax[2, col].set_ylabel("instantaneous |gain| / small-signal\n(2 µs moving average)")
        ax[2, col].set_xlabel("time (µs)")
    ax[1, 1].legend(fontsize=7, frameon=False)
    fig.suptitle("Step 12 — short vs long dwell: slow compression onset and post-excursion gain suppression appear only in M_FULL")
    fig.tight_layout()
    fig.savefig(OUT / "fig12_short_vs_long_dwell_response.png", dpi=170)
    plt.close(fig)


def fig_shuffle(sp, sw):
    powers = sorted(sp.Pavg_over_P1dB_dB.unique())
    vs = ["BLOCK_SHUFFLED_CYCLES", "BLOCK_SHUFFLED_FIXED", "SHUFFLED"]
    fig, ax = plt.subplots(2, len(powers), figsize=(3.6 * len(powers), 6.4), sharex=True)
    for j, p in enumerate(powers):
        for i, met in enumerate(("D_BLA", "AIR")):
            a = ax[i, j]
            for k, m in enumerate(("atomic", "static_mh", "static_lti_mh")):
                g = sp[(sp.metric == met) & (sp.model == m) & (sp.Pavg_over_P1dB_dB == p)].set_index("variant").reindex(vs)
                xx = np.arange(len(vs)) + (k - 1) * .22
                a.bar(xx, g.Delta_mean, width=.2, color=C[m], label=LBL[m])
                a.errorbar(xx, g.Delta_mean, yerr=[g.Delta_mean - g.Delta_ci_lo, g.Delta_ci_hi - g.Delta_mean], fmt="none", ecolor="k", capsize=2, lw=.8)
            a.axhline(0, color="k", lw=.7)
            a.set_xticks(range(len(vs)))
            a.set_xticklabels(["block (cycles)", "block (fixed)", "sample shuffle"], fontsize=7)
            if i == 0:
                a.set_title(f"Pavg/P1dB = {p:+d} dB", fontsize=9)
        ax[0, 0].set_ylabel("Δ D_BLA vs ORIGINAL")
        ax[1, 0].set_ylabel("Δ AIR vs ORIGINAL")
    ax[0, -1].legend(fontsize=7, frameon=False)
    fig.suptitle("Step 10 — temporal-order shuffle test (exact amplitude multiset); paired 95% CI over 8 realizations")
    fig.tight_layout()
    fig.savefig(OUT / "fig13_temporal_shuffle_test.png", dpi=170)
    plt.close(fig)


def fig_power(pairs, pw):
    fig, ax = plt.subplots(1, 3, figsize=(13, 3.8))
    for v, c in (("xi_0.25", "#1baf7a"), ("xi_1", "#2a78d6"), ("xi_4", "#e34948")):
        for a, (met, m, lab) in zip(ax[:2], (("D_BLA", "atomic", "Δ D_BLA full (B−A)"), ("AIR", "atomic", "Δ AIR full (B−A)"))):
            g = pairs[(pairs.metric == met) & (pairs.model == m) & (pairs.variant == v)].sort_values("Pavg_over_P1dB_dB")
            a.errorbar(g.Pavg_over_P1dB_dB, g.Delta_mean, yerr=[g.Delta_mean - g.Delta_ci_lo, g.Delta_ci_hi - g.Delta_mean], fmt="o-", capsize=2, color=c, label=v)
            a.set_ylabel(lab)
        g = pw[pw.variant == v].sort_values("Pavg_over_P1dB_dB")
        ax[2].errorbar(g.Pavg_over_P1dB_dB, g.dynexcess_D_vs_static_lti_mh_mean, yerr=[g.dynexcess_D_vs_static_lti_mh_mean - g.dynexcess_D_vs_static_lti_mh_ci_lo, g.dynexcess_D_vs_static_lti_mh_ci_hi - g.dynexcess_D_vs_static_lti_mh_mean], fmt="o-", capsize=2, color=c, label=v)
    ax[2].set_ylabel("D_BLA(full) − D_BLA(LTI+static)")
    for a in ax:
        a.axhline(0, color="k", lw=.7)
        a.set_xlabel("Pavg / P1dB (dB, converged E1dB)")
        a.legend(fontsize=7, frameon=False)
    fig.suptitle("Step 6 — power dependence of the dwell effect (B vs A = ξ 0.05); paired 95% CI")
    fig.tight_layout()
    fig.savefig(OUT / "fig14_power_dependence_of_dwell_effect.png", dpi=170)
    plt.close(fig)


# ----------------------------------------------------------------------------- numerical validation
def numerical(dw, sh):
    rows = []
    for tag in ("num_Nd8001", "num_dt0p5"):
        for design, main, ref in (("dwell", dw, "xi_0.05"), ("shuffle", sh, "ORIGINAL")):
            nv = load(tag, design)
            if nv.empty:
                continue
            for r in nv.itertuples():
                if r.model not in ("atomic",):
                    continue
                m = main[(main.realization == r.realization) & (main.variant == r.variant) & (main.Pavg_over_P1dB_dB == r.Pavg_over_P1dB_dB) & (main.model == "atomic")]
                if m.empty:
                    continue
                m = m.iloc[0]
                rows.append({"check": tag, "design": design, "realization": r.realization, "variant": r.variant, "Pavg_over_P1dB_dB": r.Pavg_over_P1dB_dB,
                             "D_BLA_main": m.D_BLA, "D_BLA_check": r.D_BLA, "D_BLA_abs_diff": abs(r.D_BLA - m.D_BLA),
                             "AIR_main": m.AIR, "AIR_check": r.AIR, "AIR_abs_diff": abs(r.AIR - m.AIR)})
    df = pd.DataFrame(rows)
    # effect-level check: Delta(B - ref) recomputed with the check runs
    eff = []
    for (tag, design, r), g in df.groupby(["check", "design", "realization"]) if not df.empty else []:
        ref = "xi_0.05" if design == "dwell" else "ORIGINAL"
        if ref not in set(g.variant):
            continue
        g0 = g[g.variant == ref].iloc[0]
        for x in g[g.variant != ref].itertuples():
            for met in ("D_BLA", "AIR"):
                dm = getattr(x, f"{met}_main") - g0[f"{met}_main"]
                dc = getattr(x, f"{met}_check") - g0[f"{met}_check"]
                eff.append({"check": tag, "design": design, "realization": r, "variant": x.variant, "metric": met,
                            "Delta_main": dm, "Delta_check": dc, "Delta_abs_diff": abs(dc - dm), "margin_ratio": abs(dm) / max(abs(dc - dm), 1e-12)})
    out = pd.concat([df.assign(level="per_waveform"), pd.DataFrame(eff).assign(level="pair_delta")], ignore_index=True)
    out.to_csv(OUT / "07_numerical_validation.csv", index=False)
    return out


def main():
    cfgdf = waveform_tables()
    dw, sh = load("main", "dwell"), load("main", "shuffle")
    pairs = paired(dw, "xi_0.05")
    pw = per_waveform_excess(dw)
    pairs["T_dwell_over_tau_B_mean"] = pairs.variant.map(cfgdf[cfgdf.design == "dwell"].groupby("variant").T_dwell_over_tau.mean())
    pairs.to_csv(OUT / "05_controlled_pair_results.csv", index=False)
    pw.to_csv(OUT / "05b_per_waveform_dynamic_excess.csv", index=False)
    if not sh.empty:
        sp = paired(sh, "ORIGINAL")
        sw = per_waveform_excess(sh)
        sp.to_csv(OUT / "06_temporal_shuffle_results.csv", index=False)
        sw.to_csv(OUT / "06b_shuffle_per_waveform_excess.csv", index=False)
        fig_shuffle(sp, sw)
    numerical(dw, sh)
    fig_nd()
    fig_tau()
    fig_waveforms(cfgdf)
    fig_pairs(pairs, pw)
    fig_loss_vs_dwell(dw, sh, cfgdf)
    fig_traces()
    fig_power(pairs, pw)


if __name__ == "__main__":
    main()
