"""P1 analysis: metric robustness of the decisive Section IV ordering contrasts (preregistered, P1).

D (self-normalized BLA residual) and D_ref (residual / linear-model BLA fit AC power) for FIR spans 2/4/8/16 us,
lag columns from the full periodic input.  The archived zero-fill 2-us D is recomputed as a pipeline check.
"""
from __future__ import annotations

import os
os.environ.setdefault("OMP_NUM_THREADS", "2")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "2")
os.environ.setdefault("MKL_NUM_THREADS", "2")

import glob
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from rev_common import REV, ROOT
import run06 as R

CP, SKIP, DEC, N = R.CP, R.SKIP, R.DEC, R.N
NPER = N // DEC
MODELS = ["atomic", "static", "static_lti", "static_mh", "static_lti_mh", "linear"]
SPANS = {2: (-10, 30), 4: (-20, 60), 8: (-40, 120), 16: (-80, 240)}
REFS = {"dwell": "xi_0.05", "shuffle": "ORIGINAL"}
P1 = REV / "p1"
ARCH = ROOT / "results" / "stage06_dwell_physics"


def design(x, w0, w1, lags, method):
    n = w1 - w0
    X = np.empty((n, len(lags) + 1), complex)
    if method == "zerofill":                      # identical to run06.fir_cols on xs = x[w0:w1]
        xs = x[w0:w1]
        for j, lag in enumerate(lags):
            c = np.zeros(n, complex)
            if lag >= 0:
                c[lag:] = xs[:n - lag]
            else:
                c[:lag] = xs[-lag:]
            X[:, j] = c
    else:                                         # full periodic context: x has period NPER after the cyclic prefix
        idx = np.arange(w0, w1)
        for j, lag in enumerate(lags):
            k = idx - lag
            k = np.where(k >= len(x), k - NPER, k)
            X[:, j] = x[k]
    X[:, -1] = 1
    return X


def job_metrics(path):
    z = np.load(path)
    x = z["x"]
    w0, w1 = (CP + SKIP) // DEC, len(x) - 2000 // DEC
    ys = {m: z[f"y_{m}"][w0:w1] for m in MODELS}
    rows = []
    for span, (lo, hi) in SPANS.items():
        for method in (("context", "zerofill") if span == 2 else ("context",)):
            Q, _ = np.linalg.qr(design(x, w0, w1, np.arange(lo, hi + 1), method))
            out = {}
            for m, y in ys.items():
                fit = Q @ (Q.conj().T @ y)
                out[m] = (float(np.mean(np.abs(y - fit) ** 2)), float(np.mean(np.abs(fit - fit.mean()) ** 2)))
            for m, (res, fac) in out.items():
                rows.append({"model": m, "span_us": span, "method": method, "D": res / fac, "D_ref": res / out["linear"][1],
                             "res_pow": res, "fit_ac_pow": fac})
    return rows


def parse(name):
    stem = Path(name).stem                      # e.g. dwell_r3_xi_0.05_p+3_Nd4001_dt1
    design, rest = stem.split("_r", 1)
    r, rest = rest.split("_", 1)
    variant, rest = rest.rsplit("_p", 1)
    pdb, nd, dt = rest.split("_")
    return {"job_id": stem, "design": design, "realization": int(r), "variant": variant, "Pavg_over_P1dB_dB": int(pdb),
            "Nd": int(nd[2:]), "dt_ns": float(dt[2:])}


def ci(v):
    v = np.asarray(v, float)
    n = len(v)
    h = stats.t.ppf(.975, n - 1) * v.std(ddof=1) / np.sqrt(n) if n > 1 else np.nan
    return float(v.mean()), float(v.mean() - h), float(v.mean() + h), n


def contrasts(df, keys=("method", "span_us", "metric")):
    out = []
    long = df.melt(id_vars=["design", "realization", "variant", "Pavg_over_P1dB_dB", "model", "span_us", "method"],
                   value_vars=["D", "D_ref"], var_name="metric")
    for (design, pdb, method, span, metric), g in long.groupby(["design", "Pavg_over_P1dB_dB", "method", "span_us", "metric"]):
        piv = g.pivot_table(index=["realization", "model"], columns="variant", values="value")
        ref = REFS[design]
        for var in [v for v in piv.columns if v != ref]:
            d = (piv[var] - piv[ref]).unstack("model").dropna()
            if len(d) < 2:
                continue
            row = {"design": design, "Pavg_over_P1dB_dB": pdb, "method": method, "span_us": span, "metric": metric, "variant": var}
            for m in MODELS:
                row[f"d_{m}"], row[f"d_{m}_lo"], row[f"d_{m}_hi"], row["n"] = ci(d[m])
            for sur in ("static_mh", "static_lti_mh"):
                row[f"exc_{sur}"], row[f"exc_{sur}_lo"], row[f"exc_{sur}_hi"], _ = ci(d["atomic"] - d[sur])
            out.append(row)
    return pd.DataFrame(out)


def numerical(df_all):
    """Numerical differences of the decisive contrasts from the archived realization-0 Nd=8001 / dt=0.5 ns traces."""
    rows = []
    for tag in ("num_Nd8001", "num_dt0p5"):
        for f in sorted(glob.glob(str(ARCH / "traces" / tag / "*.npz"))):
            meta = parse(f)
            for r in job_metrics(f):
                rows.append({**meta, "check": tag, **r})
    num = pd.DataFrame(rows)
    base = df_all[(df_all.realization == 0) & (df_all.Pavg_over_P1dB_dB == 3)]
    out = []
    for (design, var) in (("dwell", "xi_4"), ("dwell", "xi_0.25"), ("dwell", "xi_1"), ("shuffle", "BLOCK_SHUFFLED_CYCLES")):
        ref = REFS[design]
        for check, g in num.groupby("check"):
            for (method, span), gg in g.groupby(["method", "span_us"]):
                b = base[(base.method == method) & (base.span_us == span)]
                for metric in ("D", "D_ref"):
                    for m in ("atomic", "static_mh", "static_lti_mh"):
                        def val(frame, v):
                            s_ = frame[(frame.design == design) & (frame.variant == v) & (frame.model == m)][metric]
                            return float(s_.iloc[0]) if len(s_) else np.nan
                        d_chk, d_base = val(gg, var) - val(gg, ref), val(b, var) - val(b, ref)
                        out.append({"design": design, "variant": var, "check": check, "method": method, "span_us": span,
                                    "metric": metric, "model": m, "delta_base": d_base, "delta_check": d_chk,
                                    "abs_numerical_change": abs(d_chk - d_base)})
    return pd.DataFrame(out)


def classify(con):
    """Preregistered P1 rules on C1 (dwell xi_4, +3 dB) and C2 (cycle block-shuffle, +3 dB); context lag columns."""
    sign = {"C1": +1, "C2": -1}
    sel = {"C1": ("dwell", "xi_4"), "C2": ("shuffle", "BLOCK_SHUFFLED_CYCLES")}
    arch_exc_sign = {("C1", "static_lti_mh"): +1, ("C1", "static_mh"): +1, ("C2", "static_lti_mh"): -1, ("C2", "static_mh"): -1}
    rows = []
    for metric in ("D", "D_ref"):
        for span in SPANS:
            ok = {}
            for c, (design, var) in sel.items():
                r = con[(con.design == design) & (con.variant == var) & (con.Pavg_over_P1dB_dB == 3) & (con.method == "context")
                        & (con.span_us == span) & (con.metric == metric)].iloc[0]
                a = (np.sign(r.d_atomic) == sign[c]) and (r.d_atomic_lo > 0 if sign[c] > 0 else r.d_atomic_hi < 0)
                b = (np.sign(r.exc_static_lti_mh) == arch_exc_sign[(c, "static_lti_mh")]) and \
                    (r.exc_static_lti_mh_lo > 0 if arch_exc_sign[(c, "static_lti_mh")] > 0 else r.exc_static_lti_mh_hi < 0)
                cc = (np.sign(r.exc_static_mh) == arch_exc_sign[(c, "static_mh")]) and \
                    (r.exc_static_mh_lo > 0 if arch_exc_sign[(c, "static_mh")] > 0 else r.exc_static_mh_hi < 0)
                ok[c] = (bool(a), bool(b), bool(cc))
                rows.append({"metric": metric, "span_us": span, "contrast": c, "full_change": r.d_atomic,
                             "full_ci": f"[{r.d_atomic_lo:+.4f}, {r.d_atomic_hi:+.4f}]",
                             "excess_vs_LTI_static": r.exc_static_lti_mh, "excess_vs_LTI_static_ci": f"[{r.exc_static_lti_mh_lo:+.4f}, {r.exc_static_lti_mh_hi:+.4f}]",
                             "excess_vs_static": r.exc_static_mh, "excess_vs_static_ci": f"[{r.exc_static_mh_lo:+.4f}, {r.exc_static_mh_hi:+.4f}]",
                             "a_full_sign_ci": a, "b_excess_LTI_static": b, "c_excess_static": cc, "all_pass": a and b and cc})
    tab = pd.DataFrame(rows)
    ref_pass = tab[tab.metric == "D_ref"].all_pass.all()
    d2_pass = tab[(tab.metric == "D") & (tab.span_us == 2)].all_pass.all()
    all_pass = tab.all_pass.all()
    verdict = "SURVIVES" if all_pass else ("SURVIVES_WITH_QUALIFICATION" if ref_pass and d2_pass else
                                           "DOES_NOT_SURVIVE" if not ref_pass else "OTHER (D at 2 us fails)")
    return tab, verdict


def posthoc_decomposition(df):
    """POST HOC (not preregistered): split the self-normalized D into its coherent and residual parts.

    Per job and model, relative to the linear model's BLA fit AC power (2-us span, context columns):
      fit_rel = BLA fit AC power (coherent, best-linear part of the output), res_rel = residual power (= D_ref).
    Also the fitted single complex gain relative to the small-signal (linear) model, from run_job (gain_abs).
    Paired contrasts (member minus reference) with 95% CIs for the full model and the all-zone surrogates.
    """
    rows = pd.read_csv(P1 / "rows_p1.csv")
    g = rows.pivot_table(index="job_id", columns="model", values="gain_abs")
    gain_rel = (g.div(g["linear"], axis=0)).stack().rename("gain_rel").reset_index()
    d = df[(df.method == "context") & (df.span_us == 2)].copy()
    lin = d[d.model == "linear"].set_index("job_id").fit_ac_pow
    d["fit_rel"] = d.fit_ac_pow / d.job_id.map(lin)
    d["res_rel"] = d.res_pow / d.job_id.map(lin)
    d = d.merge(gain_rel, on=["job_id", "model"])
    out = []
    for (design, pdb), gg in d.groupby(["design", "Pavg_over_P1dB_dB"]):
        ref = REFS[design]
        for q in ("gain_rel", "fit_rel", "res_rel", "D"):
            piv = gg.pivot_table(index=["realization", "model"], columns="variant", values=q)
            for var in [v for v in piv.columns if v != ref]:
                dd = (piv[var] - piv[ref]).unstack("model").dropna()
                if len(dd) < 2:
                    continue
                row = {"design": design, "Pavg_over_P1dB_dB": pdb, "variant": var, "quantity": q,
                       "ref_mean_full": float(gg[(gg.variant == ref) & (gg.model == "atomic")][q].mean()),
                       "ref_mean_lti": float(gg[(gg.variant == ref) & (gg.model == "static_lti_mh")][q].mean())}
                for m in ("atomic", "static_mh", "static_lti_mh"):
                    row[f"d_{m}"], row[f"d_{m}_lo"], row[f"d_{m}_hi"], row["n"] = ci(dd[m])
                for sur in ("static_mh", "static_lti_mh"):
                    row[f"exc_{sur}"], row[f"exc_{sur}_lo"], row[f"exc_{sur}_hi"], _ = ci(dd["atomic"] - dd[sur])
                out.append(row)
    return pd.DataFrame(out)


def posthoc_gain_vs_dwell():
    """POST HOC: fitted gain (relative to the small-signal model) and fixed-reference residual versus dwell.

    QPSK carrier: the re-run P1 rows (n = 8); unmodulated carrier: the archived carrier-ablation rows (n = 4).
    The preregistered single-slow-state surrogate (tau = 2.535 us) comes from the archived *_dsh rows.
    Output: per (carrier, power, variant, model) the mean and 95% CI of the level and of the paired change vs xi = 0.05.
    """
    srcs = {"qpsk": (P1 / "rows_p1.csv", ARCH / "rows_main_dwell_dsh.csv"),
            "cw": (ARCH / "rows_fu_cwcarrier_dwell.csv", ARCH / "rows_fu_cwcarrier_dwell_dsh.csv")}
    out = []
    for car, (base, dsh) in srcs.items():
        b = pd.read_csv(base)
        b = b[b.design == "dwell"] if "design" in b else b
        d = pd.concat([b, pd.read_csv(dsh)], ignore_index=True)
        d = d[d.model.isin(["atomic", "static_mh", "static_lti_mh", "linear", "dsh_tau2.535us"])]
        lin = d[d.model == "linear"].set_index("job_id").gain_abs
        d = d[d.job_id.isin(lin.index)].copy()
        d["gain_rel"] = d.gain_abs.values / lin.loc[d.job_id].values
        for (pdb, m), g in d[d.model != "linear"].groupby(["Pavg_over_P1dB_dB", "model"]):
            piv = g.pivot_table(index="realization", columns="variant", values="gain_rel")
            if "xi_0.05" not in piv:
                continue
            for var in piv.columns:
                lv = ci(piv[var].dropna())
                dv = ci((piv[var] - piv["xi_0.05"]).dropna()) if var != "xi_0.05" else (0.0, 0.0, 0.0, lv[3])
                out.append({"carrier": car, "Pavg_over_P1dB_dB": pdb, "variant": var, "xi": float(var[3:]), "model": m,
                            "gain_rel_mean": lv[0], "gain_rel_lo": lv[1], "gain_rel_hi": lv[2],
                            "d_gain_rel_mean": dv[0], "d_gain_rel_lo": dv[1], "d_gain_rel_hi": dv[2], "n": lv[3]})
    return pd.DataFrame(out)


def main():
    files = sorted(glob.glob(str(P1 / "traces" / "p1" / "*.npz")))
    from multiprocessing import Pool
    with Pool(12) as pool:
        results = pool.map(job_metrics, files)
    rows = [{**parse(f), **r} for f, rs in zip(files, results) for r in rs]
    df = pd.DataFrame(rows)
    df.to_csv(P1 / "01_p1_metric_rows.csv", index=False)
    # pipeline check: zero-fill 2-us D reproduces the archived D_BLA (rows written by run_job, identical code)
    arch = pd.read_csv(P1 / "rows_p1.csv")[["job_id", "model", "D_BLA", "D_ref_linear"]]
    chk = df[(df.method == "zerofill")].merge(arch, on=["job_id", "model"])
    rel = np.max(np.abs(chk.D - chk.D_BLA) / np.abs(chk.D_BLA))
    rel_ref = np.max(np.abs(chk.D_ref - chk.D_ref_linear) / np.abs(chk.D_ref_linear))
    main_arch = pd.concat([pd.read_csv(ARCH / "rows_main_dwell.csv"), pd.read_csv(ARCH / "rows_main_shuffle.csv")])[["job_id", "model", "D_BLA"]]
    chk2 = df[df.method == "zerofill"].merge(main_arch, on=["job_id", "model"])
    rel_arch = np.max(np.abs(chk2.D - chk2.D_BLA) / np.abs(chk2.D_BLA))
    con = contrasts(df)
    con.to_csv(P1 / "02_p1_contrasts.csv", index=False)
    num = numerical(df)
    num.to_csv(P1 / "03_p1_numerical_differences.csv", index=False)
    tab, verdict = classify(con)
    tab.to_csv(P1 / "04_p1_decisive_table.csv", index=False)
    posthoc_decomposition(df).to_csv(P1 / "06_p1_posthoc_gain_residual_decomposition.csv", index=False)
    posthoc_gain_vs_dwell().to_csv(P1 / "07_p1_posthoc_gain_vs_dwell.csv", index=False)
    summary = {"n_jobs": int(df.job_id.nunique()), "zerofill_2us_vs_rerun_rows_max_rel": rel, "zerofill_2us_D_ref_vs_rerun_rows_max_rel": rel_ref,
               "zerofill_2us_vs_archived_campaign_max_rel": rel_arch, "n_matched_archived": int(len(chk2)), "verdict": verdict}
    (P1 / "05_p1_summary.json").write_text(json.dumps(summary, indent=1, default=float))
    print(json.dumps(summary, indent=1, default=float))
    print(tab.to_string())


if __name__ == "__main__":
    main()
