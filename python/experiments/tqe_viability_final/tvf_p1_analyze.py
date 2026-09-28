"""Part I analysis: fair cross-configuration matching (plan: results/tqe_viability_final/00_plan.md, section A).

Outputs (01_matching/): matching_results.csv (one row per configuration and matched state), matching_pairs.csv (paired tests),
matching_verdict.json, matching_summary.md, fig_matching_gain.png, fig_matching_signal_lo.png.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from tvf_common import CONFIGS, FIN, ROOT, SEEDS6, TV, git_state, now

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "revision2"))
import rev2_p2 as P2  # noqa: E402

OUT = FIN / "01_matching"
HA1_CFGS, HA2_CFGS = ("C2", "C3", "C4"), ("C5", "C6")
COLORS = {"C1": ("#2a78d6", "o"), "C2": ("#4a3aa7", "s"), "C3": ("#008300", "D"), "C4": ("#e34948", "^"), "C5": ("#eda100", "v"),
          "C6": ("#e87ba4", "P")}
SHORT = {"C1": "C1 5 MHz (ref.)", "C2": "C2 10 MHz", "C3": "C3 15 MHz", "C4": "C4 weak probe", "C5": "C5 5 MHz, LO 0.425",
         "C6": "C6 5 MHz, LO 0.35"}
METRICS = ("X", "X_coh", "X_mag", "X_rms", "lag_rad", "jitter_rad", "env_red")


def tci(v):
    v = np.asarray(v, float)
    n = len(v)
    h = stats.t.ppf(.975, n - 1) * v.std(ddof=1) / np.sqrt(n)
    return float(v.mean()), float(v.mean() - h), float(v.mean() + h)


def bci(v, n_boot=10_000):
    v = np.asarray(v, float)
    rng = np.random.default_rng(0)
    m = v[rng.integers(0, len(v), (n_boot, len(v)))].mean(axis=1)
    return float(np.percentile(m, 2.5)), float(np.percentile(m, 97.5))


def wrap(x):
    return float((x + np.pi) % (2 * np.pi) - np.pi)


def paths(p):
    d = FIN / "01_matching" / "runs"
    cw = (TV / "p1" / "runs" / f"{p['prior_tag']}_CW.npz") if p["prior_tag"] else d / f"{p['tag']}_CW.npz"
    refs = {}
    for s_ in SEEDS6:
        refs[s_] = (TV / "p1" / "runs" / f"{p['prior_tag']}_REF_s{s_}.npz") if s_ in p["seeds_prior"] else d / f"{p['tag']}_REF_s{s_}.npz"
    return cw, refs


def run_metrics(path):
    z = np.load(path)
    m = P2.window(z["t"])
    out = P2.metrics(path)
    out["mean_abs_full"] = float(np.mean(np.abs(z["z_full"][m])))
    out["median_abs_full"] = float(np.median(np.abs(z["z_full"][m])))
    out["median_abs_lin"] = float(np.median(np.abs(z["z_linear"][m])))
    out["diag"] = json.loads(str(z["diag"]))
    return out


def point(p):
    cw_path, refs = paths(p)
    c = run_metrics(cw_path)
    per = {}
    for s_, rp in refs.items():
        q = run_metrics(rp)
        per[s_] = {"X": 1 - q["g_med"] / c["g_med"], "X_coh": 1 - q["g_coh"] / c["g_coh"], "X_mag": 1 - q["g_mag"] / c["g_mag"],
                   "X_rms": 1 - q["g_rms"] / c["g_rms"], "lag_rad": wrap(q["phase_mean_rad"] - c["phase_mean_rad"]),
                   "jitter_rad": q["phase_circ_sd_rad"], "env_red": 1 - q["mean_abs_full"] / c["mean_abs_full"],
                   "trace": q["diag"]["thermal_trace_dev"], "mineig": q["diag"]["thermal_min_eig"],
                   "trace_cls": q["diag"]["trace_dev_max"], "mineig_cls": q["diag"]["min_eig_min"]}
    return c, per


def main():
    plan = json.loads((OUT / "plan.json").read_text())
    rows, perseed = [], {}
    for p in plan["points"]:
        c, per = point(p)
        key = (p["rule"], p["target"], p["config"])
        perseed[key] = per
        row = {"configuration": p["config"], "description": CONFIGS[p["config"]][5], "matching_rule": p["rule"], "target": p["target"],
               "status": "run", "signal_amplitude_Vpm": p["amp_Vpm"], "signal_to_LO": p["signal_to_LO"], "E_over_E1dB": p["E_over_E1dB"],
               "signal_gt_LO": p["signal_gt_LO"], "table_CW_gain": p["table_CW_gain"], "local_CW_slope": p["local_CW_slope"],
               "pathological": p["pathological"], "CW_gain_rel_sim": c["g_med"],
               "abs_CW_gain_per_Vpm": c["median_abs_full"] / p["amp_Vpm"], "abs_small_signal_gain_per_Vpm": c["median_abs_lin"] / p["amp_Vpm"],
               "n_seq": len(per), "seeds": " ".join(str(s_) for s_ in per), "reused_prior_seeds": " ".join(str(s_) for s_ in p["seeds_prior"])}
        for k in METRICS:
            v = [per[s_][k] for s_ in per]
            m, lo, hi = tci(v)
            blo, bhi = bci(v)
            row.update({k: m, f"{k}_ci_lo": lo, f"{k}_ci_hi": hi, f"{k}_boot_lo": blo, f"{k}_boot_hi": bhi})
        row["effect_present"] = bool(row["X"] >= .05 and row["X_ci_lo"] > 0)
        row["effect_present_Xcoh"] = bool(row["X_coh"] >= .05 and row["X_coh_ci_lo"] > 0)
        row["max_thermal_trace_dev"] = max(max(v["trace"] for v in per.values()), c["diag"]["thermal_trace_dev"])
        row["min_thermal_eig"] = min(min(v["mineig"] for v in per.values()), c["diag"]["thermal_min_eig"])
        row["max_class_trace_dev"] = max(max(v["trace_cls"] for v in per.values()), c["diag"]["trace_dev_max"])
        row["min_class_eig"] = min(min(v["mineig_cls"] for v in per.values()), c["diag"]["min_eig_min"])
        rows.append(row)
    for x in plan["infeasible"]:
        rows.append({"configuration": x["config"], "description": CONFIGS[x["config"]][5], "matching_rule": x["rule"], "target": x["target"],
                     "status": "infeasible: " + x["reason"], "signal_amplitude_Vpm": x.get("amp_Vpm"),
                     "signal_to_LO": (x["amp_Vpm"] / CONFIGS[x["config"]][1]) if x.get("amp_Vpm") else None})
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "matching_results.csv", index=False)

    # ---------------------------------------------------------------- paired tests
    patho = {(r.matching_rule, r.target, r.configuration): bool(r.pathological) for r in df[df.status == "run"].itertuples()}

    def pair_rows(metric):
        out = []
        states = sorted({(k[0], k[1]) for k in perseed})
        for rule, tgt in states:
            if (rule, tgt, "C1") not in perseed:
                continue
            ref = perseed[(rule, tgt, "C1")]
            xr = np.array([ref[s_][metric] for s_ in SEEDS6])
            mr, lr, hr = tci(xr)
            informative = bool(mr >= .05 and lr > 0)
            for cfg in HA1_CFGS + HA2_CFGS:
                if (rule, tgt, cfg) not in perseed:
                    continue
                oth = perseed[(rule, tgt, cfg)]
                xc = np.array([oth[s_][metric] for s_ in SEEDS6])
                d = xr - xc
                dm, dlo, dhi = tci(d)
                blo, bhi = bci(d)
                path_ = patho[(rule, tgt, "C1")] or patho[(rule, tgt, cfg)]
                hyp = "HA1" if cfg in HA1_CFGS else "HA2"
                holds = (bool(xc.mean() < .5 * mr and dm >= .05 and dlo > 0) if hyp == "HA1" else bool(dm >= .05 and dlo > 0))
                convention = {"matched_CW_gain": "A", "matched_signal_to_LO": "B", "same_absolute_level": "C"}[rule]
                out.append({"metric": metric, "convention": convention, "matching_rule": rule, "target": tgt, "config": cfg, "hypothesis": hyp,
                            "X_ref": mr, "X_ref_ci_lo": lr, "X_ref_ci_hi": hr, "ref_informative": informative, "X_c": float(xc.mean()),
                            "diff": dm, "diff_ci_lo": dlo, "diff_ci_hi": dhi, "diff_boot_lo": blo, "diff_boot_hi": bhi,
                            "pathological": path_, "evaluable": bool(informative and not path_ and convention in "AB"), "holds": holds,
                            "reversal": bool(-dm >= .05 and dhi < 0)})
        # same absolute level (convention C): C1 at signal/LO 0.2 (a = 0.1) and 0.358 (a = 0.179) vs C5/C6 at the same amplitude
        for a, ratio in ((0.1, 0.2), (0.179, 0.358)):
            ref = perseed[("matched_signal_to_LO", ratio, "C1")]
            xr = np.array([ref[s_][metric] for s_ in SEEDS6])
            mr, lr, hr = tci(xr)
            for cfg in HA2_CFGS:
                oth = perseed[("same_absolute_level", a, cfg)]
                xc = np.array([oth[s_][metric] for s_ in SEEDS6])
                dm, dlo, dhi = tci(xr - xc)
                out.append({"metric": metric, "convention": "C", "matching_rule": "same_absolute_level", "target": a, "config": cfg,
                            "hypothesis": "HA2 (descriptive)", "X_ref": mr, "X_ref_ci_lo": lr, "X_ref_ci_hi": hr,
                            "ref_informative": bool(mr >= .05 and lr > 0), "X_c": float(xc.mean()), "diff": dm, "diff_ci_lo": dlo, "diff_ci_hi": dhi,
                            "pathological": False, "evaluable": False, "holds": bool(dm >= .05 and dlo > 0), "reversal": bool(-dm >= .05 and dhi < 0)})
        return pd.DataFrame(out)

    def classify(pr):
        ev = pr[pr.evaluable]
        h1 = ev[ev.hypothesis == "HA1"]
        h2 = ev[ev.hypothesis == "HA2"]
        n_a1, n_b1 = int((h1.convention == "A").sum()), int((h1.convention == "B").sum())
        if len(h1) == 0 or not h1.holds.any():
            cls = "DOES NOT SURVIVE"
        elif h1.holds.all() and n_a1 > 0 and n_b1 > 0 and h2.holds.all():
            cls = "SURVIVES"
        else:
            cls = "PARTIALLY SURVIVES"
        return cls, {"HA1_evaluable_pairs": {"A": n_a1, "B": n_b1}, "HA1_pairs_holding": int(h1.holds.sum()),
                     "HA1_pairs_failing": [f"{r.convention}:{r.matching_rule}={r.target:g}:{r.config}" for r in h1[~h1.holds].itertuples()],
                     "HA1_reversals": [f"{r.convention}:{r.matching_rule}={r.target:g}:{r.config}" for r in h1[h1.reversal].itertuples()],
                     "HA2_evaluable_pairs": int(len(h2)), "HA2_pairs_holding": int(h2.holds.sum()),
                     "HA2_pairs_failing": [f"{r.convention}:{r.matching_rule}={r.target:g}:{r.config}" for r in h2[~h2.holds].itertuples()],
                     "HA2_reversals": [f"{r.convention}:{r.matching_rule}={r.target:g}:{r.config}" for r in h2[h2.reversal].itertuples()]}

    pr_x, pr_c = pair_rows("X"), pair_rows("X_coh")
    pd.concat([pr_x, pr_c]).to_csv(OUT / "matching_pairs.csv", index=False)
    cls, det = classify(pr_x)
    cls_coh, det_coh = classify(pr_c)

    # ---------------------------------------------------------------- secondary reading: operating-point claim at matched CW gain
    op = []
    for tgt in (0.8, 0.6, 0.4):
        r = pr_x[(pr_x.matching_rule == "matched_CW_gain") & (pr_x.target == tgt) & (pr_x.config == "C6")]
        if len(r):
            r = r.iloc[0]
            op.append({"target": tgt, "diff_C1_minus_C6": r["diff"], "ci": [r["diff_ci_lo"], r["diff_ci_hi"]],
                       "C1_gt_C6": bool(r["diff"] >= .05 and r["diff_ci_lo"] > 0), "C6_gt_C1": bool(-r["diff"] >= .05 and r["diff_ci_hi"] < 0)})
    n_gt, n_lt = sum(o["C1_gt_C6"] for o in op), sum(o["C6_gt_C1"] for o in op)
    op_claim = "KEEP" if n_gt >= 2 else ("REMOVE" if n_lt >= 2 else "WEAKEN")

    # ---------------------------------------------------------------- prior-pass P1 rule (seeds 1-3, prior states incl. 0.6)
    prior_states = {}
    for (rule, tgt, cfg), per in perseed.items():
        if cfg in ("C1", "C2", "C3", "C4") and all(s_ in per for s_ in SEEDS6[:3]):
            pl = [p for p in plan["points"] if p["config"] == cfg and p["rule"] == rule and p["target"] == tgt][0]
            if pl["prior_tag"]:
                prior_states.setdefault((rule, tgt), {})[cfg] = np.array([per[s_]["X"] for s_ in SEEDS6[:3]])
    conv, prior_rows = {}, []
    for (rule, tgt), d in sorted(prior_states.items()):
        if "C1" not in d or len(d) < 2:
            continue
        mr, lr, _ = tci(d["C1"])
        inf = bool(mr >= .05 and lr > 0)
        dep = False
        for cfg, xc in d.items():
            if cfg == "C1":
                continue
            dm, dlo, _ = tci(d["C1"] - xc)
            dep |= bool(xc.mean() < .5 * mr and dm >= .05 and dlo > 0)
        prior_rows.append({"rule": rule, "target": tgt, "informative": inf, "dependence": bool(inf and dep)})
        if inf:
            conv.setdefault(rule, []).append(bool(dep))
    flat = [x for xs in conv.values() for x in xs]
    prior_cls = ("CONFIGURATION_DEPENDENCE_SURVIVES" if flat and all(flat) else ("PARTIALLY_SURVIVES" if any(flat) else "DOES_NOT_SURVIVE"))

    rev, dirty = git_state()
    verdict = {"generated": now(), "code_version": rev, "code_dirty": dirty, "classification": cls, "details": det,
               "robustness_X_coh": {"classification": cls_coh, "details": det_coh},
               "prior_pass_rule_seeds1to3": {"classification": prior_cls, "states": prior_rows},
               "operating_point_claim_reading": {"decision": op_claim, "per_target": op},
               "stop_rule_S1": cls == "DOES NOT SURVIVE"}
    (OUT / "matching_verdict.json").write_text(json.dumps(verdict, indent=1, default=float))
    figures(df)
    summary(df, pr_x, pr_c, verdict)
    print(json.dumps(verdict, indent=1, default=float))


def figures(df):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    run = df[df.status == "run"]

    def eb(ax, sub, xcol, ycol, cfg, **kw):
        col, mk = COLORS[cfg]
        sub = sub.sort_values(xcol)
        yerr = np.vstack([sub[ycol] - sub[f"{ycol}_ci_lo"], sub[f"{ycol}_ci_hi"] - sub[ycol]])
        ok = ~sub.pathological.astype(bool)
        ax.errorbar(sub[xcol][ok], sub[ycol][ok], yerr=yerr[:, ok.values], color=col, marker=mk, ms=4.5, lw=1.1, elinewidth=0.8, capsize=2,
                    label=SHORT[cfg], **kw)
        if (~ok).any():
            ax.errorbar(sub[xcol][~ok], sub[ycol][~ok], yerr=yerr[:, (~ok).values], color=col, marker=mk, ms=4.5, mfc="white", lw=0,
                        elinewidth=0.8, capsize=2)

    # matched CW gain
    fig, axs = plt.subplots(1, 3, figsize=(11, 3.6))
    a = run[run.matching_rule == "matched_CW_gain"]
    for j, (ycol, lab) in enumerate((("X", "X (median-based)"), ("X_coh", r"$X_\mathrm{coh}$"), ("jitter_rad", "phase jitter (circ. SD, rad)"))):
        ax = axs[j]
        for cfg in ("C1", "C3", "C5", "C6"):
            sub = a[a.configuration == cfg]
            if len(sub):
                eb(ax, sub, "target", ycol, cfg)
        ax.axhline(0, color="#898781", lw=0.6)
        ax.set_xlabel(r"matched CW gain target $g_\mathrm{CW}/g_\mathrm{small}$")
        ax.set_ylabel(lab)
        ax.invert_xaxis()
    axs[0].legend(fontsize=7, frameon=False)
    fig.suptitle("I-A matched CW gain: reference QPSK phase case, 6 sequences, mean and 95% t-CI "
                 "(C2 and C4 cannot be CW-gain matched with signal <= LO)", fontsize=9)
    fig.tight_layout()
    fig.savefig(OUT / "fig_matching_gain.png", dpi=160)
    plt.close(fig)
    # matched signal/LO
    fig, axs = plt.subplots(1, 3, figsize=(11, 3.6))
    b = run[run.matching_rule == "matched_signal_to_LO"]
    for j, (ycol, lab) in enumerate((("X", "X (median-based)"), ("X_coh", r"$X_\mathrm{coh}$"), ("CW_gain_rel_sim", "CW gain / small-signal (simulated)"))):
        ax = axs[j]
        for cfg in CONFIGS:
            sub = b[b.configuration == cfg]
            if not len(sub):
                continue
            if ycol == "CW_gain_rel_sim":
                col, mk = COLORS[cfg]
                sub = sub.sort_values("target")
                ax.plot(sub.target, sub[ycol], color=col, marker=mk, ms=4.5, lw=1.1, label=SHORT[cfg])
            else:
                eb(ax, sub, "target", ycol, cfg)
        ax.axhline(0 if ycol != "CW_gain_rel_sim" else 1, color="#898781", lw=0.6)
        ax.set_xlabel("matched signal/LO ratio")
        ax.set_ylabel(lab)
    axs[0].legend(fontsize=7, frameon=False)
    fig.suptitle("I-B matched signal/LO: reference QPSK phase case, 6 sequences (open markers: pathological state, excluded)", fontsize=9)
    fig.tight_layout()
    fig.savefig(OUT / "fig_matching_signal_lo.png", dpi=160)
    plt.close(fig)


def summary(df, pr_x, pr_c, v):
    run = df[df.status == "run"].copy()
    L = ["# Part I — fair cross-configuration matching", "",
         f"**Classification: {v['classification']}** (rule: `00_plan.md`, section A; generated {v['generated']}, code {v['code_version'][:7]}).", "",
         "Reference QPSK phase case (1 µs symbols, 300 ns ramps, zero drift) at constant amplitude, six sequences (20260701–06) shared by "
         "all configurations; mean with 95% Student-t CI (bootstrap CIs in `matching_results.csv`).", "",
         "## Matched states", "",
         "| Config | Rule | Target | a (V/m) | a/LO | E/E1dB | CW gain (sim) | local slope | X | X_coh | X_mag | lag (rad) | jitter (rad) |",
         "|---|---|---:|---:|---:|---:|---:|---:|---|---|---|---:|---:|"]
    for r in run.sort_values(["matching_rule", "target", "configuration"]).itertuples():
        flag = " (pathological)" if r.pathological else ""
        L.append(f"| {r.configuration} | {r.matching_rule}{flag} | {r.target:g} | {r.signal_amplitude_Vpm:.4f} | {r.signal_to_LO:.3f} | "
                 f"{r.E_over_E1dB:.2f} | {r.CW_gain_rel_sim:.3f} | {r.local_CW_slope:+.2f} | {r.X:+.3f} [{r.X_ci_lo:+.3f}, {r.X_ci_hi:+.3f}] | "
                 f"{r.X_coh:+.3f} [{r.X_coh_ci_lo:+.3f}, {r.X_coh_ci_hi:+.3f}] | {r.X_mag:+.3f} [{r.X_mag_ci_lo:+.3f}, {r.X_mag_ci_hi:+.3f}] | "
                 f"{r.lag_rad:+.2f} | {r.jitter_rad:.2f} |")
    L += ["", "Infeasible (not run, not extrapolated):", ""]
    for r in df[df.status != "run"].itertuples():
        L.append(f"- {r.configuration} {r.matching_rule} {r.target:g}: {r.status}")
    L += ["", "## Paired tests (X; reference C1 minus configuration; same sequences)", "",
          "| Conv. | Rule | Target | Config | Hyp. | X_C1 | X_c | X_C1 − X_c [95% CI] | evaluable | holds |", "|---|---|---:|---|---|---:|---:|---|---|---|"]
    for r in pr_x.itertuples():
        L.append(f"| {r.convention} | {r.matching_rule} | {r.target:g} | {r.config} | {r.hypothesis} | {r.X_ref:+.3f} | {r.X_c:+.3f} | "
                 f"{r.diff:+.3f} [{r.diff_ci_lo:+.3f}, {r.diff_ci_hi:+.3f}] | {'yes' if r.evaluable else 'no'} | {'yes' if r.holds else 'no'} |")
    L += ["", "## Decision", "",
          f"- HA1 (headline; C2–C4): evaluable informative pairs A = {v['details']['HA1_evaluable_pairs']['A']}, B = {v['details']['HA1_evaluable_pairs']['B']}; "
          f"holding {v['details']['HA1_pairs_holding']}; failing {v['details']['HA1_pairs_failing']}; reversals {v['details']['HA1_reversals']}.",
          f"- HA2 (operating point; C5, C6): evaluable pairs {v['details']['HA2_evaluable_pairs']}; holding {v['details']['HA2_pairs_holding']}; "
          f"failing {v['details']['HA2_pairs_failing']}; reversals {v['details']['HA2_reversals']}.",
          f"- **Class: {v['classification']}.** With X_coh in place of X: {v['robustness_X_coh']['classification']}.",
          f"- Prior-pass P1 rule (seeds 1–3, including ratio 0.6): {v['prior_pass_rule_seeds1to3']['classification']}.",
          f"- Operating-point claim reading (C1 vs C6 at matched CW gain 0.8/0.6/0.4): {v['operating_point_claim_reading']['decision']}.",
          f"- Stop rule S1 triggered: {v['stop_rule_S1']}."]
    (OUT / "matching_summary.md").write_text("\n".join(L) + "\n")


if __name__ == "__main__":
    main()
