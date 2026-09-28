"""P1 analysis: fair cross-configuration matching (rules: results/tqe_viability/00_decision_rules.json, P1)."""
from __future__ import annotations

import json

import numpy as np
import pandas as pd
from scipy import stats

from tv_common import TV, SEEDS3, RP
import rev2_p2 as P2

A_LO = 0.5
REFC = "IF5"


def ci(v):
    v = np.asarray(v, float)
    h = stats.t.ppf(.975, len(v) - 1) * v.std(ddof=1) / np.sqrt(len(v))
    return float(v.mean()), float(v.mean() - h), float(v.mean() + h)


def point_metrics(tag):
    d = TV / "p1" / "runs"
    c = P2.metrics(d / f"{tag}_CW.npz")
    per = []
    for s in SEEDS3:
        q = P2.metrics(d / f"{tag}_REF_s{s}.npz")
        per.append({f"X_{k}": 1 - q[f"g_{k}"] / c[f"g_{k}"] for k in ("med", "coh", "mag", "rms")})
    return c, per


def main():
    plan = json.loads((TV / "p1" / "plan.json").read_text())
    rows, per_seed = [], {}
    for pt in plan["run"]:
        tag = f"{pt['config']}_{pt['rule']}_{pt['target']:g}"
        c, per = point_metrics(tag)
        per_seed[(pt["rule"], pt["target"], pt["config"])] = [p["X_med"] for p in per]
        row = {"configuration": pt["config"], "matching_rule": pt["rule"], "target": pt["target"], "status": "run",
               "signal_amplitude_Vpm": pt["amp_Vpm"], "signal_to_LO": pt["amp_Vpm"] / A_LO, "CW_gain": c["g_med"]}
        for k, name in (("X_med", "X"), ("X_coh", "X_coh"), ("X_mag", "X_mag"), ("X_rms", "X_rms")):
            m, lo, hi = ci([p[k] for p in per])
            row.update({name: m, f"{name}_lo": lo, f"{name}_hi": hi})
        row["effect_present"] = bool(row["X"] >= .05 and row["X_lo"] > 0)
        diag = json.loads(str(np.load(TV / "p1" / "runs" / f"{tag}_CW.npz")["diag"]))
        row["final_state_trace_dev_max"] = diag["trace_dev_max"]
        rows.append(row)
    for pt in plan["infeasible"]:
        a = pt["amp_Vpm"]
        rows.append({"configuration": pt["config"], "matching_rule": pt["rule"], "target": pt["target"],
                     "status": "infeasible: " + ("no crossing on the CW grid" if a is None else f"needs signal/LO = {a / A_LO:.2f} > 1"),
                     "signal_amplitude_Vpm": a, "signal_to_LO": None if a is None else a / A_LO})
    df = pd.DataFrame(rows)
    df.to_csv(TV / "01_matching_comparison.csv", index=False)
    # classification
    states, verdict_rows = {}, []
    for (rule, tgt, cfg), xs in per_seed.items():
        states.setdefault((rule, tgt), {})[cfg] = np.array(xs)
    for (rule, tgt), d in sorted(states.items()):
        if REFC not in d or len(d) < 2:
            continue
        xr = d[REFC]
        mr, lr, hr = ci(xr)
        informative = bool(mr >= .05 and lr > 0)
        comps = {}
        for cfg, xc in d.items():
            if cfg == REFC:
                continue
            dm, dlo, dhi = ci(xr - xc)
            comps[cfg] = {"X_c": float(xc.mean()), "d": dm, "d_lo": dlo, "d_hi": dhi,
                          "dependent": bool(xc.mean() < .5 * mr and dm >= .05 and dlo > 0)}
        verdict_rows.append({"rule": rule, "target": tgt, "X_ref": mr, "X_ref_lo": lr, "X_ref_hi": hr, "informative": informative,
                             "comparisons": comps, "dependence": bool(informative and any(v["dependent"] for v in comps.values()))})
    conv = {}
    for v in verdict_rows:
        if v["informative"]:
            conv.setdefault(v["rule"], []).append(v["dependence"])
    flat = [x for xs in conv.values() for x in xs]
    if flat and all(flat) and all(len(conv.get(r, [])) > 0 for r in conv):
        cls = "CONFIGURATION_DEPENDENCE_SURVIVES"
    elif any(flat):
        cls = "PARTIALLY_SURVIVES"
    else:
        cls = "DOES_NOT_SURVIVE"
    out = {"classification": cls, "informative_states_per_convention": {k: len(v) for k, v in conv.items()},
           "states": verdict_rows}
    (TV / "01_matching_verdict.json").write_text(json.dumps(out, indent=1, default=float))
    pd.set_option("display.width", 250)
    print(df[["configuration", "matching_rule", "target", "status", "signal_amplitude_Vpm", "signal_to_LO", "CW_gain", "X", "X_lo", "X_hi",
              "X_coh", "X_mag", "X_rms"]].round(3).to_string())
    print(json.dumps({k: v for k, v in out.items() if k != "states"}, indent=1))
    for v in verdict_rows:
        print(v["rule"], v["target"], "X_ref", round(v["X_ref"], 3), "informative", v["informative"], "dependence", v["dependence"],
              {k: (round(c["X_c"], 3), round(c["d"], 3), c["dependent"]) for k, c in v["comparisons"].items()})


if __name__ == "__main__":
    main()
