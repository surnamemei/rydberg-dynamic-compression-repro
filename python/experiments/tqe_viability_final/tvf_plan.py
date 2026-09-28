"""Part I state plan from the existing CW tables only (no simulation): matched amplitudes, feasibility, state descriptors.

Writes results/tqe_viability_final/01_matching/plan.json and prints the table quoted in 00_plan.md.
"""
from __future__ import annotations

import json

from tvf_common import CONFIGS, FIN, PRIOR_NAME, TV, describe, table, target_amplitude

TARGETS_A = (0.95, 0.9, 0.85, 0.8, 0.6, 0.4)          # matched CW gain (0.8/0.6/0.4 preferred; 0.95/0.9/0.85 nearest common feasible)
RATIOS_B = (0.1, 0.2, 0.358, 0.4, 0.6)                # matched signal/LO (0.1/0.2/0.4 requested; 0.358/0.6 from the prior pass)
ABS_C = (0.1, 0.179)                                  # same absolute level as C1 at signal/LO 0.2 and 0.358 (C5, C6 only)


def prior_tag(cfg, rule, target):
    if cfg not in PRIOR_NAME:
        return None
    tag = f"{PRIOR_NAME[cfg]}_{rule}_{target:g}"
    return tag if (TV / "p1" / "runs" / f"{tag}_CW.npz").exists() else None


def main():
    pts, infeasible = [], []
    for cfg, (f_if, lo, *_rest) in CONFIGS.items():
        tab = table(cfg)
        targets = TARGETS_A if cfg in ("C1", "C2", "C3", "C4") else (0.8, 0.6, 0.4)
        for tgt in targets:
            a = target_amplitude(tab, tgt)
            if a is None:
                infeasible.append({"config": cfg, "rule": "matched_CW_gain", "target": tgt, "reason": "no crossing on the CW table grid "
                                   f"(max amplitude {tab['amps'][-1]:.3f} V/m = {tab['amps'][-1] / lo:.2f} x LO)"})
                continue
            if a / lo > 1.0:
                infeasible.append({"config": cfg, "rule": "matched_CW_gain", "target": tgt, "amp_Vpm": a,
                                   "reason": f"needs signal/LO = {a / lo:.2f} > 1 (outside the superheterodyne regime)"})
                continue
            pts.append({"config": cfg, "rule": "matched_CW_gain", "target": tgt, **describe(cfg, a)})
        ratios = RATIOS_B if cfg in ("C1", "C2", "C3", "C4") else (0.1, 0.2, 0.358, 0.4)
        for r in ratios:
            pts.append({"config": cfg, "rule": "matched_signal_to_LO", "target": r, **describe(cfg, r * lo)})
        if cfg in ("C5", "C6"):
            for a in ABS_C:
                pts.append({"config": cfg, "rule": "same_absolute_level", "target": a, **describe(cfg, a)})
    for p in pts:
        p["prior_tag"] = prior_tag(p["config"], p["rule"], p["target"])
        p["seeds_prior"] = [20260701, 20260702, 20260703] if p["prior_tag"] else []
        p["tag"] = f"{p['config']}_{p['rule']}_{p['target']:g}"
    # minimum CW gain reachable inside the regime (documents why no compression target is feasible)
    reach = {}
    for cfg, (f_if, lo, *_rest) in CONFIGS.items():
        tab = table(cfg)
        m = (tab["amps"][1:] <= lo)
        g = 10 ** (tab["gain_db"][1:][m] / 20)
        reach[cfg] = {"min_CW_gain_signal_le_LO": float(g.min()), "max_CW_gain_signal_le_LO": float(g.max()), "E1dB_Vpm": tab["e1db"],
                      "E1dB_over_LO": tab["e1db"] / lo}
    out = {"points": pts, "infeasible": infeasible, "reachable_CW_gain": reach}
    (FIN / "01_matching").mkdir(parents=True, exist_ok=True)
    (FIN / "01_matching" / "plan.json").write_text(json.dumps(out, indent=1, default=float))
    print("| Config | Rule | Target | a (V/m) | a/LO | E/E1dB | CW gain (table) | local slope | pathological | prior 3-seed run |")
    print("|---|---|---:|---:|---:|---:|---:|---:|---|---|")
    for p in pts:
        print(f"| {p['config']} | {p['rule']} | {p['target']:g} | {p['amp_Vpm']:.4f} | {p['signal_to_LO']:.3f} | {p['E_over_E1dB']:.2f} | "
              f"{p['table_CW_gain']:.3f} | {p['local_CW_slope']:+.2f} | {'YES' if p['pathological'] else 'no'} | {'yes' if p['prior_tag'] else 'no'} |")
    print("\nInfeasible:")
    for x in infeasible:
        print(f"- {x['config']} {x['rule']} {x['target']:g}: {x['reason']}")
    print("\nReachable CW gain with signal <= LO:")
    for c, v in reach.items():
        print(f"- {c}: {v['min_CW_gain_signal_le_LO']:.3f} .. {v['max_CW_gain_signal_le_LO']:.3f}; E1dB {v['E1dB_Vpm']:.4f} V/m = {v['E1dB_over_LO']:.2f} x LO")
    print(f"\n{len(pts)} points")


if __name__ == "__main__":
    main()
