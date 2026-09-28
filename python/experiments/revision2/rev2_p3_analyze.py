"""P3 analysis: weak-probe generality check (rules: 00_decision_rules.json P3; 01_rule_clarifications.json)."""
from __future__ import annotations

import json

import numpy as np
import pandas as pd
from scipy import stats

from rev2_common import KERNEL_OP1, OMEGA_P_WEAK, REV2, SEEDS3, s
import rev_phase as RP


def ci(v):
    v = np.asarray(v, float)
    n = len(v)
    h = stats.t.ppf(.975, n - 1) * v.std(ddof=1) / np.sqrt(n) if n > 1 else np.nan
    return float(v.mean()), float(v.mean() - h), float(v.mean() + h), n


def probe_side(runs, dwell):
    gcw = RP.g_final(runs / "CW_H.npz")
    xs = [1 - RP.g_final(runs / f"REF_s{sd}_H.npz") / gcw for sd in SEEDS3]
    d_full, d_st, g_ref = [], [], []
    for r in range(4):
        a, b = np.load(dwell / f"dwell_r{r}_xi_0.05_p+3.npz"), np.load(dwell / f"dwell_r{r}_xi_4_p+3.npz")
        rel = {(m, k): abs(complex(z[f"g_{m}"])) / abs(complex(z["g_linear"])) for k, z in (("ref", a), ("long", b)) for m in ("full", "static")}
        d_full.append(rel[("full", "long")] - rel[("full", "ref")])
        d_st.append(rel[("static", "long")] - rel[("static", "ref")])
        g_ref.append(rel[("full", "ref")])
    return gcw, xs, d_full, d_st, g_ref


def main():
    w_tab = np.load(REV2 / "tables" / "OP1_IF5_weakprobe_ext.npz")      # deviation D3: no -1 dB crossing on the archived grid
    s_tab = np.load(REV2 / "tables" / "OP1_IF5.npz")
    kern = np.load(REV2 / "tables" / "kernel_weakprobe.npz")
    out = {}
    for probe, runs, dwell, tab in (("weak", REV2 / "p3" / "runs", REV2 / "p3" / "dwell", w_tab),
                                    ("standard", REV2 / "p1" / "IF5" / "runs", REV2 / "p1" / "IF5" / "dwell", s_tab)):
        gcw, xs, d_full, d_st, g_ref = probe_side(runs, dwell)
        mx, lo, hi, n = ci(xs)
        dm, dlo, dhi, _ = ci(d_full)
        sm, _, _, _ = ci(d_st)
        steps = {k: json.loads((REV2 / "p3" / "steps" / f"{probe}_{k}.json").read_text()) for k in ("plus5pct_at_1.0", "linear_0.05to0.10")}
        st = steps["plus5pct_at_1.0"]
        out[probe] = {"Omega_p_over_2pi_MHz": (OMEGA_P_WEAK if probe == "weak" else s.old.raqr.Omega_p) / 2 / np.pi / 1e6,
                      "E1dB_Vpm": float(tab["e1db"]), "small_signal_slope": float(tab["slope"]), "cw_settle_rel": float(tab["cw_settle_rel"]),
                      "CW_gain_at_aH_rel_small_signal": gcw, "X_REF_mean": mx, "X_REF_lo": lo, "X_REF_hi": hi, "X_REF_seeds": xs,
                      "gain_rel_full_xi0.05": float(np.mean(g_ref)), "d_gain_full": dm, "d_gain_full_lo": dlo, "d_gain_full_hi": dhi, "d_gain_static": sm,
                      "step_plus5_tau_1e_s": st.get("tau_1e_s"), "step_plus5_tau_90_s": st.get("tau_90_s"),
                      "step_plus5_fit2_tau_slow_s": st.get("fit2_tau_slow_s"), "step_plus5_fit2_a_slow": st.get("fit2_a_slow"),
                      "step_plus5_fit2_tau_fast_s": st.get("fit2_tau_fast_s"), "step_plus5_fit2_a_fast": st.get("fit2_a_fast"),
                      "step_plus5_fit2_rms_resid": st.get("fit2_rms_resid"),
                      "step_linear_tau_1e_s": steps["linear_0.05to0.10"].get("tau_1e_s"), "step_linear_tau_90_s": steps["linear_0.05to0.10"].get("tau_90_s"),
                      "step_linear_fit2_a_slow": steps["linear_0.05to0.10"].get("fit2_a_slow"),
                      "step_linear_fit2_tau_slow_s": steps["linear_0.05to0.10"].get("fit2_tau_slow_s"),
                      "memory_persists": bool(st.get("fit2_tau_slow_s", 0) >= 1e-6 and abs(st.get("fit2_a_slow", 0)) >= .1),
                      "phase_loss_persists": bool(mx >= .05 and lo > 0),
                      "dwell_gain_persists": bool(dm < 0 and dhi < 0 and abs(dm) >= 3 * abs(sm))}
    w, sd = out["weak"], out["standard"]
    ratios = {"X_REF": w["X_REF_mean"] / sd["X_REF_mean"], "d_gain_full": w["d_gain_full"] / sd["d_gain_full"],
              "tau_slow": w["step_plus5_fit2_tau_slow_s"] / sd["step_plus5_fit2_tau_slow_s"],
              "tau_1e (descriptive)": w["step_plus5_tau_1e_s"] / sd["step_plus5_tau_1e_s"]}
    verdict = {"memory_persists": w["memory_persists"], "phase_loss_persists": w["phase_loss_persists"], "dwell_gain_persists": w["dwell_gain_persists"],
               "weak_over_standard": ratios,
               "material_change": bool(any(not (.5 <= v <= 2) for k, v in ratios.items() if "descriptive" not in k)),
               "kernel_lin_check_weak": float(kern["lin_check"])}
    pd.DataFrame(out).to_csv(REV2 / "30_p3_weak_vs_standard.csv")
    (REV2 / "31_p3_verdict.json").write_text(json.dumps({"per_probe": out, **verdict}, indent=1, default=float))
    print(pd.DataFrame(out).to_string())
    print(json.dumps(verdict, indent=1, default=float))


if __name__ == "__main__":
    main()
