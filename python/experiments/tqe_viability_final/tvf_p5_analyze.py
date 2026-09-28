"""Part V: low-cost numerical and regime checks (plan: results/tqe_viability_final/00_plan.md, Part V).

V-A trace/positivity vs duration (prefix runs) + free diagnostics of every stored run; V-B E1dB slope sensitivity (tables only);
V-C RWA / intended-regime table (analysis only). The final NUMERICAL SANITY class also needs the claim matrix (Part VII); this
script writes the inputs to that decision.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from tvf_common import (CONFIGS, F_RF_CARRIER_HZ, FIN, KERNEL_OP1, R, TV, git_state, lo_rabi_hz, md_table, now, s, target_amplitude)

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "revision2"))
from rev2_common import e1db_rule  # noqa: E402

OUT = FIN / "05_numerical"


# ------------------------------------------------------------------ V-A
def va():
    rows = []
    for p in sorted((OUT / "runs").glob("VA_*.npz")):
        z = np.load(p)
        d = json.loads(str(z["diag"]))
        rows.append({"case": "phase" if "phase" in p.stem else "dwell", "T_us": float(z["T"]) * 1e6, **d, "wall_s": float(z["wall_s"])})
    df = pd.DataFrame(rows).sort_values(["case", "T_us"])
    res = {}
    for case, g in df.groupby("case"):
        ok = bool((g.thermal_trace_dev <= 1e-4).all() and (g.thermal_min_eig >= -1e-4).all() and (g.trace_dev_max <= 1e-3).all() and (g.min_eig_min >= -1e-3).all())
        first, last = g.iloc[0], g.iloc[-1]
        acc = bool(last.trace_dev_max > 3 * first.trace_dev_max and last.trace_dev_max > 1e-4)
        res[case] = {"within_thresholds": ok, "accumulation": acc, "trace_dev_max_first_last": [float(first.trace_dev_max), float(last.trace_dev_max)],
                     "min_eig_min_first_last": [float(first.min_eig_min), float(last.min_eig_min)]}
    return df, res


def free_diagnostics():
    """Final-state diagnostics stored with every run of the prior pass and the campaign."""
    groups = {"prior P1": TV / "p1" / "runs", "prior P2": TV / "p2", "prior P4": TV / "p4" / "jobs", "campaign Part I": FIN / "01_matching" / "runs",
              "campaign Part II": FIN / "02_resonance" / "runs"}
    rows = []
    for name, d in groups.items():
        vals = []
        for p in d.rglob("*.npz"):
            z = np.load(p)
            if "diag" in z.files:
                vals.append(json.loads(str(z["diag"])))
        if vals:
            rows.append({"group": name, "n_runs": len(vals), "max_thermal_trace_dev": max(v["thermal_trace_dev"] for v in vals),
                         "min_thermal_eig": min(v["thermal_min_eig"] for v in vals), "max_class_trace_dev": max(v["trace_dev_max"] for v in vals),
                         "min_class_eig": min(v["min_eig_min"] for v in vals), "max_hermiticity_dev": max(v["hermiticity_dev_max"] for v in vals)})
    return pd.DataFrame(rows)


# ------------------------------------------------------------------ V-B
def vb():
    plan = json.loads((FIN / "01_matching" / "plan.json").read_text())
    rows, amp_rows = [], []
    for cfg, (f_if, lo, *_r) in CONFIGS.items():
        z = np.load(CONFIGS[cfg][3])
        amps, harm = np.asarray(z["amps"], float), np.abs(np.asarray(z["harm"]))
        e_arch = float(z["e1db"])
        variants = {}
        for n in (2, 3, 4):
            variants[f"slope_{n}_smallest"] = float(np.mean(harm[1:1 + n] / amps[1:1 + n]))
        A = np.column_stack([np.ones(4), amps[1:5] ** 2])
        coef, *_ = np.linalg.lstsq(A, np.log(harm[1:5] / amps[1:5]), rcond=None)
        variants["local_fit_4_smallest"] = float(np.exp(coef[0]))
        for name, sl in variants.items():
            g = np.full(len(amps), np.nan)
            g[1:] = 20 * np.log10(harm[1:] / (sl * amps[1:]))
            e = e1db_rule(amps, g)
            rows.append({"config": cfg, "variant": name, "slope": sl, "E1dB_Vpm": e, "E1dB_archived_Vpm": e_arch,
                         "change_dB_power": float(20 * np.log10(e / e_arch)) if np.isfinite(e) and np.isfinite(e_arch) else np.nan,
                         "slope_rel_to_archived": sl / variants["slope_2_smallest"]})
            tab = {"amps": amps, "gain_db": g}
            for p in plan["points"]:
                if p["config"] == cfg and p["rule"] == "matched_CW_gain":
                    a = target_amplitude(tab, p["target"])
                    amp_rows.append({"config": cfg, "target": p["target"], "variant": name, "amp_Vpm": a, "amp_planned_Vpm": p["amp_Vpm"],
                                     "change_pct": None if a is None else 100 * (a / p["amp_Vpm"] - 1)})
    return pd.DataFrame(rows), pd.DataFrame(amp_rows)


# ------------------------------------------------------------------ V-C
def rabi_mhz(e_vpm):
    return lo_rabi_hz(1.0) * e_vpm / 1e6


def vc():
    rows = []

    def add(src, desc, lo, a_peak, frac_above=None, used_by=""):
        rows.append({"source": src, "state": desc, "A_LO_Vpm": lo, "peak_signal_Vpm": a_peak, "peak_signal_to_LO": a_peak / lo,
                     "fraction_of_samples_signal_gt_LO": frac_above, "peak_total_Rabi_MHz": rabi_mhz(lo + a_peak),
                     "Rabi_to_carrier": rabi_mhz(lo + a_peak) * 1e6 / F_RF_CARRIER_HZ,
                     "OUTSIDE_INTENDED_REGIME": bool(a_peak / lo > 1), "RWA_QUESTIONABLE": bool(rabi_mhz(lo + a_peak) * 1e6 / F_RF_CARRIER_HZ > .01),
                     "used_by": used_by})
    # archived revision2 configuration checks at 2.45 E1dB(config)
    for cfg, e1, lo in (("C1", 0.0732, 0.5), ("C2", 0.581, 0.5), ("C3", 0.218, 0.5), ("C4", 0.713, 0.5)):
        e1 = float(np.load(CONFIGS[cfg][3])["e1db"])
        add("revision2 (archived)", f"{cfg} reference phase case at 2.4468 E1dB", lo, 2.4468 * e1, None, "manuscript Sec. 7.4-7.5 (C4 claim)")
    for op, lo, e1 in (("OP2", 0.35, float(np.load(CONFIGS["C6"][3])["e1db"])), ("OP3", 0.425, float(np.load(CONFIGS["C5"][3])["e1db"]))):
        add("revision (archived)", f"{op} reference phase case at 2.4468 E1dB", lo, 2.4468 * e1, None, "manuscript Sec. 7.3")
    # campaign states
    plan = json.loads((FIN / "01_matching" / "plan.json").read_text())
    for p in plan["points"]:
        add("campaign Part I", f"{p['config']} {p['rule']} {p['target']:g}", CONFIGS[p["config"]][1], p["amp_Vpm"], None, "Part I")
    for lo in (0.35, 0.425, 0.5):
        add("campaign Part II", f"scan, A_LO {lo:g}, signal/LO 0.358", lo, 0.358 * lo, None, "Part II")
    # fair waveforms (OP1, peak envelope per format and power)
    e1 = float(np.load(KERNEL_OP1)["e1db"])
    for mod in ("CE", "QPSK", "16QAM", "OFDM"):
        seeds = sorted({int(p.stem.split("_s")[1].split("_")[0]) for p in (TV / "p4" / "jobs").glob(f"regen_{mod}_s*.npz")})
        pk = []
        for sd in seeds:
            tx, _, _ = s.waveform(mod, sd)
            pk.append(np.abs(tx))
        mag = np.concatenate(pk)
        for pdb in (-6, 0, 3, 6):
            a = e1 * 10 ** (pdb / 20) * mag
            add("fair waveforms (archived, Sec. 3)", f"{mod} at {pdb:+d} dB re P1dB (peak over {len(seeds)} realizations)", 0.5, float(a.max()),
                float(np.mean(a > 0.5)), "Sections 3 and Part III/IV")
    # revision2 configuration checks of the dwell contrast (xi 4 vs 0.05, +3 dB re each configuration's E1dB) and the step test
    amps_d, *_ = R.realization("dwell", 0)
    mx_d = max(float(np.max(np.abs(amps_d[k]))) for k in ("xi_4", "xi_0.05"))
    for cfg in ("C2", "C3", "C4"):
        e1c = float(np.load(CONFIGS[cfg][3])["e1db"])
        add("revision2 (archived)", f"{cfg} dwell contrast, +3 dB re E1dB({cfg}) (peak, r0)", 0.5, e1c * 10 ** (3 / 20) * mx_d, None,
            "manuscript Sec. 7.4-7.5 (dwell configuration checks)")
    add("revision2 (archived)", "C4 (weaker probe) amplitude step x1.05 at E1dB (upper level)", 0.5, float(np.load(CONFIGS["C4"][3])["e1db"]) * 1.05,
        None, "manuscript Sec. 7.5 (memory persists with the weaker probe)")
    add("revision2 (archived)", "C1 amplitude step x1.05 at E1dB (upper level)", 0.5, float(np.load(CONFIGS["C1"][3])["e1db"]) * 1.05, None,
        "manuscript Sec. 7.5 (standard-probe comparison)")
    # dwell / shuffle (+3 and +6 dB)
    for design in ("dwell", "shuffle"):
        allamps = [R.realization(design, r)[0] for r in range(8)]
        mx = max(float(np.max(np.abs(v))) for amps in allamps for v in amps.values())
        for pdb in (0, 3, 6):
            a = e1 * 10 ** (pdb / 20) * mx
            frac = max(float(np.mean(e1 * 10 ** (pdb / 20) * np.abs(v) > 0.5)) for amps in allamps for v in amps.values())
            add("dwell/shuffle (archived, Sec. 4)", f"{design} r0-r7, max over variants, {pdb:+d} dB", 0.5, a, frac, "Sections 4-6, Part III")
    return pd.DataFrame(rows)


SIM_NOTE = (
    "Simulator inspection (python/utils/transient_quantum.py, unmodified upstream): the RF input is passed as a complex field "
    "envelope E(t) = A_LO + a(t) e^{j2pi f_IF t} in the frame rotating at the LO carrier; Omega_RF(t) = mu_MW E(t)/hbar enters the "
    "four-level ladder Hamiltonian as Omega_RF/2 and conj(Omega_RF)/2 (rotating-wave form; counter-rotating terms at twice the carrier "
    "are absent by construction). The transient solver integrates the full Lindblad equation for each velocity class without "
    "linearizing about the LO, so a signal larger than the LO is represented without formal change. The rotating-wave treatment "
    "requires |Omega_RF(t)| and the envelope's instantaneous frequency to be small compared with the carrier (the 6.946 GHz "
    "transition of the model's source configuration); this holds for every tabulated state (ratio column). What changes when the "
    "signal exceeds the LO is the receiver regime: the superheterodyne picture (a weak signal beating against a dominant LO, and "
    "the small-signal transfer function linearized about the LO steady state) no longer describes the operating point, which is "
    "therefore outside the intended regime of the inherited receiver model. Not quantified here: the four-level truncation "
    "(neighbouring Rydberg transitions), whose importance grows with the total RF Rabi frequency.")


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    va_df, va_res = va()
    fd = free_diagnostics()
    vb_df, vb_amp = vb()
    vc_df = vc()
    va_df.to_csv(OUT / "va_trace_positivity.csv", index=False)
    fd.to_csv(OUT / "free_diagnostics.csv", index=False)
    vb_df.to_csv(OUT / "vb_e1db_slope_sensitivity.csv", index=False)
    vb_amp.to_csv(OUT / "vb_matched_amplitude_sensitivity.csv", index=False)
    vc_df.to_csv(OUT / "vc_regime_table.csv", index=False)
    c1 = vb_df[vb_df.config == "C1"]
    material_c1 = bool((c1.change_dB_power.abs() > 0.5).any())
    material_amp = bool((vb_amp.change_pct.abs() > 5).any())
    free_ok = bool((fd.max_thermal_trace_dev <= 1e-4).all() and (fd.min_thermal_eig >= -1e-4).all() and (fd.max_class_trace_dev <= 1e-3).all()
                   and (fd.min_class_eig >= -1e-3).all())
    va_ok = all(v["within_thresholds"] and not v["accumulation"] for v in va_res.values())
    camp = vc_df[vc_df.source.str.startswith("campaign")]
    s2 = bool((camp.Rabi_to_carrier > .05).any())
    s4 = bool((fd.max_thermal_trace_dev > 1e-2).any() or (fd.min_thermal_eig < -1e-2).any()
              or any(v["trace_dev_max_first_last"][1] > 1e-2 for v in va_res.values()))
    flagged = vc_df[vc_df.OUTSIDE_INTENDED_REGIME | vc_df.RWA_QUESTIONABLE]
    rev, dirty = git_state()
    verdict = {"generated": now(), "code_version": rev, "code_dirty": dirty, "VA": va_res, "VA_pass": va_ok, "free_diagnostics_pass": free_ok,
               "VB_material_C1_E1dB": material_c1, "VB_material_matched_amplitude": material_amp, "VB_max_abs_change_dB": {c: float(g.change_dB_power.abs().max()) for c, g in vb_df.groupby("config")},
               "VC_flagged_states": flagged[["source", "state", "peak_signal_to_LO", "Rabi_to_carrier", "OUTSIDE_INTENDED_REGIME", "RWA_QUESTIONABLE", "used_by"]].to_dict("records"),
               "VC_max_Rabi_to_carrier": float(vc_df.Rabi_to_carrier.max()), "stop_rule_S2": s2, "stop_rule_S4": s4,
               "numerical_checks_pass_excluding_claim_matrix": bool(va_ok and free_ok and not material_c1 and not material_amp),
               "simulator_note": SIM_NOTE}
    (OUT / "numerical_verdict.json").write_text(json.dumps(verdict, indent=1, default=float))
    L = ["# Part V — numerical and regime checks", "", f"Generated {verdict['generated']} (code {rev[:7]}). Rule: `00_plan.md`, Part V.", "",
         "## V-A trace / positivity vs duration (prefix runs)", "", md_table(va_df, ".3g"), "",
         f"Result: {json.dumps(va_res, default=float)}", "", "## Free diagnostics (final states of every stored run)", "", md_table(fd, ".3g"), "",
         "## V-B E1dB slope sensitivity", "", md_table(vb_df, ".4g"), "",
         f"Material change of C1's E1dB (> 0.5 dB): {material_c1}; of a matched I-A amplitude (> 5%): {material_amp} "
         f"(max |change| {vb_amp.change_pct.abs().max():.2f}%).", "",
         "## V-C rotating-wave / intended-regime check", "", SIM_NOTE, "",
         md_table(vc_df[["source", "state", "peak_signal_to_LO", "fraction_of_samples_signal_gt_LO", "peak_total_Rabi_MHz", "Rabi_to_carrier", "OUTSIDE_INTENDED_REGIME",
                "RWA_QUESTIONABLE", "used_by"]], ".3g"), "",
         f"Stop rules: S2 {s2}; S4 {s4}. Numerical checks pass (excluding the claim-matrix condition): {verdict['numerical_checks_pass_excluding_claim_matrix']}."]
    (OUT / "numerical_summary.md").write_text("\n".join(L) + "\n")
    print(json.dumps({k: v for k, v in verdict.items() if k not in ("VC_flagged_states", "simulator_note")}, indent=1, default=float))
    print(flagged[["source", "state", "peak_signal_to_LO", "used_by"]].to_string())


if __name__ == "__main__":
    main()
