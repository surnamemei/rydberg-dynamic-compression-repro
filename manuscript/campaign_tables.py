"""Supplementary tables S11-S15 of the final adversarial-validation campaign, rendered from the stored result files.

Used by build_numbers_ledger.py (each table's rendering is a ledger row, so every number in it is traceable) and by the
supplement update, so the two can never disagree. No simulation, no new analysis.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
TVF = ROOT / "results" / "tqe_viability_final"


def _m(x, n=2):
    return f"{x:.{n}f}".replace("-", "−")


def table_s11():
    mr = pd.read_csv(TVF / "01_matching" / "matching_results.csv")
    run = mr[mr.status == "run"]
    cols = [("C1", "5 MHz"), ("C2", "10 MHz"), ("C3", "15 MHz"), ("C4", "5 MHz, weaker probe"), ("C5", "5 MHz, A_LO 0.425"), ("C6", "5 MHz, A_LO 0.35")]
    rows = [("matched_CW_gain", t_, f"CW gain {t_:g}") for t_ in (0.95, 0.9, 0.85, 0.8, 0.6, 0.4)]
    rows += [("matched_signal_to_LO", t_, f"signal/LO {t_:g}" + (" (excluded)" if t_ == 0.6 else "")) for t_ in (0.1, 0.2, 0.358, 0.4, 0.6)]
    rows += [("same_absolute_level", t_, f"a = {t_:g} V/m") for t_ in (0.1, 0.179)]
    out = ["| State | " + " | ".join(c for _, c in cols) + " |", "|---|" + "---:|" * len(cols)]
    for rule, t_, lab in rows:
        cells = []
        for c, _ in cols:
            r = run[(run.configuration == c) & (run.matching_rule == rule) & (abs(run.target - t_) < 1e-9)]
            if len(r):
                r = r.iloc[0]
                cells.append(f"{_m(r.X, 3)} [{_m(r.X_ci_lo, 3)}, {_m(r.X_ci_hi, 3)}]")
            else:
                inf = mr[(mr.configuration == c) & (mr.matching_rule == rule) & (abs(mr.target - t_) < 1e-9)]
                cells.append("infeasible" if len(inf) else "—")
        out.append(f"| {lab} | " + " | ".join(cells) + " |")
    return "\n".join(out)


def table_s12():
    v = json.loads((TVF / "02_resonance" / "resonance_verdict.json").read_text())
    out = ["| A_LO (V/m) | Ω_LO/2π (MHz) | CW-gain minimum: IF (MHz), r₂ | selectivity maximum: IF, r₂ | X_coh maximum: IF, r₂, value |", "|---|---:|---|---|---|"]
    for lo in ("0.35", "0.425", "0.5"):
        f = v["features"][lo]
        rabi = 2 * f["g_min"]["f"] / f["g_min"]["r2"]
        out.append(f"| {lo} | {rabi:.2f} | {f['g_min']['f']:.2f}, {f['g_min']['r2']:.3f} | {f['S_max']['f']:.2f}, {f['S_max']['r2']:.3f} | "
                   f"{f['Xcoh_max']['f']:.2f}, {f['Xcoh_max']['r2']:.3f}, {f['Xcoh_max']['value']:.2f} |")
    out += ["", "| Normalization | features aligned (of 3) | collapse ratio, ln CW gain | collapse ratio, X_coh |", "|---|---:|---:|---:|"]
    names = {"r1": "f/(Ω_LO/2π)", "r2": "f/(Ω_LO/4π)", "d1": "f − Ω_LO/2π", "d2": "f − Ω_LO/4π", "s3": "f/(Ω₃/4π)", "d3": "f − Ω₃/4π"}
    for n_, lab in names.items():
        a = sum(bool(x["aligned"]) for x in v["alignment"][n_].values())
        c = v["collapse_ratio_vs_absolute_f"][n_]
        out.append(f"| {lab} | {a} | {c['ln_g_CW']:.2f} | {c['X_coh']:.2f} |")
    return "\n".join(out)


def table_s13():
    v = json.loads((TVF / "03_memory_model" / "gmp_verdict.json").read_text())
    ms = [("full", "full model"), ("gmp", "GMP"), ("static_mh", "all-zone static"), ("static_lti_mh", "all-zone LTI+static"), ("dsh", "single-slow-state")]
    out = ["| Quantity (test set) | " + " | ".join(l_ for _, l_ in ms) + " |", "|---|" + "---:|" * len(ms)]
    for fam, lab in (("dwell", "median NMSE, dwell"), ("shuffle", "median NMSE, shuffle"), ("phase", "median NMSE, phase"), ("fair", "median NMSE, fair waveforms")):
        out.append(f"| {lab} | — | " + " | ".join(f"{v['nmse'][fam][m]:.3f}" for m, _ in ms[1:]) + " |")
    for k, lab, n in (("dwell_gain", "dwell gain contrast (+3 dB)", 3), ("decluster_Dref", "declustering D_ref contrast (+3 dB)", 4),
                      ("X_REF", "X, reference case", 2), ("X_coh_REF", "X_coh, reference case", 2), ("X_mag_REF", "X_mag, reference case", 2),
                      ("X_FAST", "X, fast rate", 2), ("X_SLOW", "X, slow rate", 2), ("X_JUMP", "X, steps", 2), ("X_LONGRAMP", "X, 900 ns ramps", 2)):
        out.append(f"| {lab} | " + " | ".join(_m(v["contrasts"][k][m], n) for m, _ in ms) + " |")
    for mod, lab in (("QPSK", "QPSK"), ("16QAM", "16-QAM"), ("OFDM", "OFDM"), ("CE", "CE")):
        out.append(f"| mean abs. AIR error at +3 dB, {lab} (bit) | — | " + " | ".join(f"{v['air_err'][mod][m]:.3f}" for m, _ in ms[1:]) + " |")
    return "\n".join(out)


def table_s14():
    v = json.loads((TVF / "04_noise" / "noise_verdict.json").read_text())
    t = pd.read_csv(TVF / "04_noise" / "noise_sensitivity.csv")
    snr = t[(t.table == "effective_SNR_dB") & (t.Pavg_over_P1dB_dB == 3)].groupby("noise_scale").eff_SNR_dB.mean()
    out = ["| Noise scale | effective SNR at +3 dB (dB, mean over formats) | loss at +3 dB: CE / QPSK / 16-QAM / OFDM (bit) | P5: 16-QAM / OFDM (dB) | "
           "surrogate failure (a) and (b) | broad ordering unchanged |", "|---:|---:|---|---|---|---|"]
    for c in ("0.25", "0.5", "1.0", "2.0", "4.0"):
        p = v["per_scale"][c]
        loss = " / ".join(f"{p['loss_+3dB_bit'][m]:.2f}" for m in ("CE", "QPSK", "16QAM", "OFDM"))
        p5 = f"{_m(p['P5_dB']['16QAM'])} / {_m(p['P5_dB']['OFDM'])}"
        out.append(f"| {float(c):g}× | {snr[float(c)]:.1f} | {loss} | {p5} | {'yes' if p['surrogate_failure_holds'] else 'no'} | "
                   f"{'yes' if p.get('broad_ordering_unchanged') else 'no'} |")
    return "\n".join(out)


def table_s15():
    va = pd.read_csv(TVF / "05_numerical" / "va_trace_positivity.csv")
    vb = pd.read_csv(TVF / "05_numerical" / "vb_e1db_slope_sensitivity.csv")
    vc = pd.read_csv(TVF / "05_numerical" / "vc_regime_table.csv")
    out = ["| Check | Result |", "|---|---|"]
    for case, g in va.groupby("case"):
        sci = lambda x: f"{x:.1e}".replace("-", "−", 1) if x < 0 else f"{x:.1e}"  # noqa: E731
        out.append(f"| V-A {case} case, {g.T_us.min():.0f}–{g.T_us.max():.0f} µs: largest per-class trace deviation; smallest eigenvalue | "
                   f"{sci(g.trace_dev_max.max())}; {sci(g.min_eig_min.min())} |")
    names = {"C1": "reference", "C2": "10 MHz", "C3": "15 MHz", "C4": "weaker probe", "C5": "A_LO 0.425", "C6": "A_LO 0.35"}
    for c, g in vb.groupby("config"):
        g = g.set_index("variant")
        out.append(f"| V-B E1dB, {names[c]} (V/m): 2 / 3 / 4 smallest amplitudes / local fit | "
                   f"{g.loc['slope_2_smallest', 'E1dB_Vpm']:.4f} / {g.loc['slope_3_smallest', 'E1dB_Vpm']:.4f} / {g.loc['slope_4_smallest', 'E1dB_Vpm']:.4f} / "
                   f"{g.loc['local_fit_4_smallest', 'E1dB_Vpm']:.4f} |")
    label = {"C2 reference phase case at 2.4468 E1dB": "10 MHz, reference phase case at 2.45 E1dB (archived)",
             "C3 reference phase case at 2.4468 E1dB": "15 MHz, reference phase case at 2.45 E1dB (archived)",
             "C4 reference phase case at 2.4468 E1dB": "weaker probe, reference phase case at 2.45 E1dB (archived)",
             "C2 dwell contrast, +3 dB re E1dB(C2) (peak, r0)": "10 MHz, dwell high level at +3 dB (archived)",
             "C3 dwell contrast, +3 dB re E1dB(C3) (peak, r0)": "15 MHz, dwell high level at +3 dB (archived)",
             "C4 dwell contrast, +3 dB re E1dB(C4) (peak, r0)": "weaker probe, dwell high level at +3 dB (archived)",
             "C4 (weaker probe) amplitude step x1.05 at E1dB (upper level)": "weaker probe, step to 1.05 E1dB (archived)"}
    for _, r in vc[vc.OUTSIDE_INTENDED_REGIME].iterrows():
        frac = "" if pd.isna(r.fraction_of_samples_signal_gt_LO) else f"; {100 * r.fraction_of_samples_signal_gt_LO:.2f}% of samples above the LO"
        lab = label.get(r.state, r.state.replace("shuffle r0-r7, max over variants,", "shuffle family (8 realizations, all variants),"))
        out.append(f"| V-C above the LO: {lab} | peak signal/LO {r.peak_signal_to_LO:.2f}{frac} |")
    out.append(f"| V-C largest peak RF Rabi frequency / carrier, all tabulated states | {vc.Rabi_to_carrier.max():.4f} |")
    return "\n".join(out)


TABLES = {"S11": table_s11, "S12": table_s12, "S13": table_s13, "S14": table_s14, "S15": table_s15}

if __name__ == "__main__":
    for k, f in TABLES.items():
        print(f"== {k}\n{f()}\n")
