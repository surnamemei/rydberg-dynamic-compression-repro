"""Numbers ledger for the manuscript draft.

Re-extracts every headline number quoted in manuscript/main_draft.md directly from the existing
result files (no new simulation, no new analysis beyond reproducing already-reported values).
Output: manuscript/numbers_ledger.csv with id, description, value, text (exact rendering used in
the draft), source.  The consistency audit checks the draft against this file.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python"))
R06 = ROOT / "results" / "stage06_dwell_physics"
PM = R06 / "phase_mechanism"
FV = ROOT / "results" / "final_validation"
PH = ROOT / "manuscript" / "posthoc"
S05 = ROOT / "results" / "p1db_waveform_stage05"
L = []


def add(i, desc, value, text, source):
    L.append({"id": i, "description": desc, "value": value, "text": text, "source": str(Path(source).relative_to(ROOT)) if isinstance(source, Path) else source})


def f(x, n=3):
    return f"{x:.{n}f}"


def pm(x, n=2):
    return f"{x:+.{n}f}".replace("-", "−")


def neg(x, n=3):
    return f"{x:.{n}f}".replace("-", "−")


def ci_paired(x):
    x = np.asarray(x, float)
    h = stats.t.ppf(.975, len(x) - 1) * x.std(ddof=1) / np.sqrt(len(x))
    return x.mean(), x.mean() - h, x.mean() + h


def model_and_numerics():
    from utils.transient_quantum import configure_raqr, get_normal_quadrature
    r = configure_raqr("Transit")
    src = ROOT / "python" / "utils" / "transient_quantum.py"
    add("M_LO", "LO field amplitude (original point)", r.A_LO, "0.5 V/m", src)
    add("M_OmegaP", "probe Rabi frequency /2pi", r.Omega_p / 2 / np.pi, f"{r.Omega_p / 2 / np.pi / 1e6:.2f} MHz", src)
    add("M_OmegaC", "coupling Rabi frequency /2pi", r.Omega_c / 2 / np.pi, f"{r.Omega_c / 2 / np.pi / 1e6:.2f} MHz", src)
    add("M_gamma2", "intermediate-state decay /2pi", r.gamma_2 / 2 / np.pi, f"{r.gamma_2 / 2 / np.pi / 1e6:.1f} MHz", src)
    add("M_gamma3", "Rydberg decay gamma_3 /2pi", r.gamma_3 / 2 / np.pi, f"{r.gamma_3 / 2 / np.pi / 1e3:.1f} kHz", src)
    add("M_gamma4", "Rydberg decay gamma_4 /2pi", r.gamma_4 / 2 / np.pi, f"{r.gamma_4 / 2 / np.pi / 1e3:.1f} kHz", src)
    add("M_transit", "transit relaxation rate (s^-1) and 1/rate", r.gamma, f"{r.gamma / 1e5:.2f}×10⁵ s⁻¹", src)
    add("M_transit_inv", "1/transit rate (us)", 1e6 / r.gamma, f"{1e6 / r.gamma:.1f} µs", src)
    add("M_wp", "probe beam width", r.w_p, f"{r.w_p * 1e3:.2f} mm", src)
    add("M_d", "cell length", r.d, f"{r.d * 100:.0f} cm", src)
    add("M_N0", "atomic density", r.N0, f"{r.N0 / 1e6:.2e} cm⁻³".replace("e+10", "×10¹⁰"), src)
    hbar = 6.626e-34 / (2 * np.pi)
    mu = 1443.45 * 1.6e-19 * 5.2918e-11
    add("M_LORabi", "LO Rabi frequency /2pi at 0.5 V/m", mu * 0.5 / hbar / 2 / np.pi, f"{mu * 0.5 / hbar / 2 / np.pi / 1e6:.2f} MHz", src)
    add("M_LORabi2", "LO Rabi frequency /2pi at 0.35 V/m", mu * 0.35 / hbar / 2 / np.pi, f"{mu * 0.35 / hbar / 2 / np.pi / 1e6:.2f} MHz", src)
    for nd in (1501, 4001, 8001):
        x, _ = get_normal_quadrature(nd)
        add(f"M_Nd{nd}_kept", f"velocity classes retained within ±3 sigma for Nd={nd}", int(np.count_nonzero(np.abs(x) <= 3)), f"{int(np.count_nonzero(np.abs(x) <= 3))}", src)
    sig = float(np.load(ROOT / "artifacts" / "stage1_5_frozen.npz")["noise_sigma"])
    add("M_sigma", "intensity noise sigma (receiver)", sig, f"{sig:.2e}", ROOT / "artifacts" / "stage1_5_frozen.npz")

    e = {nd: float(np.load(R06 / "artifacts" / f"controls_Nd{nd}.npz")["e1db"]) for nd in (1501, 4001, 8001)}
    for nd, v in e.items():
        add(f"N_E1_{nd}", f"CW E1dB at Nd={nd}", v, f"{v:.4f} V/m", R06 / "artifacts" / f"controls_Nd{nd}.npz")
    e_arch = float(np.load(ROOT / "artifacts" / "p1db_reference.npz")["e1db"])
    add("N_E1_arch", "archived Stage-0.5 E1dB (Nd=1501, 40 us CW)", e_arch, f"{e_arch:.4f} V/m", ROOT / "artifacts" / "p1db_reference.npz")
    add("N_shift", "power-axis shift archived->converged (dB)", 20 * np.log10(e_arch / e[4001]), f"{20 * np.log10(e_arch / e[4001]):.2f} dB", ROOT / "artifacts" / "p1db_reference.npz")
    c = {nd: np.load(R06 / "artifacts" / f"controls_Nd{nd}.npz") for nd in (1501, 8001)}
    g = (abs(c[1501]["harm"][1]) / c[1501]["amps"][1]) / (abs(c[8001]["harm"][1]) / c[8001]["amps"][1])
    add("N_ssgain1501", "small-signal CW gain Nd1501 / Nd8001", g, f"{100 * (1 - g):.1f}%", R06 / "artifacts")
    settle = float(np.load(R06 / "artifacts" / "controls_Nd4001.npz")["cw_settle_rel"])
    add("N_cwsettle", "CW-table settling check (rel change between last two 10 us windows)", settle, f"{settle:.1e}", R06 / "artifacts" / "controls_Nd4001.npz")

    nd = pd.read_csv(R06 / "02_nd_convergence.csv")
    ab = nd[nd.case != "E_stage05_dwell_pair"]
    for k in (1501, 3001, 4001):
        s = ab[ab.Nd == k]
        add(f"N_nd{k}_AIR", f"max |AIR - AIR(8001)| at Nd={k}", s.AIR_diff_vs_top.abs().max(), f"{s.AIR_diff_vs_top.abs().max():.4f}", R06 / "02_nd_convergence.csv")
        allk = nd[nd.Nd == k]
        add(f"N_nd{k}_NMSE", f"max baseband NMSE vs Nd=8001 at Nd={k}", allk.baseband_NMSE_vs_top.max(), f"{allk.baseband_NMSE_vs_top.max():.1e}", R06 / "02_nd_convergence.csv")
        add(f"N_nd{k}_RMS", f"max |output RMS rel diff| at Nd={k} (%)", 100 * s.out_rms_rel_diff_vs_top.abs().max(), f"{100 * s.out_rms_rel_diff_vs_top.abs().max():.2g}%", R06 / "02_nd_convergence.csv")
    nv = pd.read_csv(R06 / "07_numerical_validation.csv")
    pw = nv[nv.level == "per_waveform"]
    add("N_val_D_Nd8001", "max |D_BLA| change Nd4001->8001 (per waveform)", pw[pw.check == "num_Nd8001"].D_BLA_abs_diff.max(), f"{pw[pw.check == 'num_Nd8001'].D_BLA_abs_diff.max():.1e}", R06 / "07_numerical_validation.csv")
    add("N_val_D_dt", "max |D_BLA| change dt 1->0.5 ns (per waveform)", pw[pw.check == "num_dt0p5"].D_BLA_abs_diff.max(), f"{pw[pw.check == 'num_dt0p5'].D_BLA_abs_diff.max():.1e}", R06 / "07_numerical_validation.csv")
    pdl = nv[(nv.level == "pair_delta") & (nv.metric == "D_BLA")]
    add("N_val_pair", "max |pair-effect change| over Nd/dt checks", pdl.Delta_abs_diff.max(), f"{pdl.Delta_abs_diff.max():.1e}", R06 / "07_numerical_validation.csv")
    nn = pdl[pdl.Delta_main.abs() > .05]
    add("N_val_margin_min", "min margin ratio (non-null effects)", nn.margin_ratio.min(), f"{nn.margin_ratio.min():.0f}", R06 / "07_numerical_validation.csv")
    add("N_val_margin_max", "max margin ratio (non-null effects)", nn.margin_ratio.max(), f"{nn.margin_ratio.max():.0f}", R06 / "07_numerical_validation.csv")
    cv = pd.read_csv(R06 / "07b_control_cw_validation.csv")
    mh = cv[cv.model.isin(["static_mh", "static_lti_mh"])].rel_rms_err_vs_full
    f1 = cv[cv.model.isin(["static_F1", "static_lti_F1"])].rel_rms_err_vs_full
    add("N_ctl_mh", "all-zone controls: CW rel. error range", (mh.min(), mh.max()), f"{100 * mh.min():.2f}–{100 * mh.max():.1f}%", R06 / "07b_control_cw_validation.csv")
    add("N_ctl_f1", "fundamental-only controls: CW rel. error range", (f1.min(), f1.max()), f"{100 * f1.min():.0f}–{100 * f1.max():.0f}%", R06 / "07b_control_cw_validation.csv")
    M = np.load(R06 / "artifacts" / "controls_MH_Nd4001.npz")
    amps, C = M["amps"], M["C"]
    E1 = e[4001]
    sel = (amps >= .75 * E1) & (amps <= 3.2 * E1)
    ratio = np.abs(C[sel, 0]) / np.abs(C[sel, 1])
    near = np.argmin(np.abs(amps - E1))
    add("N_c0c1_nearE1", "|c0/c1| at the grid amplitude nearest E1dB", abs(C[near, 0]) / abs(C[near, 1]), f"{abs(C[near, 0]) / abs(C[near, 1]):.2f}", R06 / "artifacts" / "controls_MH_Nd4001.npz")
    add("N_c0c1_max", "max |c0/c1| for 0.75-3.2 E1dB", ratio.max(), f"{ratio.max():.1f}", R06 / "artifacts" / "controls_MH_Nd4001.npz")
    osc = pd.read_csv(PH / "posthoc_oscillation_vs_Nd.csv")
    on = osc[osc.experiment == "IF_plateau_0.5to1.5E1_on"].osc_period_us
    add("N_osc_on", "plateau-onset oscillation period at Nd 1501/4001/8001 (us)", tuple(on), f"{on.min():.3f} µs", PH / "posthoc_oscillation_vs_Nd.csv")


def regen():
    Lo = pd.read_csv(FV / "regen_stage05" / "02_losses.csv")
    for r in Lo[Lo.Pavg_over_P1dB_dB.isin([0, 3, 6])].itertuples():
        add(f"R_loss_{r.modulation}_{r.Pavg_over_P1dB_dB:+d}", f"full-model AIR loss vs -6 dB, {r.modulation}, {r.Pavg_over_P1dB_dB:+d} dB", r.loss_mean,
            f"{neg(r.loss_mean, 2)} [{neg(r.loss_ci_lo, 2)}, {neg(r.loss_ci_hi, 2)}]", FV / "regen_stage05" / "02_losses.csv")
    E = pd.read_csv(FV / "regen_stage05" / "03_excess_over_controls.csv")
    for r in E[E.Pavg_over_P1dB_dB.isin([0, 3, 6]) & E.control.isin(["static_mh", "static_lti_mh"])].itertuples():
        add(f"R_exc_{r.modulation}_{r.Pavg_over_P1dB_dB:+d}_{r.control}", f"AIR(control)-AIR(full) {r.control} {r.modulation} {r.Pavg_over_P1dB_dB:+d} dB", r.excess_mean,
            f"{pm(r.excess_mean)} [{pm(r.ci_lo)}, {pm(r.ci_hi)}]", FV / "regen_stage05" / "03_excess_over_controls.csv")
    T = pd.read_csv(FV / "regen_stage05" / "04_thresholds.csv")
    for r in T.itertuples():
        add(f"R_thr_{r.modulation}_{r.limit_pct}", f"P{r.limit_pct} - P1dB {r.modulation}", r.P_minus_P1dB_dB, f"{neg(r.P_minus_P1dB_dB, 2)} [{neg(r.boot_lo, 2)}, {neg(r.boot_hi, 2)}]", FV / "regen_stage05" / "04_thresholds.csv")
    rows = pd.read_csv(FV / "regen_stage05" / "rows.csv").drop_duplicates(["job_id", "model"])
    add("R_nruns", "number of regenerated full-atomic runs", rows.job_id.nunique(), f"{rows.job_id.nunique()}", FV / "regen_stage05" / "rows.csv")
    W = pd.read_csv(S05 / "03_waveform_statistics_extended.csv")
    W = W[W.config_version.isin(["fair_v2", "fair_v4_randomized_balanced_ofdm_training"]) & ~((W.modulation == "OFDM") & (W.config_version == "fair_v2"))]
    b = W.groupby("modulation").B_occ_99_Hz.mean() / 1e6
    add("R_B99", "mean 99% bandwidth per format (MHz)", dict(b.round(3)), f"{b.min():.2f}–{b.max():.2f} MHz", S05 / "03_waveform_statistics_extended.csv")
    add("R_B99_mis", "max relative B99 mismatch across formats (%)", 100 * (b.max() - b.min()) / b.mean(), f"{100 * (b.max() - b.min()) / b.mean():.2f}%", S05 / "03_waveform_statistics_extended.csv")


def dwell_and_shuffle():
    cfg = pd.read_csv(R06 / "04_controlled_pair_configs.csv")
    d = cfg[cfg.design == "dwell"]
    add("D_sorted_identical", "sorted samples identical in every realization (dwell + shuffle)", bool(cfg.sorted_samples_identical.all()), "identical", R06 / "04_controlled_pair_configs.csv")
    add("D_B99_mis", "dwell family B99 mismatch range (%)", (d.B99_complex_env_Hz_mismatch_pct.min(), d.B99_complex_env_Hz_mismatch_pct.max()),
        f"{neg(d.groupby('variant').B99_complex_env_Hz_mismatch_pct.mean().min(), 1)}% to +{d.groupby('variant').B99_complex_env_Hz_mismatch_pct.mean().max():.1f}%", R06 / "04_controlled_pair_configs.csv")
    ib = d.groupby("variant").inband_pm3MHz_power_fraction_mismatch_pct.mean()
    add("D_inband", "dwell in-band (±3 MHz) power fraction mismatch (%)", (ib[ib > 0].min(), ib.max()), f"+{ib[ib > 0].min():.1f}% to +{ib.max():.1f}%", R06 / "04_controlled_pair_configs.csv")
    sp = d.groupby("variant").meta_spike_share_of_high_time.mean()
    add("D_spikes", "spike share of high-field time (mean over realizations)", (sp[sp > 0].min(), sp.max()), f"{100 * sp[sp > 0].min():.1f}–{100 * sp.max():.1f}%", R06 / "04_controlled_pair_configs.csv")
    tw = d.groupby("variant").T_dwell_over_tau.mean()
    add("D_Tdw", "time-weighted T_dwell/tau range", (tw.min(), tw.max()), f"{tw.min():.2f}–{tw.max():.1f}", R06 / "04_controlled_pair_configs.csv")
    s = cfg[(cfg.design == "shuffle") & (cfg.variant == "BLOCK_SHUFFLED_CYCLES")]
    add("D_cyc_B99", "cycle block-shuffle B99 mismatch (%)", s.B99_complex_env_Hz_mismatch_pct.mean(), f"{abs(s.B99_complex_env_Hz_mismatch_pct.mean()):.2f}%", R06 / "04_controlled_pair_configs.csv")
    P = pd.read_csv(R06 / "05_controlled_pair_results.csv")
    DB = P[P.metric == "D_BLA"]
    for p in (-6, 0, 3, 6):
        for xi in ("xi_0.1", "xi_1", "xi_4"):
            r = DB[(DB.Pavg_over_P1dB_dB == p) & (DB.variant == xi) & (DB.model == "atomic")].iloc[0]
            add(f"D_full_{p:+d}_{xi}", f"dwell dD_BLA full, {p:+d} dB, {xi}", r.Delta_mean, f"{pm(r.Delta_mean, 3)} [{pm(r.Delta_ci_lo, 3)}, {pm(r.Delta_ci_hi, 3)}]", R06 / "05_controlled_pair_results.csv")
            if p in (0, 3):
                add(f"D_exc_{p:+d}_{xi}", f"dwell full - LTI_MH, {p:+d} dB, {xi}", r.full_minus_static_lti_mh_mean,
                    f"{pm(r.full_minus_static_lti_mh_mean, 3)} [{pm(r.full_minus_static_lti_mh_ci_lo, 3)}, {pm(r.full_minus_static_lti_mh_ci_hi, 3)}]", R06 / "05_controlled_pair_results.csv")
        c = DB[(DB.Pavg_over_P1dB_dB == p) & (DB.variant == "xi_4") & (DB.model == "static_lti_mh")].iloc[0]
        add(f"D_lti_{p:+d}_xi4", f"dwell dD_BLA LTI_MH at xi=4, {p:+d} dB", c.Delta_mean, f"{pm(c.Delta_mean, 3)}", R06 / "05_controlled_pair_results.csv")
        if p in (0, 3):
            c = DB[(DB.Pavg_over_P1dB_dB == p) & (DB.variant == "xi_4") & (DB.model == "static_mh")].iloc[0]
            add(f"D_smh_{p:+d}_xi4", f"dwell dD_BLA all-zone static at xi=4, {p:+d} dB", c.Delta_mean, f"{pm(c.Delta_mean, 2)}", R06 / "05_controlled_pair_results.csv")
            r = DB[(DB.Pavg_over_P1dB_dB == p) & (DB.variant == "xi_4") & (DB.model == "atomic")].iloc[0]
            add(f"D_excsmh_{p:+d}_xi_4", f"dwell full - all-zone static, {p:+d} dB, xi_4", r.full_minus_static_mh_mean,
                f"{pm(r.full_minus_static_mh_mean, 2)} [{pm(r.full_minus_static_mh_ci_lo, 2)}, {pm(r.full_minus_static_mh_ci_hi, 2)}]", R06 / "05_controlled_pair_results.csv")
    lin = DB[DB.model == "linear"].Delta_mean.abs().max()
    add("D_lin_max", "max |dD_BLA| of the linear model (dwell)", lin, f"{lin:.3f}", R06 / "05_controlled_pair_results.csv")
    la = P[(P.metric == "AIR") & (P.model == "linear")].Delta_mean
    add("D_lin_AIR", "linear-model dAIR range across dwell pairs (bit)", (la.min(), la.max()), f"{la.min():.2f}–{la.max():.2f} bit", R06 / "05_controlled_pair_results.csv")
    S = pd.read_csv(R06 / "06_temporal_shuffle_results.csv")
    SB = S[(S.metric == "D_BLA") & (S.variant == "BLOCK_SHUFFLED_CYCLES")]
    for p in (-6, 0, 3, 6):
        r = SB[(SB.Pavg_over_P1dB_dB == p) & (SB.model == "atomic")].iloc[0]
        add(f"S_full_{p:+d}", f"cycle-shuffle dD full {p:+d} dB", r.Delta_mean, f"{pm(r.Delta_mean, 3)} [{pm(r.Delta_ci_lo, 3)}, {pm(r.Delta_ci_hi, 3)}]", R06 / "06_temporal_shuffle_results.csv")
        add(f"S_exc_{p:+d}", f"cycle-shuffle full - LTI_MH {p:+d} dB", r.full_minus_static_lti_mh_mean,
            f"{pm(r.full_minus_static_lti_mh_mean, 3)} [{pm(r.full_minus_static_lti_mh_ci_lo, 3)}, {pm(r.full_minus_static_lti_mh_ci_hi, 3)}]", R06 / "06_temporal_shuffle_results.csv")
        m = SB[(SB.Pavg_over_P1dB_dB == p) & (SB.model == "static_mh")].iloc[0]
        add(f"S_static_{p:+d}", f"cycle-shuffle dD all-zone static {p:+d} dB", m.Delta_mean, f"{pm(m.Delta_mean, 3)} [{pm(m.Delta_ci_lo, 3)}, {pm(m.Delta_ci_hi, 3)}]", R06 / "06_temporal_shuffle_results.csv")
    W = pd.read_csv(R06 / "06b_shuffle_per_waveform_excess.csv")
    for p in (3, 6):
        o = W[(W.Pavg_over_P1dB_dB == p) & (W.variant == "ORIGINAL")].D_BLA_atomic_mean.iloc[0]
        cc = W[(W.Pavg_over_P1dB_dB == p) & (W.variant == "BLOCK_SHUFFLED_CYCLES")].D_BLA_atomic_mean.iloc[0]
        add(f"S_relchg_{p:+d}", f"cycle-shuffle relative change of mean full D, {p:+d} dB", cc / o - 1, f"{100 * (1 - cc / o):.0f}%", R06 / "06b_shuffle_per_waveform_excess.csv")
    rows = pd.read_csv(R06 / "rows_main_dwell.csv").drop_duplicates(["job_id", "model"], keep="last")
    lin = rows[rows.model == "linear"].set_index(["realization", "variant", "Pavg_over_P1dB_dB"]).gain_abs
    for mdl in ("atomic", "static_lti_mh"):
        r3 = rows[(rows.model == mdl) & (rows.Pavg_over_P1dB_dB == 3)].copy()
        r3["g"] = r3.gain_abs.values / lin.reindex(pd.MultiIndex.from_frame(r3[["realization", "variant", "Pavg_over_P1dB_dB"]])).values
        gm = r3.groupby("variant").g.mean()
        add(f"D_gain_{mdl}_+3", f"fitted gain rel. linear, {mdl}, +3 dB, xi 0.05 and 4", (gm["xi_0.05"], gm["xi_4"]), f"{gm['xi_0.05']:.2f} → {gm['xi_4']:.2f}", R06 / "rows_main_dwell.csv")


def timescale():
    t = pd.read_csv(R06 / "03_tau_atom_results.csv")
    t4 = t[t.Nd == 4001].set_index("experiment")
    add("T_LO_1e", "LO step tau_1e", t4.loc["LO_step_+1pct", "tau_1e_s"], f"{t4.loc['LO_step_+1pct', 'tau_1e_s'] * 1e6:.2f} µs", R06 / "03_tau_atom_results.csv")
    add("T_LO_dom", "LO step tail tau", t4.loc["LO_step_+1pct", "tau_dom_s"], f"{t4.loc['LO_step_+1pct', 'tau_dom_s'] * 1e6:.2f} µs", R06 / "03_tau_atom_results.csv")
    add("T_lin_1e", "linear-point IF step tau_1e", t4.loc["IF_step_small_0.05to0.10E1", "tau_1e_s"], f"{t4.loc['IF_step_small_0.05to0.10E1', 'tau_1e_s'] * 1e9:.0f} ns", R06 / "03_tau_atom_results.csv")
    add("T_lin_90", "linear-point IF step tau_90", t4.loc["IF_step_small_0.05to0.10E1", "tau_90_s"], f"{t4.loc['IF_step_small_0.05to0.10E1', 'tau_90_s'] * 1e9:.0f} ns", R06 / "03_tau_atom_results.csv")
    add("T_lin_slowweight", "linear-point IF step slow-component weight", abs(t4.loc["IF_step_small_0.05to0.10E1", "fit2_a_slow"]), f"{100 * abs(t4.loc['IF_step_small_0.05to0.10E1', 'fit2_a_slow']):.0f}%", R06 / "03_tau_atom_results.csv")
    r = t4.loc["IF_step_P1dB_1.00to1.05E1"]
    add("T_atom", "tau_atom (tau_1e, IF step at 0.0877 V/m)", r.tau_1e_s, f"{r.tau_1e_s * 1e6:.3f} µs", R06 / "03_tau_atom_results.csv")
    add("T_atom_90", "tau_90 same step", r.tau_90_s, f"{r.tau_90_s * 1e6:.2f} µs", R06 / "03_tau_atom_results.csv")
    add("T_rec_90", "post-plateau recovery tau_90 (Step 2)", t4.loc["IF_plateau_0.5to1.5E1_off_recovery", "tau_90_s"], f"{t4.loc['IF_plateau_0.5to1.5E1_off_recovery', 'tau_90_s'] * 1e6:.1f} µs", R06 / "03_tau_atom_results.csv")
    stable = t[t.experiment == "IF_step_P1dB_1.00to1.05E1"].tau_1e_s
    add("T_atom_Nd", "tau_atom across Nd 1501/4001/8001 (us)", tuple(stable * 1e6), f"{stable.min() * 1e6:.3f}–{stable.max() * 1e6:.3f} µs", R06 / "03_tau_atom_results.csv")
    s = pd.read_csv(R06 / "03b_tau_atom_sensitivity.csv").set_index("experiment")
    add("T_slowpole", "slow-pole constants at 0.8/1.0/1.2 E1dB", tuple(s.fit2_tau_slow_s * 1e6), " / ".join(f"{v * 1e6:.2f}" for v in s.fit2_tau_slow_s) + " µs", R06 / "03b_tau_atom_sensitivity.csv")
    add("T_1e_sens", "tau_1e at 0.8/1.0/1.2 E1dB", tuple(s.tau_1e_s * 1e6), " / ".join(f"{v * 1e6:.2f}" for v in s.tau_1e_s) + " µs", R06 / "03b_tau_atom_sensitivity.csv")
    add("T_slowweight", "slow-component weights at 0.8/1.0/1.2 E1dB", tuple(abs(s.fit2_a_slow)), " / ".join(f"{abs(v):.2f}" for v in s.fit2_a_slow), R06 / "03b_tau_atom_sensitivity.csv")
    pl = pd.read_csv(R06 / "12_plateau_steady_state.csv").set_index("experiment")
    for p in ("+0", "+3"):
        r = pl.loc[f"plateau_p{p}_cw"]
        add(f"P_onset_{p}", f"CW plateau gain first us / final ({p} dB level)", (r.gain_plateau_first_us, r.gain_plateau_last10us), f"{r.gain_plateau_first_us:.2f} → {r.gain_plateau_last10us:.2f}", R06 / "12_plateau_steady_state.csv")
        add(f"P_rec_{p}", f"post-plateau gain first us / 1-5 us / final ({p} dB level)", (r.gain_low_first_us_after, r.gain_low_1_to_5us_after, r.gain_low_end),
            f"{r.gain_low_first_us_after:.2f}, {r.gain_low_1_to_5us_after:.2f}, {r.gain_low_end:.2f}", R06 / "12_plateau_steady_state.csv")
        add(f"P_QS_{p}", f"CW plateau vs CW-static ({p} dB level)", r.plateau_vs_cw_static_rel, f"{pm(100 * r.plateau_vs_cw_static_rel)}%", R06 / "12_plateau_steady_state.csv")
        q = pl.loc[f"plateau_p{p}_qpsk"]
        add(f"P_carrier_{p}", f"QPSK plateau vs CW-static ({p} dB level)", q.plateau_vs_cw_static_rel, f"{neg(100 * q.plateau_vs_cw_static_rel, 1)}%", R06 / "12_plateau_steady_state.csv")
    ds = pd.read_csv(R06 / "16_followup_single_slow_state.csv")
    for p in (3, 6):
        r = ds[(ds.set.str.startswith("shuffle")) & (ds.Pavg_over_P1dB_dB == p)].iloc[0]
        add(f"E_dsh_clu_{p:+d}", f"slow-state vs full, clustering contrast {p:+d} dB", (r.Delta_DSH, r.Delta_full), f"{neg(r.Delta_DSH)} vs {neg(r.Delta_full)}", R06 / "16_followup_single_slow_state.csv")
    dw = ds[ds.contrast.str.startswith("xi")]
    add("E_dsh_dwell", "dwell contrasts within 25% (slow-state model)", (int(dw.within_25pct.sum()), len(dw)), f"{int(dw.within_25pct.sum())} of {len(dw)}", R06 / "16_followup_single_slow_state.csv")
    nm = pd.read_csv(R06 / "16b_followup_dsh_trace_nmse.csv")
    fr = nm.groupby("set").DSH_better.mean()
    add("E_dsh_trace", "fraction of realization-0 waveforms where slow-state beats LTI_MH (by set)", dict(fr.round(3)), f"{100 * fr.min():.0f}–{100 * fr.max():.0f}%", R06 / "16b_followup_dsh_trace_nmse.csv")
    add("E_dsh_trace_sets", "sets meeting the pre-registered >=75% rule", (int((fr >= .75).sum()), len(fr)), f"{int((fr >= .75).sum())} of {len(fr)}", R06 / "16b_followup_dsh_trace_nmse.csv")
    med = nm.groupby("set")[["NMSE_vs_full_trace", "NMSE_vs_full_trace_static_lti_mh"]].median()
    for st, key in (("fu_cwcarrier/dwell", "cwcarrier"), ("main/dwell", "main_dwell"), ("fu_long_qpsk/dwell", "long_qpsk")):
        add(f"E_dsh_nmse_{key}", f"median trace NMSE slow-state vs LTI_MH, {st}", tuple(med.loc[st]), f"{med.loc[st].iloc[0]:.3f} vs {med.loc[st].iloc[1]:.3f}", R06 / "16b_followup_dsh_trace_nmse.csv")


def quasi_static():
    q = pd.read_csv(R06 / "15_followup_long_dwell_per_waveform.csv")
    b = q[(q.metric == "D_BLA") & (q.carrier == "cw") & (q.model == "full_minus_static_mh")].set_index("xi")["mean"]
    add("Q_DBLA_cw", "pre-registered D_BLA full-static gap, CW carrier, xi 4 and 32", (b[4.0], b[32.0]), f"{neg(b[4.0])} → {neg(b[32.0])}", R06 / "15_followup_long_dwell_per_waveform.csv")
    add("Q_DBLA_frac", "remaining fraction of the xi=4 gap at xi=32 (D_BLA)", b[32.0] / b[4.0], f"{100 * b[32.0] / b[4.0]:.0f}%", R06 / "15_followup_long_dwell_per_waveform.csv")
    st = q[(q.metric == "D_BLA") & (q.carrier == "cw") & (q.model == "static_mh") & (q.xi == .05)]["mean"].iloc[0]
    add("Q_DBLA_patho", "static-model D_BLA at xi=0.05, CW carrier", st, f"{st:.0f}", R06 / "15_followup_long_dwell_per_waveform.csv")
    for car in ("cw", "qpsk"):
        rows = pd.read_csv(R06 / f"rows_fu_long_{car}_dwell.csv").drop_duplicates(["job_id", "model"], keep="last")
        w = rows.pivot_table(index=["variant", "realization"], columns="model", values="D_ref_linear")
        linr = rows[rows.model == "linear"].set_index(["variant", "realization"]).gain_abs
        g0 = w.loc["xi_0.05"]
        add(f"Q_ref_{car}_xi_0.05", f"D_ref full vs static, {car}, xi_0.05", (g0.atomic.mean(), g0.static_mh.mean()), f"{g0.atomic.mean():.3f} vs {g0.static_mh.mean():.2f}", R06 / f"rows_fu_long_{car}_dwell.csv")
        for xi in ("xi_4", "xi_32"):
            g = w.loc[xi]
            m, lo, hi = ci_paired(g.atomic - g.static_mh)
            base = g.static_mh.mean()
            add(f"Q_ref_{car}_{xi}", f"D_ref full vs static, {car}, {xi}", (g.atomic.mean(), base), f"{g.atomic.mean():.4f} vs {base:.4f}", R06 / f"rows_fu_long_{car}_dwell.csv")
            add(f"Q_refgap_{car}_{xi}", f"relative full-static D_ref gap, {car}, {xi}", m / base, f"+{100 * m / base:.1f}% [+{100 * lo / base:.1f}%, +{100 * hi / base:.1f}%]" if xi == "xi_32" and car == "cw" else f"+{100 * m / base:.0f}%",
                R06 / f"rows_fu_long_{car}_dwell.csv")
        if car == "cw":
            add("Q_ref_lti_cw_xi_32", "D_ref all-zone LTI+static, cw, xi_32", w.loc["xi_32"].static_lti_mh.mean(), f"{w.loc['xi_32'].static_lti_mh.mean():.4f}", R06 / "rows_fu_long_cw_dwell.csv")
        for xi in ("xi_0.05", "xi_32"):
            ga = rows[(rows.model == "atomic") & (rows.variant == xi)].set_index("realization").gain_abs / linr.loc[xi]
            gs = rows[(rows.model == "static_mh") & (rows.variant == xi)].set_index("realization").gain_abs / linr.loc[xi]
            add(f"Q_gain_{car}_{xi}", f"fitted gain rel. linear full vs static, {car}, {xi}", (ga.mean(), gs.mean()), f"{ga.mean():.2f} vs {gs.mean():.2f}", R06 / f"rows_fu_long_{car}_dwell.csv")


def carrier_ablation():
    c = pd.read_csv(R06 / "14_followup_carrier_ablation_pairs.csv")
    cw = c[(c.carrier == "cw") & (c.metric == "D_BLA") & (c.Pavg_over_P1dB_dB == 3) & (c.xi == 4)].set_index("model")
    add("C_fuc_exc", "FU-C: full - LTI_MH at +3 dB, xi=4, unmodulated carrier (pre-registered test)", cw.loc["full_minus_static_lti_mh", "Delta_mean"],
        f"{pm(cw.loc['full_minus_static_lti_mh', 'Delta_mean'], 3)} [{pm(cw.loc['full_minus_static_lti_mh', 'Delta_ci_lo'], 3)}, {pm(cw.loc['full_minus_static_lti_mh', 'Delta_ci_hi'], 3)}]",
        R06 / "14_followup_carrier_ablation_pairs.csv")
    add("C_fuc_full", "FU-C: full-model dD at +3 dB, xi=4, unmodulated carrier", cw.loc["atomic", "Delta_mean"],
        f"{pm(cw.loc['atomic', 'Delta_mean'], 3)} [{pm(cw.loc['atomic', 'Delta_ci_lo'], 3)}, {pm(cw.loc['atomic', 'Delta_ci_hi'], 3)}]", R06 / "14_followup_carrier_ablation_pairs.csv")
    add("C_fuc_lti", "FU-C: LTI_MH dD at +3 dB, xi=4, unmodulated carrier", cw.loc["static_lti_mh", "Delta_mean"], f"{pm(cw.loc['static_lti_mh', 'Delta_mean'], 3)}", R06 / "14_followup_carrier_ablation_pairs.csv")


def phase():
    cfg = pd.read_csv(PM / "02_phase_mechanism_configs.csv").set_index("case")
    add("F_REF_ntr", "REF nonzero transitions / phase excursion / mean |df|", (cfg.loc["QPSK_REF", "n_transitions_nonzero"], cfg.loc["QPSK_REF", "total_phase_excursion_rad"], cfg.loc["QPSK_REF", "mean_abs_freq_offset_Hz"]),
        f"{cfg.loc['QPSK_REF', 'n_transitions_nonzero']} transitions, {cfg.loc['QPSK_REF', 'total_phase_excursion_rad']:.1f} rad, {cfg.loc['QPSK_REF', 'mean_abs_freq_offset_Hz'] / 1e3:.0f} kHz", PM / "02_phase_mechanism_configs.csv")
    pk = cfg.loc[["QPSK_JUMP", "QPSK_REF", "QPSK_LONGRAMP"], "peak_abs_freq_offset_Hz"]
    add("F_peak", "peak |f_inst - f_IF| JUMP/REF/LONGRAMP (MHz)", tuple(pk / 1e6), f"{pk.iloc[0] / 1e6:.0f} / {pk.iloc[1] / 1e6:.2f} / {pk.iloc[2] / 1e6:.2f} MHz", PM / "02_phase_mechanism_configs.csv")
    add("F_peak_ratio", "ratio of peak excursion JUMP / LONGRAMP", pk.iloc[0] / pk.iloc[2], f"≈{pk.iloc[0] / pk.iloc[2]:.0f}", PM / "02_phase_mechanism_configs.csv")
    rate = cfg.loc[["QPSK_SLOW", "QPSK_REF", "QPSK_FAST"], "transition_rate_per_us"]
    add("F_rates", "transition rates SLOW/REF/FAST per us", tuple(rate), f"{rate.iloc[0]:.2f} / {rate.iloc[1]:.2f} / {rate.iloc[2]:.2f}", PM / "02_phase_mechanism_configs.csv")
    add("F_drift", "mean signed frequency offset SLOW / REF / FAST, +pi convention (kHz)", tuple(cfg.loc[["QPSK_SLOW", "QPSK_REF", "QPSK_FAST"], "mean_signed_freq_offset_Hz"] / 1e3),
        " / ".join(f"{v / 1e3:.0f}" for v in cfg.loc[["QPSK_SLOW", "QPSK_REF", "QPSK_FAST"], "mean_signed_freq_offset_Hz"]) + " kHz", PM / "02_phase_mechanism_configs.csv")
    res = pd.read_csv(PM / "03_phase_mechanism_results.csv")
    m = res[(res.Nd == 4001) & (res.dt_ns == 1)].set_index("case")
    add("F_gCW", "full-model CW gain at a_H", m.loc["CW", "g_final"], f"{m.loc['CW', 'g_final']:.3f}", PM / "03_phase_mechanism_results.csv")
    add("F_gstatic", "CW-static prediction at a_H", m.loc["CW", "static_prediction_g"], f"{m.loc['CW', 'static_prediction_g']:.3f}", PM / "03_phase_mechanism_results.csv")
    add("F_gREF", "full-model QPSK_REF gain", m.loc["QPSK_REF", "g_final"], f"{m.loc['QPSK_REF', 'g_final']:.3f}", PM / "03_phase_mechanism_results.csv")
    pp = m[["g_final_static_lti_mh", "g_final_dsh"]].stack()
    add("F_gmodels", "pipeline gain of Hammerstein & slow-state models over all cases", (pp.min(), pp.max()), f"{pp.min():.3f}–{pp.max():.3f}", PM / "03_phase_mechanism_results.csv")
    for c in m.index:
        add(f"F_X_{c}", f"X {c}", m.loc[c, "X"], neg(m.loc[c, "X"]), PM / "03_phase_mechanism_results.csv")
    num = []
    for c in ("CW", "QPSK_REF", "QPSK_JUMP", "OFF_+1.31"):
        for nd, dt in ((8001, 1), (4001, .5)):
            num.append(abs(float(res[(res.case == c) & (res.Nd == nd) & (res.dt_ns == dt)].X.iloc[0]) - m.loc[c, "X"]))
    add("F_numerics", "max |dX| over Nd/dt checks", max(num), f"{max(num):.4f}", PM / "03_phase_mechanism_results.csv")
    cw = m.loc["CW"]
    for c in ("OFF_+0.125", "OFF_+1.31", "OFF_-1.31", "OFF_+2.62", "OFF_-2.62"):
        add(f"F_abs_{c}", f"large/small-signal output vs CW, {c}", (m.loc[c, "abs_out_last20"] / cw.abs_out_last20, m.loc[c, "abs_lin_last20"] / cw.abs_lin_last20),
            f"{m.loc[c, 'abs_out_last20'] / cw.abs_out_last20:.2f}× / {m.loc[c, 'abs_lin_last20'] / cw.abs_lin_last20:.2f}×", PM / "03_phase_mechanism_results.csv")
    ph = [c for c in m.index if c.startswith(("QPSK", "TOGGLE"))]
    ry = (m.loc[ph, "pop_rydberg_3"] + m.loc[ph, "pop_rydberg_4"]) / (cw.pop_rydberg_3 + cw.pop_rydberg_4) - 1
    add("F_ryd", "Rydberg population change vs CW, phase cases (%)", (ry.min(), ry.max()), f"{pm(100 * ry.min(), 1)}% to {pm(100 * ry.max(), 1)}%", PM / "03_phase_mechanism_results.csv")
    pb = (m.loc[ph, "probe_db_last20us"] - cw.probe_db_last20us).abs().max()
    add("F_probe", "max |probe-baseline shift| vs CW, phase cases (dB)", pb, f"{pb:.2f} dB", PM / "03_phase_mechanism_results.csv")
    ob = (m.loc[[c for c in m.index if c.startswith("OFF")], "probe_db_last20us"] - cw.probe_db_last20us).abs().max()
    add("F_probe_off", "max |probe-baseline shift| vs CW, offsets (dB)", ob, f"{ob:.2f} dB", PM / "03_phase_mechanism_results.csv")
    t63 = m.loc[["QPSK_REF", "QPSK_REF_s2", "QPSK_JUMP", "QPSK_LONGRAMP", "TOGGLE"], "t63_us"]
    add("F_t63", "build-up t63 (REF, s2, JUMP, LONGRAMP, TOGGLE)", (t63.min(), t63.max()), f"{t63.min():.1f}–{t63.max():.1f} µs", PM / "03_phase_mechanism_results.csv")
    it = pd.read_csv(PH / "posthoc_isolated_transitions.csv")
    h = it[np.isclose(it.dphi_rad.abs(), np.pi / 2, atol=.1)]
    p_ = it[np.isclose(it.dphi_rad.abs(), np.pi, atol=.1)]
    add("H_pi2_mag", "pi/2 steps: min |z_full| rel. pre (first 0.5 us)", (h.min_absfull_rel_first_0p5us.min(), h.min_absfull_rel_first_0p5us.max()),
        f"{h.min_absfull_rel_first_0p5us.min():.2f}–{h.min_absfull_rel_first_0p5us.max():.2f}", PH / "posthoc_isolated_transitions.csv")
    add("H_pi2_lin", "pi/2 steps: min |z_lin| rel. pre", (h.min_abslin_rel_first_0p5us.min(), h.min_abslin_rel_first_0p5us.max()),
        f"{h.min_abslin_rel_first_0p5us.min():.2f}–{h.min_abslin_rel_first_0p5us.max():.2f}", PH / "posthoc_isolated_transitions.csv")
    well = h[h.min_absfull_rel_first_0p5us > .4]
    add("H_pi2_phase", "pi/2 steps: max phase lag (rad), excluding the case with near-vanishing magnitude", (well.max_abs_dtheta_rad.min(), well.max_abs_dtheta_rad.max()),
        f"{well.max_abs_dtheta_rad.min():.1f}–{well.max_abs_dtheta_rad.max():.1f} rad", PH / "posthoc_isolated_transitions.csv")
    add("H_pi2_time", "pi/2 steps: time of max phase lag (us)", (h.time_max_dtheta_us.min(), h.time_max_dtheta_us.max()), f"{h.time_max_dtheta_us.min():.1f}–{h.time_max_dtheta_us.max():.1f} µs", PH / "posthoc_isolated_transitions.csv")
    add("H_pi_mag", "pi steps: max |z_full| rel. pre", (p_.max_absfull_rel_first_0p5us.min(), p_.max_absfull_rel_first_0p5us.max()),
        f"{p_.max_absfull_rel_first_0p5us.min():.1f}–{p_.max_absfull_rel_first_0p5us.max():.1f}×", PH / "posthoc_isolated_transitions.csv")
    add("H_pi_lin", "pi steps: min |z_lin| rel. pre", p_.min_abslin_rel_first_0p5us.min(), f"{p_.min_abslin_rel_first_0p5us.min():.2f}", PH / "posthoc_isolated_transitions.csv")
    add("H_n", "number of isolated pi/2 and pi transitions", (len(h), len(p_)), f"{len(h)} and {len(p_)}", PH / "posthoc_isolated_transitions.csv")
    sr = pd.read_csv(PH / "posthoc_steady_relative_phase.csv").set_index("case")
    phc = [c for c in sr.index if c.startswith(("QPSK", "TOGGLE"))]
    add("H_jitter", "steady relative-phase circular SD, phase cases (rad)", (sr.loc[phc, "circ_std_final_rad"].min(), sr.loc[phc, "circ_std_final_rad"].max()),
        f"{sr.loc[phc, 'circ_std_final_rad'].min():.2f}–{sr.loc[phc, 'circ_std_final_rad'].max():.2f} rad", PH / "posthoc_steady_relative_phase.csv")
    offc = [c for c in sr.index if c.startswith("OFF")]
    add("H_jitter_off", "steady relative-phase circular SD, offsets (max)", sr.loc[offc, "circ_std_final_rad"].max(), f"{sr.loc[offc, 'circ_std_final_rad'].max():.2f} rad", PH / "posthoc_steady_relative_phase.csv")
    ts = pd.read_csv(PH / "posthoc_timescales.csv")
    pw_ = ts[ts.transient.str.startswith("power")].osc_period_us
    ph_ = ts[ts.transient.str.startswith("phase")].osc_period_us
    add("H_per_power", "power-step oscillation period at a_H (us)", tuple(pw_), f"{pw_.min():.2f}–{pw_.max():.2f} µs", PH / "posthoc_timescales.csv")
    add("H_per_phase", "phase-step oscillation period at a_H (us)", (ph_.min(), ph_.max()), f"{ph_.min():.2f}–{ph_.max():.2f} µs", PH / "posthoc_timescales.csv")


def replication():
    cfg = json.loads((FV / "replication" / "02_run_configs.json").read_text())
    add("V_E1p", "E1dB at LO 0.35 V/m", cfg["E1dB_prime_Vpm"], f"{cfg['E1dB_prime_Vpm']:.4f} V/m", FV / "replication" / "02_run_configs.json")
    add("V_ap", "operating amplitude a'", cfg["a_prime_Vpm"], f"{cfg['a_prime_Vpm']:.4f} V/m", FV / "replication" / "02_run_configs.json")
    add("V_aH", "operating amplitude a_H (original)", 3 * .5774 * 10 ** (3 / 20) * 0.07321886766113006, "0.179 V/m", R06 / "phase_mechanism" / "00_preregistration.json")
    r = pd.read_csv(FV / "replication" / "03_replication_results.csv").set_index("case")
    add("V_gCW", "full-model CW gain at a' / static prediction", (r.loc["LO2_CW", "g_final"], r.loc["LO2_CW", "g_static_prediction"]),
        f"{r.loc['LO2_CW', 'g_final']:.3f} / {r.loc['LO2_CW', 'g_static_prediction']:.3f}", FV / "replication" / "03_replication_results.csv")
    for c in r.index:
        add(f"V_X_{c}", f"X {c}", r.loc[c, "X"], neg(r.loc[c, "X"]), FV / "replication" / "03_replication_results.csv")
    xr = (r.loc["LO2_QPSK_REF_zd", "X"] + r.loc["LO2_QPSK_REF_zd_s2", "X"]) / 2
    add("V_Xref_mean", "mean X_REF' over two seeds", xr, f"{xr:.3f}", FV / "replication" / "03_replication_results.csv")
    add("V_drift", "zero-drift mean signed offset REF / FAST (kHz)", (r.loc["LO1_QPSK_REF_zd", "mean_signed_freq_offset_Hz"] / 1e3, r.loc["LO1_QPSK_FAST_zd", "mean_signed_freq_offset_Hz"] / 1e3),
        f"{r.loc['LO1_QPSK_REF_zd', 'mean_signed_freq_offset_Hz'] / 1e3:.0f} / {r.loc['LO1_QPSK_FAST_zd', 'mean_signed_freq_offset_Hz'] / 1e3:.0f} kHz", FV / "replication" / "03_replication_results.csv")
    v = pd.read_csv(FV / "replication" / "04_replication_verdicts.csv").set_index("item")
    add("V_class", "preregistered classification", v.loc["CLASSIFICATION", "outcome"], "B", FV / "replication" / "04_replication_verdicts.csv")
    num = v.loc["numerics (Nd=8001 closure)", "value"]
    add("V_numerics", "Nd=8001 closure checks", num, "−0.1730 → −0.1746; 0.04%", FV / "replication" / "04_replication_verdicts.csv")
    wb = json.loads((FV / "replication" / "05_posthoc_C9a_wideband.json").read_text())
    add("V_c9a", "post-hoc wide-band periods at LO' power / phase median (us)", (wb["LO2_power_step_period_us"], wb["LO2_phase_step_median_us"]),
        f"{wb['LO2_power_step_period_us']:.3f} vs {wb['LO2_phase_step_median_us']:.3f} µs", FV / "replication" / "05_posthoc_C9a_wideband.json")
    add("V_c9a_range", "post-hoc phase-step period range at LO' (us)", (min(wb["LO2_phase_step_periods_us"]), max(wb["LO2_phase_step_periods_us"])),
        f"{min(wb['LO2_phase_step_periods_us']):.3f}–{max(wb['LO2_phase_step_periods_us']):.3f} µs", FV / "replication" / "05_posthoc_C9a_wideband.json")
    add("V_orig_wide", "wide-band estimator at the original point (us)", wb["original_power_step_period_wideband_us"], f"{wb['original_power_step_period_wideband_us']:.2f} µs", FV / "replication" / "05_posthoc_C9a_wideband.json")


def revision():
    """Pre-submission revision (preregistered in results/revision/00_revision_preregistration.json)."""
    RV = ROOT / "results" / "revision"
    P1R = RV / "p1"
    ci3 = lambda m, lo, hi, n=3: f"{pm(m, n)} [{pm(lo, n)}, {pm(hi, n)}]"  # noqa: E731
    ci3u = lambda m, lo, hi, n=3: f"{m:.{n}f} [{lo:.{n}f}, {hi:.{n}f}]"  # noqa: E731
    # P1: preregistered robustness table (context lag columns)
    t = pd.read_csv(P1R / "04_p1_decisive_table.csv")
    s1 = json.loads((P1R / "05_p1_summary.json").read_text())
    add("RV_p1_verdict", "P1 preregistered verdict", s1["verdict"], "does not survive", P1R / "05_p1_summary.json")
    add("RV_p1_repro", "P1 zero-fill 2-us D vs archived campaign (max rel.)", s1["zerofill_2us_vs_archived_campaign_max_rel"], "7e-13", P1R / "05_p1_summary.json")
    d1 = t[(t.metric == "D") & (t.contrast == "C1")]
    add("RV_p1_D_C1_range", "C1 full-model change in D over FIR spans 2-16 us", (d1.full_change.min(), d1.full_change.max()),
        f"{d1.full_change.min():.3f}–{d1.full_change.max():.3f}", P1R / "04_p1_decisive_table.csv")
    d2 = t[(t.metric == "D") & (t.contrast == "C2")]
    add("RV_p1_D_C2_range", "C2 full-model change in D over FIR spans 2-16 us", (d2.full_change.min(), d2.full_change.max()),
        f"{neg(d2.full_change.max())} to {neg(d2.full_change.min())}", P1R / "04_p1_decisive_table.csv")
    c = pd.read_csv(P1R / "02_p1_contrasts.csv")
    def crow(design, var, pdb, metric, span=2):
        return c[(c.design == design) & (c.variant == var) & (c.Pavg_over_P1dB_dB == pdb) & (c.method == "context") & (c.metric == metric) & (c.span_us == span)].iloc[0]
    r = crow("dwell", "xi_4", 3, "D_ref")
    add("RV_p1_Dref_C1_full", "C1 full-model change in D_ref (2 us)", r.d_atomic, ci3(r.d_atomic, r.d_atomic_lo, r.d_atomic_hi), P1R / "02_p1_contrasts.csv")
    add("RV_p1_Dref_C1_lti", "C1 LTI+static change in D_ref (2 us)", r.d_static_lti_mh, f"{pm(r.d_static_lti_mh, 3)}", P1R / "02_p1_contrasts.csv")
    add("RV_p1_Dref_C1_exc", "C1 excess over LTI+static in D_ref (2 us)", r.exc_static_lti_mh, ci3(r.exc_static_lti_mh, r.exc_static_lti_mh_lo, r.exc_static_lti_mh_hi), P1R / "02_p1_contrasts.csv")
    add("RV_p1_Dref_C1_exc_static", "C1 excess over static in D_ref (2 us)", r.exc_static_mh, f"{pm(r.exc_static_mh, 2)}", P1R / "02_p1_contrasts.csv")
    r0 = crow("dwell", "xi_4", 0, "D_ref")
    add("RV_p1_Dref_C1_0dB", "C1 at 0 dB: full vs LTI+static change in D_ref (2 us)", (r0.d_atomic, r0.d_static_lti_mh), f"{pm(r0.d_atomic, 3)} vs {pm(r0.d_static_lti_mh, 3)}", P1R / "02_p1_contrasts.csv")
    r = crow("shuffle", "BLOCK_SHUFFLED_CYCLES", 3, "D_ref")
    add("RV_p1_Dref_C2_full", "C2 full-model change in D_ref (2 us)", r.d_atomic, ci3(r.d_atomic, r.d_atomic_lo, r.d_atomic_hi, 4), P1R / "02_p1_contrasts.csv")
    add("RV_p1_Dref_C2_exc", "C2 excess over LTI+static in D_ref (2 us)", r.exc_static_lti_mh, ci3(r.exc_static_lti_mh, r.exc_static_lti_mh_lo, r.exc_static_lti_mh_hi, 4), P1R / "02_p1_contrasts.csv")
    ex = t[t.metric == "D_ref"]
    add("RV_p1_Dref_C1_exc_range", "C1 D_ref excess over LTI+static, spans 2-16 us", (ex[ex.contrast == "C1"].excess_vs_LTI_static.min(), ex[ex.contrast == "C1"].excess_vs_LTI_static.max()),
        f"{neg(ex[ex.contrast == 'C1'].excess_vs_LTI_static.min(), 4)} to {neg(ex[ex.contrast == 'C1'].excess_vs_LTI_static.max(), 4)}", P1R / "04_p1_decisive_table.csv")
    # Table S7: per metric x contrast, ranges over the four FIR spans and whether every criterion holds
    for metric in ("D", "D_ref"):
        for cc in ("C1", "C2"):
            q = t[(t.metric == metric) & (t.contrast == cc)]
            nd = 3 if metric == "D" else 4
            rng = lambda col: (pm(q[col].min(), nd) if pm(q[col].min(), nd) == pm(q[col].max(), nd) else f"{pm(q[col].min(), nd)} to {pm(q[col].max(), nd)}")  # noqa: E731
            add(f"RV_T7_{metric}_{cc}_full", f"Table S7 {metric} {cc}: full-model change, range over spans", (q.full_change.min(), q.full_change.max()), rng("full_change"), P1R / "04_p1_decisive_table.csv")
            add(f"RV_T7_{metric}_{cc}_lti", f"Table S7 {metric} {cc}: excess over LTI+static, range over spans", (q.excess_vs_LTI_static.min(), q.excess_vs_LTI_static.max()), rng("excess_vs_LTI_static"), P1R / "04_p1_decisive_table.csv")
            add(f"RV_T7_{metric}_{cc}_static", f"Table S7 {metric} {cc}: excess over static, range over spans", (q.excess_vs_static.min(), q.excess_vs_static.max()), rng("excess_vs_static"), P1R / "04_p1_decisive_table.csv")
            add(f"RV_T7_{metric}_{cc}_pass", f"Table S7 {metric} {cc}: all criteria hold at all spans", bool(q.all_pass.all()), "yes" if q.all_pass.all() else "no", P1R / "04_p1_decisive_table.csv")
    # P1 post hoc: fitted gain decomposition
    g = pd.read_csv(P1R / "07_p1_posthoc_gain_vs_dwell.csv")
    def gv(car, pdb, xi, m, col):
        return float(g[(g.carrier == car) & (g.Pavg_over_P1dB_dB == pdb) & np.isclose(g.xi, xi) & (g.model == m)][col].iloc[0])
    add("RV_g_qpsk_ref_xi4", "fitted gain QPSK +3 dB: reference vs xi=4 (full)", (gv("qpsk", 3, .05, "atomic", "gain_rel_mean"), gv("qpsk", 3, 4, "atomic", "gain_rel_mean")),
        f"{gv('qpsk', 3, .05, 'atomic', 'gain_rel_mean'):.3f} to {gv('qpsk', 3, 4, 'atomic', 'gain_rel_mean'):.3f}", P1R / "07_p1_posthoc_gain_vs_dwell.csv")
    add("RV_g_qpsk_d_xi4", "paired gain change xi=4, QPSK +3 dB (full)", gv("qpsk", 3, 4, "atomic", "d_gain_rel_mean"),
        ci3(gv("qpsk", 3, 4, "atomic", "d_gain_rel_mean"), gv("qpsk", 3, 4, "atomic", "d_gain_rel_lo"), gv("qpsk", 3, 4, "atomic", "d_gain_rel_hi")), P1R / "07_p1_posthoc_gain_vs_dwell.csv")
    sur = max(abs(gv("qpsk", 3, 4, m, "d_gain_rel_mean")) for m in ("static_mh", "static_lti_mh"))
    add("RV_g_qpsk_d_sur", "max |paired gain change| of the all-zone surrogates, xi=4, QPSK +3 dB", sur, f"{sur:.3f}", P1R / "07_p1_posthoc_gain_vs_dwell.csv")
    add("RV_g_qpsk_d_xi4_0dB", "paired gain change xi=4, QPSK 0 dB (full)", gv("qpsk", 0, 4, "atomic", "d_gain_rel_mean"),
        ci3(gv("qpsk", 0, 4, "atomic", "d_gain_rel_mean"), gv("qpsk", 0, 4, "atomic", "d_gain_rel_lo"), gv("qpsk", 0, 4, "atomic", "d_gain_rel_hi")), P1R / "07_p1_posthoc_gain_vs_dwell.csv")
    add("RV_g_qpsk_xi025", "fitted gain QPSK +3 dB at xi=0.25 (full)", gv("qpsk", 3, .25, "atomic", "gain_rel_mean"), f"{gv('qpsk', 3, .25, 'atomic', 'gain_rel_mean'):.3f}", P1R / "07_p1_posthoc_gain_vs_dwell.csv")
    add("RV_g_cw_ref_xi4", "fitted gain unmodulated +3 dB: reference, xi=0.25, xi=4 (full)",
        (gv("cw", 3, .05, "atomic", "gain_rel_mean"), gv("cw", 3, .25, "atomic", "gain_rel_mean"), gv("cw", 3, 4, "atomic", "gain_rel_mean")),
        f"{gv('cw', 3, .05, 'atomic', 'gain_rel_mean'):.2f}, {gv('cw', 3, .25, 'atomic', 'gain_rel_mean'):.2f}, {gv('cw', 3, 4, 'atomic', 'gain_rel_mean'):.2f}", P1R / "07_p1_posthoc_gain_vs_dwell.csv")
    add("RV_g_cw_d_xi4", "paired gain change xi=4, unmodulated +3 dB (full)", gv("cw", 3, 4, "atomic", "d_gain_rel_mean"),
        ci3(gv("cw", 3, 4, "atomic", "d_gain_rel_mean"), gv("cw", 3, 4, "atomic", "d_gain_rel_lo"), gv("cw", 3, 4, "atomic", "d_gain_rel_hi"), 2), P1R / "07_p1_posthoc_gain_vs_dwell.csv")
    add("RV_g_dsh_xi4", "single-slow-state vs full paired gain change at xi=4 (+3 dB): QPSK; unmodulated",
        (gv("qpsk", 3, 4, "dsh_tau2.535us", "d_gain_rel_mean"), gv("qpsk", 3, 4, "atomic", "d_gain_rel_mean"), gv("cw", 3, 4, "dsh_tau2.535us", "d_gain_rel_mean"), gv("cw", 3, 4, "atomic", "d_gain_rel_mean")),
        f"{neg(gv('qpsk', 3, 4, 'dsh_tau2.535us', 'd_gain_rel_mean'))} vs {neg(gv('qpsk', 3, 4, 'atomic', 'd_gain_rel_mean'))}; {neg(gv('cw', 3, 4, 'dsh_tau2.535us', 'd_gain_rel_mean'))} vs {neg(gv('cw', 3, 4, 'atomic', 'd_gain_rel_mean'))}",
        P1R / "07_p1_posthoc_gain_vs_dwell.csv")
    add("RV_g_dsh_xi025", "single-slow-state vs full paired gain change at xi=0.25 (+3 dB): QPSK; unmodulated",
        (gv("qpsk", 3, .25, "dsh_tau2.535us", "d_gain_rel_mean"), gv("qpsk", 3, .25, "atomic", "d_gain_rel_mean"), gv("cw", 3, .25, "dsh_tau2.535us", "d_gain_rel_mean"), gv("cw", 3, .25, "atomic", "d_gain_rel_mean")),
        f"{neg(gv('qpsk', 3, .25, 'dsh_tau2.535us', 'd_gain_rel_mean'))} vs {neg(gv('qpsk', 3, .25, 'atomic', 'd_gain_rel_mean'))}; {neg(gv('cw', 3, .25, 'dsh_tau2.535us', 'd_gain_rel_mean'))} vs {neg(gv('cw', 3, .25, 'atomic', 'd_gain_rel_mean'))}",
        P1R / "07_p1_posthoc_gain_vs_dwell.csv")
    # P2: six zero-drift sequences
    s2 = pd.read_csv(RV / "11_p2_summary.csv").set_index(["op", "case"])
    for (op, case), row in s2.iterrows():
        add(f"RV_X_{op}_{case}", f"X {case} {op}, six zero-drift sequences (mean [95% CI])", row.X_mean, ci3u(row.X_mean, row.ci_lo, row.ci_hi), RV / "11_p2_summary.csv")
        add(f"RV_Xm_{op}_{case}", f"X {case} {op}, six-sequence mean (2 d.p.)", row.X_mean, f"{row.X_mean:.2f}", RV / "11_p2_summary.csv")
    for op in ("OP1", "OP2"):
        add(f"RV_sd_{op}", f"sequence-to-sequence SD of X_REF at {op}", s2.loc[(op, "REF"), "sd"], f"{s2.loc[(op, 'REF'), 'sd']:.2f}", RV / "11_p2_summary.csv")
        add(f"RV_rng_{op}", f"range of X_REF over six sequences at {op}", (s2.loc[(op, "REF"), "min"], s2.loc[(op, "REF"), "max"]),
            f"{s2.loc[(op, 'REF'), 'min']:.2f}–{s2.loc[(op, 'REF'), 'max']:.2f}", RV / "11_p2_summary.csv")
    c2 = pd.read_csv(RV / "12_p2_contrasts.csv")
    for _, row in c2.iterrows():
        key = row.contrast.replace("OP1-OP2 ", "OPdiff_").replace("-REF", "")
        add(f"RV_c_{row.op}_{key}", f"paired contrast {row.contrast} ({row.op}), six sequences", row["mean"], ci3(row["mean"], row.ci_lo, row.ci_hi, 2), RV / "12_p2_contrasts.csv")
    per = pd.read_csv(RV / "10_p2_per_sequence.csv")
    for case in ("LONGRAMP", "FAST"):
        a = per[(per.op == "OP1") & (per.case == case)].set_index("seed").X - per[(per.op == "OP1") & (per.case == "REF")].set_index("seed").X
        b = per[(per.op == "OP2") & (per.case == case)].set_index("seed").X - per[(per.op == "OP2") & (per.case == "REF")].set_index("seed").X
        m, lo, hi = ci_paired((a - b).values)
        add(f"RV_opdiff_{case}", f"({case}-REF) at OP1 minus at OP2, paired by sequence", m, ci3(m, lo, hi, 2), RV / "10_p2_per_sequence.csv")
    gref = per[(per.op == "OP1") & (per.case == "REF")]
    gg = (gref.g_CW * (1 - gref.X)).mean()
    add("RV_gREF_OP1", "OP1 gain under REF, six-sequence mean (relative to small-signal)", gg, f"{gg:.2f}", RV / "10_p2_per_sequence.csv")
    det = pd.read_csv(RV / "13_p2_determinism_check.csv")
    add("RV_p2_identical", "archived zero-drift runs reproduced bit-identically", int(det.z_full_identical.sum()), f"{int(det.z_full_identical.sum())}", RV / "13_p2_determinism_check.csv")
    # P3
    s3 = pd.read_csv(RV / "21_p3_summary.csv").set_index(["op", "level_dB"])
    for (op, lv), row in s3.iterrows():
        add(f"RV_P3_{op}_{lv:+d}", f"X_REF at {lv:+d} dB re P1dB, {op} (3 sequences)", row.X_mean, f"{row.X_mean:.2f}", RV / "21_p3_summary.csv")
    sm = json.loads((RV / "50_p2_p6_summary.json").read_text())
    add("RV_P3_onset", "onset level (dB re P1dB) OP1 / OP2", (sm["P3_onset_dB_re_P1dB"]["OP1"], sm["P3_onset_dB_re_P1dB"]["OP2"]), "0 dB / −3 dB", RV / "50_p2_p6_summary.json")
    # P4
    sw = pd.read_csv(RV / "30_p4_offset_sweep.csv")
    Rv = lambda op, f: float(sw[(sw.op == op) & np.isclose(sw.offset_MHz, f)].R.iloc[0])  # noqa: E731
    add("RV_R_OP1_near", "OP1 R at -0.25 and +0.125 MHz", (Rv("OP1", -.25), Rv("OP1", .125)), f"{Rv('OP1', -.25):.2f} and {Rv('OP1', .125):.2f}", RV / "30_p4_offset_sweep.csv")
    near2 = sw[(sw.op == "OP2") & (sw.offset_MHz.abs() <= .75)].R
    add("RV_R_OP2_near", "OP2 R range for |df| <= 0.75 MHz", (near2.min(), near2.max()), f"{near2.min():.2f}–{near2.max():.2f}", RV / "30_p4_offset_sweep.csv")
    add("RV_R_OP1_max", "OP1 R at +3 MHz", Rv("OP1", 3.0), f"{Rv('OP1', 3.0):.1f}", RV / "30_p4_offset_sweep.csv")
    o3 = pd.read_csv(RV / "31_p4_op3_offsets.csv")
    add("RV_R_OP3_near", "OP3 R at +0.125 MHz", float(o3[np.isclose(o3.offset_MHz, .125)].R.iloc[0]), f"{float(o3[np.isclose(o3.offset_MHz, .125)].R.iloc[0]):.2f}", RV / "31_p4_op3_offsets.csv")
    ops = pd.read_csv(RV / "32_p4_operating_points.csv").set_index("op")
    add("RV_LO_OP3", "LO field at OP3", ops.loc["OP3", "A_LO_Vpm"], f"{ops.loc['OP3', 'A_LO_Vpm']:.3f} V/m", RV / "32_p4_operating_points.csv")
    add("RV_E1_OP3", "E1dB at OP3 (A_LO = 0.425 V/m)", ops.loc["OP3", "E1dB_Vpm"], f"{ops.loc['OP3', 'E1dB_Vpm']:.4f} V/m", RV / "32_p4_operating_points.csv")
    add("RV_S5", "selectivity S5 at OP2 / OP3 / OP1", tuple(ops.loc[o, "S5"] for o in ("OP2", "OP3", "OP1")),
        " / ".join(f"{ops.loc[o, 'S5']:.2f}" for o in ("OP2", "OP3", "OP1")), RV / "32_p4_operating_points.csv")
    add("RV_XREF3", "X_REF (3 sequences) at OP2 / OP3 / OP1", tuple(ops.loc[o, "X_REF_mean_3seeds"] for o in ("OP2", "OP3", "OP1")),
        " / ".join(f"{ops.loc[o, 'X_REF_mean_3seeds']:.2f}" for o in ("OP2", "OP3", "OP1")), RV / "32_p4_operating_points.csv")
    add("RV_XREF_OP3_ci", "X_REF at OP3, 3 sequences, 95% CI", ops.loc["OP3", "X_REF_mean_3seeds"],
        ci3u(ops.loc["OP3", "X_REF_mean_3seeds"], ops.loc["OP3", "X_REF_ci_lo"], ops.loc["OP3", "X_REF_ci_hi"], 2), RV / "32_p4_operating_points.csv")
    add("RV_P4_class", "P4 preregistered classification", sm["P4_classification"], "MONOTONIC_IN_LO", RV / "50_p2_p6_summary.json")
    # P5
    t1 = pd.read_csv(RV / "41_p5_T1.csv").set_index("op")
    add("RV_P5_T1", "g(a,f) surrogate X_REF at OP1 / OP2 (preregistered median)", (t1.loc["OP1", "X_sur_REF"], t1.loc["OP2", "X_sur_REF"]),
        f"{t1.loc['OP1', 'X_sur_REF']:.3f} / {neg(t1.loc['OP2', 'X_sur_REF'])}", RV / "41_p5_T1.csv")
    mb = pd.read_csv(RV / "45_p5_posthoc_mean_based.csv")
    add("RV_P5_mean", "g(a,f) surrogate mean-based X range over cases and OPs (post hoc)", (mb.X_sur_meanbased.min(), mb.X_sur_meanbased.max()),
        f"{neg(mb.X_sur_meanbased.min(), 2)} to {neg(mb.X_sur_meanbased.max(), 3)}", RV / "45_p5_posthoc_mean_based.csv")
    t3 = pd.read_csv(RV / "43_p5_T3.csv")
    add("RV_P5_T3max", "max |X_sur| over the amplitude sweep", float(t3.X_sur_mean.abs().max()), f"{float(t3.X_sur_mean.abs().max()):.3f}", RV / "43_p5_T3.csv")
    add("RV_P5_verdict", "P5 preregistered verdict", sm["P5_verdict"], "FAILS", RV / "50_p2_p6_summary.json")
    # P6
    p6 = sm["P6"]
    add("RV_P6_X", "X_REF FP64 vs FP32 (OP1, seed 20260701)", (p6["phase_fp64"]["X_fp64"], p6["phase_fp64"]["X_fp32_gpu"], p6["phase_fp64"]["abs_dX"]),
        f"{p6['phase_fp64']['X_fp64']:.4f} vs {p6['phase_fp64']['X_fp32_gpu']:.4f}", RV / "50_p2_p6_summary.json")
    add("RV_P6_long", "long-dwell D_ref full FP64 vs FP32; gap change (pp)", (p6["long_dwell_fp64"]["D_ref_full_fp64"], p6["long_dwell_fp64"]["D_ref_full_fp32"], p6["long_dwell_fp64"]["gap_change_pp"]),
        f"{p6['long_dwell_fp64']['D_ref_full_fp64']:.5f} vs {p6['long_dwell_fp64']['D_ref_full_fp32']:.5f}; {p6['long_dwell_fp64']['gap_change_pp']:.3f} percentage points", RV / "50_p2_p6_summary.json")
    w = pd.DataFrame(p6["extended_segment"]["windows"])
    add("RV_P6_ext", "extended 200-us segment: X per 20-us window min-max; first; last", (w.X.min(), w.X.max(), w.X.iloc[0], w.X.iloc[-1]),
        f"{w.X.min():.2f}–{w.X.max():.2f}; {w.X.iloc[0]:.2f} → {w.X.iloc[-1]:.2f}", RV / "50_p2_p6_summary.json")
    tm = p6["timing"]
    add("RV_P6_time", "wall time: GPU 110-us phase run, 840-us dwell job, 2.4-ms long dwell; CPU FP64 110 us, 2.4 ms (s)",
        (tm["gpu_phase_run_110us_s"], tm["gpu_dwell_job_840us_median_s"], tm["gpu_long_dwell_2.44ms_archived_runtime_s"], tm["cpu_fp64_phase_run_110us_s"], tm["cpu_fp64_long_dwell_2.44ms_s"]),
        f"{tm['gpu_phase_run_110us_s']:.0f} s, {tm['gpu_dwell_job_840us_median_s']:.0f} s, {tm['gpu_long_dwell_2.44ms_archived_runtime_s']:.0f} s; {tm['cpu_fp64_phase_run_110us_s']:.0f} s, {tm['cpu_fp64_long_dwell_2.44ms_s']:.0f} s",
        RV / "50_p2_p6_summary.json")


def revision2():
    """Final receiver-physics validation pass (decision rules: results/revision2/00_decision_rules.json, 01_*, 02_deviations.json)."""
    R2 = ROOT / "results" / "revision2"
    P1R = ROOT / "results" / "revision" / "p1"
    ci3 = lambda m, lo, hi, n=3: f"{pm(m, n)} [{pm(lo, n)}, {pm(hi, n)}]"  # noqa: E731
    # ---- archived revision P1 declustering contrast at +6 dB and over spans (Section 4.3)
    c = pd.read_csv(P1R / "02_p1_contrasts.csv")
    sh = c[(c.design == "shuffle") & (c.method == "context") & (c.metric == "D_ref")]
    r6 = sh[(sh.Pavg_over_P1dB_dB == 6) & (sh.span_us == 2)].iloc[0]
    add("R2_C2_Dref_full_6dB", "C2 full-model change in D_ref at +6 dB (2 us)", r6.d_atomic, ci3(r6.d_atomic, r6.d_atomic_lo, r6.d_atomic_hi, 4), P1R / "02_p1_contrasts.csv")
    add("R2_C2_Dref_exc_6dB", "C2 excess over LTI+static in D_ref at +6 dB (2 us)", r6.exc_static_lti_mh,
        ci3(r6.exc_static_lti_mh, r6.exc_static_lti_mh_lo, r6.exc_static_lti_mh_hi, 4), P1R / "02_p1_contrasts.csv")
    r3 = sh[sh.Pavg_over_P1dB_dB == 3]
    add("R2_C2_Dref_span_range", "C2 full-model change in D_ref at +3 dB over FIR spans 2-16 us (min, max)", (r3.d_atomic.min(), r3.d_atomic.max()),
        f"{neg(r3.d_atomic.min(), 4)} and {neg(r3.d_atomic.max(), 4)}", P1R / "02_p1_contrasts.csv")
    # ---- P1: IF dependence
    s1 = pd.read_csv(R2 / "11_p1_if_summary.csv").set_index("IF_MHz")
    v1 = json.loads((R2 / "15_p1_verdict.json").read_text())
    src = R2 / "11_p1_if_summary.csv"
    add("R2_IFs", "IFs of the P1 check (MHz)", (5, 10, 15), "5, 10 and 15 MHz", R2 / "00_decision_rules.json")
    add("R2_E1_IF", "E1dB at IF 5/10/15 MHz (V/m)", tuple(s1.E1dB_Vpm), f"{s1.E1dB_Vpm[5]:.4f}, {s1.E1dB_Vpm[10]:.3f} and {s1.E1dB_Vpm[15]:.3f} V/m", src)
    add("R2_E1_ratio10", "E1dB(10 MHz)/E1dB(5 MHz)", s1.E1dB_Vpm[10] / s1.E1dB_Vpm[5], f"{s1.E1dB_Vpm[10] / s1.E1dB_Vpm[5]:.1f}", src)
    add("R2_aH_IF10", "matched phase level at 10 MHz (2.4468 E1dB, V/m)", 2.4468 * s1.E1dB_Vpm[10], f"{2.4468 * s1.E1dB_Vpm[10]:.2f} V/m", src)
    add("R2_gCW_IF", "CW gain at the matched level rel. small-signal, 5/10/15 MHz", tuple(s1.CW_gain_at_aH_rel_small_signal),
        f"{s1.CW_gain_at_aH_rel_small_signal[5]:.2f}, {s1.CW_gain_at_aH_rel_small_signal[10]:.2f} and {s1.CW_gain_at_aH_rel_small_signal[15]:.2f}", src)
    for f_ in (5, 10, 15):
        r = s1.loc[f_]
        add(f"R2_XREF_IF{f_}", f"X_REF at {f_} MHz (3 sequences, mean [95% CI])", r.X_REF_mean, ci3(r.X_REF_mean, r.X_REF_lo, r.X_REF_hi), src)
        add(f"R2_dg_IF{f_}", f"paired fitted-gain change xi=4 minus 0.05 at {f_} MHz: full [CI]; static", (r.d_gain_full, r.d_gain_full_lo, r.d_gain_full_hi, r.d_gain_static),
            f"{ci3(r.d_gain_full, r.d_gain_full_lo, r.d_gain_full_hi)}; static {pm(r.d_gain_static, 3)}", src)
    add("R2_dg_full_IF", "paired gain change (full) at 5/10/15 MHz", tuple(s1.d_gain_full), f"{neg(s1.d_gain_full[5])}, {neg(s1.d_gain_full[10])} and {neg(s1.d_gain_full[15])}", src)
    add("R2_dg_static_IF", "paired gain change (static surrogate) at 5/10/15 MHz", tuple(s1.d_gain_static),
        f"{neg(s1.d_gain_static[5])}, {neg(s1.d_gain_static[10])} and {neg(s1.d_gain_static[15])}", src)
    add("R2_S5_IF", "frequency selectivity S5 at 5/10/15 MHz", tuple(s1.S5), f"{s1.S5[5]:.2f}, {s1.S5[10]:.2f} and {s1.S5[15]:.2f}", src)
    add("R2_R0125_IF", "R(+0.125 MHz) at 5/10/15 MHz", tuple(s1["R_+0.125"]), f"{s1['R_+0.125'][5]:.2f}, {s1['R_+0.125'][10]:.2f} and {s1['R_+0.125'][15]:.2f}", src)
    add("R2_P1_class", "P1 prespecified classification", v1["classification"], "does not survive", R2 / "15_p1_verdict.json")
    add("R2_P1_ratio", "15/5 MHz ratio of the dwell gain change", v1["d_gain_ratio_15_over_5"], f"{100 * v1['d_gain_ratio_15_over_5']:.0f}%", R2 / "15_p1_verdict.json")
    val = pd.read_csv(R2 / "10_p1_validation_5MHz.csv")
    add("R2_val5", "5 MHz validation: phase runs bit-identical; max rel. diff of dwell gains", (bool(v1["validation_5MHz_all_phase_runs_identical"]), v1["validation_5MHz_max_dwell_gain_rel_diff"]),
        "bit-identical; 1e-14", R2 / "10_p1_validation_5MHz.csv")
    t10 = np.load(R2 / "tables" / "OP1_IF10_ext.npz")
    t15 = np.load(R2 / "tables" / "OP1_IF15.npz")
    add("R2_exp_IF10", "10 MHz: max CW fundamental gain re weak-tone slope (dB) and its amplitude", (float(np.nanmax(t10["gain_db"])), float(t10["amps"][np.nanargmax(t10["gain_db"])])),
        f"+{np.nanmax(t10['gain_db']):.1f} dB near {t10['amps'][np.nanargmax(t10['gain_db'])]:.1f} V/m", R2 / "tables" / "OP1_IF10_ext.npz")
    add("R2_cw15_07", "15 MHz: CW fundamental gain at 0.7 V/m (dB)", float(t15["gain_db"][-1]), f"{neg(t15['gain_db'][-1], 1)} dB", R2 / "tables" / "OP1_IF15.npz")
    dw = pd.read_csv(R2 / "12_p1_dwell_per_realization.csv")
    dw["fms"] = dw.gain_rel_full - dw.gain_rel_static
    fm = dw.groupby(["IF_MHz", "variant"]).fms.mean()
    add("R2_fms", "full minus static fitted gain (xi 0.05; xi 4) at 5/10/15 MHz", tuple(fm.values),
        f"{pm(fm[(5, 'xi_0.05')])}/{pm(fm[(5, 'xi_4')])}, {pm(fm[(10, 'xi_0.05')], 3)}/{pm(fm[(10, 'xi_4')], 3)}, {pm(fm[(15, 'xi_0.05')], 3)}/{pm(fm[(15, 'xi_4')], 3)}",
        R2 / "12_p1_dwell_per_realization.csv")
    lk = pd.read_csv(R2 / "13_p1_zone_leakage.csv").groupby(["IF_MHz", "variant"]).mean(numeric_only=True)
    add("R2_leak0", "zone-0 in-band power re zone 1 (dB), static model, xi 0.05 at 5/10/15 MHz", tuple(lk.xs("xi_0.05", level=1).zone0_rel_zone1_dB),
        ", ".join(f"{pm(v, 1)}" for v in lk.xs("xi_0.05", level=1).zone0_rel_zone1_dB) + " dB", R2 / "13_p1_zone_leakage.csv")
    add("R2_leak0_xi4", "zone-0 in-band power re zone 1 (dB), static model, xi 4 at 5/10/15 MHz", tuple(lk.xs("xi_4", level=1).zone0_rel_zone1_dB),
        ", ".join(f"{pm(v, 1)}" for v in lk.xs("xi_4", level=1).zone0_rel_zone1_dB) + " dB", R2 / "13_p1_zone_leakage.csv")
    add("R2_leak23", "max in-band power of zones 2-3 re zone 1 (dB), any IF/variant", float(lk[["zone2_rel_zone1_dB", "zone3_rel_zone1_dB"]].max().max()),
        f"{lk[['zone2_rel_zone1_dB', 'zone3_rel_zone1_dB']].max().max():.0f} dB".replace("-", "−"), R2 / "13_p1_zone_leakage.csv")
    de = pd.read_csv(R2 / "14_p1_ref_decomposition.csv").groupby("IF_MHz").mean(numeric_only=True)
    add("R2_Xcoh_IF", "X_coh of the reference case (3 sequences) at 5/10/15 MHz", tuple(de.X_coh),
        f"{pm(de.X_coh[5], 3)}, {pm(de.X_coh[10], 3)} and {pm(de.X_coh[15], 3)}", R2 / "14_p1_ref_decomposition.csv")
    add("R2_criteria", "prespecified criterion constants quoted in the evidence table", None,
        "0.8 X_REF; 0.25; 0.10; 0.20; 0.05; 3x; slow tau >= 1 us, weight >= 0.1; ratio outside [0.5, 2]; dip >= 20%, within 0.3 us, r >= 0.8; 25%; 2%",
        "results/stage06_dwell_physics/phase_mechanism/00_preregistration.json; results/final_validation/replication/00_preregistration.json; results/revision/00_revision_preregistration.json; results/revision2/00_decision_rules.json; results/stage06_dwell_physics/11_followup_preregistration.json")
    add("R2_sd_IF", "circular SD of the response phase (rad), reference case, 5/10/15 MHz", tuple(de.phase_circ_sd_rad),
        f"{de.phase_circ_sd_rad[5]:.2f}, {de.phase_circ_sd_rad[10]:.2f} and {de.phase_circ_sd_rad[15]:.2f} rad", R2 / "14_p1_ref_decomposition.csv")
    # ---- P2: coherent vs magnitude decomposition
    s2 = pd.read_csv(R2 / "21_p2_summary.csv")
    v2 = json.loads((R2 / "22_p2_verdict.json").read_text())
    ref = s2[s2.case == "REF"].set_index("op")
    for op in ("OP1", "OP2", "OP3"):
        r = ref.loc[op]
        add(f"R2_dec_{op}", f"{op} reference: X, X_mag, X_rms, X_coh; magnitude share; phase share; circ. SD; lag (rad)",
            (r.X, r.X_mag, r.X_rms, r.X_coh, r.share_magnitude, r.share_phase, r.phase_circ_sd_rad, r.phase_lag_vs_CW_rad),
            f"X {r.X:.3f}, X_mag {r.X_mag:.3f}, X_rms {r.X_rms:.3f}, X_coh {r.X_coh:.3f}; {100 * r.share_magnitude:.0f}% / {100 * r.share_phase:.0f}%; "
            f"{r.phase_circ_sd_rad:.2f} rad; {r.phase_lag_vs_CW_rad:.2f} rad", R2 / "21_p2_summary.csv")
    add("R2_gcoh_rel_OP1", "OP1 reference: coherent gain relative to CW, 1 - X_coh (six sequences)", 1 - ref.loc["OP1"].X_coh, f"{1 - ref.loc['OP1'].X_coh:.2f}", R2 / "21_p2_summary.csv")
    add("R2_share_mag", "magnitude share of ln(1-X_coh), reference, OP1/OP2/OP3", tuple(ref.share_magnitude[["OP1", "OP2", "OP3"]]),
        f"{100 * ref.share_magnitude['OP1']:.0f}%, {100 * ref.share_magnitude['OP2']:.0f}% and {100 * ref.share_magnitude['OP3']:.0f}%", R2 / "21_p2_summary.csv")
    add("R2_share_ph", "phase-misalignment share, reference, OP1/OP2/OP3", tuple(ref.share_phase[["OP1", "OP2", "OP3"]]),
        f"{100 * ref.share_phase['OP1']:.0f}%, {100 * ref.share_phase['OP2']:.0f}% and {100 * ref.share_phase['OP3']:.0f}%", R2 / "21_p2_summary.csv")
    fa = s2[(s2.op == "OP1") & (s2.case == "FAST")].iloc[0]
    add("R2_share_ph_fast", "phase-misalignment share, FAST at OP1", fa.share_phase, f"{100 * fa.share_phase:.0f}%", R2 / "21_p2_summary.csv")
    dx = v2["X_vs_X_mag_max_abs_diff_by_op_all_cases"]
    add("R2_XvsXmag", "max |X - X_mag| over cases and sequences at OP1/OP2/OP3", (dx["OP1"], dx["OP2"], dx["OP3"]),
        f"{dx['OP1']:.2f} at OP1 ({dx['OP2']:.2f} and {dx['OP3']:.2f} at OP2 and OP3)", R2 / "22_p2_verdict.json")
    lr = s2[(s2.op == "OP2") & (s2.case == "LONGRAMP")].iloc[0]
    add("R2_share_ph_lr", "phase-misalignment share, 900-ns ramps at OP2 (X small)", lr.share_phase, f"{100 * lr.share_phase:.0f}%", R2 / "21_p2_summary.csv")
    fx = s2[(s2.op == "OP1") & (s2.case == "FAST")].iloc[0]
    add("R2_fast_X_Xmag", "OP1 FAST: mean X vs X_mag", (fx.X, fx.X_mag), f"{fx.X:.3f} vs {fx.X_mag:.3f}", R2 / "21_p2_summary.csv")
    # ---- P3: weaker probe
    v3 = json.loads((R2 / "31_p3_verdict.json").read_text())
    w, sd = v3["per_probe"]["weak"], v3["per_probe"]["standard"]
    src3 = R2 / "31_p3_verdict.json"
    add("R2_wp_Op", "weak-probe Rabi frequency /2pi (MHz)", w["Omega_p_over_2pi_MHz"], f"{w['Omega_p_over_2pi_MHz']:.2f} MHz", src3)
    add("R2_wp_E1", "weak-probe E1dB (V/m) and ratio to standard", (w["E1dB_Vpm"], w["E1dB_Vpm"] / sd["E1dB_Vpm"]), f"{w['E1dB_Vpm']:.3f} V/m", src3)
    add("R2_wp_aH", "weak-probe matched phase level (V/m)", 2.4468 * w["E1dB_Vpm"], f"{2.4468 * w['E1dB_Vpm']:.2f} V/m", src3)
    add("R2_wp_gCW", "weak-probe CW gain at the matched level rel. small-signal", w["CW_gain_at_aH_rel_small_signal"], f"{w['CW_gain_at_aH_rel_small_signal']:.3f}", src3)
    add("R2_wp_X", "weak-probe X_REF (3 sequences)", w["X_REF_mean"], ci3(w["X_REF_mean"], w["X_REF_lo"], w["X_REF_hi"]), src3)
    add("R2_wp_dg", "weak-probe dwell gain change: full [CI]; static", (w["d_gain_full"], w["d_gain_full_lo"], w["d_gain_full_hi"], w["d_gain_static"]),
        f"{ci3(w['d_gain_full'], w['d_gain_full_lo'], w['d_gain_full_hi'])}; static {pm(w['d_gain_static'], 3)}", src3)
    add("R2_wp_tau", "slow component of the +5% step at 1.0 E1dB: weak (tau us, weight); standard (tau us, weight)",
        (w["step_plus5_fit2_tau_slow_s"], w["step_plus5_fit2_a_slow"], sd["step_plus5_fit2_tau_slow_s"], sd["step_plus5_fit2_a_slow"]),
        f"{w['step_plus5_fit2_tau_slow_s'] * 1e6:.2f} µs, weight {w['step_plus5_fit2_a_slow']:.2f}; {sd['step_plus5_fit2_tau_slow_s'] * 1e6:.2f} µs, weight {abs(sd['step_plus5_fit2_a_slow']):.2f}", src3)
    add("R2_wp_tau1e", "1/e time of the +5% step at 1.0 E1dB: weak; standard (us)", (w["step_plus5_tau_1e_s"], sd["step_plus5_tau_1e_s"]),
        f"{w['step_plus5_tau_1e_s'] * 1e6:.2f} µs vs {sd['step_plus5_tau_1e_s'] * 1e6:.2f} µs", src3)
    add("R2_wp_verdict", "P3 verdicts: memory, phase loss, dwell gain persist; material change", (v3["memory_persists"], v3["phase_loss_persists"], v3["dwell_gain_persists"], v3["material_change"]),
        "memory persists; phase loss and dwell gain excess do not; material change", src3)
    wf = {}
    for v in ("xi_0.05", "xi_4"):
        fr, st = [], []
        for r in range(4):
            z = np.load(R2 / "p3" / "dwell" / f"dwell_r{r}_{v}_p+3.npz")
            gl = abs(complex(z["g_linear"]))
            fr.append(abs(complex(z["g_full"])) / gl)
            st.append(abs(complex(z["g_static"])) / gl)
        wf[v] = float(np.mean(fr) - np.mean(st))
    add("R2_wp_fms", "weak probe: full minus static fitted gain (xi 0.05; xi 4)", (wf["xi_0.05"], wf["xi_4"]),
        f"{pm(wf['xi_0.05'])}/{pm(wf['xi_4'])}", R2 / "p3" / "dwell")
    e1r = (s1.E1dB_Vpm[15] / s1.E1dB_Vpm[5], s1.E1dB_Vpm[10] / s1.E1dB_Vpm[5], w["E1dB_Vpm"] / sd["E1dB_Vpm"])
    add("R2_E1_ratios", "E1dB relative to the reference configuration: 15 MHz, 10 MHz, weak probe", e1r,
        f"{e1r[0]:.1f}, {e1r[1]:.1f} and {e1r[2]:.1f}", src3)
    add("R2_wp_lin", "small-signal kernel linearity check: weak; standard", (v3["kernel_lin_check_weak"], 0.07652518991777806), f"{v3['kernel_lin_check_weak']:.3f} vs 0.077", src3)
    tw = np.load(R2 / "tables" / "OP1_IF5_weakprobe_ext.npz")
    add("R2_wp_exp", "weak probe: max CW gain re weak-tone slope (dB) and amplitude", (float(np.nanmax(tw["gain_db"])), float(tw["amps"][np.nanargmax(tw["gain_db"])])),
        f"+{np.nanmax(tw['gain_db']):.1f} dB near {tw['amps'][np.nanargmax(tw['gain_db'])]:.1f} V/m", R2 / "tables" / "OP1_IF5_weakprobe_ext.npz")
    # ---- P4: internal traces
    v4 = json.loads((R2 / "41_p4_verdict.json").read_text())
    mA, mS = v4["meta"]
    src4 = R2 / "41_p4_verdict.json"
    add("R2_p4_window", "P4 recorded window and step centre (us)", (49, 55, 50), "from 49 to 55 µs (step centred at 50 µs)", R2 / "00_decision_rules.json")
    add("R2_level", "matched level relative to E1dB (a_H/E1dB)", 2.4468, "2.45", ROOT / "python" / "experiments" / "revision" / "rev_phase.py")
    add("R2_p4_chain", "chained vs continuous probe output, max deviation (a_H; small-signal)", (mA["chain_dev"], mS["chain_dev"]), f"{mA['chain_dev']:.0e} and {mS['chain_dev']:.0e}", src4)
    add("R2_p4_out", "a_H step: output IF magnitude minimum (fraction of pre-step), time after step (us); max phase deviation (rad)",
        (1 - mA["output_dip"], mA["output_t_min_us"] - 50.0, mA["output_phase_max_abs_dev_rad"]),
        f"{1 - mA['output_dip']:.2f} at {mA['output_t_min_us'] - 50.0:.1f} µs; {mA['output_phase_max_abs_dev_rad']:.2f} rad", src4)
    add("R2_p4_out_ss", "small-signal step: output magnitude minimum (fraction), max phase deviation (rad)", (1 - mS["output_dip"], mS["output_phase_max_abs_dev_rad"]),
        f"{1 - mS['output_dip']:.2f}; {mS['output_phase_max_abs_dev_rad']:.2f} rad", src4)
    d4 = pd.read_csv(R2 / "40_p4_internal_tracking.csv")
    th = d4[(d4["class"] == "thermal") & (d4.case == "step_aH")].set_index(["element", "component"])
    r43 = th.loc[("rho43 (RF)", "+IF")]
    add("R2_p4_rho43", "rho43 +IF: minimum (fraction of pre-step); time vs output minimum (us); r", (1 - r43.dip, r43.dt_min_vs_output_us, r43.r_vs_output),
        f"{1 - r43.dip:.2f}; {r43.dt_min_vs_output_us:+.2f} µs; r = {r43.r_vs_output:.2f}".replace("-", "−"), R2 / "40_p4_internal_tracking.csv")
    r41 = th.loc[("rho41 (three-photon)", "-IF")]
    add("R2_p4_rho41", "rho41 -IF: dip; time vs output minimum; r", (r41.dip, r41.dt_min_vs_output_us, r41.r_vs_output),
        f"{100 * r41.dip:.0f}%; {r41.dt_min_vs_output_us:+.2f} µs; r = {r41.r_vs_output:.2f}".replace("-", "−"), R2 / "40_p4_internal_tracking.csv")
    ims = th.loc[("rho21 (probe)", "Im part, IF")]
    add("R2_p4_cons", "probe observable Im(rho21) IF component vs output: r", ims.r_vs_output, f"{ims.r_vs_output:.3f}", R2 / "40_p4_internal_tracking.csv")
    pops = [abs(th.loc[(e, "slow")].dip) for e in ("rho22", "rho33", "rho44")]
    add("R2_p4_pops", "max |relative change| of thermally averaged populations (slow part)", max(pops), f"{100 * max(pops):.1f}%", R2 / "40_p4_internal_tracking.csv")
    tr, trs = np.load(R2 / "p4" / "step_aH_traces.npz"), np.load(R2 / "p4" / "step_smallsignal_traces.npz")
    ok = lambda z: (z["t"] >= 49.2e-6) & (z["t"] < 54.75e-6)  # noqa: E731
    after = lambda z: (z["t"] >= 49.8e-6) & (z["t"] < 54.75e-6)  # noqa: E731
    ph = tr["thermal|rho43 (RF)|+IF|phase"][ok(tr)]
    ph_s = trs["thermal|rho43 (RF)|+IF|phase"][ok(trs)]
    add("R2_p4_drift", "rho43 +IF phase relative to the drive at 54.74 us (rad): a_H; small-signal max |dev|", (float(ph[-1]), float(np.abs(ph_s).max())),
        f"{abs(ph[-1]):.1f} rad within {tr['t'][ok(tr)][-1] * 1e6 - 50:.1f} µs; {np.abs(ph_s).max():.2f} rad", R2 / "p4" / "step_aH_traces.npz")
    sl, sls = tr["thermal|rho43 (RF)|slow"][after(tr)], trs["thermal|rho43 (RF)|slow"][after(trs)]
    add("R2_p4_slow43", "|slow rho43| relative to pre-step: a_H range; small-signal max |change|", (float(sl.min()), float(sl.max()), float(np.abs(sls - 1).max())),
        f"{sl.min():.2f}–{sl.max():.2f}; {100 * np.abs(sls - 1).max():.1f}%", R2 / "p4" / "step_aH_traces.npz")
    ssmax = d4[(d4["class"] == "thermal") & (d4.case == "step_smallsignal") & d4.component.isin(["+IF", "-IF"]) & (d4.element == "rho43 (RF)") & (d4.component == "+IF")].dip.iloc[0]
    add("R2_p4_ss43", "small-signal: rho43 +IF dip", ssmax, f"{100 * ssmax:.0f}%", R2 / "40_p4_internal_tracking.csv")
    add("R2_p4_verdict", "any internal IF component meets the prespecified tracking criterion", v4["any_internal_IF_component_tracks"], "none", src4)


def campaign():
    """Final adversarial-validation campaign (plan: results/tqe_viability_final/00_plan.md; deviations 00_deviations.md)."""
    TVF = ROOT / "results" / "tqe_viability_final"
    M1, M2, M3, M4, M5 = (TVF / d for d in ("01_matching", "02_resonance", "03_memory_model", "04_noise", "05_numerical"))
    ci2 = lambda m, lo, hi, n=2: f"{pm(m, n)} [{pm(lo, n)}, {pm(hi, n)}]"  # noqa: E731
    ci2u = lambda m, lo, hi, n=2: f"{f(m, n)} [{f(lo, n)}, {f(hi, n)}]"  # noqa: E731
    mr = pd.read_csv(M1 / "matching_results.csv")
    run = mr[mr.status == "run"]
    pr = pd.read_csv(M1 / "matching_pairs.csv")
    px = pr[pr.metric == "X"]
    v1 = json.loads((M1 / "matching_verdict.json").read_text())
    plan = json.loads((M1 / "plan.json").read_text())
    src1 = M1 / "matching_results.csv"

    def x(cfg, rule, tgt, col="X"):
        return run[(run.configuration == cfg) & (run.matching_rule == rule) & (abs(run.target - tgt) < 1e-9)].iloc[0][col]
    sig = (0.2, 0.358, 0.4)
    add("TV_plan", "campaign plan registered before any campaign run (commit time)", "2026-09-28T02:37:07+10:00", "2026-09-28", TVF / "00_plan.md")
    add("TV_nseq", "sequences per matched state (REF seeds 20260701-06)", 6, "six", TVF / "00_plan.md")
    add("TV_ratios", "matched signal/LO ratios entering the class", sig, "0.2, 0.358 and 0.4", TVF / "00_plan.md")
    for cfg, lab in (("C1", "5 MHz"), ("C2", "10 MHz"), ("C3", "15 MHz"), ("C4", "weaker probe")):
        vals = [x(cfg, "matched_signal_to_LO", r) for r in sig]
        add(f"TV_Xsig_{cfg}", f"X at matched signal/LO 0.2/0.358/0.4, {lab}", tuple(vals), ", ".join(pm(v, 2) if v < 0 else f(v, 2) for v in vals[:-1]) + " and " + (pm(vals[-1], 2) if vals[-1] < 0 else f(vals[-1], 2)), src1)
    b = px[(px.convention == "B") & px.evaluable & (px.hypothesis == "HA1")]
    add("TV_dsig_range", "paired X difference 5 MHz minus other configuration, evaluable HA1 pairs of convention B (range)", (b["diff"].min(), b["diff"].max()),
        f"{pm(b['diff'].min(), 2)} to {pm(b['diff'].max(), 2)}", M1 / "matching_pairs.csv")
    add("TV_dsig_cilo", "smallest lower CI bound of those differences", b.diff_ci_lo.min(), f"{pm(b.diff_ci_lo.min(), 2)}", M1 / "matching_pairs.csv")
    ratio = float((b.X_c / b.X_ref).max())
    add("TV_ratio_max", "largest X_c / X_ref over convention-B HA1 pairs", ratio, "a fifth", M1 / "matching_pairs.csv")
    a = px[(px.convention == "A") & px.evaluable & (px.hypothesis == "HA1")].sort_values("target", ascending=False)
    add("TV_Xcw_C1", "X at matched CW gain 0.9 / 0.85, 5 MHz", (x("C1", "matched_CW_gain", .9), x("C1", "matched_CW_gain", .85)),
        f"{f(x('C1', 'matched_CW_gain', .9), 2)} and {f(x('C1', 'matched_CW_gain', .85), 2)}", src1)
    add("TV_Xcw_C3", "X at matched CW gain 0.9 / 0.85, 15 MHz", (x("C3", "matched_CW_gain", .9), x("C3", "matched_CW_gain", .85)),
        f"{pm(x('C3', 'matched_CW_gain', .9), 2)} and {pm(x('C3', 'matched_CW_gain', .85), 2)}", src1)
    add("TV_dcw", "paired X difference 5 minus 15 MHz at CW gain 0.9 / 0.85 [95% CI]", tuple(a["diff"]),
        " and ".join(ci2(r["diff"], r.diff_ci_lo, r.diff_ci_hi) for _, r in a.iterrows()), M1 / "matching_pairs.csv")
    add("TV_amp_cw_C3", "15 MHz amplitudes at CW gain 0.9 / 0.85 (V/m) and 5 MHz amplitudes", (x("C3", "matched_CW_gain", .9, "signal_amplitude_Vpm"), x("C3", "matched_CW_gain", .85, "signal_amplitude_Vpm")),
        f"{x('C1', 'matched_CW_gain', .9, 'signal_amplitude_Vpm'):.3f} and {x('C1', 'matched_CW_gain', .85, 'signal_amplitude_Vpm'):.3f} V/m at 5 MHz; "
        f"{x('C3', 'matched_CW_gain', .9, 'signal_amplitude_Vpm'):.3f} and {x('C3', 'matched_CW_gain', .85, 'signal_amplitude_Vpm'):.3f} V/m at 15 MHz", src1)
    reach = plan["reachable_CW_gain"]
    add("TV_reach", "CW gain range with signal <= LO: 10 MHz, weaker probe; 15 MHz minimum", (reach["C2"]["min_CW_gain_signal_le_LO"], reach["C2"]["max_CW_gain_signal_le_LO"],
        reach["C4"]["min_CW_gain_signal_le_LO"], reach["C4"]["max_CW_gain_signal_le_LO"], reach["C3"]["min_CW_gain_signal_le_LO"]),
        f"{reach['C2']['min_CW_gain_signal_le_LO']:.2f}–{reach['C2']['max_CW_gain_signal_le_LO']:.2f}; {reach['C4']['min_CW_gain_signal_le_LO']:.2f}–{reach['C4']['max_CW_gain_signal_le_LO']:.2f}; "
        f"{reach['C3']['min_CW_gain_signal_le_LO']:.2f}", M1 / "plan.json")
    add("TV_E1_over_LO", "E1dB / A_LO at 10 MHz and with the weaker probe", (reach["C2"]["E1dB_over_LO"], reach["C4"]["E1dB_over_LO"]),
        f"{reach['C2']['E1dB_over_LO']:.2f} and {reach['C4']['E1dB_over_LO']:.2f}", M1 / "plan.json")
    inf = pd.DataFrame(plan["infeasible"])
    need = lambda c: inf[(inf.config == c) & inf.amp_Vpm.notna()].amp_Vpm / 0.5  # noqa: E731
    add("TV_need", "signal/LO needed for CW-gain targets 0.95-0.4: 10 MHz; weaker probe", (need("C2").min(), need("C2").max(), need("C4").min(), need("C4").max()),
        f"{need('C2').min():.2f}–{need('C2').max():.2f}; {need('C4').min():.2f}–{need('C4').max():.2f}", M1 / "plan.json")
    r8 = run[(run.configuration == "C1") & (run.matching_rule == "matched_CW_gain") & (run.target == 0.8)].iloc[0]
    r6 = run[(run.configuration == "C1") & (run.matching_rule == "matched_CW_gain") & (run.target == 0.6)].iloc[0]
    rH = run[(run.configuration == "C1") & (run.matching_rule == "matched_signal_to_LO") & (run.target == 0.358)].iloc[0]
    add("TV_N5_08", "5 MHz, CW gain 0.8: X [CI]; local CW slope", (r8.X, r8.local_CW_slope), f"{ci2u(r8.X, r8.X_ci_lo, r8.X_ci_hi)}; {pm(r8.local_CW_slope, 2)}", src1)
    add("TV_N5_06", "5 MHz, CW gain 0.6: X; local CW slope", (r6.X, r6.local_CW_slope), f"{f(r6.X, 2)}; {pm(r6.local_CW_slope, 2)}", src1)
    add("TV_slope_ref", "local CW slope at the reference level a_H (signal/LO 0.358)", rH.local_CW_slope, pm(rH.local_CW_slope, 2), src1)
    add("TV_ref_Xcoh", "reference state (signal/LO 0.358, 5 MHz): X_coh [CI], six sequences", rH.X_coh, ci2u(rH.X_coh, rH.X_coh_ci_lo, rH.X_coh_ci_hi, 3), src1)
    op = {c: [x(c, "matched_CW_gain", t_) for t_ in (0.8, 0.6, 0.4)] for c in ("C1", "C5", "C6")}
    for c, lab in (("C1", "OP1"), ("C5", "OP3"), ("C6", "OP2")):
        add(f"TV_op_{lab}", f"X at matched CW gain 0.8/0.6/0.4, {lab}", tuple(op[c]), " / ".join(f(v, 2) for v in op[c]), src1)
    for c, lab in (("C6", "OP2"), ("C5", "OP3")):
        rr = px[(px.convention == "A") & (px.config == c)].sort_values("target", ascending=False)
        add(f"TV_opd_{lab}", f"paired X difference OP1 minus {lab} at CW gain 0.8/0.6/0.4 [CI]", tuple(rr["diff"]), ", ".join(ci2(r["diff"], r.diff_ci_lo, r.diff_ci_hi) for _, r in rr.iterrows()),
            M1 / "matching_pairs.csv")
    add("TV_HA1", "HA1 pairs holding / evaluable", (v1["details"]["HA1_pairs_holding"], sum(v1["details"]["HA1_evaluable_pairs"].values())),
        f"{v1['details']['HA1_pairs_holding']} of {sum(v1['details']['HA1_evaluable_pairs'].values())}", M1 / "matching_verdict.json")
    add("TV_HA2", "HA2 pairs holding / evaluable", (v1["details"]["HA2_pairs_holding"], v1["details"]["HA2_evaluable_pairs"]),
        f"{v1['details']['HA2_pairs_holding']} of {v1['details']['HA2_evaluable_pairs']}", M1 / "matching_verdict.json")
    add("TV_P1_class", "Part I class (X; X_coh robustness)", (v1["classification"], v1["robustness_X_coh"]["classification"]), "partially survives; survives with X_coh", M1 / "matching_verdict.json")
    add("TV_patho", "pathological matched states (ratio 0.6 at 5 MHz and weaker probe)", 0.6, "0.6", M1 / "plan.json")
    add("TV_targets", "CW-gain targets of the prespecified list", (0.95, 0.9, 0.85, 0.8, 0.6, 0.4), "0.95–0.4", TVF / "00_plan.md")
    add("TV_cwtargets", "CW-gain targets matchable between 5 and 15 MHz with an informative reference state", (0.9, 0.85), "0.9 and 0.85", M1 / "matching_pairs.csv")
    absl = {c: x(c, "same_absolute_level", .179) for c in ("C5", "C6")}
    add("TV_abs", "same absolute amplitude 0.179 V/m: X at OP3, OP2", tuple(absl.values()), f"{f(absl['C5'], 2)} and {f(absl['C6'], 2)}", src1)
    # ---- Part II
    v2 = json.loads((M2 / "resonance_verdict.json").read_text())
    sc = pd.read_csv(M2 / "resonance_scan.csv")
    sc = sc[sc.status == "ok"]
    add("TV_scan", "IF scan: IFs per LO field (min-max), range (MHz)", (int(sc.groupby("A_LO_Vpm").size().min()), int(sc.groupby("A_LO_Vpm").size().max()), sc.IF_MHz.min(), sc.IF_MHz.max()),
        f"{sc.IF_MHz.min():.1f} to {sc.IF_MHz.max():.1f} MHz", M2 / "resonance_scan.csv")
    fl = [lo_ for lo_ in ("0.35", "0.425", "0.5")]
    rab = [v2["features"][lo_]["g_min"]["f"] / v2["features"][lo_]["g_min"]["r2"] * 2 for lo_ in fl]
    add("TV_rabi", "LO Rabi frequencies (MHz) at 0.35/0.425/0.5 V/m", tuple(rab), f"{rab[0]:.2f}, {rab[1]:.2f} and {rab[2]:.2f} MHz", M2 / "resonance_verdict.json")
    g2 = [v2["features"][lo_]["g_min"]["r2"] for lo_ in fl]
    s2 = [v2["features"][lo_]["S_max"]["r2"] for lo_ in fl]
    x2 = [v2["features"][lo_]["Xcoh_max"]["r2"] for lo_ in fl]
    add("TV_r2_gmin", "r2 of the CW-gain minimum and selectivity maximum (range over LO fields)", (min(g2 + s2), max(g2 + s2)), f"{min(g2 + s2):.2f}–{max(g2 + s2):.2f}", M2 / "resonance_verdict.json")
    add("TV_r2_xcoh", "r2 of the X_coh maximum (range over LO fields; identical to 2 decimals)", (min(x2), max(x2)),
        f"{min(x2):.2f}" if f"{min(x2):.2f}" == f"{max(x2):.2f}" else f"{min(x2):.2f}–{max(x2):.2f}", M2 / "resonance_verdict.json")
    fg = [v2["features"][lo_]["g_min"]["f"] for lo_ in fl]
    add("TV_spread_f", "spread (max/min) of the CW-gain-minimum IF across LO fields", max(fg) / min(fg), f"{max(fg) / min(fg):.2f}", M2 / "resonance_verdict.json")
    fx = [v2["features"][lo_]["Xcoh_max"]["f"] for lo_ in fl]
    add("TV_f_xcoh", "IF of the X_coh maximum at 0.35/0.425/0.5 V/m (MHz)", tuple(fx), f"{fx[0]:.2f}, {fx[1]:.2f} and {fx[2]:.2f} MHz", M2 / "resonance_verdict.json")
    vx = [v2["features"][lo_]["Xcoh_max"]["value"] for lo_ in fl]
    add("TV_xcoh_max", "X_coh maximum at 0.35/0.425/0.5 V/m", tuple(vx), f"{vx[0]:.2f}, {vx[1]:.2f} and {vx[2]:.2f}", M2 / "resonance_verdict.json")
    cr = [vv for n_, c_ in v2["collapse_ratio_vs_absolute_f"].items() if n_ != "d1" for vv in c_.values()]
    add("TV_collapse", "collapse ratios (normalizations with aligned features)", (min(cr), max(cr)), f"{min(cr):.2f}–{max(cr):.2f}", M2 / "resonance_verdict.json")
    add("TV_II_class", "Part II class", v2["classification"], "suggestive only", M2 / "resonance_verdict.json")
    xg = []
    for lo_ in fl:
        f0 = v2["features"][lo_]["g_min"]["f"]
        r_ = sc[(abs(sc.A_LO_Vpm - float(lo_)) < 1e-9) & (abs(sc.IF_MHz - f0) < 1e-6)].iloc[0]
        xg.append(float(r_.X_coh))
    add("TV_xcoh_gmin", "X_coh at the CW-gain minimum (min, max over LO fields; post hoc)", (min(xg), max(xg)), f"{neg(min(xg), 2)} to {neg(max(xg), 3)}",
        M2 / "resonance_scan.csv")
    sd = v2["sequence_sd_at_5MHz_ratio0.358"]
    add("TV_seqsd", "sequence SD of X_coh at 5 MHz, signal/LO 0.358 (min, max over LO fields)", (min(s_["X_coh"] for s_ in sd.values()), max(s_["X_coh"] for s_ in sd.values())),
        f"{min(s_['X_coh'] for s_ in sd.values()):.2f}–{max(s_['X_coh'] for s_ in sd.values()):.2f}", M2 / "resonance_verdict.json")
    near = []
    for lo_, g in sc.groupby("A_LO_Vpm"):
        gg = g[(g.r1 > .9) & (g.r1 < 1.0)]
        k = gg.X_coh.idxmax()
        near.append((float(gg.loc[k, "r1"]), float(gg.loc[k, "X_coh"])))
    add("TV_N3b", "second X_coh maximum just below f = Omega_LO/2pi: r1 range; X_coh at 0.35/0.425/0.5", near,
        f"{min(n_[0] for n_ in near):.2f}–{max(n_[0] for n_ in near):.2f}; {near[0][1]:.2f}, {near[1][1]:.2f} and {near[2][1]:.2f}", M2 / "resonance_scan.csv")
    r2c = {"OP1": 5 / (rab[2] / 2), "OP3": 5 / (rab[1] / 2), "OP2": 5 / (rab[0] / 2)}
    add("TV_r2_5MHz", "r2 of the 5 MHz IF at OP1/OP3/OP2", tuple(r2c.values()), f"{r2c['OP1']:.2f}, {r2c['OP3']:.2f} and {r2c['OP2']:.2f}", M2 / "resonance_verdict.json")
    add("TV_r2_1015", "r2 of 10 and 15 MHz at OP1", (10 / (rab[2] / 2), 15 / (rab[2] / 2)), f"{10 / (rab[2] / 2):.2f} and {15 / (rab[2] / 2):.2f}", M2 / "resonance_verdict.json")
    osc = pd.read_csv(M2 / "oscillation_comparison.csv")
    o2 = osc[(osc.A_LO_Vpm == .35) & (osc.candidate == "|f_IF - Omega_LO/4pi|")].iloc[0]
    o1 = osc[(osc.A_LO_Vpm == .5) & (osc.candidate == "|f_IF - Omega_LO/4pi|")].iloc[0]
    add("TV_osc", "|f_IF - Omega_LO/4pi| (MHz) at OP2 and OP1 vs observed oscillation frequencies", (o2.candidate_MHz, o1.candidate_MHz),
        f"{o2.candidate_MHz:.2f} MHz against 1.79 MHz at OP2; {o1.candidate_MHz:.2f} MHz against 0.58–0.63 MHz at OP1", M2 / "oscillation_comparison.csv")
    # ---- Part III
    v3 = json.loads((M3 / "gmp_verdict.json").read_text())
    cfg3 = pd.read_csv(M3 / "gmp_config.csv")
    sel = v3["selected"]
    ncoef = (sel["K"] + 9) * sel["M"] + 2
    add("TV_gmp_sel", "selected GMP: order, memory depth (us), taps, complex coefficients", (sel["K"], sel["D"], sel["M"], ncoef),
        f"order {sel['K']}, {sel['D']:.0f} µs ({sel['M']} taps), {ncoef} complex coefficients", M3 / "gmp_config.csv")
    add("TV_gmp_grid", "GMP grid size", len(cfg3), f"{len(cfg3)}", M3 / "gmp_config.csv")
    dep = cfg3.groupby("D").val_nmse_family_mean.min()
    add("TV_gmp_val", "validation family-mean NMSE at 0.25 us and 4 us (best over order and ridge)", (dep[0.25], dep[4.0]), f"{dep[0.25]:.3f} to {dep[4.0]:.3f}", M3 / "gmp_config.csv")
    nm = v3["nmse"]
    best = {fam: min(nm[fam][m] for m in ("static_mh", "static_lti_mh", "dsh")) for fam in nm}
    add("TV_gmp_nmse", "median test NMSE vs full model: GMP / best CW-derived surrogate (dwell, shuffle, phase, fair)",
        {fam: (nm[fam]["gmp"], best[fam]) for fam in nm},
        f"{nm['dwell']['gmp']:.2f} vs {best['dwell']:.2f} (dwell), {nm['shuffle']['gmp']:.2f} vs {best['shuffle']:.2f} (declustering), "
        f"{nm['phase']['gmp']:.2f} vs {best['phase']:.2f} (phase), {nm['fair']['gmp']:.3f} vs {best['fair']:.3f} (fair waveforms)", M3 / "gmp_verdict.json")
    fac = [best[fam] / nm[fam]["gmp"] for fam in ("dwell", "shuffle", "phase")]
    add("TV_gmp_fac", "NMSE improvement factor (dwell, shuffle, phase)", tuple(fac), f"{min(fac):.1f}–{max(fac):.1f}", M3 / "gmp_verdict.json")
    ae = v3["air_err"]
    add("TV_gmp_air", "mean |AIR error| at +3 dB: GMP for QPSK, 16-QAM, OFDM, CE", tuple(ae[m]["gmp"] for m in ("QPSK", "16QAM", "OFDM", "CE")),
        f"{ae['QPSK']['gmp']:.3f}, {ae['16QAM']['gmp']:.3f}, {ae['OFDM']['gmp']:.3f} and {ae['CE']['gmp']:.2f} bit", M3 / "gmp_verdict.json")
    surr = [ae[m][c] for m in ("QPSK", "16QAM") for c in ("static_mh", "static_lti_mh", "dsh")]
    add("TV_surr_air", "mean |AIR error| at +3 dB of the three surrogates for QPSK and 16-QAM (range)", (min(surr), max(surr)), f"{min(surr):.2f}–{max(surr):.2f}", M3 / "gmp_verdict.json")
    add("TV_surr_air_ofdm_ce", "OFDM: LTI+static and single-slow-state errors; CE: LTI+static", (ae["OFDM"]["static_lti_mh"], ae["OFDM"]["dsh"], ae["CE"]["static_lti_mh"]),
        f"{ae['OFDM']['static_lti_mh']:.3f} and {ae['OFDM']['dsh']:.3f}; {ae['CE']['static_lti_mh']:.2f}", M3 / "gmp_verdict.json")
    cc = v3["contrasts"]
    add("TV_gmp_dwell", "test dwell gain contrast: full, GMP, static, LTI+static, single-slow-state", tuple(cc["dwell_gain"][m] for m in ("full", "gmp", "static_mh", "static_lti_mh", "dsh")),
        f"{neg(cc['dwell_gain']['full'])}; {neg(cc['dwell_gain']['gmp'])}; {neg(cc['dwell_gain']['static_mh'])}; {neg(cc['dwell_gain']['static_lti_mh'])}; {neg(cc['dwell_gain']['dsh'])}",
        M3 / "gmp_verdict.json")
    add("TV_gmp_dwell_pct", "GMP overestimate of the dwell contrast (%)", 100 * (cc["dwell_gain"]["gmp"] / cc["dwell_gain"]["full"] - 1),
        f"{100 * (cc['dwell_gain']['gmp'] / cc['dwell_gain']['full'] - 1):.0f}%", M3 / "gmp_verdict.json")
    add("TV_gmp_decl", "test declustering D_ref contrast: full, GMP", (cc["decluster_Dref"]["full"], cc["decluster_Dref"]["gmp"]),
        f"{neg(cc['decluster_Dref']['full'], 4)} and {neg(cc['decluster_Dref']['gmp'], 4)}", M3 / "gmp_verdict.json")
    add("TV_gmp_phase", "test REF case: X, X_coh, X_mag of full vs GMP", tuple((cc[k]["full"], cc[k]["gmp"]) for k in ("X_REF", "X_coh_REF", "X_mag_REF")),
        f"X {cc['X_REF']['full']:.2f} vs {cc['X_REF']['gmp']:.2f}; X_coh {cc['X_coh_REF']['full']:.2f} vs {cc['X_coh_REF']['gmp']:.2f}; X_mag {cc['X_mag_REF']['full']:.2f} vs {cc['X_mag_REF']['gmp']:.2f}",
        M3 / "gmp_verdict.json")
    sm = max(abs(cc[k][m]) for k in ("X_REF", "X_coh_REF") for m in ("static_mh", "static_lti_mh", "dsh"))
    add("TV_surr_phase", "largest |X| or |X_coh| of the surrogates in the REF case", sm, f"{sm:.2f}", M3 / "gmp_verdict.json")
    lv = v3["levels"]["6"]
    add("TV_gmp_lev6", "held-out +6 dB level: X full vs GMP", (lv["full"]["X"], lv["gmp"]["X"]), f"{lv['full']['X']:.2f} vs {lv['gmp']['X']:.2f}", M3 / "gmp_verdict.json")
    add("TV_III_class", "Part III class", v3["classification"], "partial success", M3 / "gmp_verdict.json")
    g = cfg3.groupby(["K", "D"]).val_nmse_family_mean.min()
    add("TV_gmp_edge", "validation family-mean NMSE: order 7 at 2 -> 4 us; 4 us at order 5 -> 7", (g[(7, 2.0)], g[(7, 4.0)], g[(5, 4.0)]),
        f"{g[(7, 2.0)]:.4f} to {g[(7, 4.0)]:.4f}; {g[(5, 4.0)]:.4f} to {g[(7, 4.0)]:.4f}", M3 / "gmp_config.csv")
    qam = max(ae["QPSK"]["gmp"], ae["16QAM"]["gmp"])
    add("TV_gmp_air_qam", "largest GMP |AIR error| at +3 dB for QPSK and 16-QAM (rounded up to 0.01 bit)", qam, f"{np.ceil(qam * 100) / 100:.2f} bit", M3 / "gmp_verdict.json")
    fr = (cc["X_REF"]["gmp"] / cc["X_REF"]["full"], cc["X_coh_REF"]["gmp"] / cc["X_coh_REF"]["full"])
    add("TV_gmp_quarter", "fraction of X and X_coh captured by the GMP (test sequences)", fr, "a quarter", M3 / "gmp_verdict.json")
    # ---- Part IV
    v4 = json.loads((M4 / "noise_verdict.json").read_text())
    t4 = pd.read_csv(M4 / "noise_sensitivity.csv")
    snr = t4[(t4.table == "effective_SNR_dB") & (t4.Pavg_over_P1dB_dB == 3)].groupby("noise_scale").eff_SNR_dB.agg(["min", "max"])
    add("TV_snr", "effective SNR at +3 dB (dB) at c = 0.25, 1, 4 (range over formats)", {c: tuple(snr.loc[c]) for c in (0.25, 1.0, 4.0)},
        f"{snr.loc[0.25, 'min']:.0f}–{snr.loc[0.25, 'max']:.0f} dB at 0.25×, {snr.loc[1.0, 'min']:.0f}–{snr.loc[1.0, 'max']:.0f} dB at 1× and {snr.loc[4.0, 'min']:.0f}–{snr.loc[4.0, 'max']:.0f} dB at 4×",
        M4 / "noise_sensitivity.csv")
    ex = t4[t4.table == "surrogate_minus_full"]
    q = ex[(ex.control == "static_mh") & (ex.modulation == "16QAM") & (ex.Pavg_over_P1dB_dB == 0)].set_index("noise_scale").excess_mean
    sall = t4[(t4.table == "effective_SNR_dB") & (t4.Pavg_over_P1dB_dB == 3)].eff_SNR_dB
    add("TV_snr_all", "effective SNR range at +3 dB over all noise scales and formats (dB)", (sall.min(), sall.max()), f"{sall.min():.0f}–{sall.max():.0f} dB",
        M4 / "noise_sensitivity.csv")
    add("TV_noise_16qam0", "all-zone static, 16-QAM, 0 dB: AIR error at c = 0.25, 1, 4", (q[0.25], q[1.0], q[4.0]), f"{pm(q[0.25])}, {pm(q[1.0])} and {pm(q[4.0])} bit", M4 / "noise_sensitivity.csv")
    q = ex[(ex.control == "static_mh") & (ex.modulation == "OFDM") & (ex.Pavg_over_P1dB_dB == 3)].set_index("noise_scale").excess_mean
    add("TV_noise_ofdm3", "all-zone static, OFDM, +3 dB: AIR error over c (range)", (q.min(), q.max()), f"{pm(q.min())} to {pm(q.max())} bit", M4 / "noise_sensitivity.csv")
    lf = t4[(t4.table == "loss_full") & (t4.Pavg_over_P1dB_dB == 3) & (t4.noise_scale == 4.0)].set_index("modulation")
    add("TV_noise_loss4", "full-model loss at +3 dB, c = 4: CE and OFDM [CI]", (lf.loc["CE", "loss_mean"], lf.loc["OFDM", "loss_mean"]),
        f"{ci2(lf.loc['CE', 'loss_mean'], lf.loc['CE', 'loss_ci_lo'], lf.loc['CE', 'loss_ci_hi'])} and {ci2(lf.loc['OFDM', 'loss_mean'], lf.loc['OFDM', 'loss_ci_lo'], lf.loc['OFDM', 'loss_ci_hi'])} bit",
        M4 / "noise_sensitivity.csv")
    ps = v4["per_scale"]
    p5 = {c: (ps[c]["P5_dB"]["16QAM"], ps[c]["P5_dB"]["OFDM"]) for c in ("0.25", "1.0", "4.0")}
    add("TV_noise_P5", "P5 (dB re P1dB) of 16-QAM and OFDM at c = 0.25, 1, 4", p5,
        f"16-QAM {pm(p5['0.25'][0])}, {pm(p5['1.0'][0])} and {pm(p5['4.0'][0])} dB; OFDM {pm(p5['0.25'][1])}, {pm(p5['1.0'][1])} and {pm(p5['4.0'][1])} dB", M4 / "noise_verdict.json")
    add("TV_IV_class", "Part IV class", v4["classification"], "partially robust", M4 / "noise_verdict.json")
    add("TV_IV_valid", "max |AIR difference| at c = 1 vs archived (all models)", max(v4["validation_c1"]["max_abs_surrogate_and_linear"], v4["validation_c1"]["max_abs_atomic"]),
        "4e-13", M4 / "noise_verdict.json")
    # ---- Part V
    v5 = json.loads((M5 / "numerical_verdict.json").read_text())
    va = pd.read_csv(M5 / "va_trace_positivity.csv")
    add("TV_va", "V-A: max per-class |Tr-1|; min eigenvalue; longest duration (us)", (va.trace_dev_max.max(), va.min_eig_min.min(), va.T_us.max()),
        f"{va.trace_dev_max.max():.1e}; {va.min_eig_min.min():.1e}; {va.T_us.max():.0f} µs", M5 / "va_trace_positivity.csv")
    vb = pd.read_csv(M5 / "vb_e1db_slope_sensitivity.csv")
    c1 = vb[vb.config == "C1"].change_dB_power
    add("TV_vb_C1", "V-B: change of the reference E1dB under the slope variants (dB, min..max)", (c1.min(), c1.max()), f"{pm(c1.min())} to {pm(c1.max())} dB", M5 / "vb_e1db_slope_sensitivity.csv")
    c56 = vb[vb.config.isin(["C5", "C6"])].change_dB_power.abs().max()
    add("TV_vb_C56", "V-B: largest change at the lower LO fields (dB)", c56, f"{c56:.2f} dB", M5 / "vb_e1db_slope_sensitivity.csv")
    vba = pd.read_csv(M5 / "vb_matched_amplitude_sensitivity.csv")
    add("TV_vb_amp", "V-B: largest shift of a matched CW-gain amplitude (%)", vba.change_pct.abs().max(), f"{vba.change_pct.abs().max():.1f}%", M5 / "vb_matched_amplitude_sensitivity.csv")
    vc = pd.read_csv(M5 / "vc_regime_table.csv")
    add("TV_rwa", "V-C: largest peak RF Rabi frequency / carrier", vc.Rabi_to_carrier.max(), f"{vc.Rabi_to_carrier.max():.3f}", M5 / "vc_regime_table.csv")
    arch = vc[vc.state.str.contains("reference phase case at 2.4468") & vc.state.str.startswith(("C2", "C3", "C4"))].sort_values("state")
    add("TV_vc_arch", "V-C: signal/LO of the archived E1dB-matched phase checks at 10 MHz, 15 MHz, weaker probe", tuple(arch.peak_signal_to_LO),
        f"{arch.peak_signal_to_LO.iloc[0]:.2f}, {arch.peak_signal_to_LO.iloc[1]:.2f} and {arch.peak_signal_to_LO.iloc[2]:.2f}", M5 / "vc_regime_table.csv")
    dw = vc[vc.state.str.contains("dwell contrast")].sort_values("state")
    add("TV_vc_dwell", "V-C: dwell high level / LO in the archived configuration checks (10, 15 MHz, weaker probe)", tuple(dw.peak_signal_to_LO),
        f"{dw.peak_signal_to_LO.iloc[0]:.2f}, {dw.peak_signal_to_LO.iloc[1]:.2f} and {dw.peak_signal_to_LO.iloc[2]:.2f}", M5 / "vc_regime_table.csv")
    st = vc[vc.state.str.startswith("C4 (weaker probe) amplitude step")].iloc[0]
    add("TV_vc_step", "V-C: weaker-probe step level / LO", st.peak_signal_to_LO, f"{st.peak_signal_to_LO:.2f}", M5 / "vc_regime_table.csv")
    sh3 = vc[vc.state.str.startswith("shuffle") & vc.state.str.endswith("+3 dB")].iloc[0]
    sh6 = vc[vc.state.str.startswith("shuffle") & vc.state.str.endswith("+6 dB")].iloc[0]
    add("TV_vc_shuffle", "V-C: shuffle family peak / LO at +3 dB; fraction of samples above the LO at +3 and +6 dB",
        (sh3.peak_signal_to_LO, sh3.fraction_of_samples_signal_gt_LO, sh6.fraction_of_samples_signal_gt_LO),
        f"{sh3.peak_signal_to_LO:.2f}; {100 * sh3.fraction_of_samples_signal_gt_LO:.2f}% and {100 * sh6.fraction_of_samples_signal_gt_LO:.2f}%", M5 / "vc_regime_table.csv")
    add("TV_va_T", "V-A durations: phase case and long-dwell case (us)", (float(va[va.case == "phase"].T_us.max()), float(va[va.case == "dwell"].T_us.max())),
        f"{va[va.case == 'phase'].T_us.max():.0f} µs phase case and an {va[va.case == 'dwell'].T_us.max():.0f} µs long-dwell case", M5 / "va_trace_positivity.csv")
    sys.path.insert(0, str(ROOT / "manuscript"))
    import campaign_tables as CT
    for k, fn in CT.TABLES.items():
        add(f"TV_tab{k}", f"Supplementary Table {k} (full rendering; campaign_tables.py)", "table", fn(), TVF)
    add("TV_carrier", "RF carrier of the model's source configuration (Cs 47D5/2 -> 48P3/2), used only for the rotating-wave check", 6.946e9, "6.946 GHz",
        ROOT / "python" / "experiments" / "tqe_viability_final" / "tvf_common.py")
    add("TV_scan_grid", "IF-scan grid: prior-pass IFs, campaign IFs per LO field, IF definition", (29, 31), "29 IFs; 31 IFs (1000/P MHz for integer periods P from 500 to 80 ns)",
        M2 / "resonance_scan.csv")
    fw = vc[vc.source.str.startswith("fair")]
    add("TV_vc_fair", "V-C: largest peak signal/LO of the fair waveforms (+6 dB)", fw.peak_signal_to_LO.max(), f"{fw.peak_signal_to_LO.max():.2f}", M5 / "vc_regime_table.csv")


if __name__ == "__main__":
    model_and_numerics()
    regen()
    dwell_and_shuffle()
    timescale()
    quasi_static()
    carrier_ablation()
    phase()
    replication()
    revision()
    revision2()
    campaign()
    df = pd.DataFrame(L)
    df.to_csv(Path(__file__).resolve().parent / "numbers_ledger.csv", index=False)
    pd.set_option("display.width", 250)
    pd.set_option("display.max_colwidth", 70)
    pd.set_option("display.max_rows", 500)
    print(df[["id", "text", "source"]].to_string())
