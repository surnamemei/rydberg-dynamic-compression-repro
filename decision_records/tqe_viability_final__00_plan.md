# Final adversarial-validation campaign — pre-registered plan

- **Written:** 2026-09-28, before any run of this campaign (the in-file time of the commit that adds this file is the registration time; internal time-stamped project record, not an external registration).
- **Branch:** `tqe-viability-final` (from `tqe-submission` at 07a3cf0). Nothing is pushed, published or submitted.
- **Scope:** the four risks named in the request (fair matching, LO/dressed-state resonance, strength of the surrogate baselines, noise dependence of the AIR conclusions), plus low-cost numerical checks. No new waveform families, no neural networks, no Liouvillian eigenanalysis, no other extensions.
- **Simulator and numerics:** the unmodified upstream simulator at the converged settings (Nd = 4001, dt = 1 ns, 300 K, GPU FP32 kernel). They change only where Part V says so.
- **Class names:** used exactly as given in the request. The operational definitions below are fixed now and are not altered later.

## 0. Starting state, reused data and what has been seen

**Archive** (Step 0, `00_archive/`):
- Both manuscript PDFs, copied, with SHA256:
  - `main.pdf`: 1daf7471…
  - `supplement.pdf`: d0538868…
- The release-package manifest `SHA256SUMS` (itself 92d36e1c…), plus a recomputation of all 1,894 file hashes (0 mismatches).

**Pre-existing uncommitted changes.** 109 tracked files differ from HEAD only in line endings (zero diff with `--ignore-cr-at-eol`). They predate this work, are not committed and are not touched.

**Prior records, unchanged and still binding for their own classifications:**
- `results/tqe_viability/00_decision_rules.json` (fixed 2026-09-28T01:55:49+10:00, commit e941f0f);
- `results/tqe_viability/01_deviations.json` (V1, n_box fix);
- `results/tqe_viability/p1/plan.json`.

Where this plan's rules differ from those of the prior pass, both results are reported. This plan's rules decide the campaign classes.

**Data produced by the prior-pass runner (`tv_run.py p1 p2 p4`, launched before this campaign; code in `python/experiments/tqe_viability/`, committed 07a3cf0).** All of it is reused as input:

| Data | Location | State at registration | Inspected? |
|---|---|---|---|
| P1: 21 matched points × (CW + REF seeds 20260701–03) | `results/tqe_viability/p1/runs/` | complete (84 files) | **one point only**: C1 (IF5), CW-gain target 0.95. g_CW = 0.9502; X = 0.038 / 0.027 / 0.030 for seeds 1–3 (a check of the pipeline). Nothing else. |
| P2: scan at signal/LO 0.358, 3 LO fields × 29 IFs (2.5–12.5 MHz), CW + REF seed 20260701 | `results/tqe_viability/p2/LO*/` | complete (174 files) | no |
| P4: the 192 regen_stage05 fair-waveform jobs rerun, storing every model's noiseless receiver baseband and the baseband of the archived noise draw | `results/tqe_viability/p4/jobs/` | running (47/192 at registration) | no |

**Execution.**
- The prior runner is left to finish its P4 stage in its own directory; stopping it early was not permitted.
- The campaign runner waits for that process to exit before using the GPU, so there are never two concurrent GPU jobs.
- Every campaign job then logs per Step 2.
- New outputs go only to `results/tqe_viability_final/`.

**Archived inputs, read-only:**
- CW tables and kernels: `results/revision2/tables/`, `results/final_validation/artifacts/`, `results/revision/artifacts/`, `results/stage06_dwell_physics/artifacts/controls_Nd4001.npz` and `controls_MH_Nd4001.npz`;
- dwell/shuffle traces: `results/revision/p1/traces/p1/`;
- phase runs: `results/revision/p2/` (seeds 20260701–06 at a_H) and `results/revision/p3/` (OP1 REF, seeds 1–3, at −6 … +8 dB re P1dB);
- archived AIR: `results/final_validation/regen_stage05/jobs/*.json` and `02_…04_*.csv`;
- the revision2 decision records.

## Common definitions

**Phase-case runs.** The phase-case runs use the prior-pass `phase_run`:
- 0–50 µs CW at amplitude a, then 60 µs of phase modulation at constant amplitude a.
- Demodulation: boxcar of n_box samples (the smallest boxcar of at least 200 samples spanning an integer number of IF periods), with base removal.
- The linear reference is the small-signal kernel of the configuration.

**Per-run metrics** (`rev2_p2.metrics` on the window [89.5, 109.5) µs):
- g_med, g_mag, g_rms, g_coh;
- phase mean and circular SD of arg(z_full · conj(z_lin)).

**Point metrics** (each relative to the CW run at the same point):
- X_k = 1 − g_k(REF)/g_k(CW) for k ∈ {med, coh, mag, rms}; X ≡ X_med.
- **Mean phase lag** = circular mean phase (REF) − circular mean phase (CW), wrapped to (−π, π].
- **Phase jitter** = circular SD (REF).
- **Envelope-magnitude reduction** = 1 − mean|z_full|_REF / mean|z_full|_CW, without linear normalization.
- **Absolute CW gain** = median |z_full,CW| / a (probe transmission per V/m). The small-signal gain is median |z_lin,CW| / a.
- **Relative CW gain** g_CW = g_med(CW).

**Table descriptors** (from the configuration's CW table; `tvf_common.describe`):
- E/E1dB;
- the table CW gain;
- the **local CW slope** d ln|output| / d ln a (1 = linear, 0 = saturated, < 0 = the output falls with amplitude);
- signal > LO.

**Uncertainty.**
- Six sequences, REF seeds 20260701–06, shared by all configurations, so that comparisons are paired.
- Mean with a two-sided 95% Student-t CI over sequences; paired differences are formed per sequence.
- A percentile bootstrap 95% CI (10,000 resamples, generator seed 0) is reported alongside.
- Decisions use the t-CI.

**Effect present** (unchanged from the prior pass): mean X ≥ 0.05 and the 95% CI excludes zero.

**Regime.** Matched points require signal/LO ≤ 1 (the superheterodyne regime) and an amplitude inside the calibrated CW table. No extrapolation.

## A. Fair-matching hypotheses (Part I)

**Configurations:**

| Config | Setting | Table / kernel |
|---|---|---|
| C1 | 5 MHz, standard probe (Ω_p/2π = 8.08 MHz), A_LO = 0.5 V/m (reference) | `OP1_IF5` / OP1 kernel |
| C2 | 10 MHz, standard probe, A_LO = 0.5 V/m | `OP1_IF10_ext` / OP1 kernel |
| C3 | 15 MHz, standard probe, A_LO = 0.5 V/m | `OP1_IF15_ext` / OP1 kernel |
| C4 | 5 MHz, weaker probe (2.02 MHz), A_LO = 0.5 V/m | `OP1_IF5_weakprobe_ext` / weak-probe kernel |
| C5 | 5 MHz, standard probe, A_LO = 0.425 V/m (OP3) | `controls_LO0.425_Nd4001` |
| C6 | 5 MHz, standard probe, A_LO = 0.35 V/m (OP2) | `controls_LO0.35_Nd4001` |

C5 and C6 are cheap (about 5 s per run) and are included.

**Hypotheses:**
- **HA1 (headline, IF/probe; the title's "configuration dependence").** At fairly matched states where the reference configuration C1 shows the phase effect, C2, C3 and C4 do not show it. The criterion is X_c < 0.5 X_C1 **and** a paired difference X_C1 − X_c ≥ 0.05 whose 95% CI excludes zero; this is called *dependence* at the pair (state, c).
- **HA2 (operating point).** At fairly matched states where C1 shows the effect, it is larger at the higher LO field. The criterion is X_C1 − X_c ≥ 0.05, 95% CI excluding zero, for c ∈ {C5, C6}. This reproduces the manuscript's OP1 > OP3 > OP2 ordering, which was obtained at E1dB-matched levels.

**Target list, fixed before any QPSK run of this campaign** (computed from the CW tables only by `tvf_plan.py`; full table in `01_matching/plan.json`):
- **I-A, matched CW gain** g_CW/g_small (first crossing on the table, interpolated in dB):
  - Preferred targets 0.8 / 0.6 / 0.4.
  - Nearest common feasible targets with C3: 0.95 / 0.9 / 0.85.
  - Feasible points:
    - C1: 0.95, 0.9, 0.85, 0.8, 0.6, 0.4 (a = 0.0488 … 0.1856 V/m);
    - C3: 0.95, 0.9, 0.85 (a = 0.133, 0.202, 0.381 V/m);
    - C5: 0.8, 0.6, 0.4 (0.086, 0.150, 0.232 V/m);
    - C6: 0.8, 0.6, 0.4 (0.059, 0.111, 0.224 V/m).
  - **Infeasible, documented, not extrapolated:**
    - C2: every target needs signal/LO 1.11–1.63. Its CW gain with signal ≤ LO spans 0.968–1.265 (gain expansion); E1dB = 1.16 × LO.
    - C4: every target needs signal/LO 1.41–1.62. Its CW gain spans 0.999–1.872; E1dB = 1.43 × LO.
    - C3: 0.8 needs ratio 1.05; 0.6 and 0.4 have no crossing up to 0.7 V/m (1.4 × LO). Its minimum reachable CW gain is 0.808.
  - Consequently C2 and C4 cannot be matched in CW gain inside the regime at all. The matched-CW-gain comparison for HA1 therefore rests on C3 (0.95 / 0.9 / 0.85), and for HA2 on C5 and C6 (0.8 / 0.6 / 0.4).
- **I-B, matched signal/LO:**
  - Ratios 0.1, 0.2 and 0.4 as requested, plus the prior-pass ratios 0.358 (the OP1 reference state a_H/A_LO) and 0.6.
  - Configurations: C1–C4 at all five ratios; C5 and C6 at 0.1 / 0.2 / 0.358 / 0.4.
- **Pathology rule** (fixed before running; I-B asks to replace 0.4 if it is pathological):
  - A matched amplitude is *pathological* if the CW gain curve has a local extremum within [0.85 a, 1.15 a]. The test uses the grid points in the window plus the nearest point outside on each side; steps below 0.02 dB count as flat.
  - Ratio 0.4 is **not** pathological in any configuration, so it is kept.
  - Ratio 0.6 is pathological in C1 (the CW-gain minimum near 0.26–0.30 V/m) and in C4 (the gain maximum near 0.3 V/m). Ratio-0.6 points are run and reported but **excluded from the classification**. The prior-pass classification, which includes 0.6, is reported as well.
  - This rule was defined knowing the CW tables (pre-existing calibration data), but not any matched-state X value except the one point listed in Section 0.
  - Note, not a pathology: the local CW slope is negative at C1's reference state and nearby (−0.93 at ratio 0.358, −1.19 at 0.4, −1.01 at CW gain 0.4). The archived reference level therefore lies beyond the CW output maximum. This is reported and carried into the claim matrix and the reviewer attack.
- **I-C, same absolute amplitude:**
  - Among C1–C4 (all at A_LO = 0.5 V/m), same absolute amplitude ≡ same signal/LO, so I-B already covers it.
  - The one additional comparison is C1 vs C5 vs C6 at a = 0.179 V/m (C1's reference level; ratios 0.42 and 0.51, neither pathological).
  - The point a = 0.1 V/m (ratios 0.235 and 0.286) is also run, because it costs about 1 min. Both are descriptive and do not enter the class.
- **Runs:**
  - The 21 prior points get seeds 20260704–06 (63 runs).
  - The 26 new points get CW plus six seeds each (182 runs).
  - Total: 245 runs, about 5 s each.
  - Outputs: `01_matching/runs/{Ck}_{rule}_{target}_{CW|REF_s<seed>}.npz`. The prior 3-seed files are read in place.
- **Pairs.**
  - A state s is *informative* if C1 shows the effect at s (six sequences).
  - A pair (s, c) is *evaluable* if c was run at s and s is not pathological.
  - E_A and E_B are the evaluable informative pairs under I-A and I-B.
  - I-C pairs are descriptive.
- **FAIR MATCHING class** (the request's definitions, operationalized):
  - **SURVIVES:**
    - HA1 dependence holds at **every** pair of E_A ∪ E_B with c ∈ {C2, C3, C4}, with at least one such pair in each of E_A and E_B;
    - **and** HA2 holds at every pair of E_A ∪ E_B with c ∈ {C5, C6}.
    - The qualitative configuration differences then remain under both conventions.
  - **DOES NOT SURVIVE:** HA1 dependence holds at **no** pair of E_A ∪ E_B with c ∈ {C2, C3, C4}, so the headline ordering disappears or reverses. This includes the case where no such pair is informative.
  - **PARTIALLY SURVIVES:** every other case, for example:
    - dependence at some but not all HA1 pairs;
    - no evaluable informative HA1 pair under one convention;
    - HA2 failing at some pair (a material change in magnitude or order).
- **Robustness reports** (they do not change the class):
  - the class with X_coh in place of X;
  - the prior-pass rule on the seeds 1–3 data, including ratio 0.6.
- **I-D rule.** If DOES NOT SURVIVE, the configuration-dependent headline is marked invalid immediately (stop rule S1). No experiment is run to rescue it.
- **Outputs:** `01_matching/matching_results.csv`, `matching_summary.md`, `fig_matching_gain.png`, `fig_matching_signal_lo.png`.

## B. Resonance-scaling hypotheses (Part II)

**HB (the request's question).** The 5 MHz phase effect is organized by an LO-related frequency scale. The alternative is generic configuration dependence.

**Matched state.** Matched signal/LO = 0.358, the same rule at every IF and LO field.
- Reason: it is the only Part I convention feasible at every IF and LO field inside the regime. The Part I plan shows that matched CW gain is infeasible for most configurations, and E1dB exceeds the LO at 10 MHz.
- It is not 2.45 E1dB: at A_LO = 0.35 and 0.425 V/m the amplitudes (0.125 and 0.152 V/m) differ from 2.45 E1dB.
- This choice is fixed by feasibility, not by any Part I result.

**Grid.**
- LO fields: A_LO = 0.35, 0.425 and 0.5 V/m. LO Rabi frequencies from the model's dipole: Ω_LO/2π = 6.456, 7.839 and 9.222 MHz.
- Coarse IFs: the 29 prior-pass IFs (1000/P MHz, P = 400 … 80 ns; 2.5–12.5 MHz; spacing 0.17–0.74 MHz) plus two new ones, P = 500 and 450 ns (2.0 and 2.22 MHz), so the scan covers 2–12.5 MHz.
- At each point: CW plus one REF sequence (seed 20260701), as in the prior pass.
- New runs: 12, in `02_resonance/runs/LO{lo}/`.

**Refinement, once and only if:**
- It is triggered for an LO field whose X_coh maximum satisfies all of:
  - X_coh,max ≥ 0.2;
  - it is not at a scan edge;
  - an adjacent grid spacing exceeds 0.2 MHz.
- Then up to 4 extra IFs are added (the integer-ns periods nearest to ±0.1 and ±0.2 MHz from the maximum), CW plus the same sequence.
- No other refinement is done.

**Uncertainty.** The sequence-to-sequence SD of X, X_coh and X_mag at 5 MHz and ratio 0.358 comes from the six-sequence Part I data (C1, C5, C6). It is used as the approximate uncertainty of single-sequence scan values.

**Per-point measures:**
- the large-signal CW gain g_CW;
- local frequency selectivity S(f) = max over |δ| ≤ 1.31 MHz (inside the scan) of |ln g_CW(f+δ) − ln g_CW(f)|, from the CW scan;
- X, X_coh, X_mag, X_rms;
- mean phase lag and phase jitter (as in Part I).

**Normalizations (no fitted parameters):**
- r1 = f/(Ω_LO/2π);
- r2 = f/(Ω_LO/4π);
- d1 = f − Ω_LO/2π;
- d2 = f − Ω_LO/4π;
- one further Hamiltonian scale: the dressed splitting of the coupling + RF ladder, Ω₃ = (Ω_c² + Ω_LO²)^½, as s3 = f/(Ω₃/4π) and d3 = f − Ω₃/4π.

No other scales are tried.

**Features** per LO field, located on the raw grid:
- the minimum of g_CW;
- the maximum of S;
- the maximum of X_coh (the phase-effect feature).

Features at the first or last scan point are not eligible.

**Alignment and collapse:**
- *Aligned* in a normalization N:
  - ratio-type N (r1, r2, s3): the across-LO spread max/min of the location is ≤ 1.08 while the spread in absolute f is ≥ 1.25;
  - difference-type N (d1, d2, d3): the across-LO range of the location is ≤ 0.4 MHz while the range in absolute f is ≥ 1.0 MHz.
- *Collapse ratio* of a quantity q ∈ {ln g_CW, X_coh}: q is interpolated on a common abscissa over the overlapping range of the three LO curves. The ratio is the RMS across-LO SD in N divided by the same in absolute f.

**RESONANCE class:**
- **CLEAR SCALING:** some N has the X_coh-maximum aligned **and** at least one CW feature (g_CW minimum or S maximum) aligned **and** collapse ratios ≤ 0.5 for both ln g_CW and X_coh. The LO curves then align around the same normalized scale and the phase effect tracks that feature.
- **SUGGESTIVE ONLY:** not CLEAR, and some N has at least one aligned feature or a collapse ratio ≤ 0.75 for ln g_CW or X_coh.
- **NO SIMPLE SCALING:** otherwise. The search then stops; no further scales are tried.
- The prior-pass P2 rule (r1/r2, X) is also reported.

**II-C.** The archived transient oscillation frequencies are compared descriptively with a set of candidate frequencies:
- Archived frequencies: OP1, periods 1.60–1.62 µs (power steps) and 1.63–1.74 µs (phase steps); OP2, 0.560 and 0.558 µs.
- Candidates:
  - Ω_LO/2π and Ω_LO/4π;
  - the two-level dressed splittings of each drive from the Hamiltonian parameters (Ω_LO, Ω_c, Ω_p over 2π);
  - Ω₃/2π and Ω₃/4π;
  - the IF detunings |f_IF − Ω_LO/2π|, |f_IF − Ω_LO/4π| and |2f_IF − Ω_LO/2π|.
- A match is a candidate within 10% of the observed frequency.
- It is called *consistent* only if the same candidate matches at both LO fields with observations. No match is called causal.

**Stated caveat.** The RF detuning of the model is Δ_l/2π = 10 Hz (resonant LO). All scales are Hamiltonian (coherent) scales; decay rates are not used.

**Outputs:** `02_resonance/resonance_scan.csv`, `resonance_summary.md`, `fig_if_scan.png`, `fig_scaled_if.png`, `fig_phase_effect_vs_scaled_if.png`.

## C. Memory-model evaluation protocol (Part III)

**Model.** Exactly one standard model: the generalized memory polynomial (GMP; Morgan et al. 2006) on the complex baseband envelope, in the reference configuration (OP1, 5 MHz). The structure is that of the prior pass (`tv_p3_gmp.py`):
- **Aligned terms** x(n−m)|x(n−m)|^k, for k = 0 … K−1 and m = 0 … M−1.
- **Lagging cross terms** x(n)|x(n−l)|^k, for k = 1, 2 and l = 1 … M.
- **Harmonic-zone terms** for the zones that fold into the 5 MHz receive band, for m = 0 … M−1:
  - zone 0: |x(n−m)|^k e^{−jωt}, k = 1 … 4;
  - zone 2: x(n−m)²|x(n−m)|^k e^{+jωt}, k = 0 … 2.
- **Constant terms:** a constant and the zone-0 tone.
- **Input:** the 1 GHz envelope through an 8 MHz zero-phase 8th-order Butterworth, decimated to 20 MHz and normalized by E1dB.
- **Pipeline:** every basis column passes through the receiver low-pass (H20 FIR) before a weighted least-squares fit (normal equations, column scaling, relative ridge).
- **Parameter count:** (K + 9)·M + 2 complex coefficients, for example 62 at K = 3, M = 5 and 1282 at K = 7, M = 80. The count is reported for every grid point, and twice that for real parameters.

**Grid.** Nonlinear order K ∈ {3, 5, 7}; memory depth D ∈ {0.25, 0.5, 1, 2, 4} µs, i.e. M = 20 D taps ∈ {5, 10, 20, 40, 80}; relative ridge λ ∈ {1e-8, 1e-5}. That makes 30 fits. The prior pass's P3 grid is superseded and not run.

**Split, fixed before fitting.** No realization or sequence appears in two splits:

| Split | Fair-waveform set | Dwell / shuffle | Phase cases |
|---|---|---|---|
| Train | r0–r1 (all formats and powers; prior-pass P4 basebands) | r0–r1 (stored revision P1 traces) | seeds 20260705–06, cases REF, FAST, SLOW, JUMP, LONGRAMP at a_H, plus the CW run (revision P2) |
| Validation | r2–r3 | r2–r3 | seed 20260704 |
| Test | r4–r7 | r4–r7 | seeds 20260701–03 |

- **Held-out drive levels.** This is the request's "operating point not used in training", interpreted as a drive level at OP1: revision P3 OP1 REF, seeds 1–3, at −6, −3, 0, +3 and +6 dB re P1dB, with the CW run at each level. The +8 dB level (2.51 E1dB) is almost the training level a_H and is reported separately. Transfer to another LO field would need re-identification, so it is not a behavioral-model test and is not done.
- **Weighting.** Each family (fair, dwell, shuffle, phase) has equal weight in the fit.
- **Selection.** The minimum of the family-averaged validation NMSE. The selected model is trained on the training split only. The test and held-out-level sets are never used for fitting or selection.

**Comparators** (on the same receiver-chain baseband): all-zone static, all-zone LTI+static and single-slow-state (τ_ref = 2.535 µs). The linear model serves as reference.

**Test metrics:**
- NMSE vs the full model;
- fitted-gain error, via the dwell gain contrast (ξ = 4 minus 0.05, +3 dB);
- D_ref error, via the declustering contrast (block-shuffled minus original, +3 dB);
- X, X_coh and X_mag errors on the baseband. The window is [89.5, 109.5) µs, relative to each model's own CW output. The cases are REF, the rate variants FAST and SLOW, and the ramp variants JUMP and LONGRAMP.
- AIR error for the fair waveforms at +3 dB, archived noise draw and level.

**MEMORY MODEL class** (the request's definitions, operationalized):
- A key contrast is *reproduced* by a model under the following criteria:
  - dwell gain contrast and declustering D_ref contrast: same sign and |model − full| ≤ 0.25 |full| (means over the test realizations);
  - phase REF: |X_model − X_full| ≤ 0.25 X_full **and** |X_coh,model − X_coh,full| ≤ 0.25 X_coh,full;
  - fair AIR at +3 dB: mean |AIR_model − AIR_full| ≤ 0.05 bit for every format.
- **SUCCEEDS:** the GMP reproduces all four key contrasts **and** its median test NMSE is ≤ 0.5 × the best comparator's in every family. It then captures the major held-out effects within a practically small error and substantially outperforms all surrogates.
- **FAILS:** no material improvement, i.e. neither of the following holds:
  - (i) median test NMSE ≤ 0.8 × the best comparator's in at least two families;
  - (ii) the GMP reproduces at least one key contrast that no comparator reproduces.
- **PARTIAL SUCCESS:** otherwise. The summary names exactly which key contrasts and families remain uncaptured.
- The prior-pass P3 rule, applied to the same fits, is reported as well.

**Outputs:** `03_memory_model/gmp_config.csv` (grid, parameter counts, validation NMSE, selection), `gmp_results.csv`, `gmp_summary.md`, `fig_model_comparison.png`.

## D. Noise-sensitivity criteria (Part IV)

**Method.**
- The receiver is linear up to the modulation-specific step. The AIR at noise scale c is therefore that of the stored noiseless baseband + c × the baseband of the archived noise draw (`rng(seed + 9991)`).
- The draw is identical across models and scales (paired).
- The symbol partition is unchanged: the same `split_eval` with the archived CE timing offset.
- Scales: c ∈ {0.25, 0.5, 1, 2, 4}. All five are used because the computation is trivial.
- **Validation.** At c = 1, every model's AIR must reproduce the archived AIR_native:
  - to ≤ 1e-9 bit for the surrogate and linear models;
  - to ≤ 1e-3 bit for the full model (GPU rerun).
  - The archived loss, excess and threshold tables must also be reproduced.
  - Any exceedance is recorded as a deviation.

**Reported per scale:**
- AIR by format and power;
- loss relative to each realization's −6 dB AIR at the same c, in bit and relative;
- surrogate − full AIR for all-zone static, all-zone LTI+static and single-slow-state at 0, +3 and +6 dB (paired t-CIs over 8 realizations);
- the format ordering at +3 dB;
- P5 and P10, where the loss crosses inside the power grid (bootstrap CIs);
- effective SNR = 10 log10(mean|bb_linear|² / mean|c·bb_noise|²) over the analysis window, per format and power;
- the baseline AIR, i.e. the full-model AIR at −6 dB.

**Informative formats.** A format is *near ceiling/floor* at scale c if its full-model AIR at both −6 and +3 dB lies within 0.05·log2(M) of log2(M), or below 0.05·log2(M). Orderings use only the other formats. A scale with fewer than two informative formats is excluded from the ordering test (documented).

**Tests at each scale:**
- **Surrogate-failure conclusion** (the central Section III conclusion; it holds when (a) and (b) both hold):
  - (a) each CW-matched surrogate (all-zone static; all-zone LTI+static) mispredicts at least one format-power cell (0, +3, +6 dB) by ≥ 0.1 bit with a 95% CI excluding zero;
  - (b) errors of both signs occur among cells whose CIs exclude zero.
- **Broad format ordering:** among informative formats, the most-affected and least-affected formats by full-model loss (bit) at +3 dB are the same as at c = 1 (16-QAM and OFDM in the archived table).
- **Detailed ranking:** the full order of the informative formats.

**NOISE class:**
- **ROBUST:** (a) and (b) hold at all five scales and the broad ordering is unchanged at every scale where it is testable.
- **PARTIALLY ROBUST:** (a) and (b) hold at 0.5, 1 and 2, but the broad ordering changes at some scale, or (a)/(b) fails only at 0.25 or 4. Detailed-ranking changes alone are reported and do not demote ROBUST.
- **NOISE DEPENDENT:** (a) or (b) fails at 0.5 or 2.
- The prior-pass P4 rule (scales 0.25–2) is also reported.

**Outputs:** `04_noise/noise_sensitivity.csv`, `noise_summary.md`, `fig_air_vs_noise.png`, `fig_surrogate_error_vs_noise.png`.

## Part V. Numerical checks (low cost)

**V-A, trace and positivity.** Duration dependence uses prefix runs: the final state of a prefix run is the state at that time of the full run, because the field prefix is identical and the simulator is deterministic.
- **Phase case:** C1 REF seed 20260701 at a_H. CW for 50 µs, then modulation; total durations 55, 110, 165 and 220 µs (4 runs).
- **Long-dwell case:** dwell ξ = 4, realization r0, +3 dB (840 µs with cyclic prefix), truncated at 105, 210, 420 and 840 µs (4 runs).
- **Reported:** max |Tr ρ − 1| and the minimum Hermitian eigenvalue, per velocity class and for the thermal average.
- **Free diagnostics:** the same quantities for every campaign run and every prior-pass P1/P2/P4 run.
- **PASS, all at every duration:**
  - thermal-average |Tr − 1| ≤ 1e-4 and thermal-average minimum eigenvalue ≥ −1e-4;
  - per-class |Tr − 1| ≤ 1e-3 and per-class minimum eigenvalue ≥ −1e-3;
  - no accumulation, where *accumulation* means the per-class max trace deviation at the longest duration is > 3× that at the shortest and > 1e-4.

**V-B, E1dB slope sensitivity** (tables only, all six configurations):
- The archived rule takes the small-signal slope from the 2 smallest nonzero amplitudes. Variants use the 3 and 4 smallest, and a local fit: ln(|harm|/a) = ln s₀ + κa² by least squares over the 4 smallest, with s₀ the slope.
- Reported: E1dB in V/m, and the change in dB (power, 20 log10 of the amplitude ratio); the implied shift of the I-A matched amplitudes.
- *Material:* a change of more than 0.5 dB in C1's E1dB, or of more than 5% in any matched I-A amplitude used in a classification. The axis is not redefined unless the change is material. A material change is a CAVEAT and is reported; it does not rerun anything.

**V-C, RWA regime check** (analysis only):
- Inspect how the simulator enters the RF field: a complex envelope in the frame rotating at the LO carrier.
- For every state behind a manuscript or campaign comparison, tabulate:
  - signal/LO (peak |envelope|/A_LO for the fair waveforms);
  - the peak total RF Rabi frequency (A_LO + a_peak)·μ/h;
  - its ratio to the model's RF carrier. The carrier is the Cs 47D5/2 → 48P3/2 transition, 6.946 GHz, from the model's source configuration; it is not a simulator input.
- Flags:
  - **OUTSIDE INTENDED REGIME** where signal/LO > 1: the superheterodyne receiver regime of the inherited model, which assumes the LO dominates;
  - **RWA QUESTIONABLE** where the Rabi-to-carrier ratio > 0.01.
- Not quantified: the four-level truncation (other Rydberg transitions). It is stated as a limitation, with no invented numbers.

**NUMERICAL SANITY class:**
- **PASS** requires all of:
  - V-A passes;
  - V-B is not material;
  - no comparison used by a KEEP or WEAKEN claim of the final claim matrix involves a flagged state;
  - all free diagnostics meet the V-A thresholds.
- **CAVEAT** otherwise, with the caveats listed.

## E. Global decision rules

**Stop rules** (the request's four, operationalized). Stop immediately and summarize if:
- **S1.** The Part I class is DOES NOT SURVIVE.
- **S2.** The inherited simulator is clearly outside its valid regime for a required comparison, meaning RWA QUESTIONABLE with a ratio > 0.05 for any state entering a Part I–IV classification.
  - signal/LO > 1 alone is outside the receiver's intended regime but not outside the simulator's validity. It is flagged, not a stop.
  - Part I–II states are ≤ 1 by construction.
- **S3.** The Part III class is SUCCEEDS **and** the GMP reproduces X and X_coh within 25% at every held-out drive level where X_full ≥ 0.05. The strong memory model then completely explains the dynamics.
- **S4.** Any run entering a classification has thermal-average |Tr ρ − 1| > 1e-2 or minimum eigenvalue < −1e-2, or V-A shows accumulation reaching these values.

**On a stop:**
- no further simulation;
- the triggering analysis is completed;
- the claim matrix (Part VII) is written with the invalidated claims marked;
- the manuscript is **not** revised (a rescope needs the author's approval);
- the final report carries all Part X headings, with unexecuted parts marked NOT RUN (stop rule).

**Order of execution:**
1. Part I runs, then the Part I class (S1 check).
2. Part II runs (12 + possible refinement), then the Part II class. Part III fitting runs on the CPU concurrently, after the prior-pass P4 data are complete (S3 check).
3. Part IV analysis.
4. Part V runs and analyses (S2, S4 checks).
5. Parts VI–X.

The analysis of each part runs only after that part's data are complete.

**Story selection (Part VI), from the classes only:**
- **A (resonance-conditioned):** Fair = SURVIVES and Resonance = CLEAR SCALING.
- **B (engineering memory modeling):** not A, and all of:
  - Fair ∈ {SURVIVES, PARTIALLY SURVIVES} (matching remains meaningful);
  - Memory ∈ {SUCCEEDS, PARTIAL SUCCESS};
  - Resonance ≠ CLEAR SCALING.
- **C (cautionary / model transfer):** Fair ∈ {PARTIALLY SURVIVES, DOES NOT SURVIVE} and Resonance ≠ CLEAR SCALING, and the claim matrix's KEEP claims are convergence, temporal-order and model-transfer results.
- If no story's conditions all hold, choose the story with the largest fraction of satisfied conditions; ties go to C.

**Claim matrix (Part VII)**, for every candidate claim:
- Columns: exact wording; supporting result; fair-matching status; noise robustness; baseline robustness (GMP); generalization status; prespecified or post hoc; DIRECT / INFERRED / INTERPRETIVE; decision.
- **KEEP:** an explicit comparison with a 95% CI (or a prespecified classification) that survives every campaign test relevant to it, inside the intended regime.
- **WEAKEN:** supported only in a named, restricted domain. Examples:
  - PARTIALLY SURVIVES, SUGGESTIVE ONLY, PARTIALLY ROBUST or PARTIAL SUCCESS limits;
  - reliance on archived comparisons outside the intended regime.
- **REMOVE:** contradicted by a campaign test, or resting only on comparisons the campaign invalidated.
- Only KEEP claims may appear in the abstract. Failed prespecified tests stay reported.
- **Secondary Part I reading for the manuscript's operating-point claim** ("X smaller at lower LO fields"):
  - KEEP if X_C1 − X_C6 ≥ 0.05 (CI excluding zero) at ≥ 2 of the 3 matched-CW-gain targets;
  - REMOVE if X_C6 − X_C1 ≥ 0.05 (CI excluding zero) at ≥ 2 of them;
  - WEAKEN otherwise.

**Manuscript (Part VIII).** Revised only after the story is selected and only if no stop rule fired. Scope: title, abstract, contributions, Sections III, IV and VII, Discussion, Limitations, Conclusion, the evidence table, and figures/captions as needed. Unaffected methods and results are not rewritten, and failed tests are preserved.

**Reviewer attack (Part IX).** The six questions of the request, no requests for new experiments. Each comment is classed as:
- **BLOCKING:** it invalidates a retained claim or makes the paper unfit for the venue as written, and cannot be fixed by wording;
- **IMPORTANT:** it must be addressed in the text before submission, but the retained claims stand;
- **MINOR:** presentation.

**FINAL PAPER class:**
- **TQE STRONG SUBMIT:** no stop rule; Fair = SURVIVES; Noise = ROBUST; Numerical = PASS; story A or B; no BLOCKING comment.
- **TQE BORDERLINE SUBMIT:** no stop rule; Fair ∈ {SURVIVES, PARTIALLY SURVIVES}; Noise ∈ {ROBUST, PARTIALLY ROBUST}; story A or B; no BLOCKING comment.
- **BETTER TO RETARGET:** otherwise, including story C, any stop rule, or any BLOCKING comment.

## Step 2. Progress reporting and reproducibility

**Progress.** After every job, the runner (`python/experiments/tqe_viability_final/tvf_run.py`):
- rewrites `00_progress.txt` with [done/total] %, stage, case, elapsed time, the rolling s/job (last 10 jobs), ETA and the latest result;
- appends to `00_progress.csv` the columns timestamp, stage, job_id, config, matching_rule, parameter_value, elapsed_s, avg_job_s, eta_s, status.

A watchdog logs a job that exceeds 5× its stage's median runtime, together with process and GPU state. It never kills anything.

**Reproducibility.**
- Every new output gets a row in `00_run_manifest.csv`: config, seed, command, git commit and dirty flag, runtime, source path, output path and SHA256.
- Analyses record their inputs and code version in their outputs.
- The reproducibility package is updated at the end (SHA256 manifest, clean-room reproduction, figure/number/citation audits) and is not published.
- Deviations from this plan go to `00_deviations.md` with a time stamp and the state at the time of writing.
