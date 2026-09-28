# Stage-0.6 Decision

**CONDITIONAL GO.**

Two waveforms with exactly the same amplitude samples can produce different nonlinear distortion in the full thermal atomic receiver. Only the time order of the samples differs; the RMS, peak, high-field occupancy and amplitude histogram are identical by construction. The difference is reproducible and strongly significant. It is 100–13,000× larger than the numerical error, weak in the linear regime, and not reproduced by any tested static or LTI+static control. That answers the central question **yes**.

The effect is organized on the atomic timescale in one specific sense. The full-model effect keeps growing across T_dwell/τ_atom ≈ 0.25–4, while every control saturates at T_dwell/τ_atom ≈ 0.5.

It is **not** a clean dwell law, for four reasons:
- there is no interior maximum near T_dwell/τ_atom ~ 1;
- the quasi-static limit is not reached by ξ = 4;
- the pattern changes character at +6 dB;
- the symbol-level AIR metric is confounded by the linear receiver chain.

That places the result under CONDITIONAL GO, not GO.

Three process changes relative to the brief:
1. **Nd.** Nd = 1501 failed the pre-registered convergence test. All decisive runs use **Nd = 4001**, which passed, with 8001 as the check.
2. **Controls.** The Stage-0.5 control models were incomplete. Stage-0.6 adds complete all-zone versions (see "Static / LTI Controls").
3. **Stage-0.5 impact.** Both findings materially change how Stage-0.5 should be read (last section).

---

# Nd Convergence

The tolerances were fixed before any results (`00_preregistered_criteria.json`): AIR < 0.005 bit/native symbol, baseband NMSE < 10⁻³, output RMS < 0.5%. The test used four Stage-0.5 cases plus the Stage-0.5 Stage-K pair (`02_nd_convergence.csv`, `fig01`):
- A: QPSK at −10 dB
- B: OFDM at 0 dB
- C: 16-QAM at +3 dB
- D: QPSK at +6 dB, the strongest full-minus-control case
- E: the Stage-0.5 Stage-K dwell pair

Results against Nd = 8001 (worst case over A–D):

| Nd | worst \|ΔAIR\| | worst NMSE | worst ΔRMS | Verdict |
|---:|---:|---:|---:|---|
| 1501 | 0.019 | 8.9e-3 | 0.57% | **FAIL** (also fails against the pre-registered Nd = 3001 reference) |
| 2001 | 0.012 | 2.2e-3 | 0.08% | fail |
| 3001 | 0.0039 | 1.8e-4 | 0.02% | pass (small AIR margin) |
| **4001** | **0.0015** | **1.8e-5** | **0.015%** | **pass; used for all decisive runs** |

- The NMSE falls roughly as Nd⁻⁴, which is regular quadrature convergence. The drift is systematic, not noise: the Stage-K pair contrast moves monotonically, 1.106 → 1.078 → 1.057 → 1.052 ×10⁻³.
- The pre-registered rule said to STOP if 1501 failed. I recorded the failure and did not use Nd = 1501 for any decisive claim. Because a converged Nd was cheap on this GPU (~30 s per 800 µs run), I continued at Nd = 4001 rather than halting. This is a documented deviation from "use Nd = 1501".
- **CW quantities are far more biased than waveform runs at Nd = 1501.** The converged E1dB is **0.0732 V/m** (8001: 0.0727), against **0.0878** at 1501. That shifts the power axis by ≈1.6 dB. The small-signal CW gain at 1501 is 9.5% low, and the compression-curve shape is off by up to ±10%. A plausible cause is that the 1501-grid two-photon Doppler comb recurs every ≈1.98 µs, which is ≈10 IF periods, so a steady 5 MHz tone accumulates the aliasing error. This is a hypothesis, not tested.
- Stage-0.6 therefore rebuilt E1dB, the CW compression table and the LTI kernel at Nd = 4001. All power levels here are relative to the converged E1dB.

# Atomic Timescale

Full thermal step responses were run at Nd = 1501, 4001 and 8001 (`03_tau_atom_results.csv`, `fig02`). The definitions are explicit in `tau_atom.py`.

| Experiment | τ_1e | τ_90 | tail τ_dom | Notes |
|---|---:|---:|---:|---|
| IF-envelope step +5% at E = 0.0877 V/m ("P1dB") | **2.535 µs** | 7.04 µs | ≈4.0 µs | multi-exponential and oscillatory; slow part ~50% weight |
| same, step down | 2.530 µs | 7.02 µs | 4.2 µs | symmetric |
| small LO step ±1% | 1.94 µs | 5.47 µs | 2.77 µs | clean single-exponential tail (log-residual 0.008) |
| IF-envelope step at the linear point (0.05 → 0.10·E1) | 0.047 µs | 0.074 µs | – | linear transconductance is wideband |
| recovery after a 30 µs plateau at 1.5·E1 | 0.94 µs | 6.3 µs | 4.6 µs | undershoot, see below |

- **τ_atom ≡ 2.535 µs**, the τ_1e of the P1dB IF-envelope step at Nd = 4001. It was fixed in `stage06_config.json` before any dwell run and is stable to 0.2% across Nd.
- At Nd = 1501 the linear-point step shows a spurious 16% slow tail (τ_90 = 1.66 µs) that disappears at Nd ≥ 4001. This is another 1501 artifact. The inherited Stage-1.5 value of 1.727 µs is not reproduced by any definition here.
- **Caveat.** The "P1dB" step amplitude was the Stage-0.5 E1dB (0.0877 V/m), which is 1.20× the converged E1dB. The planned sensitivity run at the converged E1dB was cancelled under the "no further sweeps" instruction.
- **Clean dynamic signatures.** These come from the unmodulated CW plateau, IF-demodulated with a boxcar, and are stable across Nd:
  - **slow compression onset:** the output overshoots its final compressed level by up to 44% and settles over ≈2–5 µs;
  - **post-excursion gain suppression:** after the field drops back, the output is 31% below its eventual low-level value at 1 µs, 16% below at 5 µs and 3% below at 10 µs.

# Controlled Pair Quality

Configuration is in `04_controlled_pair_configs.csv`; see `fig03`–`fig06`.

**Dwell family.** Each of 8 realizations contains 7 members:
- the reference **A** with 127 ns plateaus (T_dwell/τ = 0.05);
- members **B** with plateau lengths ξ·τ for ξ = 0.1, 0.25, 0.5, 1, 2, 4.

Construction:
- two levels in a 1:3 ratio, 25% occupancy, 10 ns raised-cosine edges;
- a constant-envelope QPSK phase carrier at 1 Msym/s, shared by every member;
- 800 µs per waveform plus a 40 µs cyclic prefix.

Every member has the same number of rise and fall edges. The surplus edge pairs in the long-dwell members become brief spikes that carry 3.6–7.2% of the high-field time (mean over realizations; 10_verification_report.txt). The consequence is that **every member has bitwise-identical sorted samples**. The spike count is nearly constant from ξ = 0.5 to ξ = 4 (1313 → 1441, +10%) while the full-model ΔD_BLA rises from 0.108 to 0.389 at +3 dB, so the spikes do not drive the long-dwell effect.

**Shuffle set.** The original is a clustered envelope, |fast Gaussian| × a log-normal slow modulation (correlation ≈ 2τ). Three variants reuse its exact samples:
- a permutation of whole excursion cycles, cut at median up-crossings;
- a fixed-length block permutation;
- a sample-level permutation.

| Quantity | Target | Dwell family | Cycle block-shuffle |
|---|---|---|---|
| RMS, peak, energy, occupancy mismatch | <0.5%, <0.5%, –, <1 pt | **0 (exact)** | **0 (exact)** |
| KS distance / sorted samples identical | – | **0 / yes** | **0 / yes** |
| B99 mismatch vs reference | <5% | −1.0 … +5.6% (ξ = 2, 4 slightly over) | −0.15% |
| in-band (±3 MHz) power fraction | – | +6 … +9% | +0.01% |
| time-weighted T_dwell/τ | – | 0.05 → 3.8 | 0.97 → 0.97 (unchanged) |
| C_env = τ_env/τ | – | 0.03 → 2.0 | 1.20 → 0.57 |

- In the dwell family, B99 (≈34–36 MHz) is set by the 10 ns edges, far outside the receiver and atomic bandwidths. Because the ξ = 2 and 4 members exceed the 5% target, the LTI control and the purely linear model are included. The linear model's ΔD_BLA is ≤ 0.006 at every ξ, so the physics metric is insensitive to these spectral differences.
- The fixed-block and sample shuffles change the bandwidth strongly (+200% and +19,000%). They are reported but not used as causal evidence.

# Temporal Ordering Effect

**Metrics.**
- **Primary (physics):** D_BLA, the noise-free in-band distortion left after the best linear FIR equalizer (41 taps, 2 µs) maps the ideal-receiver input to the output.
- **Secondary (communication-like):** Gaussian-auxiliary AIR, EVM and SER on the embedded QPSK carrier, with the Stage-0.5 noise σ.
- Statistics are paired over 8 realizations with 95% t-intervals. A and B share the samples, the carrier and the noise draw.

**Dwell family, full model, ΔD_BLA = D(B) − D(A)** (`05_controlled_pair_results.csv`, `fig07`):

| ξ of B | −6 dB | 0 dB | +3 dB | +6 dB |
|---:|---:|---:|---:|---:|
| 0.1 | −0.024 [−0.026, −0.022] | −0.040 [−0.048, −0.032] | −0.048 [−0.061, −0.034] | −0.157 |
| 0.25 | −0.027 | −0.014 | −0.003 [−0.015, +0.009] | −0.157 |
| 0.5 | −0.032 | +0.010 | +0.108 [+0.090, +0.126] | −0.034 |
| 1 | −0.032 | +0.074 [+0.057, +0.090] | +0.234 [+0.199, +0.269] | +0.056 |
| 2 | −0.029 | +0.139 | +0.306 | +0.168 |
| 4 | −0.026 [−0.030, −0.022] | **+0.178 [+0.151, +0.206]** | **+0.389 [+0.368, +0.410]** | +0.283 [+0.260, +0.307] |

**Shuffle test, cycle block-shuffle versus original** (`06_temporal_shuffle_results.csv`, `fig13`). This keeps the samples, every excursion and B99, and removes only the clustering beyond about τ.

| | ΔD_BLA full | full − LTI+static (all-zone) |
|---|---:|---:|
| −6 dB | −0.034 [−0.079, +0.011] n.s. | n.s. |
| 0 dB | −0.015 n.s. | −0.033 [−0.068, +0.002] |
| +3 dB | **−0.075 [−0.106, −0.044]** | **−0.102 [−0.145, −0.060]** |
| +6 dB | **−0.124 [−0.173, −0.076]** | **−0.149 [−0.210, −0.087]** |

Declustering lowers full-atomic distortion by 10–15% in the nonlinear regime. The memoryless controls predict no change: for the all-zone static model ΔD = +0.004 (n.s.) at +3 dB.

**Traces** (`fig12`). At long dwell, only the full model shows the two dynamic signatures:
- **slow compression onset:** each plateau starts near the linear gain and compresses over microseconds;
- **post-excursion suppression:** after the plateau the gain drops well below its steady low-level value and recovers over ≈10–15 µs.

At short dwell, the full model does not follow individual 127 ns excursions. It holds an averaged, depressed gain (≈0.47 of linear at +3 dB), while the static models follow each excursion.

# Static / LTI Controls

**Critical finding.** The Stage-0.5 controls (labeled F1) model only the IF fundamental zone. With a 5 MHz IF and a MHz-wide envelope, two other intensity components land inside the ±3 MHz receive band:
- the atom's rectified baseline c₀(|E|), which is **as large as or larger than the fundamental** (|c₀/c₁| ≈ 0.5 near E1dB and up to 2.8 at 3·E1dB);
- the 2nd-harmonic zone.

So I built complete controls (`models06.py`) from all-zone CW tables (k = 0…8) at Nd = 4001:
- **M_STATIC (all-zone):** memoryless; reproduces the full CW periodic steady state.
- **M_LTI_STATIC (all-zone):** the Stage-0.5 static→LTI (Hammerstein) convention for the fundamental; static c₀ followed by the small-signal baseline kernel from the LO step; static k ≥ 2 zones.

Both reproduce off-grid CW tones from the full model to 0.06–1.7% (`07b`); the F1 controls miss 4–67%. All four controls are reported; the all-zone versions are decisive.

**Dwell family.** Every control's ΔD_BLA(ξ) **saturates by ξ ≈ 0.5**. At 0 / +3 dB the LTI+static values are +0.016 / +0.094 and the F1 values +0.018 / +0.12 (`fig08`, `fig09`). The full model keeps rising to +0.18 / +0.39.

Full minus LTI+static (all-zone), paired:

| ξ | 0 dB | +3 dB |
|---:|---:|---:|
| 1 | +0.058 [+0.041, +0.074] | +0.142 [+0.107, +0.176] |
| 4 | **+0.162 [+0.135, +0.190]** | **+0.294 [+0.274, +0.315]** |

Against the all-zone static model the full-minus-control difference is even larger (+1.6 to +1.8 at +3 dB), because that control moves the opposite way: its rectified-zone leakage is largest for fast envelopes.

**+6 dB is different.** The controls predict *more* dwell sensitivity than the full model (LTI+static ΔD +0.50 at ξ = 4 against full +0.28). The full model still differs from every control, but in the opposite direction.

# Dwell / Tau Scaling

See `fig10`, `fig11` and `05b_per_waveform_dynamic_excess.csv`.

- **Organized by τ_atom.** At 0 and +3 dB the full-model ΔD changes most steeply between ξ ≈ 0.25 and 2. That is T_dwell ≈ 0.6–5 µs, which brackets τ_1e = 2.5 µs and τ_90 = 7 µs. The controls stop changing at ξ ≈ 0.5, i.e. 1.3 µs. So the extra dwell sensitivity lives on the atomic relaxation scale, not on the LTI (~50 ns) or receiver-filter scale.
- **No interior maximum at ξ ~ 1.** The per-waveform dynamic excess in D (full minus LTI+static) *increases monotonically* through ξ = 4 at 0 and +3 dB. At +6 dB it is non-monotone: +0.37 at ξ = 0.05, −0.11 at ξ = 0.25–0.5, +0.15 at ξ = 4.
- **C_env.** In the shuffle set, reducing C_env from 1.20 to 0.57 at fixed T_dwell reduces distortion only in the full model. Clustering beyond the excursion scale matters, which is consistent with a slow atomic state integrating recent high-field exposure.
- **AIR is not a clean dwell diagnostic here.** Even the purely linear model's AIR changes by +0.24 to +0.45 bit between A and B, because of within-symbol envelope variation through the receiver filter. The per-waveform AIR excess falls with ξ, the opposite of D. I checked one explanation, lower average gain costing SNR at short dwell, and it is rejected: the full model's fitted gain falls with ξ (0.47 → 0.36 of linear at +3 dB). The AIR trend is left unexplained and is not used for mechanism claims.
- No fitted law is proposed.

# Quasi-Static Limit (Step 12)

- **Short dwell (ξ ≪ 1): supported.** The atom averages rapid excursions and holds a roughly constant, mean-power-set gain.
- **Long dwell (ξ ≫ 1): not reached.** At ξ = 4 the full model is *further* from the CW-static prediction. Plateau-interior gain relative to the linear fit (realization 0, +3 dB, `10_verification_report.txt`): full 0.57 in the first µs, **0.21** in the second half of the plateau, and **0.29** in the low level 0.25–2 µs after it; all-zone static 0.48 / 0.44 / 0.93. (The whole-waveform fitted gains are 0.36 for full and 0.56 for static.) In the CW-carrier Step-2 experiment, by contrast, the plateau output does converge to the CW value.
- Two things remain open:
  - whether ξ = 4 is simply short of the limit (a 10 µs plateau is only ≈1.4 τ_90);
  - whether the QPSK phase carrier's transitions deepen compression beyond CW, which would be a separate dynamic effect.
- The planned CW-carrier ablation that would separate these was **not run** (cancelled pending your review).

# Power Dependence

See `fig14`.
- **−6 dB:** full ΔD is small (|ΔD| ≤ 0.03) with no ξ trend, and the shuffle effect is n.s. The effect is weak in the near-linear regime, which is the expected behavior for a nonlinear mechanism.
- **0 and +3 dB:** a strong monotone dwell effect that grows with power.
- **+6 dB:** the character changes. At a_H ≈ 3.5·E1dB the CW response enters its non-monotone, collapsing-fundamental region, and the controls then overshoot the full model's dwell sensitivity.

# Numerical Stability

`07_numerical_validation.csv` checks realizations 0–1 at +3 dB with the decisive dwell and shuffle variants.

| Check | Max change in per-waveform D_BLA | Max change in AIR | Max change in pair ΔD |
|---|---:|---:|---:|
| Nd 4001 → 8001 | ≤ 5.2e-4 | ≤ 2.7e-3 | ≤ 9e-4 |
| dt 1 → 0.5 ns | ≤ 3.7e-4 | ≤ 5e-4 | ≤ 4e-4 |

Every non-null effect exceeds the numerical difference by **116× to 13,800×**. The Stage-0.5 results also reproduce on the new Linux/CUDA-13 build (full ΔAIR ≤ 7.5e-7; `00_stage05_reproducibility_check.csv`).

# What Is Supported

- Temporal ordering alone (exact same amplitude samples) changes the full thermal atomic receiver's in-band nonlinear distortion, reproducibly across 8 realizations and robust to Nd and dt.
- The effect is absent or weak near-linearly and present from P1dB upward.
- No tested static (fundamental-only or all-zone) or LTI+static (fundamental-only or all-zone) control reproduces it.
- The full model's extra dwell sensitivity accumulates over T_dwell ≈ 0.25–4 τ_atom, where all controls have saturated.
- Two dynamic atomic signatures are directly observed and are Nd-stable: slow compression onset (µs) and post-excursion gain suppression (~10 µs).

# What Is NOT Supported

- A dwell *law*: there is no ξ ~ 1 maximum, no quasi-static convergence by ξ = 4, and the pattern differs at +6 dB.
- Any claim about communication AIR ordering. At symbol level the AIR is confounded by linear receiver effects.
- That the phase carrier plays no role in the long-dwell compression (untested).
- Universality beyond this Transit configuration, LO = 0.5 V/m, IF = 5 MHz and this receiver chain; also any hardware relevance.
- τ_atom at the exact converged E1dB (the definition used 1.20·E1dB); sensitivity not run.
- That M_LTI_STATIC is the unique right control. Its baseline-zone kernel is a modeling choice. The all-zone static and all-zone LTI+static controls bracket that choice, and the conclusions hold against both.

# Scientific Interpretation

In this model, steady-state (CW) compression metrics do not characterize dynamic distortion. The atom carries a slow state (µs relaxation, with post-excursion suppression) that integrates recent high-field exposure. The nonlinear response therefore depends on how high-field time is organized relative to τ_atom, not only on the amplitude statistics. The resulting distortion cannot be captured by a memoryless map or a short-memory Hammerstein model, including one with all harmonic zones.

**Implications for Stage-0.5.**
1. Its P1dB axis was ≈1.6 dB off, because E1dB at Nd = 1501 was biased.
2. Its "full-minus-static" excess was partly an artifact of the incomplete fundamental-only controls. At +3 dB with all-zone static controls (`09_stage05_recheck.csv`), the QPSK and 16-QAM excess falls from ≈0.70 / 2.27 (F1 static) to ≈0.23 / 0.14 bit, OFDM reverses (all-zone static 0.54 vs full 1.77), and CE remains large (≈0.95).
3. The excess over the LTI+static control is barely changed by the all-zone version: 16-QAM 0.87, CE 0.65, OFDM 0.13, QPSK 0.63 bit. The full-atomic AIRs themselves change little at Nd = 4001 (e.g., QPSK 1.297 vs 1.277).

# Recommended Next Step

Held for your review; nothing has been run beyond the planned campaign. In priority order:
1. **CW-carrier ablation** of the dwell family (already scripted as `--cw-carrier`). It separates amplitude-dwell dynamics from phase-transition dynamics and tests the quasi-static limit.
2. **Longer dwell (ξ = 8–16)** to see whether the full model converges toward static.
3. **τ sensitivity at the converged E1dB** (`tau_sensitivity.py`, about 1 minute of GPU time).
4. Only then, a minimal low-order dynamic control, for example a Hammerstein model with a slow gain-state ODE whose time constant is taken from Step 2 rather than fitted. That would test whether one slow state explains the effect.
5. Revise the Stage-0.5 summary to the converged E1dB and the all-zone controls before any paper use.


---

# Follow-up (FU-A…FU-E): Results and Revised Reading

Pre-registered in `11_followup_preregistration.json` before any follow-up run; verdicts in `14_followup_verdicts.csv`; figures `fig15`–`fig18`. 111 GPU jobs plus CPU-only control evaluation, no errors.

| Item | Pre-registered outcome |
|---|---|
| FU-A τ at converged E1dB | **Rule triggered.** τ_1e = 0 / 1.13 / 2.54 µs at 0.8 / 1.0 / 1.2 × E1dB. τ_1e is ill-conditioned because the slow component's weight (0.24 / 0.38 / 0.56) crosses 1/e at P1dB. The robust slow-pole constant is 2.49 / 2.70 / 3.11 µs, consistent with τ_atom = 2.535 µs. With τ_1e(1.0·E1dB), every ξ axis would be multiplied by 2.25. |
| FU-B CW plateau → CW-static gain | **Supported** (+0.2% at 0 dB, +0.3% at +3 dB) |
| FU-B QPSK phase carrier deepens compression | **Supported**: plateau gain −60% (0 dB) and −55% (+3 dB) relative to CW-static |
| FU-C dwell effect without phase modulation | **Formally supported** (full − LTI+static = +0.350 [+0.265, +0.434] at +3 dB, ξ=4), but carried by the control's change. See below. |
| FU-D quasi-static convergence (D_BLA rule) | **Not supported for either carrier.** CW: the gap shrinks to 37% of its ξ=4 value (threshold 25%). QPSK: the gap grows. |
| FU-E single slow state (|E|²-driven, τ = 2.535 µs, not fitted) | **Explains the clustering effect** (2/2 contrasts, median error 8%). **Does not explain the dwell contrasts** (1/12). Better waveform-level predictor than LTI+static on 88–100% of waveforms, except the long-dwell QPSK family (67%). |

**Metric correction.** D_BLA divides each model's residual by that model's *own* compressed signal power. With an unmodulated carrier, compression flattens the envelope, so this denominator collapses; the all-zone static model reaches D_BLA ≈ 58. The long-dwell runs therefore also record residual power relative to a **fixed** reference, the linear model's signal power (exact for all realizations). Paired results at +3 dB, 4 realizations:

| | full − static at ξ=4 | full − static at ξ=32 |
|---|---:|---:|
| unmodulated carrier | +45% | **+4.9% [+2.9%, +7.0%]** (gain 0.26 vs 0.24) |
| QPSK phase carrier | +68% | **+79%** (gain 0.35 vs 0.56) |

**Revised physical reading.** There are two separable dynamic mechanisms:
1. **An |E|²-driven slow atomic state (τ ≈ 2.5–3 µs).** It causes averaging of fast excursions, slow compression onset, post-excursion gain suppression (≈10 µs), and the clustering effect. With an unmodulated carrier, the full model shows the canonical memory behavior. At short dwell it averages (gain 0.72 of linear at +3 dB, against 0.18 for static). At long dwell it converges to the CW-static response, reaching 4.9% at ξ = 32 on the fixed reference. The Stage-0.6 "no quasi-static limit" statement was an artifact of stopping at ξ = 4 and of the QPSK carrier.
2. **Phase-modulation-driven extra compression.** A communication-like phase carrier at high field drives the atom into a *different* steady state, 55–60% more compressed than the CW curve at the same amplitude. It is reached within ξ ≈ 4 and stays flat to ξ = 32. It is why long high-field dwell is more damaging with the QPSK carrier. Neither a CW-matched static model nor a single |E|²-driven slow state captures it.

**Revised decision: still CONDITIONAL GO, with a sharper claim.** Temporal organization matters beyond static and LTI+static controls, it survives every numerical check, and the amplitude-dwell part is now organized by τ_atom with a demonstrated quasi-static limit (unmodulated carrier). What is still missing for GO is a mechanistic account of mechanism 2. Its dependence on symbol rate and phase-transition shape is untested.

**Recommended next step (not run).** Test mechanism 2 directly: constant-amplitude plateaus at a_H with the QPSK carrier at 2–3 symbol rates and 2 phase-ramp durations, plus a pure frequency offset (no phase jumps). This would show whether the extra compression tracks the transition rate, the instantaneous frequency excursion, or the IF detuning. It costs ~10 short GPU runs.
