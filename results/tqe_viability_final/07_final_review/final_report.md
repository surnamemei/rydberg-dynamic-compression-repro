# Final adversarial-validation campaign — final report (Part X)

**Setup.**
- Branch `tqe-viability-final`. The plan `00_plan.md` was committed at 02:37 (e11049c), before any campaign run. The analysis scripts were committed at 02:47, before any result.
- No stop rule fired (S1–S4 all false).
- Nothing was pushed, published or submitted.
- Classes and rules are exactly those of the plan. Deviations: D1 and D2 in `00_deviations.md`, both implementation fixes that changed no rule and no result.

## Fair Matching

**PARTIALLY SURVIVES** (`01_matching/`).
- **Design.** Six configurations and six sequences per state, shared across configurations and paired. Two conventions, both with signal ≤ LO: matched signal/LO, and matched CW gain where reachable.
- **Headline ordering (IF and probe): holds at 11 of 11 evaluable comparisons.** At signal/LO 0.2 / 0.358 / 0.4:
  - 5 MHz: X = 0.32 / 0.49 / 0.46;
  - 10 MHz: 0.06 / 0.08 / 0.09;
  - 15 MHz: −0.01 / −0.03 / −0.04;
  - weaker probe: 0.01 / 0.04 / 0.05;
  - paired differences +0.27 to +0.52, every lower bound ≥ +0.20.
- **At matched CW gain 0.9 / 0.85** (the only compression-matchable pair): 5 MHz 0.12 / 0.23 vs 15 MHz −0.04 / −0.01.
- **Caveat.** 10 MHz and the weaker probe expand rather than compress below the LO field (E1dB = 1.16 and 1.43 × LO), so no compression-matched comparison with them is possible.
- **Why "partially".** The operating-point ordering holds at 11 of 12 comparisons. It fails at matched CW gain 0.4 for A_LO = 0.425 V/m (+0.02 [−0.05, +0.10]).
- **Robustness.** With X_coh in place of X: SURVIVES. The prior-pass rule (three sequences, including ratio 0.6) gives the same headline outcome.
- **Not a pathology artifact.** The archived reference level lies beyond the CW output maximum (local slope −0.93), but the effect is already X = 0.33 [0.26, 0.41] at CW gain 0.8 (slope +0.53).

## Resonance Hypothesis

**SUGGESTIVE ONLY** (`02_resonance/`).
- **Scan.** IF 2.0–12.5 MHz at A_LO 0.35 / 0.425 / 0.5 V/m, signal/LO 0.358, one sequence per IF, plus a prespecified refinement.
- **Locations align with the LO Rabi frequency.** The CW-gain dip and the selectivity maximum sit at r₂ = f_IF/(Ω_LO/4π) = 1.02–1.03. The X_coh maximum sits at r₂ = 1.11 at all three LO fields (3.57, 4.35 and 5.10 MHz), while the absolute IF moves 1.43-fold.
- **A second maximum** lies just below f_IF = Ω_LO/2π (post hoc).
- **But the curves do not collapse:** collapse ratios 0.77–1.12, against ≤ 0.5 required.
- **Limitations of the scan:**
  - ratio and difference scalings cannot be told apart over the 1.43-fold LO range;
  - the weaker probe sits at the same r₂ = 1.08 as the reference without the effect, so the LO scale is at most necessary;
  - at the CW-gain minimum itself X_coh ≈ 0.
- **Oscillation periods (II-C):** no candidate scale matches at both LO fields. No mechanism is claimed.

## Standard Memory Model

**PARTIAL SUCCESS** (`03_memory_model/`). The selected GMP has order 7, 4 µs memory and 1282 complex coefficients; the split was fixed before fitting.

**What it captures:**
- median test NMSE 3.0–5.5-fold below the best CW-derived surrogate for dwell, shuffle and phase signals (0.16 vs 0.50, 0.22 vs 0.66, 0.20 vs 1.11);
- only marginally better on the fair waveforms (0.049 vs 0.055);
- QPSK and 16-QAM AIR at +3 dB within 0.018 / 0.037 bit (surrogates 0.43–0.77).

**What remains uncaptured (all four key contrasts miss the tolerance):**
- **Phase effect:** X 0.50 vs 0.14, X_coh 0.54 vs 0.13, X_mag 0.45 vs 0.02. This is structural: under constant-envelope modulation every envelope-driven GMP term reduces to a linear filter.
- **CE link:** 0.80 bit AIR error.
- **Temporal-order contrasts:** correct sign, wrong size. Dwell −0.139 vs −0.106; declustering −0.0057 vs −0.0026.
- **Held-out drive levels:** the level dependence of the phase effect is missed.

The stop rule S3 (a model that explains everything) did not fire.

## AIR Noise Sensitivity

**PARTIALLY ROBUST** (`04_noise/`). Noise scales 0.25–4×, i.e. effective SNR 20–46 dB at +3 dB; the noise draw is paired, and at c = 1 the archived AIR is reproduced to 4e-13 bit.
- **Surrogate failure** (a miss ≥ 0.1 bit with CI, errors of both signs) holds at every scale, and the errors barely move: static 16-QAM at 0 dB is +1.51 / +1.49 / +1.26 bit at 0.25× / 1× / 4×.
- **Broad ordering** (16-QAM most, OFDM least affected at +3 dB) holds from 0.25× to 2×. At 4×, CE (0.07 [−0.18, 0.32]) falls below OFDM (0.08 [−0.01, 0.17]); neither is resolved.
- **P5 is noise dependent:** 16-QAM −2.47 dB at 1× vs −0.60 dB at 4×.

## Numerical Sanity (PASS/CAVEAT)

**CAVEAT** (`05_numerical/`).
- **Trace and positivity (V-A): pass.** Per-class trace deviation ≤ 2.4e-05 and minimum eigenvalue ≥ −7.6e-10 for a 220 µs phase case and an 840 µs long-dwell case, with no accumulation. The final states of all ~730 runs also pass.
- **Rotating-wave treatment: valid everywhere.** The peak RF Rabi frequency is at most 0.006 of the carrier.
- **Why CAVEAT:**
  - **E1dB rule (V-B):** the reference E1dB moves by only −0.05 to +0.24 dB under the slope variants. However, two matched CW-gain amplitudes move by more than 5% (up to 6.3%): one is at an uninformative state; the other is a comparison that holds with margin.
  - **Regime (V-C):** the archived E1dB-matched configuration checks put the signal above the LO (2.84 / 1.06 / 3.49 × LO), as do the archived dwell checks at other IFs, the weaker-probe step (1.50) and brief excursions of the shuffle family (0.50% of samples at +3 dB). The claims that use these states are marked WEAKEN, and flagged.

## Surviving Central Contribution

In a converged, thermally averaged transient simulation of a Rydberg superheterodyne receiver:
- **CW-derived surrogates fail.** In the 5 MHz reference configuration they mispredict communication-waveform degradation, with errors of both signs, over a 16-fold noise range.
- **A standard identified memory model recovers the amplitude-driven part** of what they miss. The generalized memory polynomial predicts QPSK and 16-QAM rates within 0.04 bit.
- **It cannot represent the phase-driven part.** A constant-amplitude, phase-modulated gain reduction (coherent gain 0.47 of CW) is:
  - present already below the CW output maximum;
  - specific to the reference configuration under fair matching;
  - located, across IF, at a fixed multiple of the LO Rabi frequency, without a quantitative collapse.
- The receiver's microsecond large-signal memory and the temporal-order effects complete the picture.

## Claims Removed

- **Title (T0),** "Configuration-Dependent Dynamic Compression…" (the only REMOVE).
- **From the abstract** (marked WEAKEN, so kept only in the body with their restriction):
  - the E1dB-matched statement "absent at 10 and 15 MHz and with a weaker probe" (A7; replaced by the fair-matched N1);
  - "the slow memory persists with the weaker probe" (A8; above the LO);
  - the declustering statement (A3; brief excursions above the LO);
  - "strongly configuration dependent" as a general statement (A9a).
- **From "what remains robust":** the dwell gain effect beyond the static prediction at 15 MHz (B5; above the LO).
- **From Section VII:** the cross-IF S5 selectivity comparison (B2; restated from the IF scan).

## Claims Strengthened

- **A1 (+ N4):** the surrogate failure holds at every tested noise level (0.25–4×; effective SNR 20–46 dB).
- **A6 (+ N5):** the phase effect is present below the CW output maximum (X = 0.33 at CW gain 0.8, slope +0.53), not only past it.
- **N1:** configuration dependence under fair matching (11 of 11 comparisons), replacing a comparison outside the intended regime.
- **B1:** the LO-field ordering at matched compression (11 of 12).
- **B7 / N2:** a standard GMP was tested. It captures amplitude-driven distortion but, structurally, not the phase effect.
- **C5:** the E1dB-rule sensitivity is small in the reference configuration (≤ 0.24 dB).

## Best Paper Story (A/B/C)

**B: engineering memory modeling.** It was selected by the prespecified rule:
- A requires fair matching SURVIVES and resonance CLEAR SCALING; neither holds.
- B's conditions all hold: fair matching PARTIALLY SURVIVES, memory model PARTIAL SUCCESS, resonance not CLEAR.

## Recommended Final Title

**"CW-Derived and Identified Models of Large-Signal Memory in a Simulated Rydberg Atomic Receiver"** (now in the manuscript).

Optional alternative from the review (R5.1, MINOR): "Testing CW-Derived and Identified Memory Models Against Transient Simulations of a Rydberg Atomic Receiver". The request's Story-B suggestion, "Modeling Configuration-Dependent Large-Signal Memory…", was not adopted. The only configuration dependence established inside the intended regime is that of the phase effect; the memory's configuration dependence rests on checks above the LO field.

## Final TQE Reviewer Attack (BLOCKING/IMPORTANT/MINOR)

Full text: `reviewer_attack.md`.

**BLOCKING: none.**

**IMPORTANT:**
- **R1.1:** at matched signal/LO, configuration and compression state are confounded.
- **R2.1:** the LO scale is supported in location only. Normalizations are not discriminated; grid and tolerance limit the test; one sequence per IF.
- **R2.2:** the weaker probe shows the LO scale is not sufficient.
- **R3.1:** the GMP is structurally unable to represent constant-envelope phase effects.
- **R3.2:** the GMP sits at the grid edge, with diminishing validation gains.
- **R4.1:** the tested SNR range (≥ 20 dB) is high.
- **R6.1:** the significance beyond "simple models fail" must be argued.

All seven were addressed by wording in the manuscript (commit c719743), with no new experiments.

**MINOR:**
- **R1.2:** exclusion of ratio 0.6 — addressed in S11.
- **R1.3:** CW-gain matching of LO fields changes signal/LO — addressed in S11.
- **R1.4:** reference level beyond the CW output maximum — already addressed.
- **R2.3:** X_coh ≈ 0 at the CW-gain minimum — addressed.
- **R3.3:** reason for the CE failure — addressed.
- **R4.2:** noise model — addressed in Limitations.
- **R5.1:** title wording — optional.

## Publication Decision

**TQE BORDERLINE SUBMIT.**
- By the prespecified rule: no stop rule fired; fair matching PARTIALLY SURVIVES; noise PARTIALLY ROBUST; story B; no BLOCKING comment.
- It is not STRONG, because fair matching did not fully survive: the operating-point ordering fails at one comparison, and the headline contrast rests on only one compression-matchable configuration. The numerical class is also CAVEAT.

## Recommended Next Step

1. **Author review** of the revised manuscript (`submission/tqe/main.pdf`, 21 pages; `supplement.pdf`, 13 pages). Check in particular the new Section VIII, the rebuilt Section VII and the title. Then decide whether to merge `tqe-viability-final` into `tqe-submission`; nothing has been merged or pushed.
2. **Author metadata, then publication.** Fill the remaining author-input placeholders (repository URL, DOI). Publish the reproducibility archive (built locally, 2,715 files, 403.5 MB; not published) only when the author decides.
3. **Do not add simulations before submission.** The campaign answered the four questions it was set; further scans would be the "rescue" the rules forbid.
   - The natural future work is experimental: test the constant-envelope phase effect and its IF dependence near f_IF ≈ 1.1 Ω_LO/4π on a real superheterodyne receiver.
   - The natural modeling follow-up is a phase-sensitive identified model (a separate study).

## Reproducibility

- **Run manifest.** Every new result has a manifest row (`00_run_manifest.csv`: configuration, seed, command, code version, runtime, source, output, SHA-256). Analyses record their code version in their verdict files.
- **Audits.** The draft check (all numbers in the ledger, all IDs resolved, all citations resolved), the figure audit (242 checks, ALL PASS) and the citation audit (49 keys, 0 BibTeX warnings) pass.
- **Package.** The reproducibility package was rebuilt with the campaign records, code, results and raw runs (2,715 files, about 404 MB), and its SHA256SUMS was regenerated. It is not published.
- **Clean-room reproduction.** A copy of the package was overlaid with the upstream simulator at 4095724 (from `git archive`), and `verify_reproduction.py` was run. As shipped: **PASS**, 97 of 97 regenerated files reproduced; the identified-memory-model fit was skipped, because its 128 revision P1 traces are not shipped, as before. With those traces linked in: **PASS**, 100 of 100, the GMP tables within 1.7e-8 relative, which is BLAS reduction order in an ill-conditioned fit (compared at 1e-6).
- **Process notes.**
  - The prior-pass runner could not be stopped early (the action was not permitted). Its P2 scan and P4 reruns therefore completed under that runner and were declared as reused inputs in the plan.
  - A local scratch path in `tvf_common.py`, used only to wait for that runner, was replaced by an optional environment variable before release.
