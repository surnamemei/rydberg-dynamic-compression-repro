# Part IX — final hostile TQE reviewer attack

**Object of review.** The Part VIII manuscript: `submission/tqe/main.pdf` (21 pages) and `supplement.pdf` (13 pages), commit 24a8b84, branch `tqe-viability-final`.

**Ground rules.**
- The review asks only the six questions of the request and requests no additional experiments.
- Classes follow `00_plan.md`, section E:
  - **BLOCKING:** invalidates a retained claim, or makes the paper unfit for the venue as written, and cannot be fixed by wording.
  - **IMPORTANT:** must be addressed in the text before submission; the retained claims stand.
  - **MINOR:** presentation.

## 1. Are the compared operating states genuinely fair?

**R1.1 (IMPORTANT). Configuration and compression state are still confounded at matched signal/LO.**
- At signal/LO 0.2–0.4 the reference configuration is compressed (CW gain 0.34–0.81 of small-signal). At the same ratios, 10 MHz and the weaker probe expand (1.07–1.64) and 15 MHz barely compresses (0.90–0.97).
- The one compression-matched comparison is 5 vs 15 MHz at mild compression (CW gain 0.9 and 0.85). It pairs very different field ratios: 0.14–0.17 × LO at 5 MHz against 0.40–0.76 × LO at 15 MHz.
- "Specific to the reference configuration among those tested" is therefore, in part, "specific to the only configuration that compresses inside the superheterodyne regime".
- The paper states the facts (Section VII-D) but does not draw this consequence where it states the claim (abstract, C4, Conclusion).
- *Fix by wording:* say explicitly that at 10 MHz and with the weaker probe the absence of the effect coincides with the absence of CW compression below the LO field, and that only 15 MHz provides a compression-matched contrast.

**R1.2 (MINOR). The ratio 0.6 was excluded by a pathology rule defined after the CW tables were known.**
- The rule was fixed before any matched-state X value was seen, and the prior-pass rule that includes 0.6 gives the same headline outcome.
- A sentence saying so (in S11) would pre-empt a forking-paths objection.

**R1.3 (MINOR). The matched-CW-gain comparison of LO fields also changes signal/LO** (e.g., OP2 at CW gain 0.4 needs 0.64 × LO). The paper reports the amplitudes (Table S11) but not this confound.

**R1.4 (MINOR). The reference level lies beyond the CW output maximum.** This is addressed: the effect is already 0.33 [0.26, 0.41] at CW gain 0.8, where the local slope is +0.53. No action needed beyond keeping this paragraph.

## 2. Is the strongest claimed organizing variable supported?

**R2.1 (IMPORTANT). The LO-scale organization is supported in location only, and only weakly discriminated.**
- The feature locations align at r₂ = 1.02–1.03 and 1.11 across three LO fields spanning a 1.43-fold range. Over so narrow a range, ratio and difference normalizations (r₂, d₂, d₃) all "align" (Table S12), so the scan cannot say which scale organizes the effect.
- The features were located on grids with 0.2–0.3 MHz spacing near the maxima, comparable to the 8% alignment tolerance.
- One sequence per IF is used, with a sequence SD of X_coh up to 0.08 against feature contrasts of 0.3–0.5.
- The paper classifies this as "suggestive only" and claims no mechanism, so the claim as written is supported. It should, however, add the non-discrimination between normalizations and the grid/tolerance caveat where it reports the alignment (Section VII-E).

**R2.2 (IMPORTANT). The weaker probe shares r₂ = 1.08 with the reference configuration but shows no effect, so the LO scale is at most necessary, not sufficient.**
- The paper says so (Section VII-E). The Discussion and Conclusion should not speak of the maximum "following the LO Rabi frequency" without this qualifier, and should name the probe Rabi frequency as an uncontrolled second variable.

**R2.3 (MINOR). The co-location of the CW-gain dip (r₂ ≈ 1.03) with a sharp drop of X_coh to about zero at the same IF** (visible in Fig. 9(e)) argues against the simple reading "the phase effect follows the local compression depth". The paper could state this observation (post hoc) in one sentence.

## 3. Is the engineering baseline strong enough?

**R3.1 (IMPORTANT). For constant-envelope phase modulation the GMP class is structurally unable to represent an envelope-magnitude reduction.**
- With |x| constant, every aligned and lagging term reduces exactly to a linear filter of x. Only the second-zone terms, through x², carry phase nonlinearity (checked algebraically for the implemented basis).
- The GMP therefore cannot express the X_mag part of the phase effect except through its linear memory, and it captured X_mag 0.02 of 0.45.
- Its "miss" of the phase effect thus demonstrates a limitation of envelope-magnitude-driven behavioural models, the standard pre-distortion class, rather than a failure of identification.
- This is the more useful engineering statement, and the paper should make it (Section VIII, Discussion). The reader should not infer that a larger GMP would help: validation gains were already diminishing (0.1565 → 0.1531 from 2 to 4 µs; 0.1535 → 0.1531 from order 5 to 7).

**R3.2 (IMPORTANT). The selected model lies at the edge of the prespecified grid** (largest order, deepest memory). The paper should report this, with the diminishing validation gains above, so that the partial-success class is not read as a capacity artifact.

**R3.3 (MINOR). The GMP's CE-link failure (0.80 bit) goes without a stated reason.** The discriminator receiver is sensitive to phase errors, which the GMP models only linearly (R3.1); one clause would connect the two results.

## 4. Does the communication conclusion survive noise assumptions?

**R4.1 (IMPORTANT). The tested SNR range is high.**
- The nominal noise corresponds to an effective SNR of 32–34 dB at +3 dB. The tested range, 0.25–4×, spans 44–46 down to 20–21 dB.
- Many wireless links operate at 10–20 dB. The surrogate failure is shown to be noise-insensitive only down to about 20 dB, and the format ordering already changes at 4×.
- The paper gives the range (Section III-C). The abstract's "persist when the receiver noise is scaled from 0.25 to 4 times its nominal level" is correct but should carry the SNR range, or the SNR range should appear in C1 and the Conclusion.

**R4.2 (MINOR). The noise model is additive white Gaussian noise after the probe detector.** Atomic projection noise and photon shot noise enter only through this aggregate. A reader would benefit from one sentence in the Limitations.

## 5. Does the title exactly match the surviving contribution?

**R5.1 (MINOR).** "CW-Derived and Identified Models of Large-Signal Memory in a Simulated Rydberg Atomic Receiver" matches the modeling contribution (CW-derived surrogates fail; one identified model partially recovers the behaviour).
- It does not signal the most distinctive result, a phase-modulated, configuration-specific gain reduction that neither model class captures.
- It also reads as if models were proposed rather than tested.
- A slightly more exact alternative: "Testing CW-Derived and Identified Memory Models Against Transient Simulations of a Rydberg Atomic Receiver". This is optional.

## 6. Is the significance more than "simple models fail"?

**R6.1 (IMPORTANT).** The paper offers four things beyond "simple models fail", and it should say so in its contribution and conclusion text:
- **A reproducible benchmark:** rate- and bandwidth-matched formats, amplitude-preserving temporal-order signals and prespecified tests.
- **A quantified division of labor between model classes.** CW-derived models fail across a 16-fold noise range. A standard envelope-driven memory model recovers amplitude-driven distortion (QPSK/16-QAM rates within 0.04 bit) but, by construction, not the phase-driven effect (R3.1).
- **A fair-matched configuration dependence** of the phase effect.
- **A suggestive LO-scale organization,** which tells experiments where to look (f_IF ≈ 1.1 Ω_LO/4π).

Two things limit the significance and must stay in view:
- The work is simulation-only, with one inherited model.
- The central phase effect has no mechanism.

This is a judgment about venue fit rather than a defect that invalidates a claim. The retained claims stand, and the significance argument can be made in the text.

## Summary

| Class | Comments |
|---|---|
| **BLOCKING** | none |
| **IMPORTANT** | R1.1, R2.1, R2.2, R3.1, R3.2, R4.1, R6.1 |
| **MINOR** | R1.2, R1.3, R1.4, R2.3, R3.3, R4.2, R5.1 |

**No BLOCKING comment.** Every IMPORTANT comment is fixable by wording; none requests an experiment or invalidates a KEEP claim. The IMPORTANT comments are addressed by wording in the manuscript after this review (see `final_report.md`); the classification above refers to the Part VIII text.
