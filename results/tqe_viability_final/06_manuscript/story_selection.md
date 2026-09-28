# Part VI — story selection (rule: `00_plan.md`, section E)

Selection was made after Parts I–V were complete. No stop rule fired: S1 (fair matching), S2 (regime), S3 (memory model) and S4 (numerics) were all false.

## Classes

| Question | Class | Source |
|---|---|---|
| Fair matching | **PARTIALLY SURVIVES** | `01_matching/matching_verdict.json` |
| Resonance | **SUGGESTIVE ONLY** | `02_resonance/resonance_verdict.json` |
| Memory model | **PARTIAL SUCCESS** | `03_memory_model/gmp_verdict.json` |
| Noise | **PARTIALLY ROBUST** | `04_noise/noise_verdict.json` |
| Numerical sanity | **CAVEAT** | `05_numerical/numerical_verdict.json`; claim matrix |

The detailed results behind each class:

- **Fair matching.** The headline IF/probe dependence (HA1) holds at all 11 evaluable pairs under both conventions. The operating-point ordering (HA2) fails at one pair.
- **Resonance.** All three features align under LO normalization (r1, r2, d2, s3, d3), but the curves do not collapse (collapse ratios 0.77–1.12).
- **Memory model.** The GMP beats the surrogates' NMSE by ≥ 20% in 3 of 4 families. It reproduces no key contrast within tolerance.
- **Noise.** Surrogate failure holds at all five scales. The broad ordering changes at 4×.
- **Numerical sanity.**
  - V-A and the free diagnostics pass.
  - V-B shows a > 5% shift in two matched amplitudes.
  - V-C: retained (WEAKEN) claims use states with signal > LO.

## Conditions of each story

| Story | Condition | Holds? |
|---|---|---|
| A | Fair = SURVIVES | no (PARTIALLY SURVIVES) |
| A | Resonance = CLEAR SCALING | no (SUGGESTIVE ONLY) |
| B | not A | yes |
| B | Fair ∈ {SURVIVES, PARTIALLY SURVIVES} | yes |
| B | Memory ∈ {SUCCEEDS, PARTIAL SUCCESS} | yes |
| B | Resonance ≠ CLEAR SCALING | yes |
| C | Fair ∈ {PARTIALLY SURVIVES, DOES NOT SURVIVE} | yes |
| C | Resonance ≠ CLEAR SCALING | yes |
| C | KEEP claims are only convergence, temporal-order and model-transfer results | no. The KEEP list includes the fair-matched configuration dependence (N1), the GMP result (N2) and the noise-robust surrogate failure (N4). |

**Selected: Story B (engineering memory modeling).** All of B's conditions hold, and the rule evaluates A → B → C in order.

- **Main contribution under B:** CW-derived surrogates fail, and a standard identified memory model partially recovers the dynamics. The manuscript must state exactly which dynamics remain uncaptured (the PARTIAL SUCCESS instruction).
- **Headline suggested by the request:** "Modeling Configuration-Dependent Large-Signal Memory in a Simulated Rydberg Atomic Receiver". The final title is checked against Part IX, question 5 (the title must match the surviving contribution); see `07_final_review/`.
