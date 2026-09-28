# Second-Operating-Point Replication: Result

**Classification (preregistered rule): B. The effect weakens, but the qualitative effect is preserved.**

Preregistration: `00_preregistration.json`, committed in `8d555b6` before any run. All 49 planned runs were executed and nothing was added.

## Operating points

| | Original | Second (LO′) |
|---|---|---|
| LO | 0.5 V/m (Rabi 2π × 9.2 MHz) | 0.35 V/m (Rabi 2π × 6.46 MHz) |
| E1dB (Nd = 4001, same CW protocol) | 0.0732 V/m | 0.0403 V/m |
| Operating amplitude (2.4468 × E1dB) | a_H = 0.179 V/m | a′ = 0.0987 V/m |
| CW gain at that amplitude (full model; CW-static prediction) | 0.43 (0.43) | 0.633 (0.635) |

The same relative level (2.45 × E1dB) is **less compressed** at LO′: CW gain 0.63 vs 0.43. The E1dB-based level matching does not equalize the degree of nonlinearity, and this cannot be separated from the operating-point dependence without further runs, which the preregistration excludes.

## Preregistered tests (`04_replication_verdicts.csv`)

| Test | Result | Outcome |
|---|---|---|
| H_gen: X_REF′ ≥ 0.20 | X_REF′ = 0.159 (seeds 0.188 / 0.129); original point 0.520 | **FAIL** |
| H_shape: \|X_JUMP′ − X_LONGRAMP′\| ≤ 0.10 | JUMP′ 0.194, REF′ (300 ns) 0.188, LONGRAMP′ (900 ns) 0.038 | **FAIL**: 900 ns ramps largely suppress the effect at LO′ |
| H_detune: all offsets X′ < 0.5 X_REF′ | +0.125 MHz: +0.007; +1.31: −0.562; −1.31: −0.178 | **PASS** |
| H_drift (original point): zero-drift X_FAST ≥ X_REF − 0.03 | zero-drift SLOW / REF / FAST = 0.178 / 0.520 / 0.406 (+π convention: 0.159 / 0.520 / 0.438) | **FAIL**: the turnover is real, not a drift artifact |
| C9a: \|power-step period − phase-step period\| ≤ 0.15 µs at LO′ | 4.00 vs 3.68 µs, both at the edge of the estimator band (0.8–4 µs) | **FAIL** as evaluated (see post-hoc note) |
| Numerics (Nd = 8001 closure) | OFF_+0.125: X −0.1730 → −0.1746; ξ = 32 long-dwell D_ref 0.03266 → 0.03268 (0.04%) | **PASS** |
| Classification | X_REF′ = 0.159 ≥ 0.05 but not all of H_gen / H_shape / H_detune | **B** |

## Descriptive results at LO′ (not separately preregistered)

- **Rate triplet** (zero-drift SLOW / REF / FAST): X′ = 0.041 / 0.188 / 0.197. The effect rises and saturates; there is no turnover at this point.
- **Frequency selectivity:** the sharp selectivity near the IF seen at the original point is **absent** at LO′ (a +125 kHz offset gives X′ = +0.007, against −0.173 at the original point). The phase-driven excess persists where the special feature is absent, so it is not a byproduct of that feature.
- **Atomic-state proxies:**
  - probe-baseline shifts are −0.02 to −0.15 dB for phase cases (≤ 0.03 dB at the original point) and up to −0.44 dB for offsets;
  - thermally averaged Rydberg population changes by ≤ +1.2%.
- **Build-up time:** t63 = 1.8–6.0 µs.

## Post-hoc note on C9a (`05_posthoc_C9a_wideband.json`; labeled post hoc)

The transients at LO′ oscillate at ≈ 0.56 µs, outside the 0.8–4 µs band of the estimator carried over from the original point. Re-estimated with a 0.3–4 µs band:
- power step: **0.560 µs**;
- nine isolated phase steps: **0.556–0.564 µs** (median 0.558);
- difference: 0.002 µs.

The same estimator gives 1.62 µs at the original point. Power-step and phase-step transients therefore share one oscillation frequency at *both* operating points, and that frequency changes with the LO (≈ 0.62 → 1.79 MHz). This strengthens the "common slow mode" observation (C9a) as a post-hoc result. The preregistered verdict remains FAIL as evaluated.

## Consequences for the paper claims

- **C5** (phase modulation adds compression beyond CW-derived models): **replicated qualitatively** at a second LO.
  - The extra compression is smaller: X = 0.13–0.19 vs 0.52, at a less compressed operating point.
  - Static, Hammerstein and slow-state models still predict X = 0.
- **C7:**
  - "Not explained by detuning": **replicated**.
  - "Independent of transition shape": **not replicated**. At LO′ abrupt and 300 ns transitions are equivalent, but 900 ns ramps suppress the effect. Shape independence becomes an original-point result, and the ramp duration that matters is operating-point dependent.
- **C8:** the rate dependence is operating-point dependent (turnover at the original point, now confirmed not to be drift; saturation at LO′). "No single scalar explains all cases" is strengthened.
- **C9a:** strengthened post hoc (the mode match holds at both points, at LO-dependent frequencies).
- **C12:** the sharp frequency selectivity is **specific to the original operating point**. Rescope it as an observation about that point.
- **Numerics:** the two open Nd = 8001 items are closed.

Per the instruction, the sweep is not broadened.
