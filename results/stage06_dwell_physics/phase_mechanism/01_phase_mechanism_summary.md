# Phase-Mechanism Decision

**CONDITIONAL GO**, reached by the pre-registered decision logic (`00_preregistration.json`, `04_phase_mechanism_verdicts.csv`).

The phase-modulation excess compression is **reproducible and numerically exact** at constant amplitude. Three of the four candidate explanations are **rejected**:
- instantaneous-frequency swing (B);
- detuning (C);
- discrete jumps (D).

The surviving explanation is **transition-driven** (A): each change of the signal phase knocks the atomic response out of its beat-locked steady state, and it takes ~τ_atom to re-lock. However, no pre-registered single variable organizes all cases cleanly:
- the compression saturates and turns over at the highest transition rate;
- the step pattern shifts it by ±25%.

So the mechanism is identified qualitatively, but not as a one-variable law.

# Experiment Design

- **Operating point:** constant |E| = a_H = 0.17915 V/m (the +3 dB dwell-family high level, 2.45 × converged E1dB), LO 0.5 V/m, IF 5 MHz, full thermal transient model.
- **Timeline:** 0–50 µs is an unmodulated CW tone at a_H, giving every case the same settled starting state. At 50 µs the case-specific phase or frequency modulation switches on **at the same amplitude**. Runs end at 110 µs, so the returned density matrices describe the end-of-plateau atomic state.
- **Numerics:** Nd = 4001, dt = 1 ns. Checks at Nd = 8001 and at dt = 0.5 ns used waveforms regenerated analytically at 2 GHz.
- **Metric:** the gain g(t) = |z_full| / |z_linear|. Both are the IF response demodulated with the known input phase; z_linear is the small-signal model driven by the identical input, which removes frequency-dependent small-signal and demodulator effects. The primary metric is X = 1 − g_final / g_CW, with g_final the median over the last 20 µs.
- **Cases (13):**
  - CW;
  - constant offsets of +0.125, ±1.31 and ±2.62 MHz (the QPSK mean drift, and the peak excursions of π/2 and π transitions with 300 ns ramps);
  - random QPSK: REF (1 µs symbols, 300 ns raised-cosine ramps), a second seed of REF, SLOW (4 µs), FAST (0.5 µs), JUMP (instantaneous steps), LONGRAMP (900 ns ramps);
  - TOGGLE, a deterministic period-4 cycle with the same transition statistics.

  JUMP and LONGRAMP use REF's exact symbol sequence.
- **Scale:** 21 runs in total (13 cases plus 8 numerical checks), about 2.5 GPU-minutes.

# Matching Quality

See `02_phase_mechanism_configs.csv` and fig_pm01/02.
- **Amplitude:** exactly constant in every case.
- **REF / JUMP / LONGRAMP** are identical in the number of nonzero transitions (44), total phase excursion (94.2 rad), mean |Δf| (250 kHz) and mean signed drift (142 kHz). They differ only in the peak frequency swing (2.62 MHz / 500 MHz, i.e. a step / 0.87 MHz) and the RMS swing.
- **SLOW / REF / FAST** share the ramp shape and peak swing (2.62 MHz), with 0.15 / 0.73 / 1.45 transitions per µs.
- **TOGGLE vs REF:** 44 vs 44 transitions, 91.1 vs 94.2 rad, mean |Δf| 242 vs 250 kHz.
- **REF seed 2:** 42 transitions, 86.4 rad.
- **Offsets** have zero transitions and constant |Δf| = |δ|.

# CW Baseline

g_CW = 0.4325. The CW-static (all-zone table) prediction is 0.4308, and the pre-modulation level in every case is 0.4329. The same number holds at Nd = 8001 (0.4317) and at dt = 0.5 ns (0.4324).

The static, Hammerstein (LTI+static) and single-slow-state models all predict g ≈ 0.43 for **every** case, because |E| is identical and constant. Through the same pipeline they give 0.434–0.436 everywhere.

# Transition-Rate Test

See fig_pm04.

| Transitions per µs (per τ) | Case | X |
|---|---|---:|
| 0 | CW | 0 |
| 0.15 (0.38) | QPSK_SLOW | 0.159 |
| 0.73 (1.86) | QPSK_REF | **0.520** (second seed 0.523) |
| 1.45 (3.7) | QPSK_FAST | 0.438 |

- X rises steeply up to ~2 transitions per τ_atom and then **turns over**. The pre-registered H1 required a monotone rise through FAST, so it is **not supported** as stated.
- The rate effect itself is large and organized on the τ_atom scale.
- TOGGLE (same rate and statistics, deterministic order) gives X = 0.651. That is +0.131 against a pre-registered limit of 0.130, so the step pattern matters at about the 25% level.

# Phase-Ramp Test

See fig_pm05.

| Ramp | X |
|---|---:|
| instantaneous steps (JUMP) | 0.523 |
| 300 ns (REF) | 0.520 |
| 900 ns (LONGRAMP) | 0.503 |

The peak instantaneous-frequency swing varies by a factor of about 575 (0.87 MHz to a step) while X changes by ≤ 0.02. The time to reach 63% of the final compression is 2.8 / 2.9 / 3.5 µs.
- **H2 (frequency swing): rejected.**
- **H4 (discrete jumps): rejected.**

How the phase gets from one value to the next does not matter, as long as it arrives within ≲1 µs (< τ_atom).

# Frequency-Offset Test

See fig_pm06.

| Offset | +0.125 MHz | +1.31 | −1.31 | +2.62 | −2.62 |
|---|---:|---:|---:|---:|---:|
| X | −0.17 | −0.83 | −0.50 | −2.02 | −0.29 |
| large-signal output vs CW | 1.16× | 1.52× | 1.74× | 2.34× | 1.21× |
| small-signal output vs CW | 0.99× | 0.83× | 1.16× | 0.78× | 0.91× |

- Every constant offset **reduces** compression. The large-signal response is sharply frequency-selective near exactly 5 MHz: a 125 kHz shift already lowers compression by 17%, while the small-signal response moves by only 1%.
- **H3: rejected**, in both the strong form (mean drift) and the partial form (instantaneous detuning). Offsets have the wrong sign.
- If the atom tracked the instantaneous frequency quasi-statically, phase ramps would *reduce* compression. They increase it.

# Atomic-State Evidence

See fig_pm07 and fig_pm08.

**Isolated transitions** (the SLOW case, 4 µs apart):
- **π steps:** a brief spike *up* (gain 0.75–0.9 for ~0.2 µs), then a shallow dip.
- **±π/2 steps:** a small bump, then a **delayed, deep dip** (gain ≈ 0.03–0.1 about 1 µs after the step) that recovers over ≥ 2–3 µs and is not complete before the next transition.

**Onset:** under repeated transitions the compression develops with t63 ≈ 2.5–3.5 µs, i.e. on the τ_atom scale.

**Populations and probe baseline at the plateau end:**
- The phase-modulated cases change the thermally averaged Rydberg population by only +0.2% to +1.3% and the probe-transmission baseline by ≤ 0.03 dB, while the coherent IF response falls by ~52%.
- The constant offsets change the Rydberg population by −1% to −7% and the baseline by up to −0.32 dB, yet they compress *less*.

The extra compression therefore sits in the **phase-coherent (beat-locked) part of the atomic response**, not in population shelving.

**Numerics:** max |ΔX| = 0.0028 over the 8 checks (Nd 4001 → 8001 and dt 1 → 0.5 ns for CW, REF, JUMP and +1.31 MHz).

# What Is Supported

- At fixed constant amplitude, QPSK-like phase transitions deepen compression by X ≈ 0.52 (the full-model gain falls from 0.43 to 0.21 of small-signal). This reproduces across seeds (0.520 / 0.523) and is exact under the Nd and dt checks.
- The excess is driven by **the occurrence of phase transitions**, and its size depends on how often they occur relative to τ_atom: 0.16 at 0.38 transitions per τ, 0.52 at 1.9 per τ.
- It is **independent of transition shape and speed** (a step, 300 ns and 900 ns give the same result within 0.02).
- It is **not a detuning effect**: constant offsets of 0.125–2.62 MHz reduce compression instead.
- It develops and relaxes on the atomic slow timescale (t63 ≈ 2.5–3.5 µs; single-transition dips last ~µs).
- It is not reproduced by static, Hammerstein, or |E|²-driven single-slow-state models, which predict X = 0 for every case.

# What Is NOT Supported

- A single-variable law. Transition rate turns over at 3.7 transitions per τ, and the step pattern (TOGGLE) and step size (π vs π/2 responses) modify X by about 25%.
- Any claim that discrete jumps or large instantaneous-frequency excursions are special.
- Any claim that the effect is ordinary detuning, or that it can be modeled as frequency-dependent static compression.
- Universality beyond this model, configuration (Transit, LO 0.5 V/m, IF 5 MHz) and amplitude (+3 dB level). Other amplitudes, IFs and LO strengths were not tested.
- Hardware relevance. This is a simulation study.

# Physical Interpretation

Under a steady LO + signal beat, the thermal ensemble settles into a periodic steady state that is **phase-locked to the beat**. That locked state has a specific, strongly frequency-selective compression. A constant frequency offset does not break the lock: after ~τ the ensemble re-locks to the new beat, and at 5 MHz ± δ that locked state happens to be *less* compressed.

A change of the signal phase relative to the LO does break the lock. Re-locking requires the slow (µs) atomic dynamics, and while it is out of lock the coherent IF response is strongly reduced: the delayed dips after each transition. Because re-locking (τ ≈ 2.5–3 µs) is much slower than any tested ramp (≤ 900 ns), the ramp shape is irrelevant. What matters is how often the lock is broken relative to τ, and the size and order of the phase steps. Beyond ~2 transitions per τ the ensemble never re-locks, and additional transitions no longer deepen the compression.

This is the same slow atomic memory found in Stage-0.6, **expressed through the phase of the drive rather than its power**.

# Implication for Stage-0.6

- Stage-0.6 "mechanism 2" is now identified as **transition-driven loss of beat-lock with τ_atom-scale re-locking**. It is not a detuning or instantaneous-frequency effect, and not specific to discrete jumps.
- Both Stage-0.6 mechanisms share the same slow atomic timescale:
  - the power-driven slow state (averaging, post-excursion suppression, clustering, quasi-static convergence for an unmodulated carrier);
  - the phase-driven de-locking (extra compression under symbol transitions).
- This explains why, with the QPSK carrier, long high-field plateaus remained more compressed than CW curves even at 32·τ. The symbol transitions keep breaking the lock during the plateau.
- A model that tracks only |E|², whether static, Hammerstein, or with a slow power state, cannot capture it. A candidate extension, not built or fitted here, is a slow state driven by changes of the complex envelope's phase relative to the LO.
- Overall Stage-0.6 status: **CONDITIONAL GO**, with both mechanisms now characterized and tied to τ_atom. A clean quantitative law for the phase-driven part is still missing.

# Recommended Paper Claim

> In a full thermal transient model of a Rydberg-atom superheterodyne receiver, nonlinear compression depends on the temporal organization of the signal relative to a microsecond-scale atomic relaxation time (τ_atom ≈ 2.5–3 µs), and is therefore not characterized by CW compression metrics such as P1dB. Two manifestations share this timescale:
> 1. Amplitude excursions shorter than τ_atom are averaged, while those much longer than τ_atom follow the CW compression curve.
> 2. At constant amplitude, phase transitions occurring more often than about once per τ_atom deepen compression by about 50% beyond the CW prediction. This excess is independent of the transition shape and is not reproduced by frequency detuning, consistent with repeated loss and slow recovery of the atomic response's phase lock to the LO–signal beat.
>
> These results are specific to the simulated configuration and are not claimed to be universal.
