# Part VII — final claim / evidence matrix

**Rules.** `00_plan.md`, section E.
- **KEEP:** an explicit comparison with a 95% CI (or a prespecified classification) that survives every campaign test relevant to it, inside the intended regime.
- **WEAKEN:** supported only in a named, restricted domain.
- **REMOVE:** contradicted, or resting only on comparisons the campaign invalidated.

Only KEEP claims may appear in the abstract.

**Scope of the matrix.** The candidate claims of the current manuscript were fixed in `claims_inventory.md` before any result was analysed. Every one of them is assessed here; none is dropped.
- The campaign's own pre-registered tests produced new results. These enter as claims **N1–N6**, each tied to its pre-registered test and class. This goes beyond the inventory's sentence that claims would be "neither added nor dropped", and is stated here explicitly.
- Inventory row A9 joins two statements with a semicolon. They are assessed separately as A9a and A9b.

**Abbreviations.**
- *Evidence*: DIR = direct measurement or comparison; INF = inferred from several results; INT = interpretive.
- *Origin*: pre = prespecified (decision rule fixed before the run); post = post hoc.
- *Fair*: the Part I result bearing on the claim.
- *Noise*: the Part IV result.
- *Baseline*: the Part III GMP result.
- *Regime*: the V-C flag.

## Decisions

| ID | Claim (current wording, condensed) | Supporting result | Fair | Noise | Baseline | Generalization / regime | Origin | Evid. | Decision |
|---|---|---|---|---|---|---|---|---|---|
| T0 | Title: "Configuration-Dependent Dynamic Compression in a Simulated Rydberg Atomic Receiver" | — | HA1 holds; HA2 partial | — | GMP partial | — | — | — | **REMOVE** (replace: the Story B title must match the surviving contribution) |
| A1 | CW-matched static and LTI+static surrogates mispredict format-dependent AIR loss (5 MHz reference); errors of either sign | Table 2 (archived); Part IV (a) and (b) hold at c = 0.25, 0.5, 1, 2, 4; e.g. static, 16-QAM at 0 dB: +1.51 / +1.49 / +1.26 bit at c = 0.25 / 1 / 4; OFDM at +3 dB: −0.88 at c ≤ 2 and −0.85 at 4× | n/a (reference configuration) | robust over 0.25–4× (errors nearly noise-invariant) | the GMP predicts QPSK/16-QAM AIR within 0.02–0.04 bit (the surrogates 0.43–0.77); the claim concerns CW-derived surrogates | fair waveforms: peak ≤ 0.92 × LO (in regime) | pre | DIR | **KEEP** (add the noise range; add the GMP comparison) |
| A1b | Losses depend strongly on format at equal power (+3 dB: 16-QAM 0.83, QPSK 0.39, CE 0.38, OFDM 0.12 bit); P5 spans five dB | Table 1; Part IV. The broad ordering (16-QAM most, OFDM least) holds at 0.25–2×; at 4×, CE's loss (0.07 [−0.18, 0.32]) falls below OFDM's (0.08 [−0.01, 0.17]); neither is resolved. The P5 span is 5.2 / 5.0 / 3.9 dB at c = 0.25 / 1 / 4 | n/a | partially robust | — | in regime | pre | DIR | **WEAKEN**: state for the archived noise level and 0.25–2× of it; the P5 values are noise dependent (16-QAM −2.47 dB at 1×, −0.60 dB at 4×) |
| A2 | Longer dwells lower the fitted gain from 0.470 to 0.360 of small-signal (paired −0.110 [−0.118, −0.102]); the surrogates do not reproduce it | Sec. 4.2 (archived). Test set r4–r7: full −0.106; all-zone static −0.007; LTI+static −0.008; single-slow-state −0.165; GMP −0.139 | n/a | n/a (gain metric) | the GMP reproduces the sign, not the size (+30%) | dwell family peak ≤ 0.51 × LO (in regime) | post (decomposition) / pre (5 MHz gain test) | DIR | **KEEP** (add: a GMP reproduces its sign but overestimates it by 30%) |
| A3 | Declustering changes the fixed-reference residual distortion (D_ref −0.0057 [−0.0098, −0.0016] at +3 dB); the surrogates do not follow | Sec. 4.3; test set: full −0.0026; GMP −0.0057 (same sign, 2.2×); all-zone surrogates of opposite sign | n/a | n/a | the GMP reproduces the sign only | shuffle family peaks reach 2.25 × LO at +3 dB; 0.50% (+3 dB) and 1.45% (+6 dB) of samples exceed the LO | pre | DIR | **WEAKEN**: state the regime caveat (brief excursions above the LO); not in the abstract |
| A4 | A prespecified distortion-based robustness test of the dwell effect failed | Sec. 4.2 record | — | — | — | — | pre | DIR | **KEEP** (record of a failed test) |
| A5 | Slow memory of 2.5–3.1 µs (slow pole at 0.8–1.2 E1dB) | Supplementary Table S1; steps at 0.059–0.088 V/m (≤ 0.18 × LO). The GMP's validation error falls with memory depth up to the grid maximum of 4 µs (family-mean NMSE 0.179 at 0.25 µs → 0.153 at 4 µs) | n/a | n/a | consistent (a memory of microseconds is selected) | in regime | pre | DIR | **KEEP** |
| A6 | At constant amplitude, phase modulation reduces the coherent large-signal gain to 0.47 of CW (X_coh = 0.525), mainly via envelope magnitude (84%) | Sec. 7.1–7.2. Part I at the same state, six sequences: X_coh 0.526 [0.493, 0.558]. Test seeds: GMP X_coh 0.13 of 0.54, X_mag 0.02 of 0.45; surrogates ≤ 0.04 | reference state; see N1, N5 | n/a | uncaptured by the GMP (about one quarter of X; magnitude part not captured) | in regime (0.358 × LO) | pre | DIR | **KEEP** (add: the GMP captures only about one quarter) |
| A7 | At levels matched to each configuration's E1dB, the phase effect is absent at 10 and 15 MHz and with the weaker probe | revision2 archived: X = −0.096 / −0.018 / −0.106 at signal/LO 2.84 / 1.06 / 3.49 | superseded by N1 (fair matching) | — | — | **outside the intended regime** (signal > LO in all three) | pre | DIR | **WEAKEN**: keep only as the archived E1dB-matched record, with the regime flag; the claim itself is carried by N1 |
| A8 | The slow memory persists with the weaker probe (4.44 µs vs 2.70 µs) | revision2 P3 step at the weak-probe E1dB (0.749 V/m = 1.50 × LO) | not re-tested | — | — | **outside the intended regime** | pre | DIR | **WEAKEN**: "in a step test at the weak-probe E1dB, which lies above the LO"; not in the abstract |
| A9a | "The dynamic large-signal behavior of this receiver is strongly configuration dependent" | N1 (phase effect); A7, A8, B5 are regime-flagged | HA1 holds; HA2 partial | — | — | only the phase effect is fair-matched | pre | INF | **WEAKEN**: restrict to the phase-modulated gain reduction (N1) |
| A9b | No microscopic mechanism is identified | Sec. 7.6; Part II SUGGESTIVE ONLY; II-C finds no consistent scale | — | — | — | — | — | DIR | **KEEP** (scope statement) |
| B1 | X is smaller at lower LO fields (OP2 / OP3 / OP1 = 0.16 / 0.29 / 0.51 at E1dB-matched levels) | Part I, matched CW gain 0.8 / 0.6 / 0.4. X = 0.334 / 0.534 / 0.482 (OP1), 0.128 / 0.273 / 0.462 (OP3), 0.172 / 0.153 / 0.238 (OP2). OP1 − OP2 = +0.16 [+0.11, +0.22], +0.38 [+0.32, +0.44], +0.24 [+0.17, +0.31]. OP1 − OP3 = +0.21, +0.26, **+0.02 [−0.05, +0.10]** | HA2 holds at 11 of 12 pairs | n/a | — | in regime | pre | DIR | **KEEP** per the pre-registered reading (OP1 > OP2 at 3 of 3 targets); state the OP3 exception at CW gain 0.4 |
| B2 | Sharp frequency selectivity near the IF at 5 MHz (S5 = 2.02 vs 0.23 / 0.20 at 10 / 15 MHz); the phase effect accompanies it | The 10 / 15 MHz S5 values come from flagged states. Part II at signal/LO 0.358: the maximum of S and the minimum of g_CW lie at r2 = 1.02–1.03, the X_coh maximum at r2 = 1.11, at all three LO fields | — | — | — | archived 10 / 15 MHz values outside the regime | pre (Part II) | DIR (co-location) / INT (organization) | **WEAKEN**: restate with the Part II scan; drop the flagged 10 / 15 MHz S5 comparison |
| B3 | Constant detuning and a frequency-dependent static model do not reproduce the phase effect | Table 3; Table 4 | — | — | — | in regime | pre | DIR | **KEEP** |
| B4 | Power and phase steps oscillate with a common, LO-dependent period (post hoc; a shared time scale only) | Sec. 7.3. II-C: no candidate scale matches at both LO fields; \|f_IF − Ω_LO/4π\| matches OP2 only (1.77 vs 1.79 MHz) | — | — | — | in regime | post | DIR | **KEEP** (body only; add the II-C negative) |
| B5 | The dwell gain change exceeds the static prediction at 5 and 15 MHz but not at 10 MHz | revision2 dwell checks. The dwell high level lies at 1.07 × LO (15 MHz), 2.85 × LO (10 MHz) and 3.50 × LO (weak probe) | not re-tested | — | — | **outside the intended regime** (10 MHz, weak probe; 15 MHz marginal) | pre | DIR | **WEAKEN**: record only, with the regime flag; remove from "what remains robust" |
| B6 | Internal-coherence observation (ρ₄₃ collapse; association only) | Sec. 7.6 | — | — | — | in regime | pre | DIR | **KEEP** (body) |
| B7 | Modeling implication: memoryless, Hammerstein, single-slow-state and frequency-dependent static structures fail; Wiener / Volterra-type structures untested | Sec. 8.6; Part III (a standard GMP, a Volterra subset, now tested: PARTIAL SUCCESS) | — | — | new evidence | — | pre (Part III) | DIR | **KEEP**, revised: a standard GMP reduces waveform error 3–5× but misses the phase effect and the sizes of the temporal-order contrasts |
| C5 | Numerical lesson: Nd = 1501 biased the CW compression point by 1.57 dB | Sec. 2.6. V-B: the slope rule changes C1's E1dB by −0.05 to +0.24 dB (C5, C6: up to 0.41 and 0.84 dB) | — | — | — | — | pre | DIR | **KEEP** (add the V-B sensitivity) |
| N1 | Fair-matched configuration dependence: at the same signal/LO ratio (0.2, 0.358, 0.4), the phase effect is present at 5 MHz and small or absent at 10 MHz, 15 MHz and with the weaker probe; at matched CW gain (feasible only at 15 MHz), present at 5 MHz and absent at 15 MHz | Part I, six sequences. X(5 MHz) = 0.32 / 0.49 / 0.46; X(10 MHz) = 0.06 / 0.08 / 0.09; X(15 MHz) = −0.01 / −0.03 / −0.04; X(weak probe) = 0.01 / 0.04 / 0.05. Paired differences +0.27 to +0.52 (all CIs exclude 0). CW gain 0.9 / 0.85: X(5 MHz) 0.12 / 0.23 vs X(15 MHz) −0.04 / −0.01; differences +0.16 [+0.11, +0.21], +0.24 [+0.16, +0.33] | HA1: 11 of 11 pairs | n/a | — | in regime by construction; 10 MHz and the weak probe cannot be CW-gain matched with signal ≤ LO (they expand: CW gain 0.97–1.26 and 1.00–1.87) | pre | DIR | **KEEP** |
| N2 | A standard identified memory model (GMP, 1282 complex coefficients, 4 µs memory, order 7; split fixed before fitting) reduces the waveform error 3–5× relative to the CW-derived surrogates and predicts QPSK / 16-QAM AIR within 0.02–0.04 bit, but it misses the phase-modulated envelope compression (captures about one quarter), CE AIR (0.80 bit error), the sizes of the dwell and declustering contrasts, and the level dependence of the phase effect | Part III test NMSE: dwell 0.16 vs 0.50; shuffle 0.22 vs 0.66; phase 0.20 vs 1.11; fair 0.049 vs 0.055 (best surrogate) | — | archived noise level | this is the baseline test | in regime (fair, dwell, phase); shuffle as A3 | pre | DIR | **KEEP** |
| N3 | The IF dependence is organized by the LO Rabi frequency: at all three LO fields, the CW-gain dip and the selectivity maximum lie at f ≈ 1.03 Ω_LO/4π and the phase-effect maximum at f ≈ 1.11 Ω_LO/4π | Part II features (spread in r2 ≤ 1.02; in absolute f 1.43); collapse ratios 0.77–1.12 (not ≤ 0.5) | — | — | — | in regime; one sequence per IF (sequence SD of X_coh 0.02–0.08) | pre | DIR (locations) / INT (organization) | **WEAKEN**: SUGGESTIVE ONLY (peak alignment without quantitative collapse; no mechanism) |
| N3b | A second maximum of X_coh lies just below f = Ω_LO/2π (r1 = 0.97–0.99; X_coh 0.30 / 0.39 / 0.31), co-located with a CW-gain maximum | Part II scan (not a pre-registered feature) | — | — | — | in regime | post | DIR | **WEAKEN** (post hoc description only; not in the abstract) |
| N4 | The surrogate-failure conclusion is insensitive to the receiver noise level over 0.25–4× (effective SNR at +3 dB from 44–46 dB down to 20–22 dB) | Part IV: (a) and (b) hold at every scale; c = 1 reproduces the archived AIR to ≤ 3e-13 bit | — | this is the noise test | — | in regime | pre | DIR | **KEEP** |
| N5 | The phase effect in the reference configuration is already present below the CW output maximum: X = 0.33 [0.26, 0.41] at CW gain 0.8 (local CW slope +0.53), rising to 0.53 at CW gain 0.6 (slope −0.26). The archived reference level lies beyond the output maximum (slope −0.93) | Part I (C1, matched CW gain; the local slope is a pre-registered descriptor) | HA1 context | n/a | — | in regime | pre (measurement) / post (framing) | DIR | **KEEP** |
| N6 | At 10 MHz and with the weaker probe, E1dB exceeds the LO field (1.16 × and 1.43 × LO), and with signal ≤ LO the CW response expands instead of compressing; at 15 MHz the CW gain stays ≥ 0.81 for signal ≤ LO | CW tables (Part I plan) | — | — | — | — | pre (plan) | DIR | **KEEP** |

## Consequences

**Abstract (KEEP claims only):** A1 (with N4), A2, A4, A5, A6 (with N5), B3 (optional), N1, N2, N6 (as context of N1), and A9b as a scope statement.
- Removed from the abstract: the declustering statement (A3, WEAKEN) and "the slow memory persists with the weaker probe" (A8, WEAKEN).
- A7's E1dB-matched wording is replaced by N1.

**Numerical sanity.** A WEAKEN claim (A3) uses a V-C-flagged state (shuffle excursions above the LO). V-B is material for two matched amplitudes (C1 at gain 0.95, not informative; C6 at 0.8, whose pair holds with margin +0.16 [+0.11, +0.22]). **Class: CAVEAT.**

**Claims removed:** T0 (the title) only.

**Claims restated from E1dB matching to fair matching:** A7 → N1; B1 (fair-matched values); B2 (Part II scan).

**Claims weakened:** A1b, A3, A7, A8, A9a, B2, B5, N3, N3b.

**Claims strengthened:** A1 (noise range, N4), A6 (present below the output maximum, N5; not captured by a standard GMP, N2), N1 (fair-matched configuration dependence replacing a comparison outside the regime), B7 (a standard GMP tested).
