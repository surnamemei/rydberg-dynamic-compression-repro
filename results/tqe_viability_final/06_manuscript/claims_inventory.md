# Candidate claims of the current manuscript (inventory for Part VII)

Written 2026-09-28, while the Part I runs were executing and before any Part I–V result was analysed. It lists every candidate claim of the rescoped manuscript (`manuscript/main_draft.md`, commit f02509d; PDFs archived in `00_archive/`) that the campaign's tests bear on, with its current wording and location. Part VII assigns KEEP / WEAKEN / REMOVE to exactly this list, by the rules of `00_plan.md` (section E); claims are not added or dropped after the results.

| ID | Current wording (condensed where long) | Where | Campaign tests that bear on it |
|---|---|---|---|
| T0 | Title: "Configuration-Dependent Dynamic Compression in a Simulated Rydberg Atomic Receiver" | title | I, II, III (title must match the surviving contribution) |
| A1 | In the 5 MHz reference configuration, CW-matched static and small-signal-plus-static surrogates mispredict how much information rate modulation formats lose at equal power (errors of either sign, depending on format and power) | abstract; Sec. 3; C1; Concl. | IV (noise), III (stronger baseline) |
| A1b | Losses depend strongly on format at equal power: at +3 dB, 16-QAM 0.83, QPSK 0.39, CE 0.38, OFDM 0.12 bit; P5 spans five dB | Sec. 3.1 | IV (ordering, P5 vs noise) |
| A2 | Reordering identical field-amplitude samples into longer dwells lowers the fitted large-signal gain from 0.470 to 0.360 of its small-signal value, which the surrogates do not reproduce | abstract; Sec. 4.2; C2 | III |
| A3 | Declustering high-field excursions changes the fixed-reference residual distortion (surrogates do not follow) | abstract; Sec. 4.3; C2 | III |
| A4 | A prespecified distortion-based robustness test of the dwell effect failed | abstract; Sec. 4.2 | none (record) |
| A5 | The response carries a slow memory of 2.5–3.1 µs (slow pole near E1dB) | abstract; Sec. 5.1; C3 | III (memory depth needed by the GMP), V-B (E1dB axis) |
| A6 | At constant amplitude, phase modulation reduces the coherent large-signal gain to 0.47 of its CW value (X_coh = 0.525), mainly through a smaller envelope magnitude (84%) | abstract; Sec. 7.1–7.2 | I (reference state and matched states), III (GMP), V (regime of the reference state) |
| A7 | At levels matched to each configuration's compression point, this phase effect is absent at IFs of 10 and 15 MHz and with a fourfold weaker probe | abstract; Sec. 7.4–7.5; C4; Concl. | I (fair matching), V-C (signal > LO at the archived matched levels) |
| A8 | The slow memory persists with the weaker probe (4.44 µs vs 2.70 µs) | abstract; Sec. 5.1, 7.5; C3 | V-C (the step test's level relative to the LO) |
| A9 | The dynamic large-signal behavior of this receiver is strongly configuration dependent; no microscopic mechanism is identified | abstract; Sec. 8.3; Concl. | I, II |
| B1 | X is smaller at lower LO fields (OP2 / OP3 / OP1: 0.16 / 0.29 / 0.51 for the same three sequences), with frequency selectivity S5 growing likewise (0.56 / 1.22 / 2.02) | Sec. 7.3 | I (HA2, matched CW gain), II (selectivity vs IF and LO) |
| B2 | The large-signal response is sharply frequency selective near the IF in the reference configuration (S5 = 2.02 vs 0.23 / 0.20 at 10 / 15 MHz), and the phase effect accompanies it | Sec. 7.1, 7.4, 8.3 | II (scan), I |
| B3 | Constant detuning and a frequency-dependent static model do not reproduce the phase effect | Sec. 7.1; Table 4 | none (record); II context |
| B4 | Power steps and isolated phase steps oscillate with a common, LO-dependent period (0.560 / 0.558 µs at OP2; 1.60–1.62 / 1.63–1.74 µs at OP1) — post hoc, a shared dynamical time scale only | Sec. 7.3 | II-C |
| B5 | The dwell gain change exceeds the static prediction at 5 and 15 MHz but not at 10 MHz ("the dwell gain contrast holds at 15 MHz only") | Sec. 4.2, 7.4, 8.4 | V-C (levels relative to the LO) |
| B6 | Internal-coherence observation: the output disruption after a +π/2 step is associated with a collapse of the IF component of ρ₄₃; no internal component met the tracking criterion; association only | Sec. 7.6 | none (record) |
| B7 | Behavioural models / CW calibrations should not be transferred between configurations without transient validation; tested structures (memoryless, Hammerstein, single-slow-state, frequency-dependent static) fail in the reference configuration; Wiener/Volterra-type structures untested | Sec. 8.6 | III (a standard GMP now tested) |
| C5 | Numerical lesson: a thermal quadrature adequate for waveform metrics (Nd = 1501) biased the CW compression point by 1.57 dB | Sec. 2.6; C5; Sec. 8.4 | V-B (E1dB rule sensitivity is a separate axis question) |
