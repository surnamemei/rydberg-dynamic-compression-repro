# Part III — one standard memory model (GMP) vs the CW-derived surrogates

**Classification: PARTIAL SUCCESS** (rule: `00_plan.md`, section C; generated 2026-09-28T03:26:40+10:00, code c1ba283).

**Selected model** (minimum family-averaged validation NMSE over 30 grid points): nonlinear order K = 7, memory depth 4.0 µs (M = 80 taps at 20 MHz; 80 lagging lags, K_b = 2), ridge 1e-08; 1282 complex coefficients (2564 real parameters). Split (fixed before fitting): train fair/dwell/shuffle r0–r1 + phase seeds 20260705–06 + CW; validation r2–r3 + seed 20260704; test r4–r7 + seeds 20260701–03; held-out drive levels: revision P3 OP1 REF seeds 1–3 at −6…+6 dB re P1dB (+8 dB ≈ training level, reported separately).

## Median test NMSE vs the full model

| Family | GMP | all-zone static | all-zone LTI+static | single-slow-state | GMP ≤ 0.5 best | GMP ≤ 0.8 best |
|---|---:|---:|---:|---:|---|---|
| dwell | 0.1600 | 0.6092 | 0.4971 | 0.6638 | True | True |
| fair | 0.0494 | 0.2938 | 0.0559 | 0.0550 | False | False |
| phase | 0.2005 | 1.1118 | 1.1300 | 1.1300 | True | True |
| shuffle | 0.2183 | 0.8336 | 0.6644 | 0.7153 | True | True |

## Held-out contrasts (test-set means)

| Quantity | Full model | GMP | all-zone static | all-zone LTI+static | single-slow-state |
|---|---:|---:|---:|---:|---:|
| dwell_gain | -0.1063 | -0.1386 | -0.0071 | -0.0079 | -0.1653 |
| decluster_Dref | -0.0026 | -0.0057 | +0.0018 | +0.0032 | -0.0071 |
| X_REF | +0.5029 | +0.1419 | +0.0033 | +0.0007 | +0.0007 |
| X_coh_REF | +0.5413 | +0.1299 | +0.0354 | +0.0003 | +0.0003 |
| X_mag_REF | +0.4544 | +0.0215 | +0.0282 | +0.0003 | +0.0003 |
| X_FAST | +0.3929 | +0.0280 | +0.0053 | +0.0004 | +0.0004 |
| X_coh_FAST | +0.4793 | +0.0290 | +0.0764 | -0.0002 | -0.0002 |
| X_mag_FAST | +0.3134 | -0.1786 | +0.0619 | -0.0003 | -0.0003 |
| X_SLOW | +0.1841 | +0.0941 | +0.0026 | -0.0001 | -0.0001 |
| X_coh_SLOW | +0.2491 | +0.1146 | +0.0111 | -0.0004 | -0.0004 |
| X_mag_SLOW | +0.2148 | +0.0850 | +0.0089 | -0.0005 | -0.0005 |
| X_JUMP | +0.5390 | +0.1685 | +0.0183 | +0.0007 | +0.0007 |
| X_coh_JUMP | +0.5680 | +0.1728 | +0.0786 | +0.0014 | +0.0014 |
| X_mag_JUMP | +0.4699 | +0.0742 | +0.0505 | -0.0004 | -0.0004 |
| X_LONGRAMP | +0.5065 | +0.1167 | -0.0008 | -0.0011 | -0.0011 |
| X_coh_LONGRAMP | +0.5640 | +0.1635 | +0.0095 | -0.0000 | -0.0000 |
| X_mag_LONGRAMP | +0.5061 | +0.0959 | +0.0091 | -0.0001 | -0.0001 |

## Mean |AIR − AIR_full| at +3 dB (bit; archived noise)

| Format | GMP | all-zone static | all-zone LTI+static | single-slow-state |
|---|---:|---:|---:|---:|
| 16QAM | 0.037 | 0.441 | 0.657 | 0.769 |
| CE | 0.804 | 0.850 | 0.544 | 0.544 |
| OFDM | 0.076 | 0.892 | 0.071 | 0.076 |
| QPSK | 0.018 | 0.432 | 0.427 | 0.428 |

## Held-out drive levels (REF case, seeds 1–3; X / X_coh)

| level (dB re P1dB) | full | GMP | all-zone static | all-zone LTI+static | single-slow-state |
|---|---|---|---|---|---|
| -6 | +0.010 / +0.014 | -0.056 / -0.072 | +0.003 / +0.035 | +0.000 / +0.000 | +0.000 / +0.000 |
| -3 | +0.037 / +0.050 | -0.022 / -0.044 | +0.003 / +0.035 | -0.000 / +0.000 | -0.000 / +0.000 |
| +0 | +0.132 / +0.158 | +0.023 / -0.005 | +0.003 / +0.035 | +0.000 / +0.000 | +0.000 / +0.000 |
| +3 | +0.319 / +0.360 | +0.068 / +0.045 | +0.002 / +0.035 | +0.000 / +0.000 | +0.000 / +0.000 |
| +6 | +0.528 / +0.566 | +0.121 / +0.101 | +0.003 / +0.036 | -0.000 / +0.001 | -0.000 / +0.001 |
| +8 | +0.491 / +0.530 | +0.143 / +0.133 | +0.003 / +0.035 | +0.000 / +0.000 | +0.000 / +0.000 |

## Decision

- Key contrasts reproduced (GMP): {'dwell gain contrast': False, 'declustering D_ref contrast': False, 'phase REF (X and X_coh)': False, 'fair AIR at +3 dB': False}.
- Key contrasts reproduced by comparators: {'static_mh': {'dwell gain contrast': False, 'declustering D_ref contrast': False, 'phase REF (X and X_coh)': False, 'fair AIR at +3 dB': False}, 'static_lti_mh': {'dwell gain contrast': False, 'declustering D_ref contrast': False, 'phase REF (X and X_coh)': False, 'fair AIR at +3 dB': False}, 'dsh': {'dwell gain contrast': False, 'declustering D_ref contrast': False, 'phase REF (X and X_coh)': False, 'fair AIR at +3 dB': False}}.
- Contrasts reproduced only by the GMP: [].
- Held-out levels reproduced (X and X_coh within 25%, levels with X_full ≥ 0.05): {0: False, 3: False, 6: False}.
- **Class: PARTIAL SUCCESS.** Stop rule S3: False. Prior-pass P3 rule on this fit: PARTIAL_SUCCESS.

Grid (validation, family-mean NMSE; lower is better): see `gmp_config.csv`.
