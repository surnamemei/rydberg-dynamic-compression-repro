# Converged Regeneration of the Fair-Waveform Communication Grid

This is a mechanical recomputation (protocol `00_protocol.json`, committed before the run), not a hypothesis test.
- 192 full-atomic runs at Nd = 4001, dt = 1 ns, with 0 failures.
- Power axis: relative to the converged E1dB = 0.0732 V/m.
- The Stage-0.5 waveforms, seeds, noise, receiver and CE-timing rule are unchanged.
- The archive `results/p1db_waveform_stage05` is untouched.

## Full-model AIR loss vs the −6 dB point

Converged axis; paired mean [95% CI]; 8 seeds (4 at −15 and −10 dB). Source: `02_losses.csv`.

| Format | 0 dB | +3 dB | +6 dB |
|---|---:|---:|---:|
| CE | −0.11 [−0.17, −0.04] | 0.38 [0.10, 0.66] | 0.77 [0.63, 0.91] |
| QPSK | 0.10 [0.09, 0.11] | 0.39 [0.36, 0.42] | 0.83 [0.76, 0.90] |
| 16QAM | 0.32 [0.28, 0.37] | 0.83 [0.72, 0.94] | 1.28 [1.20, 1.37] |
| OFDM | 0.01 [−0.00, 0.02] | 0.12 [0.05, 0.18] | 0.24 [0.18, 0.30] |

CE gains AIR up to 0 dB. Its low-power AIR is noise-limited under the simple discriminator receiver, as in Stage-0.5.

## Excess degradation over CW-derived controls

AIR(control) − AIR(full), paired. Source: `03_excess_over_controls.csv`.

| Format | Power | All-zone static | All-zone LTI+static |
|---|---|---:|---:|
| QPSK | 0 dB | +0.17 [+0.16, +0.19] | +0.12 [+0.10, +0.13] |
| | +3 dB | +0.41 [+0.37, +0.44] | +0.41 [+0.38, +0.44] |
| | +6 dB | +0.01 [−0.08, +0.11] | +0.82 [+0.75, +0.88] |
| 16QAM | 0 dB | +1.49 [+1.39, +1.59] | +0.34 [+0.29, +0.39] |
| | +3 dB | +0.61 [+0.44, +0.77] | +0.75 [+0.64, +0.85] |
| | +6 dB | −0.02 [−0.15, +0.10] | +0.89 [+0.82, +0.96] |
| CE | 0 dB | +0.11 [+0.03, +0.18] | −0.19 [−0.24, −0.14] |
| | +3 dB | +0.59 [+0.31, +0.87] | +0.29 [+0.02, +0.57] |
| | +6 dB | +0.98 [+0.84, +1.12] | +0.67 [+0.53, +0.82] |
| OFDM | 0 dB | −0.19 [−0.26, −0.11] | +0.01 [−0.00, +0.02] |
| | +3 dB | −0.88 [−0.97, −0.79] | +0.11 [+0.05, +0.17] |
| | +6 dB | −1.42 [−1.54, −1.30] | −0.11 [−0.19, −0.02] |

- Neither CW-derived control tracks the full model across power and format.
- The all-zone static model agrees with QPSK and 16-QAM only at +6 dB, which looks coincidental. It underpredicts their degradation at 0 and +3 dB, and overpredicts OFDM's at every power.
- The LTI+static model underpredicts QPSK, 16-QAM and CE degradation from +3 dB up, by 0.3–0.9 bit.

## Relative-loss thresholds on the converged axis

Highest power on the mean curve meeting the limit; paired bootstrap 95%. Source: `04_thresholds.csv`.

| Format | P5 − P1dB (dB) | P10 − P1dB (dB) | Archived P5 / P10 + 1.57 dB (prediction from the addendum) |
|---|---:|---:|---:|
| CE | 0.90 [0.51, 1.65] | 1.14 [0.70, 2.10] | −0.31 / 0.03 |
| QPSK | −0.05 [−0.34, 0.08] | 0.96 [0.86, 1.11] | −0.52 / 1.10 |
| 16QAM | −2.47 [−2.65, −2.26] | −1.07 [−1.39, −0.69] | −2.07 / −0.80 |
| OFDM | 2.56 [1.67, 3.58] | 5.04 [4.08, 6.00] | 2.32 / 3.76 |

## Consistency with the archive

The archived Nd = 1501 loss curves, shifted +1.57 dB onto the converged axis, closely overlay the regenerated ones (`fig_regen_air_loss.png`; approximate, since the archived −6 dB reference sits 1.57 dB higher). This confirms the addendum's conclusion: the archived waveform-level results were nearly right, and the dominant error was the power axis. Threshold differences are within about 0.5 dB, except for CE (receiver-limited, high variance) and OFDM P10 (a shallow curve).

## Consequence for claim C1

**Strengthened.** On the full converged grid, neither CW-derived control reproduces the full-model degradation across powers and formats. The earlier impression, from a single recheck point (+4.6 dB), that the all-zone static model "nearly predicts single-carrier losses" does not hold across power.
