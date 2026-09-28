# Stage-0.5 Decision

**CONDITIONAL GO** for a focused physical validation; no main-paper mechanism claim yet.

# Verified

- CW P1dB is E1dB = 0.08770749 V/m, from the 5 MHz CW fundamental definition. Field and power conversions reproduce the saved Stage-0 rows (max peak-relation error 2.7e-15 dB; RMS-field relation error 9e-17 V/m).
- Matched formats use 5.0 Msymbol/s and mean 99% bandwidths 16QAM 5.365 MHz, CE 5.369 MHz, OFDM 5.356 MHz, QPSK 5.372 MHz. Across four mean bandwidths, maximum mismatch is 0.18%.
- QPSK/16-QAM use identical RRC pulse shaping. OFDM has 68 active subcarriers in a 256 FFT, 16 CP samples in a 272-sample block, and 68 information-bearing symbols per block; its effective data-symbol rate is 5.0 Msymbol/s. The WOLA cyclic suffix overlaps the next prefix and adds no block duration.
- For every format, R_AIR = I_native × 5.0 Msymbol/s; eta_AIR = R_AIR/B_occ. All waveforms have unit complex-envelope RMS before scaling; P/P1dB=(E/E1dB)^2.
- Eight seeds cover -6, -3, 0, +3, +6 dB; four seeds cover -15 and -10 dB. OFDM uses an independent randomized balanced four-symbol training preamble on each carrier; all reported full-atomic runs use Nd=1501.

# Main Results

At +3 dB relative to P1dB, paired mean full-atomic loss from -6 dB reference:

| Format | AIR loss (bit/native symbol) | Relative AIR loss | Rate loss (Mbit/s) | eta loss (bit/s/Hz) |
|---|---:|---:|---:|---:|
| CE | 0.706 [0.438, 0.973] | 88.5% | 3.528 | 0.658 |
| QPSK | 0.632 [0.569, 0.695] | 33.1% | 3.159 | 0.588 |
| 16QAM | 1.037 [0.954, 1.120] | 45.2% | 5.186 | 0.967 |
| OFDM | 0.256 [0.183, 0.330] | 12.8% | 1.281 | 0.239 |

# Static-Nonlinear Control

At +3 dB, full-minus-static-plus-LTI paired excess AIR loss:
- CE: 0.723 bit/native symbol (paired 95% t interval 0.454 to 0.991; n=8).
- QPSK: 0.627 bit/native symbol (paired 95% t interval 0.566 to 0.688; n=8).
- 16QAM: 0.849 bit/native symbol (paired 95% t interval 0.779 to 0.919; n=8).
- OFDM: 0.203 bit/native symbol (paired 95% t interval 0.129 to 0.276; n=8).

Static compression plus the repository small-signal LTI response reproduces part of 16-QAM and OFDM loss, particularly at +6 dB, but does not reproduce the full-atomic moderate-power losses. OFDM converges to the static-plus-LTI loss at +6 dB, so the dynamic excess is not uniform across power.

# Peak-Power View

A common +6 dB peak-power coordinate is within the measured average-power grid for all formats. The per-seed interpolated relative-loss values are reported in the peak-normalization table below.

| Format | Relative AIR loss at common +6 dB peak |
|---|---:|
| CE | 41.9% [3.6%, 80.3%] (n=8) |
| QPSK | 8.9% [0.9%, 16.8%] (n=8) |
| 16QAM | 13.1% [1.3%, 24.8%] (n=8) |
| OFDM | 0.4% [-0.2%, 1.0%] (n=8) |

# Waveform Statistics

PAPR alone is not predictive: OFDM has the largest PAPR but smaller loss than the single-carrier formats at +3 dB. Pooled Spearman rho is -0.10 for PAPR and 0.81 for time above E1dB; power-stratified mean rho for mean dwell is 0.25. Occupancy and dwell have similar, moderate within-power associations, so Stage K was run as a controlled receiver-output diagnostic.

Stage-1.5 context reports an atomic pulse-response 1/e timescale tau_atom ≈ 1.727 µs; dwell and envelope-correlation times are also stored normalized by this inherited timescale. That pulse trace is not in the present checkout, so tau_atom is not re-estimated here.

# Temporal-Dwell Diagnostic (Stage K)

The diagnostic used one deterministic input pair, not a seed ensemble. Both patterns have Pavg/P1dB = -0.90 dB, Ppeak/P1dB = +4.56 dB, 25% occupancy above E1dB, and 99% bandwidths near 240 MHz (relative mismatch 0.002%). Mean high-field dwell is 50 ns for short excursions and 1962 ns for long plateaus; the 90th-percentile dwell is 50 ns and 2000 ns. The normalized amplitude distributions are similar but not identical (two-sample KS distance 0.071).

The normalized high/low baseband contrast changes by 0.001106 for full atomic, 0.000749 for static, and 0.000865 for static-plus-LTI. The LTI control explains most of the contrast change; the full-atomic pair leaves a small additional change of 0.000240. This indicates sensitivity to the temporal pattern in this artificial receiver-output diagnostic, but it does not establish a communication AIR effect or a causal atomic dwell law.

The artificial pair has much wider bandwidth than the communication waveforms and only one realization. Its result is a focused diagnostic, not a headline receiver operating limit.

# Constant-Envelope Control

The BT=0.36 Gaussian CPFSK waveform is constant envelope (PAPR 0 dB), with bandwidth matched within the 2% limit. Its simple discriminator plus the fixed Gaussian auxiliary AIR produces a substantially lower reference rate and high seed variability; its +3 dB interval is broad. Treat CE as a useful envelope control with a receiver limitation, not as decisive evidence about a universal modulation ranking.

# Numerical Validation

- Halving dt from 1 ns to 0.5 ns changes +3 dB AIR by at most 7.1e-05 bit/native symbol across the four spot checks.
- Nd=501 differs from Nd=1501 by 0.308 bit/native symbol at +3 dB. This screening quadrature is not converged; no result at reduced Nd is used as decisive evidence.
- A short 2 us Nd=301 CPU NumPy versus CUDA check is recorded in 05_backend_validation.csv; trace correlation was 0.999995 and maximum probe-response difference was 0.0063 dB.

# Uncertainty

Headline intervals are paired across independent waveform seeds using 95% Student-t intervals (n=8 at -6 through +6 dB; n=4 at -15 and -10 dB). They quantify seed variability in this simulated probe/IF/baseband readout, not hardware uncertainty.

# Information-Loss Thresholds

P5 and P10 are the highest powers on the mean full-atomic curve meeting 5% and 10% relative AIR-loss limits. Intervals are paired seed-bootstrap 95% ranges; the estimates apply only to this simulated receiver.

| Format | P5 - P1dB (dB) | P10 - P1dB (dB) |
|---|---:|---:|
| CE | -1.88 [-2.26, -1.09] | -1.54 [-2.04, -0.51] |
| QPSK | -2.09 [-2.19, -1.96] | -0.47 [-0.65, -0.23] |
| 16QAM | -3.64 [-3.91, -3.30] | -2.37 [-2.56, -2.12] |
| OFDM | 0.75 [0.43, 1.15] | 2.19 [1.71, 2.95] |

# What Is NOT Supported

- No universal communication-aware threshold or receiver-independent material limit.
- No proof that dwell time is causal; high-field occupancy and constellation sensitivity remain intertwined.
- No claim that every static nonlinear or linear dynamic model is excluded; the tested controls are the measured CW complex-fundamental curve and its repository small-signal LTI cascade.
- No experimental validation, quantum advantage, capacity result, or generalization beyond this operating point and receiver chain.

# Scientific Interpretation

The rate and bandwidth correction does not remove the measured waveform dependence. Full-atomic AIR losses exceed the tested static and static-plus-LTI controls for QPSK, 16-QAM, and OFDM at moderate average powers. The strongest 16-QAM point and high-power OFDM show substantial static-control contributions. The separation is therefore evidence that CW P1dB alone does not predict communication loss in this particular full-thermal transient model, but it does not identify a unique causal statistic.

The low-power -6 dB reference is a shared coordinate, not a guarantee that all instantaneous peaks are uncompressed. Peak-normalized curves remain format dependent at the common +6 dB peak coordinate, though the CE receiver uncertainty is large. Occupancy tracks the pooled trend better than PAPR; fixed-power comparisons and the single artificial dwell diagnostic limit mechanism claims.

# Next Step

Validate the matched-rate experiment at a second IF or operating point with a receiver chain calibrated against a realistic probe detector and analog bandwidth; retain the current tested controls and preamble policy.
