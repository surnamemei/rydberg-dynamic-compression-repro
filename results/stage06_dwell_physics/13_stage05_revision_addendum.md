# Stage-0.5 Revision Addendum (from Stage-0.6)

The archived Stage-0.5 bundle (`results/p1db_waveform_stage05`, commit c5d2ea0) is **not modified**. This addendum records which Stage-0.5 statements survive the Stage-0.6 findings and which must change before any paper use.

Evidence:
- `02_nd_convergence.csv`: quadrature convergence.
- `artifacts/controls_Nd{1501,4001,8001}.npz`: CW reference and E1dB at each Nd.
- `09_stage05_recheck.csv`: the same waveforms, seeds, noise, receiver and absolute fields as Stage-0.5, all 8 seeds, at the Stage-0.5 +3 dB and −6 dB points. The full model runs at Nd = 4001; the controls come from the Nd = 4001 tables.

## 1. Power axis (must change)

The CW P1dB field was computed at Nd = 1501, where CW quantities are not converged:

| Nd | E1dB (V/m) |
|---:|---:|
| 1501 | 0.0878 |
| 4001 | 0.0732 |
| 8001 | 0.0727 |

Every Stage-0.5 power coordinate is therefore **+1.57 dB** relative to the converged P1dB. Stage-0.5 "+3 dB" is +4.57 dB, and its "−6 dB" reference is −4.43 dB. The information-loss thresholds (P5, P10 − P1dB) shift by +1.57 dB, provided the Nd = 1501 full-model AIR curves are otherwise accurate (see §2). For example:
- QPSK P5: −2.09 → about −0.52 dB;
- OFDM P10: +2.19 → about +3.76 dB.

The small-signal CW gain at Nd = 1501 is 9.5% low, and the compression-curve shape is off by up to ±10%.

## 2. Full-atomic losses (survive)

Paired mean loss from the −6 dB point to the +3 dB point (Stage-0.5 coordinate), with 95% t-intervals over 8 seeds:

| Format | Full-model AIR at +3, Nd = 4001 (archived 1501) | Loss, Nd = 4001 | Loss, archived Nd = 1501 |
|---|---:|---:|---:|
| CE | 0.034 (0.092) | 0.786 [0.454, 1.118] | 0.706 |
| QPSK | 1.297 (1.277) | 0.624 [0.562, 0.686] | 0.632 |
| 16QAM | 1.275 (1.258) | 1.046 [0.964, 1.128] | 1.037 |
| OFDM | 1.767 (1.744) | 0.233 [0.160, 0.307] | 0.256 |

The waveform-level results move by ≤ 0.08 bit, well inside the Stage-0.5 intervals. The format ordering and the waveform dependence of the loss survive.

## 3. "Full-minus-control" excess (must change for the static control)

The Stage-0.5 static and static+LTI controls model only the IF fundamental zone. The atom's rectified baseline c₀(|E|) is as large as or larger than the fundamental (|c₀/c₁| ≈ 0.5 near E1dB, up to 2.8 at 3·E1dB). At a 5 MHz IF with a MHz-wide envelope, that baseline and the 2nd-harmonic zone fall inside the ±3 MHz receive band. Excess AIR at +3 dB (Stage-0.5 coordinate), paired 95% t-intervals:

| Format | over static (F1) | over static+LTI (F1) | over **all-zone static** | over **all-zone LTI+static** |
|---|---:|---:|---:|---:|
| CE | 0.955 [0.617, 1.294] | 0.668 [0.333, 1.002] | 0.951 [0.612, 1.290] | 0.650 [0.314, 0.985] |
| QPSK | 0.703 [0.641, 0.765] | 0.639 [0.580, 0.698] | **0.227 [0.149, 0.306]** | 0.634 [0.575, 0.693] |
| 16QAM | 2.266 [2.169, 2.363] | 0.887 [0.811, 0.963] | **0.142 [0.008, 0.275]** | 0.867 [0.796, 0.938] |
| OFDM | 0.154 [0.077, 0.232] | 0.158 [0.084, 0.233] | **−1.225 [−1.413, −1.036]** | 0.133 [0.053, 0.212] |

- Against a complete memoryless model, most of the single-carrier "excess over static" disappears. For OFDM the memoryless all-zone model is far *worse* than the full atom, because the real atom low-passes the baseline response that such a model passes through instantly.
- The excess over LTI+static is essentially unchanged by the all-zone extension (0.13–0.87 bit). That is the robust Stage-0.5 control result, and Stage-0.6 shows it has a genuinely dynamic, τ_atom-scale component.
- The Stage-0.5 wording "dynamic excess beyond static compression" should become: "excess beyond the fundamental-zone Hammerstein (static→LTI) model, with a static all-zone control bracketing it".

## 4. Other Stage-0.5 statements

- **τ_atom ≈ 1.727 µs (inherited from Stage-1.5):** superseded. τ_atom = 2.535 µs (τ_1e of the P1dB IF-envelope step, Nd-stable). The linear IF transconductance itself is fast (τ_1e ≈ 47 ns). At Nd = 1501 the linear-point step response contains a spurious 16% slow tail.
- **Stage-K dwell diagnostic:** its full-atomic contrast was inflated about 5% by quadrature error at Nd = 1501 (1.106 → 1.052 ×10⁻³ converged). Its reported excess over static+LTI was therefore inflated by ≈20–25% (0.240 → ≈0.186 ×10⁻³ against the same archived control value). It was also compared only against fundamental-only controls. Superseded by the Stage-0.6 exact-multiset dwell family and shuffle test.
- **"Nd = 501 is not converged":** stands, and extends to Nd = 1501 for CW-derived quantities and for any result at the 0.005-bit level.
