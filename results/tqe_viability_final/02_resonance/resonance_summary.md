# Part II — LO / dressed-state resonance hypothesis

**Classification: SUGGESTIVE ONLY** (rule: `00_plan.md`, section B; generated 2026-09-28T03:21:02+10:00, code 13fc325).

Matched state: signal/LO = 0.358 at every IF and LO field (standard probe); CW + one reference sequence (20260701) per IF; {0.35: 35, 0.425: 31, 0.5: 35} IFs per LO field. Sequence-to-sequence SD at 5 MHz (Part I, six sequences): LO 0.5: X 0.060, X_coh 0.031, X_mag 0.033; LO 0.425: X 0.072, X_coh 0.076, X_mag 0.069; LO 0.35: X 0.021, X_coh 0.023, X_mag 0.020.

## Feature locations (grid points; edge features ineligible)

| LO (V/m) | Ω_LO/2π (MHz) | feature | IF (MHz) | r1 | r2 | d1 | d2 | s3 | d3 | value |
|---|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.35 | 6.456 | g_min | 3.333 | 0.516 | 1.033 | -3.122 | +0.105 | 0.984 | -0.053 | 0.338 |
| 0.35 | 6.456 | S_max | 3.333 | 0.516 | 1.033 | -3.122 | +0.105 | 0.984 | -0.053 | 0.709 |
| 0.35 | 6.456 | Xcoh_max | 3.571 | 0.553 | 1.106 | -2.884 | +0.344 | 1.055 | +0.185 | 0.407 |
| 0.425 | 7.839 | g_min | 4.000 | 0.510 | 1.021 | -3.839 | +0.080 | 0.987 | -0.051 | 0.326 |
| 0.425 | 7.839 | S_max | 4.000 | 0.510 | 1.021 | -3.839 | +0.080 | 0.987 | -0.051 | 0.714 |
| 0.425 | 7.839 | Xcoh_max | 4.348 | 0.555 | 1.109 | -3.491 | +0.428 | 1.073 | +0.296 | 0.475 |
| 0.5 | 9.222 | g_min | 4.762 | 0.516 | 1.033 | -4.460 | +0.151 | 1.008 | +0.038 | 0.220 |
| 0.5 | 9.222 | S_max | 4.762 | 0.516 | 1.033 | -4.460 | +0.151 | 1.008 | +0.038 | 1.184 |
| 0.5 | 9.222 | Xcoh_max | 5.102 | 0.553 | 1.106 | -4.120 | +0.491 | 1.080 | +0.378 | 0.546 |

## Alignment and collapse

| normalization | g_min aligned | S_max aligned | X_coh-max aligned | collapse ratio ln g_CW | collapse ratio X_coh |
|---|---|---|---|---:|---:|
| r1 | True | True | True | 1.12 | 0.89 |
| r2 | True | True | True | 1.12 | 0.89 |
| d1 | False | False | False | 1.42 | 1.25 |
| d2 | True | True | True | 0.99 | 0.77 |
| s3 | True | True | True | 1.09 | 0.85 |
| d3 | True | True | True | 0.98 | 0.79 |

Normalizations meeting CLEAR: []; meeting SUGGESTIVE: ['r1', 'r2', 'd2', 's3', 'd3'].
Prior-pass P2 rule: SUGGESTIVE_ONLY.

## II-C transient oscillation frequencies vs candidate scales (descriptive; match = within 10%)

| LO (V/m) | observation | observed (MHz) | candidate | candidate (MHz) | rel. diff | match |
|---|---|---|---|---:|---:|---|
| 0.5 | power steps | 0.617/0.625 | Omega_LO/2pi | 9.222 | 13.76 |  |
| 0.5 | power steps | 0.617/0.625 | Omega_LO/4pi | 4.611 | 6.38 |  |
| 0.5 | power steps | 0.617/0.625 | Omega_c/2pi | 2.050 | 2.28 |  |
| 0.5 | power steps | 0.617/0.625 | Omega_p/2pi | 8.080 | 11.93 |  |
| 0.5 | power steps | 0.617/0.625 | Omega_3/2pi | 9.447 | 14.12 |  |
| 0.5 | power steps | 0.617/0.625 | Omega_3/4pi | 4.724 | 6.56 |  |
| 0.5 | power steps | 0.617/0.625 | |f_IF - Omega_LO/2pi| | 4.222 | 5.76 |  |
| 0.5 | power steps | 0.617/0.625 | |f_IF - Omega_LO/4pi| | 0.389 | 0.37 |  |
| 0.5 | power steps | 0.617/0.625 | |2 f_IF - Omega_LO/2pi| | 0.778 | 0.24 |  |
| 0.5 | phase steps | 0.575/0.613 | Omega_LO/2pi | 9.222 | 14.03 |  |
| 0.5 | phase steps | 0.575/0.613 | Omega_LO/4pi | 4.611 | 6.52 |  |
| 0.5 | phase steps | 0.575/0.613 | Omega_c/2pi | 2.050 | 2.34 |  |
| 0.5 | phase steps | 0.575/0.613 | Omega_p/2pi | 8.080 | 12.17 |  |
| 0.5 | phase steps | 0.575/0.613 | Omega_3/2pi | 9.447 | 14.40 |  |
| 0.5 | phase steps | 0.575/0.613 | Omega_3/4pi | 4.724 | 6.70 |  |
| 0.5 | phase steps | 0.575/0.613 | |f_IF - Omega_LO/2pi| | 4.222 | 5.88 |  |
| 0.5 | phase steps | 0.575/0.613 | |f_IF - Omega_LO/4pi| | 0.389 | 0.32 |  |
| 0.5 | phase steps | 0.575/0.613 | |2 f_IF - Omega_LO/2pi| | 0.778 | 0.27 |  |
| 0.35 | power steps | 1.786 | Omega_LO/2pi | 6.456 | 2.62 |  |
| 0.35 | power steps | 1.786 | Omega_LO/4pi | 3.228 | 0.81 |  |
| 0.35 | power steps | 1.786 | Omega_c/2pi | 2.050 | 0.15 |  |
| 0.35 | power steps | 1.786 | Omega_p/2pi | 8.080 | 3.52 |  |
| 0.35 | power steps | 1.786 | Omega_3/2pi | 6.773 | 2.79 |  |
| 0.35 | power steps | 1.786 | Omega_3/4pi | 3.387 | 0.90 |  |
| 0.35 | power steps | 1.786 | |f_IF - Omega_LO/2pi| | 1.456 | 0.18 |  |
| 0.35 | power steps | 1.786 | |f_IF - Omega_LO/4pi| | 1.772 | 0.01 | yes |
| 0.35 | power steps | 1.786 | |2 f_IF - Omega_LO/2pi| | 3.544 | 0.98 |  |
| 0.35 | phase steps | 1.792 | Omega_LO/2pi | 6.456 | 2.60 |  |
| 0.35 | phase steps | 1.792 | Omega_LO/4pi | 3.228 | 0.80 |  |
| 0.35 | phase steps | 1.792 | Omega_c/2pi | 2.050 | 0.14 |  |
| 0.35 | phase steps | 1.792 | Omega_p/2pi | 8.080 | 3.51 |  |
| 0.35 | phase steps | 1.792 | Omega_3/2pi | 6.773 | 2.78 |  |
| 0.35 | phase steps | 1.792 | Omega_3/4pi | 3.387 | 0.89 |  |
| 0.35 | phase steps | 1.792 | |f_IF - Omega_LO/2pi| | 1.456 | 0.19 |  |
| 0.35 | phase steps | 1.792 | |f_IF - Omega_LO/4pi| | 1.772 | 0.01 | yes |
| 0.35 | phase steps | 1.792 | |2 f_IF - Omega_LO/2pi| | 3.544 | 0.98 |  |

Candidates matching at both LO fields (consistent): none. No match is called causal.

Refinement: {"0.35": {"Xcoh_max": 0.40662716706539426, "at_IF_MHz": 3.5714285714285716, "edge": false, "max_adjacent_spacing_MHz": 0.27472527472527464, "triggered": true, "new_periods": [297, 288, 272, 265]}, "0.425": {"Xcoh_max": 0.4746561304632061, "at_IF_MHz": 4.3478260869565215, "edge": false, "max_adjacent_spacing_MHz": 0.19762845849802435, "triggered": false, "new_periods": []}, "0.5": {"Xcoh_max": 0.5300112045230763, "at_IF_MHz": 5.0, "edge": false, "max_adjacent_spacing_MHz": 0.2631578947368425, "triggered": true, "new_periods": [208, 204, 196, 192]}}.
