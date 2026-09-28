# Part V — numerical and regime checks

Generated 2026-09-28T03:29:20+10:00 (code 4122e09). Rule: `00_plan.md`, Part V.

## V-A trace / positivity vs duration (prefix runs)

| case | T_us | trace_dev_max | min_eig_min | thermal_trace_dev | thermal_min_eig | hermiticity_dev_max | wall_s |
|---|---|---|---|---|---|---|---|
| dwell | 105 | 5.31e-06 | 1.09e-10 | 4.67e-08 | 0.0122 | 6.3e-07 | 5.76 |
| dwell | 210 | 3.62e-06 | 1.09e-10 | 4.58e-08 | 0.0126 | 2.57e-07 | 6.71 |
| dwell | 420 | 4.37e-06 | 1.08e-10 | 6.52e-08 | 0.0122 | 1.18e-06 | 14.7 |
| dwell | 840 | 4.07e-06 | 1.08e-10 | 4.95e-08 | 0.0123 | 4.4e-07 | 25.9 |
| phase | 55 | 2.38e-05 | -7.62e-10 | 7.18e-08 | 0.013 | 2.6e-07 | 4.12 |
| phase | 110 | 9.32e-06 | 1.1e-10 | 1.83e-08 | 0.0127 | 1.77e-07 | 5.74 |
| phase | 165 | 1.45e-05 | 1.08e-10 | 2.38e-09 | 0.0127 | 2.8e-07 | 7.45 |
| phase | 220 | 6.78e-06 | 1.09e-10 | 2.1e-08 | 0.0127 | 9.9e-08 | 8.11 |

Result: {"dwell": {"within_thresholds": true, "accumulation": false, "trace_dev_max_first_last": [5.306977641339472e-06, 4.066489054821432e-06], "min_eig_min_first_last": [1.0887601547181915e-10, 1.0825498616423034e-10]}, "phase": {"within_thresholds": true, "accumulation": false, "trace_dev_max_first_last": [2.384073013672605e-05, 6.782873242627829e-06], "min_eig_min_first_last": [-7.618274438530764e-10, 1.0905230213565785e-10]}}

## Free diagnostics (final states of every stored run)

| group | n_runs | max_thermal_trace_dev | min_thermal_eig | max_class_trace_dev | min_class_eig | max_hermiticity_dev |
|---|---|---|---|---|---|---|
| prior P1 | 84 | 3.49e-07 | 0.00141 | 4.87e-05 | 5.35e-12 | 1.24e-06 |
| prior P2 | 174 | 3.29e-07 | 0.0117 | 7.54e-05 | 9.14e-11 | 1.23e-06 |
| prior P4 | 192 | 2.25e-07 | 0.012 | 2.53e-05 | 1.06e-10 | 1.43e-06 |
| campaign Part I | 245 | 3.03e-07 | 0.00134 | 7.52e-05 | 5.94e-12 | 1.64e-06 |
| campaign Part II | 28 | 3.7e-07 | 0.0119 | 6.06e-05 | 1.04e-10 | 6.63e-07 |

## V-B E1dB slope sensitivity

| config | variant | slope | E1dB_Vpm | E1dB_archived_Vpm | change_dB_power | slope_rel_to_archived |
|---|---|---|---|---|---|---|
| C1 | slope_2_smallest | 0.003504 | 0.07322 | 0.07322 | 0 | 1 |
| C1 | slope_3_smallest | 0.003495 | 0.074 | 0.07322 | 0.09246 | 0.9975 |
| C1 | slope_4_smallest | 0.003481 | 0.07526 | 0.07322 | 0.2388 | 0.9935 |
| C1 | local_fit_4_smallest | 0.003509 | 0.07278 | 0.07322 | -0.0518 | 1.001 |
| C2 | slope_2_smallest | 0.002927 | 0.5808 | 0.5808 | 0 | 1 |
| C2 | slope_3_smallest | 0.002904 | 0.5839 | 0.5808 | 0.04561 | 0.9923 |
| C2 | slope_4_smallest | 0.002886 | 0.5863 | 0.5808 | 0.08192 | 0.9862 |
| C2 | local_fit_4_smallest | 0.002928 | 0.5807 | 0.5808 | -0.002132 | 1 |
| C3 | slope_2_smallest | 0.001486 | 0.2175 | 0.2175 | 0 | 1 |
| C3 | slope_3_smallest | 0.001485 | 0.2181 | 0.2175 | 0.02353 | 0.9996 |
| C3 | slope_4_smallest | 0.001485 | 0.2191 | 0.2175 | 0.06153 | 0.999 |
| C3 | local_fit_4_smallest | 0.001486 | 0.2172 | 0.2175 | -0.01339 | 1 |
| C4 | slope_2_smallest | 1.885e-05 | 0.7131 | 0.7131 | 0 | 1 |
| C4 | slope_3_smallest | 1.891e-05 | 0.7128 | 0.7131 | -0.003855 | 1.003 |
| C4 | slope_4_smallest | 1.899e-05 | 0.7123 | 0.7131 | -0.009922 | 1.007 |
| C4 | local_fit_4_smallest | 1.883e-05 | 0.7133 | 0.7131 | 0.002031 | 0.9986 |
| C5 | slope_2_smallest | 0.004922 | 0.05909 | 0.05909 | 0 | 1 |
| C5 | slope_3_smallest | 0.004902 | 0.06024 | 0.05909 | 0.1685 | 0.9959 |
| C5 | slope_4_smallest | 0.004871 | 0.06193 | 0.05909 | 0.4086 | 0.9895 |
| C5 | local_fit_4_smallest | 0.004933 | 0.05844 | 0.05909 | -0.0951 | 1.002 |
| C6 | slope_2_smallest | 0.006467 | 0.04033 | 0.04033 | 0 | 1 |
| C6 | slope_3_smallest | 0.00641 | 0.04195 | 0.04033 | 0.341 | 0.9911 |
| C6 | slope_4_smallest | 0.006322 | 0.04441 | 0.04033 | 0.8366 | 0.9776 |
| C6 | local_fit_4_smallest | 0.006497 | 0.03942 | 0.04033 | -0.1978 | 1.005 |

Material change of C1's E1dB (> 0.5 dB): False; of a matched I-A amplitude (> 5%): True (max |change| 6.31%).

## V-C rotating-wave / intended-regime check

Simulator inspection (python/utils/transient_quantum.py, unmodified upstream): the RF input is passed as a complex field envelope E(t) = A_LO + a(t) e^{j2pi f_IF t} in the frame rotating at the LO carrier; Omega_RF(t) = mu_MW E(t)/hbar enters the four-level ladder Hamiltonian as Omega_RF/2 and conj(Omega_RF)/2 (rotating-wave form; counter-rotating terms at twice the carrier are absent by construction). The transient solver integrates the full Lindblad equation for each velocity class without linearizing about the LO, so a signal larger than the LO is represented without formal change. The rotating-wave treatment requires |Omega_RF(t)| and the envelope's instantaneous frequency to be small compared with the carrier (the 6.946 GHz transition of the model's source configuration); this holds for every tabulated state (ratio column). What changes when the signal exceeds the LO is the receiver regime: the superheterodyne picture (a weak signal beating against a dominant LO, and the small-signal transfer function linearized about the LO steady state) no longer describes the operating point, which is therefore outside the intended regime of the inherited receiver model. Not quantified here: the four-level truncation (neighbouring Rydberg transitions), whose importance grows with the total RF Rabi frequency.

| source | state | peak_signal_to_LO | fraction_of_samples_signal_gt_LO | peak_total_Rabi_MHz | Rabi_to_carrier | OUTSIDE_INTENDED_REGIME | RWA_QUESTIONABLE | used_by |
|---|---|---|---|---|---|---|---|---|
| revision2 (archived) | C1 reference phase case at 2.4468 E1dB | 0.358 |  | 12.5 | 0.0018 | False | False | manuscript Sec. 7.4-7.5 (C4 claim) |
| revision2 (archived) | C2 reference phase case at 2.4468 E1dB | 2.84 |  | 35.4 | 0.0051 | True | False | manuscript Sec. 7.4-7.5 (C4 claim) |
| revision2 (archived) | C3 reference phase case at 2.4468 E1dB | 1.06 |  | 19 | 0.00274 | True | False | manuscript Sec. 7.4-7.5 (C4 claim) |
| revision2 (archived) | C4 reference phase case at 2.4468 E1dB | 3.49 |  | 41.4 | 0.00596 | True | False | manuscript Sec. 7.4-7.5 (C4 claim) |
| revision (archived) | OP2 reference phase case at 2.4468 E1dB | 0.282 |  | 8.28 | 0.00119 | False | False | manuscript Sec. 7.3 |
| revision (archived) | OP3 reference phase case at 2.4468 E1dB | 0.34 |  | 10.5 | 0.00151 | False | False | manuscript Sec. 7.3 |
| campaign Part I | C1 matched_CW_gain 0.95 | 0.0976 |  | 10.1 | 0.00146 | False | False | Part I |
| campaign Part I | C1 matched_CW_gain 0.9 | 0.14 |  | 10.5 | 0.00151 | False | False | Part I |
| campaign Part I | C1 matched_CW_gain 0.85 | 0.174 |  | 10.8 | 0.00156 | False | False | Part I |
| campaign Part I | C1 matched_CW_gain 0.8 | 0.203 |  | 11.1 | 0.0016 | False | False | Part I |
| campaign Part I | C1 matched_CW_gain 0.6 | 0.293 |  | 11.9 | 0.00172 | False | False | Part I |
| campaign Part I | C1 matched_CW_gain 0.4 | 0.371 |  | 12.6 | 0.00182 | False | False | Part I |
| campaign Part I | C1 matched_signal_to_LO 0.1 | 0.1 |  | 10.1 | 0.00146 | False | False | Part I |
| campaign Part I | C1 matched_signal_to_LO 0.2 | 0.2 |  | 11.1 | 0.00159 | False | False | Part I |
| campaign Part I | C1 matched_signal_to_LO 0.358 | 0.358 |  | 12.5 | 0.0018 | False | False | Part I |
| campaign Part I | C1 matched_signal_to_LO 0.4 | 0.4 |  | 12.9 | 0.00186 | False | False | Part I |
| campaign Part I | C1 matched_signal_to_LO 0.6 | 0.6 |  | 14.8 | 0.00212 | False | False | Part I |
| campaign Part I | C2 matched_signal_to_LO 0.1 | 0.1 |  | 10.1 | 0.00146 | False | False | Part I |
| campaign Part I | C2 matched_signal_to_LO 0.2 | 0.2 |  | 11.1 | 0.00159 | False | False | Part I |
| campaign Part I | C2 matched_signal_to_LO 0.358 | 0.358 |  | 12.5 | 0.0018 | False | False | Part I |
| campaign Part I | C2 matched_signal_to_LO 0.4 | 0.4 |  | 12.9 | 0.00186 | False | False | Part I |
| campaign Part I | C2 matched_signal_to_LO 0.6 | 0.6 |  | 14.8 | 0.00212 | False | False | Part I |
| campaign Part I | C3 matched_CW_gain 0.95 | 0.266 |  | 11.7 | 0.00168 | False | False | Part I |
| campaign Part I | C3 matched_CW_gain 0.9 | 0.404 |  | 12.9 | 0.00186 | False | False | Part I |
| campaign Part I | C3 matched_CW_gain 0.85 | 0.762 |  | 16.3 | 0.00234 | False | False | Part I |
| campaign Part I | C3 matched_signal_to_LO 0.1 | 0.1 |  | 10.1 | 0.00146 | False | False | Part I |
| campaign Part I | C3 matched_signal_to_LO 0.2 | 0.2 |  | 11.1 | 0.00159 | False | False | Part I |
| campaign Part I | C3 matched_signal_to_LO 0.358 | 0.358 |  | 12.5 | 0.0018 | False | False | Part I |
| campaign Part I | C3 matched_signal_to_LO 0.4 | 0.4 |  | 12.9 | 0.00186 | False | False | Part I |
| campaign Part I | C3 matched_signal_to_LO 0.6 | 0.6 |  | 14.8 | 0.00212 | False | False | Part I |
| campaign Part I | C4 matched_signal_to_LO 0.1 | 0.1 |  | 10.1 | 0.00146 | False | False | Part I |
| campaign Part I | C4 matched_signal_to_LO 0.2 | 0.2 |  | 11.1 | 0.00159 | False | False | Part I |
| campaign Part I | C4 matched_signal_to_LO 0.358 | 0.358 |  | 12.5 | 0.0018 | False | False | Part I |
| campaign Part I | C4 matched_signal_to_LO 0.4 | 0.4 |  | 12.9 | 0.00186 | False | False | Part I |
| campaign Part I | C4 matched_signal_to_LO 0.6 | 0.6 |  | 14.8 | 0.00212 | False | False | Part I |
| campaign Part I | C5 matched_CW_gain 0.8 | 0.202 |  | 9.43 | 0.00136 | False | False | Part I |
| campaign Part I | C5 matched_CW_gain 0.6 | 0.354 |  | 10.6 | 0.00153 | False | False | Part I |
| campaign Part I | C5 matched_CW_gain 0.4 | 0.545 |  | 12.1 | 0.00174 | False | False | Part I |
| campaign Part I | C5 matched_signal_to_LO 0.1 | 0.1 |  | 8.62 | 0.00124 | False | False | Part I |
| campaign Part I | C5 matched_signal_to_LO 0.2 | 0.2 |  | 9.41 | 0.00135 | False | False | Part I |
| campaign Part I | C5 matched_signal_to_LO 0.358 | 0.358 |  | 10.6 | 0.00153 | False | False | Part I |
| campaign Part I | C5 matched_signal_to_LO 0.4 | 0.4 |  | 11 | 0.00158 | False | False | Part I |
| campaign Part I | C5 same_absolute_level 0.1 | 0.235 |  | 9.68 | 0.00139 | False | False | Part I |
| campaign Part I | C5 same_absolute_level 0.179 | 0.421 |  | 11.1 | 0.0016 | False | False | Part I |
| campaign Part I | C6 matched_CW_gain 0.8 | 0.169 |  | 7.55 | 0.00109 | False | False | Part I |
| campaign Part I | C6 matched_CW_gain 0.6 | 0.316 |  | 8.49 | 0.00122 | False | False | Part I |
| campaign Part I | C6 matched_CW_gain 0.4 | 0.641 |  | 10.6 | 0.00153 | False | False | Part I |
| campaign Part I | C6 matched_signal_to_LO 0.1 | 0.1 |  | 7.1 | 0.00102 | False | False | Part I |
| campaign Part I | C6 matched_signal_to_LO 0.2 | 0.2 |  | 7.75 | 0.00112 | False | False | Part I |
| campaign Part I | C6 matched_signal_to_LO 0.358 | 0.358 |  | 8.77 | 0.00126 | False | False | Part I |
| campaign Part I | C6 matched_signal_to_LO 0.4 | 0.4 |  | 9.04 | 0.0013 | False | False | Part I |
| campaign Part I | C6 same_absolute_level 0.1 | 0.286 |  | 8.3 | 0.00119 | False | False | Part I |
| campaign Part I | C6 same_absolute_level 0.179 | 0.511 |  | 9.76 | 0.0014 | False | False | Part I |
| campaign Part II | scan, A_LO 0.35, signal/LO 0.358 | 0.358 |  | 8.77 | 0.00126 | False | False | Part II |
| campaign Part II | scan, A_LO 0.425, signal/LO 0.358 | 0.358 |  | 10.6 | 0.00153 | False | False | Part II |
| campaign Part II | scan, A_LO 0.5, signal/LO 0.358 | 0.358 |  | 12.5 | 0.0018 | False | False | Part II |
| fair waveforms (archived, Sec. 3) | CE at -6 dB re P1dB (peak over 8 realizations) | 0.0734 | 0 | 9.9 | 0.00143 | False | False | Sections 3 and Part III/IV |
| fair waveforms (archived, Sec. 3) | CE at +0 dB re P1dB (peak over 8 realizations) | 0.146 | 0 | 10.6 | 0.00152 | False | False | Sections 3 and Part III/IV |
| fair waveforms (archived, Sec. 3) | CE at +3 dB re P1dB (peak over 8 realizations) | 0.207 | 0 | 11.1 | 0.0016 | False | False | Sections 3 and Part III/IV |
| fair waveforms (archived, Sec. 3) | CE at +6 dB re P1dB (peak over 8 realizations) | 0.292 | 0 | 11.9 | 0.00172 | False | False | Sections 3 and Part III/IV |
| fair waveforms (archived, Sec. 3) | QPSK at -6 dB re P1dB (peak over 8 realizations) | 0.135 | 0 | 10.5 | 0.00151 | False | False | Sections 3 and Part III/IV |
| fair waveforms (archived, Sec. 3) | QPSK at +0 dB re P1dB (peak over 8 realizations) | 0.268 | 0 | 11.7 | 0.00168 | False | False | Sections 3 and Part III/IV |
| fair waveforms (archived, Sec. 3) | QPSK at +3 dB re P1dB (peak over 8 realizations) | 0.379 | 0 | 12.7 | 0.00183 | False | False | Sections 3 and Part III/IV |
| fair waveforms (archived, Sec. 3) | QPSK at +6 dB re P1dB (peak over 8 realizations) | 0.536 | 0 | 14.2 | 0.00204 | False | False | Sections 3 and Part III/IV |
| fair waveforms (archived, Sec. 3) | 16QAM at -6 dB re P1dB (peak over 8 realizations) | 0.164 | 0 | 10.7 | 0.00154 | False | False | Sections 3 and Part III/IV |
| fair waveforms (archived, Sec. 3) | 16QAM at +0 dB re P1dB (peak over 8 realizations) | 0.326 | 0 | 12.2 | 0.00176 | False | False | Sections 3 and Part III/IV |
| fair waveforms (archived, Sec. 3) | 16QAM at +3 dB re P1dB (peak over 8 realizations) | 0.461 | 0 | 13.5 | 0.00194 | False | False | Sections 3 and Part III/IV |
| fair waveforms (archived, Sec. 3) | 16QAM at +6 dB re P1dB (peak over 8 realizations) | 0.651 | 0 | 15.2 | 0.00219 | False | False | Sections 3 and Part III/IV |
| fair waveforms (archived, Sec. 3) | OFDM at -6 dB re P1dB (peak over 8 realizations) | 0.232 | 0 | 11.4 | 0.00164 | False | False | Sections 3 and Part III/IV |
| fair waveforms (archived, Sec. 3) | OFDM at +0 dB re P1dB (peak over 8 realizations) | 0.464 | 0 | 13.5 | 0.00194 | False | False | Sections 3 and Part III/IV |
| fair waveforms (archived, Sec. 3) | OFDM at +3 dB re P1dB (peak over 8 realizations) | 0.655 | 0 | 15.3 | 0.0022 | False | False | Sections 3 and Part III/IV |
| fair waveforms (archived, Sec. 3) | OFDM at +6 dB re P1dB (peak over 8 realizations) | 0.925 | 0 | 17.8 | 0.00256 | False | False | Sections 3 and Part III/IV |
| revision2 (archived) | C2 dwell contrast, +3 dB re E1dB(C2) (peak, r0) | 2.85 |  | 35.5 | 0.00511 | True | False | manuscript Sec. 7.4-7.5 (dwell configuration checks) |
| revision2 (archived) | C3 dwell contrast, +3 dB re E1dB(C3) (peak, r0) | 1.07 |  | 19.1 | 0.00275 | True | False | manuscript Sec. 7.4-7.5 (dwell configuration checks) |
| revision2 (archived) | C4 dwell contrast, +3 dB re E1dB(C4) (peak, r0) | 3.5 |  | 41.5 | 0.00597 | True | False | manuscript Sec. 7.4-7.5 (dwell configuration checks) |
| revision2 (archived) | C4 (weaker probe) amplitude step x1.05 at E1dB (upper level) | 1.5 |  | 23 | 0.00332 | True | False | manuscript Sec. 7.5 (memory persists with the weaker probe) |
| revision2 (archived) | C1 amplitude step x1.05 at E1dB (upper level) | 0.154 |  | 10.6 | 0.00153 | False | False | manuscript Sec. 7.5 (standard-probe comparison) |
| dwell/shuffle (archived, Sec. 4) | dwell r0-r7, max over variants, +0 dB | 0.254 | 0 | 11.6 | 0.00167 | False | False | Sections 4-6, Part III |
| dwell/shuffle (archived, Sec. 4) | dwell r0-r7, max over variants, +3 dB | 0.359 | 0 | 12.5 | 0.0018 | False | False | Sections 4-6, Part III |
| dwell/shuffle (archived, Sec. 4) | dwell r0-r7, max over variants, +6 dB | 0.508 | 0 | 13.9 | 0.002 | False | False | Sections 4-6, Part III |
| dwell/shuffle (archived, Sec. 4) | shuffle r0-r7, max over variants, +0 dB | 1.6 | 0.00148 | 23.9 | 0.00345 | True | False | Sections 4-6, Part III |
| dwell/shuffle (archived, Sec. 4) | shuffle r0-r7, max over variants, +3 dB | 2.25 | 0.00495 | 30 | 0.00432 | True | False | Sections 4-6, Part III |
| dwell/shuffle (archived, Sec. 4) | shuffle r0-r7, max over variants, +6 dB | 3.18 | 0.0145 | 38.6 | 0.00556 | True | False | Sections 4-6, Part III |

Stop rules: S2 False; S4 False. Numerical checks pass (excluding the claim-matrix condition): False.
