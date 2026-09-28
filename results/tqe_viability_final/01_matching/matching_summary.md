# Part I — fair cross-configuration matching

**Classification: PARTIALLY SURVIVES** (rule: `00_plan.md`, section A; generated 2026-09-28T03:15:37+10:00, code d5d42f1).

Reference QPSK phase case (1 µs symbols, 300 ns ramps, zero drift) at constant amplitude, six sequences (20260701–06) shared by all configurations; mean with 95% Student-t CI (bootstrap CIs in `matching_results.csv`).

## Matched states

| Config | Rule | Target | a (V/m) | a/LO | E/E1dB | CW gain (sim) | local slope | X | X_coh | X_mag | lag (rad) | jitter (rad) |
|---|---|---:|---:|---:|---:|---:|---:|---|---|---|---:|---:|
| C1 | matched_CW_gain | 0.4 | 0.1856 | 0.371 | 2.53 | 0.401 | -1.01 | +0.482 [+0.427, +0.537] | +0.513 [+0.487, +0.540] | +0.447 [+0.416, +0.479] | +0.37 | 0.50 |
| C5 | matched_CW_gain | 0.4 | 0.2317 | 0.545 | 3.92 | 0.400 | -0.28 | +0.462 [+0.412, +0.511] | +0.449 [+0.409, +0.489] | +0.399 [+0.359, +0.439] | +0.46 | 0.42 |
| C6 | matched_CW_gain | 0.4 | 0.2244 | 0.641 | 5.56 | 0.400 | -0.00 | +0.238 [+0.202, +0.273] | +0.224 [+0.186, +0.261] | +0.193 [+0.161, +0.226] | +0.32 | 0.33 |
| C1 | matched_CW_gain | 0.6 | 0.1466 | 0.293 | 2.00 | 0.601 | -0.26 | +0.534 [+0.452, +0.616] | +0.565 [+0.512, +0.618] | +0.524 [+0.468, +0.580] | +0.30 | 0.45 |
| C5 | matched_CW_gain | 0.6 | 0.1504 | 0.354 | 2.55 | 0.599 | +0.31 | +0.273 [+0.197, +0.349] | +0.293 [+0.213, +0.373] | +0.267 [+0.194, +0.339] | +0.19 | 0.27 |
| C6 | matched_CW_gain | 0.6 | 0.1105 | 0.316 | 2.74 | 0.598 | +0.50 | +0.153 [+0.125, +0.182] | +0.163 [+0.134, +0.191] | +0.146 [+0.119, +0.174] | +0.17 | 0.20 |
| C1 | matched_CW_gain | 0.8 | 0.1016 | 0.203 | 1.39 | 0.800 | +0.53 | +0.334 [+0.259, +0.410] | +0.365 [+0.285, +0.445] | +0.347 [+0.272, +0.421] | +0.10 | 0.25 |
| C5 | matched_CW_gain | 0.8 | 0.0861 | 0.202 | 1.46 | 0.799 | +0.63 | +0.128 [+0.106, +0.150] | +0.157 [+0.132, +0.181] | +0.149 [+0.127, +0.171] | +0.06 | 0.14 |
| C6 | matched_CW_gain | 0.8 | 0.0591 | 0.169 | 1.46 | 0.798 | +0.64 | +0.172 [+0.139, +0.205] | +0.181 [+0.153, +0.209] | +0.173 [+0.146, +0.200] | +0.13 | 0.14 |
| C1 | matched_CW_gain | 0.85 | 0.0871 | 0.174 | 1.19 | 0.850 | +0.67 | +0.232 [+0.147, +0.317] | +0.256 [+0.171, +0.342] | +0.245 [+0.165, +0.326] | +0.06 | 0.18 |
| C3 | matched_CW_gain | 0.85 | 0.3810 | 0.762 | 1.75 | 0.853 | +0.84 | -0.013 [-0.017, -0.009] | -0.018 [-0.022, -0.014] | -0.018 [-0.022, -0.014] | -0.00 | 0.03 |
| C1 | matched_CW_gain | 0.9 | 0.0702 | 0.140 | 0.96 | 0.900 | +0.79 | +0.123 [+0.070, +0.175] | +0.139 [+0.083, +0.195] | +0.134 [+0.080, +0.188] | +0.03 | 0.11 |
| C3 | matched_CW_gain | 0.9 | 0.2021 | 0.404 | 0.93 | 0.903 | +0.86 | -0.038 [-0.043, -0.034] | -0.040 [-0.045, -0.035] | -0.040 [-0.045, -0.035] | -0.01 | 0.03 |
| C1 | matched_CW_gain | 0.95 | 0.0488 | 0.098 | 0.67 | 0.950 | +0.89 | +0.034 [+0.017, +0.051] | +0.042 [+0.023, +0.061] | +0.041 [+0.023, +0.059] | +0.01 | 0.05 |
| C3 | matched_CW_gain | 0.95 | 0.1332 | 0.266 | 0.61 | 0.953 | +0.91 | -0.018 [-0.022, -0.013] | -0.017 [-0.021, -0.013] | -0.017 [-0.021, -0.014] | -0.01 | 0.02 |
| C1 | matched_signal_to_LO | 0.1 | 0.0500 | 0.100 | 0.68 | 0.948 | +0.89 | +0.037 [+0.019, +0.056] | +0.046 [+0.026, +0.066] | +0.044 [+0.025, +0.064] | +0.01 | 0.05 |
| C2 | matched_signal_to_LO | 0.1 | 0.0500 | 0.100 | 0.09 | 0.983 | +1.09 | +0.002 [-0.010, +0.014] | +0.006 [-0.001, +0.014] | +0.004 [-0.003, +0.012] | +0.16 | 0.06 |
| C3 | matched_signal_to_LO | 0.1 | 0.0500 | 0.100 | 0.23 | 0.996 | +0.98 | -0.002 [-0.003, -0.001] | -0.002 [-0.003, -0.001] | -0.002 [-0.003, -0.001] | -0.00 | 0.00 |
| C4 | matched_signal_to_LO | 0.1 | 0.0500 | 0.100 | 0.07 | 1.060 | +1.10 | +0.001 [-0.000, +0.002] | +0.001 [-0.000, +0.003] | +0.001 [-0.000, +0.003] | +0.01 | 0.01 |
| C5 | matched_signal_to_LO | 0.1 | 0.0425 | 0.100 | 0.72 | 0.940 | +0.88 | +0.021 [+0.015, +0.027] | +0.028 [+0.020, +0.036] | +0.027 [+0.020, +0.035] | +0.01 | 0.04 |
| C6 | matched_signal_to_LO | 0.1 | 0.0350 | 0.100 | 0.87 | 0.913 | +0.82 | +0.116 [+0.092, +0.141] | +0.125 [+0.104, +0.146] | +0.122 [+0.102, +0.143] | +0.07 | 0.08 |
| C1 | matched_signal_to_LO | 0.2 | 0.1000 | 0.200 | 1.37 | 0.806 | +0.55 | +0.323 [+0.246, +0.400] | +0.353 [+0.271, +0.434] | +0.335 [+0.260, +0.411] | +0.09 | 0.24 |
| C2 | matched_signal_to_LO | 0.2 | 0.1000 | 0.200 | 0.17 | 1.068 | +1.12 | +0.055 [+0.042, +0.069] | +0.057 [+0.044, +0.069] | +0.054 [+0.042, +0.066] | +0.16 | 0.08 |
| C3 | matched_signal_to_LO | 0.2 | 0.1000 | 0.200 | 0.46 | 0.974 | +0.94 | -0.010 [-0.014, -0.006] | -0.010 [-0.013, -0.007] | -0.010 [-0.013, -0.007] | -0.00 | 0.01 |
| C4 | matched_signal_to_LO | 0.2 | 0.1000 | 0.200 | 0.14 | 1.207 | +1.29 | +0.009 [+0.003, +0.014] | +0.010 [+0.005, +0.015] | +0.010 [+0.005, +0.015] | +0.03 | 0.02 |
| C5 | matched_signal_to_LO | 0.2 | 0.0850 | 0.200 | 1.44 | 0.803 | +0.64 | +0.123 [+0.103, +0.143] | +0.152 [+0.129, +0.175] | +0.145 [+0.124, +0.165] | +0.06 | 0.13 |
| C6 | matched_signal_to_LO | 0.2 | 0.0700 | 0.200 | 1.74 | 0.746 | +0.57 | +0.177 [+0.143, +0.211] | +0.185 [+0.155, +0.214] | +0.175 [+0.147, +0.203] | +0.14 | 0.16 |
| C1 | matched_signal_to_LO | 0.358 | 0.1790 | 0.358 | 2.44 | 0.433 | -0.93 | +0.494 [+0.432, +0.557] | +0.526 [+0.493, +0.558] | +0.465 [+0.430, +0.500] | +0.36 | 0.49 |
| C2 | matched_signal_to_LO | 0.358 | 0.1790 | 0.358 | 0.31 | 1.178 | +1.16 | +0.083 [+0.066, +0.099] | +0.085 [+0.071, +0.100] | +0.082 [+0.068, +0.096] | +0.16 | 0.08 |
| C3 | matched_signal_to_LO | 0.358 | 0.1790 | 0.358 | 0.82 | 0.920 | +0.86 | -0.030 [-0.034, -0.025] | -0.032 [-0.036, -0.027] | -0.032 [-0.036, -0.028] | -0.01 | 0.03 |
| C4 | matched_signal_to_LO | 0.358 | 0.1790 | 0.358 | 0.25 | 1.529 | +1.53 | +0.038 [+0.026, +0.050] | +0.040 [+0.030, +0.051] | +0.039 [+0.029, +0.050] | +0.04 | 0.04 |
| C5 | matched_signal_to_LO | 0.358 | 0.1521 | 0.358 | 2.58 | 0.595 | +0.30 | +0.277 [+0.201, +0.352] | +0.297 [+0.217, +0.377] | +0.270 [+0.198, +0.343] | +0.19 | 0.28 |
| C6 | matched_signal_to_LO | 0.358 | 0.1253 | 0.358 | 3.11 | 0.562 | +0.52 | +0.153 [+0.131, +0.175] | +0.158 [+0.134, +0.183] | +0.140 [+0.119, +0.160] | +0.18 | 0.22 |
| C1 | matched_signal_to_LO | 0.4 | 0.2000 | 0.400 | 2.73 | 0.339 | -1.19 | +0.461 [+0.413, +0.510] | +0.483 [+0.454, +0.511] | +0.406 [+0.365, +0.446] | +0.39 | 0.50 |
| C2 | matched_signal_to_LO | 0.4 | 0.2000 | 0.400 | 0.34 | 1.199 | +1.16 | +0.086 [+0.079, +0.093] | +0.089 [+0.082, +0.095] | +0.085 [+0.079, +0.091] | +0.15 | 0.09 |
| C3 | matched_signal_to_LO | 0.4 | 0.2000 | 0.400 | 0.92 | 0.904 | +0.86 | -0.038 [-0.042, -0.033] | -0.039 [-0.044, -0.034] | -0.040 [-0.044, -0.035] | -0.01 | 0.03 |
| C4 | matched_signal_to_LO | 0.4 | 0.2000 | 0.400 | 0.28 | 1.637 | +1.53 | +0.052 [+0.043, +0.062] | +0.056 [+0.048, +0.065] | +0.055 [+0.047, +0.063] | +0.04 | 0.05 |
| C5 | matched_signal_to_LO | 0.4 | 0.1700 | 0.400 | 2.88 | 0.548 | +0.19 | +0.335 [+0.250, +0.419] | +0.351 [+0.268, +0.434] | +0.318 [+0.243, +0.393] | +0.26 | 0.32 |
| C6 | matched_signal_to_LO | 0.4 | 0.1400 | 0.400 | 3.47 | 0.533 | +0.53 | +0.161 [+0.139, +0.182] | +0.172 [+0.148, +0.196] | +0.150 [+0.130, +0.170] | +0.20 | 0.24 |
| C1 | matched_signal_to_LO (pathological) | 0.6 | 0.3000 | 0.600 | 4.10 | 0.212 | +1.56 | +0.397 [+0.351, +0.443] | +0.413 [+0.372, +0.454] | +0.290 [+0.232, +0.348] | +0.29 | 0.60 |
| C2 | matched_signal_to_LO | 0.6 | 0.3000 | 0.600 | 0.52 | 1.259 | +1.04 | +0.091 [+0.078, +0.103] | +0.096 [+0.082, +0.111] | +0.093 [+0.079, +0.106] | +0.14 | 0.09 |
| C3 | matched_signal_to_LO | 0.6 | 0.3000 | 0.600 | 1.38 | 0.879 | +0.94 | -0.014 [-0.019, -0.010] | -0.022 [-0.026, -0.018] | -0.022 [-0.026, -0.019] | +0.02 | 0.03 |
| C4 | matched_signal_to_LO (pathological) | 0.6 | 0.3000 | 0.600 | 0.42 | 1.879 | +1.03 | +0.032 [+0.024, +0.039] | +0.043 [+0.032, +0.055] | +0.041 [+0.031, +0.051] | +0.01 | 0.08 |
| C5 | same_absolute_level | 0.1 | 0.1000 | 0.235 | 1.69 | 0.752 | +0.56 | +0.181 [+0.127, +0.235] | +0.213 [+0.156, +0.270] | +0.200 [+0.148, +0.253] | +0.10 | 0.18 |
| C6 | same_absolute_level | 0.1 | 0.1000 | 0.286 | 2.48 | 0.629 | +0.50 | +0.165 [+0.134, +0.196] | +0.171 [+0.144, +0.199] | +0.156 [+0.129, +0.183] | +0.16 | 0.19 |
| C5 | same_absolute_level | 0.179 | 0.1790 | 0.421 | 3.03 | 0.525 | +0.12 | +0.364 [+0.278, +0.450] | +0.375 [+0.294, +0.455] | +0.340 [+0.266, +0.413] | +0.30 | 0.33 |
| C6 | same_absolute_level | 0.179 | 0.1790 | 0.511 | 4.44 | 0.472 | +0.40 | +0.188 [+0.157, +0.219] | +0.199 [+0.161, +0.237] | +0.173 [+0.140, +0.206] | +0.26 | 0.29 |

Infeasible (not run, not extrapolated):

- C2 matched_CW_gain 0.95: infeasible: needs signal/LO = 1.11 > 1 (outside the superheterodyne regime)
- C2 matched_CW_gain 0.9: infeasible: needs signal/LO = 1.15 > 1 (outside the superheterodyne regime)
- C2 matched_CW_gain 0.85: infeasible: needs signal/LO = 1.20 > 1 (outside the superheterodyne regime)
- C2 matched_CW_gain 0.8: infeasible: needs signal/LO = 1.23 > 1 (outside the superheterodyne regime)
- C2 matched_CW_gain 0.6: infeasible: needs signal/LO = 1.39 > 1 (outside the superheterodyne regime)
- C2 matched_CW_gain 0.4: infeasible: needs signal/LO = 1.63 > 1 (outside the superheterodyne regime)
- C3 matched_CW_gain 0.8: infeasible: needs signal/LO = 1.05 > 1 (outside the superheterodyne regime)
- C3 matched_CW_gain 0.6: infeasible: no crossing on the CW table grid (max amplitude 0.700 V/m = 1.40 x LO)
- C3 matched_CW_gain 0.4: infeasible: no crossing on the CW table grid (max amplitude 0.700 V/m = 1.40 x LO)
- C4 matched_CW_gain 0.95: infeasible: needs signal/LO = 1.41 > 1 (outside the superheterodyne regime)
- C4 matched_CW_gain 0.9: infeasible: needs signal/LO = 1.42 > 1 (outside the superheterodyne regime)
- C4 matched_CW_gain 0.85: infeasible: needs signal/LO = 1.44 > 1 (outside the superheterodyne regime)
- C4 matched_CW_gain 0.8: infeasible: needs signal/LO = 1.45 > 1 (outside the superheterodyne regime)
- C4 matched_CW_gain 0.6: infeasible: needs signal/LO = 1.52 > 1 (outside the superheterodyne regime)
- C4 matched_CW_gain 0.4: infeasible: needs signal/LO = 1.62 > 1 (outside the superheterodyne regime)

## Paired tests (X; reference C1 minus configuration; same sequences)

| Conv. | Rule | Target | Config | Hyp. | X_C1 | X_c | X_C1 − X_c [95% CI] | evaluable | holds |
|---|---|---:|---|---|---:|---:|---|---|---|
| A | matched_CW_gain | 0.4 | C5 | HA2 | +0.482 | +0.462 | +0.020 [-0.055, +0.095] | yes | no |
| A | matched_CW_gain | 0.4 | C6 | HA2 | +0.482 | +0.238 | +0.245 [+0.174, +0.315] | yes | yes |
| A | matched_CW_gain | 0.6 | C5 | HA2 | +0.534 | +0.273 | +0.261 [+0.158, +0.365] | yes | yes |
| A | matched_CW_gain | 0.6 | C6 | HA2 | +0.534 | +0.153 | +0.381 [+0.320, +0.441] | yes | yes |
| A | matched_CW_gain | 0.8 | C5 | HA2 | +0.334 | +0.128 | +0.207 [+0.141, +0.272] | yes | yes |
| A | matched_CW_gain | 0.8 | C6 | HA2 | +0.334 | +0.172 | +0.162 [+0.107, +0.217] | yes | yes |
| A | matched_CW_gain | 0.85 | C3 | HA1 | +0.232 | -0.013 | +0.245 [+0.159, +0.331] | yes | yes |
| A | matched_CW_gain | 0.9 | C3 | HA1 | +0.123 | -0.038 | +0.161 [+0.108, +0.214] | yes | yes |
| A | matched_CW_gain | 0.95 | C3 | HA1 | +0.034 | -0.018 | +0.052 [+0.036, +0.068] | no | yes |
| B | matched_signal_to_LO | 0.1 | C2 | HA1 | +0.037 | +0.002 | +0.035 [+0.017, +0.053] | no | no |
| B | matched_signal_to_LO | 0.1 | C3 | HA1 | +0.037 | -0.002 | +0.040 [+0.021, +0.058] | no | no |
| B | matched_signal_to_LO | 0.1 | C4 | HA1 | +0.037 | +0.001 | +0.036 [+0.018, +0.054] | no | no |
| B | matched_signal_to_LO | 0.1 | C5 | HA2 | +0.037 | +0.021 | +0.016 [-0.004, +0.036] | no | no |
| B | matched_signal_to_LO | 0.1 | C6 | HA2 | +0.037 | +0.116 | -0.079 [-0.100, -0.057] | no | no |
| B | matched_signal_to_LO | 0.2 | C2 | HA1 | +0.323 | +0.055 | +0.268 [+0.195, +0.340] | yes | yes |
| B | matched_signal_to_LO | 0.2 | C3 | HA1 | +0.323 | -0.010 | +0.333 [+0.257, +0.409] | yes | yes |
| B | matched_signal_to_LO | 0.2 | C4 | HA1 | +0.323 | +0.009 | +0.315 [+0.243, +0.387] | yes | yes |
| B | matched_signal_to_LO | 0.2 | C5 | HA2 | +0.323 | +0.123 | +0.200 [+0.134, +0.266] | yes | yes |
| B | matched_signal_to_LO | 0.2 | C6 | HA2 | +0.323 | +0.177 | +0.146 [+0.072, +0.220] | yes | yes |
| B | matched_signal_to_LO | 0.358 | C2 | HA1 | +0.494 | +0.083 | +0.412 [+0.340, +0.484] | yes | yes |
| B | matched_signal_to_LO | 0.358 | C3 | HA1 | +0.494 | -0.030 | +0.524 [+0.460, +0.589] | yes | yes |
| B | matched_signal_to_LO | 0.358 | C4 | HA1 | +0.494 | +0.038 | +0.456 [+0.395, +0.517] | yes | yes |
| B | matched_signal_to_LO | 0.358 | C5 | HA2 | +0.494 | +0.277 | +0.218 [+0.115, +0.321] | yes | yes |
| B | matched_signal_to_LO | 0.358 | C6 | HA2 | +0.494 | +0.153 | +0.341 [+0.292, +0.391] | yes | yes |
| B | matched_signal_to_LO | 0.4 | C2 | HA1 | +0.461 | +0.086 | +0.376 [+0.326, +0.425] | yes | yes |
| B | matched_signal_to_LO | 0.4 | C3 | HA1 | +0.461 | -0.038 | +0.499 [+0.448, +0.550] | yes | yes |
| B | matched_signal_to_LO | 0.4 | C4 | HA1 | +0.461 | +0.052 | +0.409 [+0.354, +0.464] | yes | yes |
| B | matched_signal_to_LO | 0.4 | C5 | HA2 | +0.461 | +0.335 | +0.127 [+0.011, +0.243] | yes | yes |
| B | matched_signal_to_LO | 0.4 | C6 | HA2 | +0.461 | +0.161 | +0.301 [+0.247, +0.354] | yes | yes |
| B | matched_signal_to_LO | 0.6 | C2 | HA1 | +0.397 | +0.091 | +0.306 [+0.258, +0.354] | no | yes |
| B | matched_signal_to_LO | 0.6 | C3 | HA1 | +0.397 | -0.014 | +0.411 [+0.367, +0.455] | no | yes |
| B | matched_signal_to_LO | 0.6 | C4 | HA1 | +0.397 | +0.032 | +0.365 [+0.320, +0.411] | no | yes |
| C | same_absolute_level | 0.1 | C5 | HA2 (descriptive) | +0.323 | +0.181 | +0.142 [+0.053, +0.230] | no | yes |
| C | same_absolute_level | 0.1 | C6 | HA2 (descriptive) | +0.323 | +0.165 | +0.158 [+0.102, +0.215] | no | yes |
| C | same_absolute_level | 0.179 | C5 | HA2 (descriptive) | +0.494 | +0.364 | +0.131 [+0.036, +0.225] | no | yes |
| C | same_absolute_level | 0.179 | C6 | HA2 (descriptive) | +0.494 | +0.188 | +0.307 [+0.234, +0.380] | no | yes |

## Decision

- HA1 (headline; C2–C4): evaluable informative pairs A = 2, B = 9; holding 11; failing []; reversals [].
- HA2 (operating point; C5, C6): evaluable pairs 12; holding 11; failing ['A:matched_CW_gain=0.4:C5']; reversals [].
- **Class: PARTIALLY SURVIVES.** With X_coh in place of X: SURVIVES.
- Prior-pass P1 rule (seeds 1–3, including ratio 0.6): CONFIGURATION_DEPENDENCE_SURVIVES.
- Operating-point claim reading (C1 vs C6 at matched CW gain 0.8/0.6/0.4): KEEP.
- Stop rule S1 triggered: False.
