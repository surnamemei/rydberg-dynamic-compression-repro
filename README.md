# Reproducibility Package: Large-Signal Memory and Model Limits in a Simulated Rydberg Atomic Receiver

**Version 1.0.0**

This repository holds the code, configurations, random seeds, time-stamped decision-rule records and processed results behind the article

> J. Mei, "Large-Signal Memory and Model Limits in a Simulated Rydberg Atomic Receiver," manuscript, 2026.

**Archive:** https://doi.org/<DOI>. **Repository:** <REPOSITORY_URL>. **License:** MIT (see `LICENSE`); the licence covers the study-specific files, not the upstream simulator.

**Quick check:** `python3 fetch_upstream.py && python3 verify_reproduction.py` (see `REPRODUCE.md`). **Provenance:** see `PROVENANCE.md`.

## What this repository is (and is not)

**The simulator is not included.** All simulations used the open-source transient density-matrix simulator of J. Zhu and L. Dai, unmodified:
- Repository: **RydbergComms**, https://github.com/johnzja/RydbergComms
- Commit: **`4095724`**
- License: MIT
- Accompanying paper: J. Zhu and L. Dai, *IEEE Trans. Wireless Commun.*, vol. 25, pp. 8292–8307, 2026, doi:10.1109/TWC.2025.3637029.

The files in this repository sit inside the upstream directory tree. `fetch_upstream.py` adds the upstream files at commit `4095724` next to them. The original project was built the same way: it only added files to the upstream tree and modified none.

**What is included** (original relative paths are kept, so the scripts run unchanged):

| Path | Contents |
|---|---|
| `python/experiments/p1db_waveform/` | Fair-waveform study: waveform generation, receivers, AIR estimation, analysis (`stage05.py`, `analyze_stage05.py`) |
| `python/experiments/stage06_dwell/` | Controlled test signals (`waveforms06.py`, `make_config.py`), surrogate/control models (`models06.py`, `build_controls06.py`, `build_controls_mh.py`, `dsh06.py`), run drivers (`run06.py`, `run_*.sh`), convergence and time-scale studies (`nd_convergence.py`, `tau_atom.py`, `fu_cw_steps.py`), analysis (`analyze06.py`, `analyze_followup.py`) |
| `python/experiments/stage06_dwell/phase_mechanism/` | Constant-amplitude phase-modulation experiments (`pm_run.py`, `pm_analyze.py`) |
| `python/experiments/final_validation/` | Converged regeneration of the fair-waveform study (`regen_stage05.py`, `analyze_regen.py`) and the second-operating-point replication (`fv_controls.py`, `fv_replication.py`, `fv_analyze_replication.py`) |
| `results/` | Processed results (CSV/JSON/NPZ), per-run job configurations with **random seeds**, run logs and diagnostic plots, as archived. Includes the thirteen **decision-rule records** (the earlier ones keep their historical file names, such as `*preregistration*.json`); copies and their provenance are in `decision_records/`. Three run manifests are sanitized (see `provenance/SANITIZATION.md`) |
| `manuscript/posthoc/` | Post-hoc analyses reported as such in the article, with their outputs |
| `python/experiments/revision/`, `python/experiments/revision2/` | Pre-submission revision and final receiver-physics checks (see the sections below) |
| `python/experiments/tqe_viability/`, `python/experiments/tqe_viability_final/` | Final adversarial-validation campaign: fair cross-configuration matching, IF scan, identified memory model (GMP), noise sensitivity, numerical and regime checks (see the section below) |
| `manuscript/campaign_tables.py` | Renders Supplementary Tables S11–S15 from the campaign results (used by the ledger and the supplement) |
| `manuscript/figures/src/` | Figure scripts (`make_all.py` builds all 20 figures and runs `audit_figures.py`) |
| `manuscript/figures/final/`, `manuscript/figures/supplement/` | Reference PDFs of the published figures |
| `manuscript/figures/audit/*_plotted_values.csv` | Reference values plotted in each figure |
| `manuscript/build_numbers_ledger.py`, `manuscript/numbers_ledger.csv` | Every number quoted in the article's text and tables, with the file it comes from |
| `p1db_comm_stage0.py`, `stage1_qpsk.py`, `p1db_reference.py`, `p1db_linear_kernel.py`, `p1db_ofdm_fair.py`, `artifacts/*.npz`, root `*.csv` | Screening-stage (N<sub>d</sub> = 1501) modules imported by the later stages, and the fixed inputs they load (CW reference, small-signal kernel, receiver noise level) |
| `decision_records/` | Byte-identical copies of the thirteen decision-rule records, and their provenance |
| `provenance/` | Sanitization record and archival checksums |
| `verify_reproduction.py`, `REPRODUCE.md`, `PROVENANCE.md` | One-command reproduction check, reproduction instructions, provenance note |
| `SHA256SUMS` | SHA-256 of every file (`sha256sum -c SHA256SUMS`) |

**Not included:**
- the upstream simulator;
- raw per-velocity-class or full-length time traces, beyond the stored traces the figures use;
- the private project history.

The receiver noise level in `artifacts/stage1_5_frozen.npz` was fixed in an early pilot stage. The pilot script is not included, because it depends on a module that is not part of this project. The value is used here as a fixed input.

## Setup

```bash
git clone <REPOSITORY_URL> repro && cd repro
python3 fetch_upstream.py        # adds RydbergComms @ 4095724 (git archive); never overwrites a file of this repository
python3 -m pip install -r requirements.txt
```

`fetch_upstream.py` needs `git` and network access. It refuses to continue if the commit cannot be resolved.

## Reproduce all figures and tables (no simulation, no GPU)

```bash
cd manuscript/figures/src && python3 make_all.py     # rebuilds Figs. 1-10 and S1-S10, then audits them
cd ../.. && python3 build_numbers_ledger.py          # rebuilds numbers_ledger.csv
```

**Expected output:**
- `make_all.py` ends with `ALL PASS`: 242 checks of plotted values against the ledger and the stored results, and of clipping, fonts and wording.
- The rebuilt `numbers_ledger.csv` is identical to the archived one.
- The rebuilt figures match the reference PDFs. Fig. 8 is pixel-identical but may differ in internal PDF object naming.
- `figure_audit_results.json` lists every check.

**Tables in the article.** Each table's values are rows of `numbers_ledger.csv`, which records each value's source:

| Table | Source files |
|---|---|
| 2 (AIR loss) | `results/final_validation/regen_stage05/02_losses.csv` |
| 3 (surrogate excess) | `results/final_validation/regen_stage05/03_excess_over_controls.csv` |
| 4 (constant offsets) | `results/stage06_dwell_physics/phase_mechanism/03_phase_mechanism_results.csv` |
| 5 (evidence summary) | the result files of the tests it lists (ledger rows `N_nd*`, `RV_p1_*`, `Q_*`, `P_QS_*`, `F_X_OFF_*`, `RV_P5_*`, `V_*`, `R2_*`) and the decision-rule records in `decision_records/` |
| S2 (thresholds) | `results/final_validation/regen_stage05/04_thresholds.csv` |
| S3 (convergence) | `results/stage06_dwell_physics/02_nd_convergence.csv` |
| S4 (shuffle) | `results/stage06_dwell_physics/06_temporal_shuffle_results.csv` |
| S5 (time scales) | `results/stage06_dwell_physics/03_tau_atom_results.csv`, `manuscript/posthoc/posthoc_timescales.csv` |
| S6 (second operating point) | `results/final_validation/replication/03_replication_results.csv`, `04_replication_verdicts.csv` |
| S7 (metric robustness) | `results/revision/p1/04_p1_decisive_table.csv` |
| S8 (six sequences) | `results/revision/11_p2_summary.csv`, `12_p2_contrasts.csv` |
| S9 (configuration checks) | `results/revision2/11_p1_if_summary.csv`, `12_p1_dwell_per_realization.csv`, `31_p3_verdict.json` |
| S10 (amplitude vs phase) | `results/revision2/21_p2_summary.csv`, `22_p2_verdict.json` |
| 6 (evidence summary, continued) | the campaign verdicts `results/tqe_viability_final/0[1-5]_*/*verdict.json` (ledger rows `TV_*`) |
| S11 (fair matching) | `results/tqe_viability_final/01_matching/matching_results.csv` |
| S12 (IF scan) | `results/tqe_viability_final/02_resonance/resonance_verdict.json` |
| S13 (identified memory model) | `results/tqe_viability_final/03_memory_model/gmp_verdict.json` |
| S14 (noise sensitivity) | `results/tqe_viability_final/04_noise/noise_verdict.json`, `noise_sensitivity.csv` |
| S15 (numerical and regime checks) | `results/tqe_viability_final/05_numerical/*.csv` |

Tables 1 and S1 describe the waveform construction and contain no computed results.

**Figures.** Each figure script names its inputs at the top. `make_all.py` needs only the files in this repository, plus two upstream modules:
- `python/utils/sim_SingleCarrier.py`, for the Fig. S2 spectra;
- `python/utils/transient_quantum.py`, for the model parameters quoted in the ledger.

## Re-running the simulations (GPU required)

**Build first.** Build the upstream CUDA/C++ extensions as described in the upstream README (`make`). The runs used an NVIDIA GPU. The drivers below write to `results/`, taking configurations and seeds from the files stored there.

| Step | Driver (run from its own directory) | Writes to |
|---|---|---|
| Fair-waveform screening study (N<sub>d</sub> = 1501) | `python/experiments/p1db_waveform/stage05.py`, `analyze_stage05.py` | `results/p1db_waveform_stage05/` |
| Test signals, CW calibration, surrogate controls, time scales, convergence | `stage06_dwell/make_config.py`, `run_step1b_2.sh`, `build_controls_mh.py`, `nd_convergence.py`, `tau_sensitivity.py`, `validate_controls_cw.py` | `results/stage06_dwell_physics/` |
| Dwell and shuffle campaign | `stage06_dwell/run_campaign.sh`, then `analyze06.py` | `results/stage06_dwell_physics/` |
| Follow-ups (plateau, carrier ablation, long dwell, single slow state) | `stage06_dwell/run_followup.sh`, then `analyze_followup.py` | `results/stage06_dwell_physics/` |
| Phase-modulation experiments | `stage06_dwell/phase_mechanism/pm_run.py`, then `pm_analyze.py` | `results/stage06_dwell_physics/phase_mechanism/` |
| Converged regeneration and second-operating-point replication | `final_validation/run_final_validation.sh`, then `analyze_regen.py`, `fv_analyze_replication.py` | `results/final_validation/` |
| Post-hoc analyses | `manuscript/posthoc/*.py` | `manuscript/posthoc/` |

**Before re-running:** re-running overwrites the archived files. Work on a copy, then compare it with the archive using `SHA256SUMS` or the stored CSVs.

**Seeds** are stored with the runs:
- `results/final_validation/regen_stage05/jobs/*.json` (fair-waveform seeds);
- `results/stage06_dwell_physics/stage06_config.json` and `04_controlled_pair_configs.csv` (test-signal seeds);
- `results/stage06_dwell_physics/phase_mechanism/02_phase_mechanism_configs.csv`;
- `results/final_validation/replication/02_run_configs.json`.

**Reruns were not repeated for this package.** Small differences from the archive are possible, from GPU floating-point nondeterminism and single-precision arithmetic. The article's numerical checks (N<sub>d</sub> = 8001 and dt = 0.5 ns) show how sensitive the reported effects are to numerics.

## Pre-submission revision (`results/revision/`)

A focused revision addressing reviewer-style objections followed decision rules fixed in `results/revision/00_revision_preregistration.json`, which was committed before any revision run.

| Item | Script (`python/experiments/revision/`) | Hardware | Output |
|---|---|---|---|
| P1: metric robustness of the Section IV contrasts (D and D_ref, FIR spans 2–16 µs) | `rev_p1_rerun.py`, then `rev_p1_analyze.py` | GPU, about 1 h | `results/revision/p1/` (`01_`–`05_*`) |
| P2–P6: six-sequence replication, amplitude sweep, frequency-offset grid, third LO operating point, 200 µs segment | `rev_gpu_plan.py`, then `rev_p2_p6_analyze.py` | GPU, about 45 min | `results/revision/1*_`–`5*_*`, run files in `p2/`, `p3/`, `p4/`, `p45/`, `p6/` |
| P6: double-precision closure runs | `rev_p6_cpu.py` | CPU (upstream C++ RK4 path) | `results/revision/p6/` |

**Building the double-precision path.** The CPU path needs upstream's C++ integrator compiled, using the same commands as upstream's `python/utils/build_rk4.sh`:

```bash
cd python/utils
g++ -O3 -std=c++17 -march=native -o bin/rk4_single_velocity rk4_single_velocity.cpp
g++ -O3 -std=c++17 -fPIC -shared -o bin/rk4_single_velocity.so rk4_single_velocity.cpp
```

**P1 traces are not shipped.** The 128 P1 baseband traces (about 320 MB) are regenerated bit-identically by `rev_p1_rerun.py`, and their per-job metrics are included.

## Final receiver-physics checks (`results/revision2/`)

A last pass tested whether the effects depend on the IF (5, 10 and 15 MHz), how much of the phase effect is envelope magnitude versus phase misalignment, whether the effects persist with a fourfold weaker probe, and what the internal density-matrix coherences do after a phase step. Its decision rules (`00_decision_rules.json`) and clarifications (`01_rule_clarifications.json`) were committed before any of its runs; `02_deviations.json` records the deviations made during execution, each with the state of the results when it was written.

| Item | Script (`python/experiments/revision2/`) | Hardware | Output |
|---|---|---|---|
| Runs: IF dependence (P1), weaker probe (P3), internal traces (P4) | `rev2_run.py p1 p4 p3` | GPU, about 45 min | `results/revision2/tables/`, `p1/`, `p3/`, `p4/` |
| P1 analysis (including the 5 MHz bit-identity validation and zone leakage) | `rev2_p1_analyze.py` | CPU | `results/revision2/10_`–`15_*` |
| P2 amplitude/phase decomposition of all saved phase runs | `rev2_p2.py` | CPU | `results/revision2/20_`–`22_*` |
| P3 and P4 analyses | `rev2_p3_analyze.py`, `rev2_p4_analyze.py` | CPU | `results/revision2/30_`–`31_*`, `40_`–`41_*`, `p4/*_traces.npz` |

The analyses regenerate every table from the shipped run files; `verify_reproduction.py` runs them.

## Final adversarial-validation campaign (`results/tqe_viability_final/`)

The campaign followed `results/tqe_viability_final/00_plan.md`, committed before any of its runs; deviations (two implementation fixes) are in `00_deviations.md`, and every run has a row in `00_run_manifest.csv` (configuration, seed, command, code version, runtime, output SHA-256). It reused data produced by the preceding TQE-viability pass with the same code (`results/tqe_viability/`, rules in `00_decision_rules.json`).

| Item | Script (`python/experiments/tqe_viability_final/`) | Hardware | Output |
|---|---|---|---|
| Runs: fair matching (Part I), IF-scan completion and refinement (Part II), trace/positivity prefix runs (Part V) | `tvf_run.py I II IIr V` | GPU, about 30 min | `01_matching/runs/`, `02_resonance/runs/`, `05_numerical/runs/` |
| Prior-pass runs reused: 3-sequence matched states, IF scan at 29 IFs, fair-waveform noiseless basebands | `../tqe_viability/tv_run.py p1 p2 p4` | GPU, about 1 h | `results/tqe_viability/p1/`, `p2/`, `p4/` |
| Part I analysis | `tvf_p1_analyze.py` | CPU | `01_matching/` |
| Part II analysis | `tvf_p2_analyze.py` | CPU | `02_resonance/` |
| Part III (GMP fit and held-out test) | `tvf_p3_gmp.py` | CPU, about 6 min | `03_memory_model/` |
| Part IV analysis | `tvf_p4_analyze.py` | CPU, about 3 min | `04_noise/` |
| Part V analysis | `tvf_p5_analyze.py` | CPU | `05_numerical/` |

The raw run files of the campaign and of the reused prior-pass runs (about 190 MB of NPZ files) are shipped although they are not in the project's version control; `verify_reproduction.py` re-derives every campaign table from them, except the identified-memory-model fit, which also needs the revision P1 traces (not shipped; see "Pre-submission revision") and is otherwise reported as skipped.

## Decision-rule records

The decision rules were written to time-stamped project records before the corresponding runs: five for the original study, one for the pre-submission revision (`results/revision/00_revision_preregistration.json`), three for the final checks (`results/revision2/00_`–`02_*.json`), two for the TQE-viability pass (`results/tqe_viability/00_`–`01_*.json`) and two for the final campaign (`results/tqe_viability_final/00_plan.md`, `00_deviations.md`). They are internal records: they were not deposited with an external registry before the runs. This archive makes them available so that the rules and their timing can be audited after the fact. For the files, their timestamps and what the record can and cannot show, see `decision_records/PROVENANCE.md`. The earlier files keep their historical names (`*preregistration*.json`); the article calls the tests they govern *prespecified*. Internal labels used in the files and scripts (stage numbers, case codes and test codes such as `H_gen`) are kept as archived. The article uses descriptive names for them.

## Citation

Please cite this archive (https://doi.org/<DOI>), the article, and the simulator: J. Zhu and L. Dai (2026), and the RydbergComms repository. See `CITATION.cff`.
