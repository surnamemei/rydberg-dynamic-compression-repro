# Reproduction Instructions (v1.0.0)

## 1. Set up (no GPU needed)

```bash
git clone <REPOSITORY_URL> rydberg-dynamic-compression-repro
cd rydberg-dynamic-compression-repro
python3 fetch_upstream.py            # adds RydbergComms @ 4095724 next to these files; never overwrites a file here
python3 -m pip install -r requirements.txt
```

Tested with Python 3.12.3 and the package versions pinned in `requirements.txt`. `verify_reproduction.py` also uses poppler's `pdftoppm` and `pdftotext`, if installed, to compare PDFs pixel by pixel.

## 2. Verify everything in one command

```bash
python3 verify_reproduction.py
```

**What it does:**
1. checks every shipped file against `SHA256SUMS`;
2. rebuilds from the shipped processed results, with no simulation:
   - all manuscript figures (Figs. 1–10) and supplementary figures (Figs. S1–S10), with the figure audit;
   - `manuscript/numbers_ledger.csv`;
   - the revision tables (`results/revision/`);
   - the final-check tables (`results/revision2/`), from the shipped run files;
   - the final-campaign tables and verdicts (`results/tqe_viability_final/01_`–`05_*`), from the shipped run files (the GMP fit, about 6 min on a CPU, only when the revision P1 traces are present; see "Expected differences");
3. compares each regenerated file with the shipped one;
4. restores the shipped files.

It ends with `REPRODUCTION: PASS`.

**Where the article's tables come from:**
- **Numbers:** every number in the article's text and tables is a row of `numbers_ledger.csv`, with its source file. The "Tables" section of `README.md` maps each table to its source.
- **Main-text tables:** Tables 2–4 are reproduced through the ledger.
- **Supplementary tables:** Tables S2–S15 are reproduced through the ledger, `results/revision/`, `results/revision2/` and `results/tqe_viability_final/`.

**Expected differences.** Outputs are byte-identical, with three tolerated exceptions:
- **Revision P1 tables** (`02_`, `04_`, `06_*.csv`): these are re-derived from the shipped per-job metrics CSV. They may differ from the shipped tables in the last floating-point digits (relative differences around 1e-11). They are accepted when the numbers agree to 1e-9 relative and every text cell, including the verdicts, is identical.
- **Final-check tables** (`results/revision2/`): CSV and JSON outputs are accepted when every number agrees to 1e-9 relative and every text value is identical; the internal-trace files (`p4/*_traces.npz`) when every array agrees to 1e-9 relative.
- **Final-campaign tables** (`results/tqe_viability_final/`): CSV and JSON outputs are accepted when every number agrees to 1e-9 relative and every text value is identical, ignoring only the generation metadata (`generated`, `code_version`, `code_dirty`, `runtime_s`); the regenerated summaries and diagnostic PNGs are restored from the shipped copies without comparison. The identified-memory-model fit (`tvf_p3_gmp.py`) also needs the 128 revision P1 dwell/shuffle traces, which are not shipped; without them it is reported as SKIPPED and its outputs are not compared. To include it, regenerate the traces with `python/experiments/revision/rev_p1_rerun.py` (GPU). Its tables are compared at 1e-6 relative, because the fit's ill-conditioned normal equations make multithreaded BLAS results differ at about 1e-8 relative.
- **A PDF whose bytes differ:** it is accepted only if it is pixel- and text-identical.

## 3. Individual steps

```bash
cd manuscript/figures/src && python3 make_all.py                      # figures + audit ("ALL PASS")
cd manuscript && python3 build_numbers_ledger.py                       # numbers ledger
cd python/experiments/revision && python3 rev_p2_p6_analyze.py         # revision tables P2-P6
cd python/experiments/revision2 && python3 rev2_p2.py && python3 rev2_p1_analyze.py && python3 rev2_p3_analyze.py && python3 rev2_p4_analyze.py   # final checks
cd python/experiments/tqe_viability_final && python3 tvf_p1_analyze.py && python3 tvf_p2_analyze.py && python3 tvf_p3_gmp.py && python3 tvf_p4_analyze.py && python3 tvf_p5_analyze.py   # final campaign
```

## 4. Re-running the simulations (optional; GPU)

Re-running the simulations needs a CUDA GPU and the upstream extensions built as described in upstream's README. The drivers and their outputs are listed in `README.md` ("Re-running the simulations", "Pre-submission revision", "Final receiver-physics checks" and "Final adversarial-validation campaign"). `tvf_run.py` can wait for another GPU job to finish before starting: set `TVF_PRIOR_PID_FILE` to a file holding that job's process ID.

- **Overwrites:** the drivers overwrite the shipped results, so work on a copy.
- **Seeds:** the random seeds are stored with the configurations listed there.
- **Precision:** the double-precision checks use upstream's CPU C++ path; `README.md` gives the build commands.

## 5. Windows notes

A few screening-stage scripts register the CUDA runtime DLL directory on Windows.
- It is taken from the environment variable `RYDBERG_CUDA_DLL_DIR`, if set.
- Otherwise it is `<repository>/.cuda-env/Library/bin`.
- On Linux the runtime is found through the rpath of the upstream `cu_ryd.so`, and nothing needs to be set.
