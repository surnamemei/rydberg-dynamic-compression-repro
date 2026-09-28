# Provenance

## Simulator
All simulations use the transient density-matrix simulator of J. Zhu and L. Dai: RydbergComms, https://github.com/johnzja/RydbergComms, commit `4095724`, MIT License. The accompanying paper is *IEEE Trans. Wireless Commun.*, vol. 25, pp. 8292–8307, 2026, doi:10.1109/TWC.2025.3637029.

It was used **without modification**. The study's project history only *added* files to the upstream tree. This repository contains those added files. `fetch_upstream.py` supplies the upstream files at the pinned commit. Upstream code is not redistributed here.

## Study-specific material
Everything in this repository was written for this study, except the upstream simulator fetched by `fetch_upstream.py`. That covers the experiment drivers, analysis code, surrogate/control models, figure and ledger scripts, configurations, seeds, decision-rule records and processed results. It is released under the MIT License (see `LICENSE`).

**Processed results** are the archived outputs of the original runs. They are shipped unchanged, apart from the sanitized manifests described below.

**The pre-submission revision** (`results/revision/`, `python/experiments/revision/`) and **the final receiver-physics checks** (`results/revision2/`, `python/experiments/revision2/`) followed decision rules committed before any of their runs; see `decision_records/`.

**The final adversarial-validation campaign** (`results/tqe_viability_final/`, `python/experiments/tqe_viability_final/`) followed a plan committed before any of its runs (`00_plan.md`). It reused data produced earlier by the same code in the TQE-viability pass (`results/tqe_viability/`), declared in the plan; one of those outputs had been inspected when the plan was written, as the plan states. Its two deviations are implementation fixes (`00_deviations.md`). After the runs, the location of an optional process-ID file used only to wait for another GPU job was made configurable (`TVF_PRIOR_PID_FILE`); this does not affect any result. The raw run files of the campaign and of the reused runs are shipped; the campaign's run manifest (`00_run_manifest.csv`) lists the SHA-256 of every campaign output.

## Decision-rule records
`decision_records/` holds byte-identical copies of the thirteen decision-rule records, together with `decision_records/PROVENANCE.md`. The records are internal, time-stamped project files; they were not deposited with an external registry before the runs, and this archive makes them available for audit after the fact. The note lists, for each file:
- its original path;
- its in-file timestamp;
- when it was first committed;
- what the record can and cannot establish.

The originals remain at their paths under `results/`.

## Sanitized files
Three archived run manifests recorded the absolute path of a local Windows CUDA directory.
- **Public versions:** in the public copies, that local root is replaced by `<LOCAL_PROJECT_ROOT>`. Nothing else is changed.
- **Documentation:** `provenance/SANITIZATION.md` documents the change.
- **Archival checksums:** `provenance/archival_sha256.txt` records the SHA-256 of the unmodified archival files, alongside the checksums of the public versions.
- **Code:** the corresponding code paths were made configurable (see `REPRODUCE.md`, "Windows notes").

## Not included
- the upstream simulator;
- raw per-velocity-class traces;
- the 128 revision P1 baseband traces (about 320 MB), which `python/experiments/revision/rev_p1_rerun.py` regenerates bit-identically on a GPU; their per-job metrics are included;
- the private project history.
