# Sanitization of archived run manifests

Three archived run manifests (written by `python/experiments/p1db_waveform/stage05.py` at run time) recorded the absolute path of the local Windows CUDA runtime directory under the author's local project root. In the public copies, **only that local root** (`D:\Research\<local project folder>`) was replaced by `<LOCAL_PROJECT_ROOT>`. The rest of the path (`\.cuda-env\Library\bin`) is kept.

**Nothing else was changed:**
- no numerical value, configuration, seed or timestamp;
- no commit, result identity or key order.

The build script parses each file before and after the substitution and verifies that the single changed field is exactly this string replacement.

| File | Changed field | SHA-256, archival original | SHA-256, public version |
|---|---|---|---|
| `results/p1db_waveform_stage05/04_job_manifest_fair_v2_high.json` | `/cuda_dll_dir` | `bd2152138c91ee22…` | `b3aa378e98ca497a…` |
| `results/p1db_waveform_stage05/screening_initial/04_job_manifest.json` | `/cuda_dll_dir` | `8f7b67eb8ae44cd3…` | `d213d0a12fdef3bf…` |
| `results/p1db_waveform_stage05/screening_v3_correlated_preamble/04_job_manifest.json` | `/cuda_dll_dir` | `592aedf06f413cfa…` | `2c54e3e865307188…` |

Full checksums are in `provenance/archival_sha256.txt`. The archival originals are kept in the private project record. The code that wrote or used this path now reads it from the environment variable `RYDBERG_CUDA_DLL_DIR` (`REPRODUCE.md`, "Windows notes").
