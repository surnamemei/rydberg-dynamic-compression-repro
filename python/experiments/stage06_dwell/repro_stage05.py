"""Step 0: recompute archived Stage-0.5 +3 dB jobs on the current backend."""
from __future__ import annotations

import json
import time

import numpy as np
import pandas as pd

from common import OUT, STAGE05, s

CASES = [("CE", 20262002), ("QPSK", 20262003), ("16QAM", 20262004), ("OFDM", 20262005)]


def main():
    rows = []
    for mod, seed in CASES:
        saved = json.loads((STAGE05 / "jobs" / f"{s.job_id(mod, seed, 3, 1501)}.json").read_text())
        ref = {r["model"]: r for r in saved["rows"]}
        tx, labels, alphabet = s.waveform(mod, seed)
        env = s.E1 * 10 ** (3 / 20) * tx
        noise = np.random.default_rng(seed + 9991).normal(0, s.SIGMA, len(tx))
        t0 = time.perf_counter()
        outs = {"atomic": s.atomic_intensity(env, nd=1501), "static": s.static_intensity(env), "static_lti": s.static_lti_intensity(env)}
        for model, y in outs.items():
            raw = s.extract(y + noise, mod)
            timing = ref[model]["ce_timing_offset"]
            air, ser, ber, evm, var, _ = s.split_eval(raw, labels, alphabet, mod, timing)
            rows.append({"modulation": mod, "seed": seed, "Pavg_over_P1dB_dB": 3, "model": model,
                         "AIR_archived": ref[model]["AIR_native"], "AIR_recomputed": air,
                         "AIR_abs_diff": abs(air - ref[model]["AIR_native"]),
                         "EVM_archived": ref[model]["EVM"], "EVM_recomputed": evm,
                         "backend": "linux_cuda13_sm120", "runtime_s": time.perf_counter() - t0})
            print(rows[-1], flush=True)
    pd.DataFrame(rows).to_csv(OUT / "00_stage05_reproducibility_check.csv", index=False)


if __name__ == "__main__":
    main()
