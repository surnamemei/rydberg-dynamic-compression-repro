"""P1: re-run the decisive Section IV jobs with stored baseband traces (archived seeds and code path).

Uses stage06_dwell/run06.run_job unchanged; only the output root is redirected to results/revision/p1.
"""
from __future__ import annotations

import time

import pandas as pd

from rev_common import REV
import run06 as R

MATRIX = [("dwell", 3, ["xi_0.05", "xi_0.1", "xi_0.25", "xi_0.5", "xi_1", "xi_2", "xi_4"]),
          ("dwell", 0, ["xi_0.05", "xi_0.1", "xi_4"]),
          ("dwell", 6, ["xi_0.05", "xi_4"]),
          ("shuffle", 3, ["ORIGINAL", "BLOCK_SHUFFLED_CYCLES"]),
          ("shuffle", 6, ["ORIGINAL", "BLOCK_SHUFFLED_CYCLES"])]


def main():
    R.controls()                    # control tables from the archived Stage-0.6 artifacts (read-only)
    R.OUT = REV / "p1"              # jobs/ and traces/ are written here
    out_csv = REV / "p1" / "rows_p1.csv"
    plan = [(d, r, v, p) for r in range(8) for d, p, vs in MATRIX for v in vs]
    cache, t_start = {}, time.perf_counter()
    for i, (d, r, v, p) in enumerate(plan, 1):
        rows = R.run_job(d, r, v, p, 4001, 1e-9, "p1", cache, True)
        pd.DataFrame(rows).to_csv(out_csv, mode="a", header=not out_csv.exists(), index=False)
        el = time.perf_counter() - t_start
        line = f"[{i}/{len(plan)}] {d} r{r} {v} {p:+d} dB | elapsed {el / 60:.1f} min | ETA {el / i * (len(plan) - i) / 60:.1f} min"
        print(line, flush=True)
        (REV / "p1" / "progress.txt").write_text(line + "\n")


if __name__ == "__main__":
    main()
