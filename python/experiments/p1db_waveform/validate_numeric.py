from __future__ import annotations

import csv
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.signal import resample_poly

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "python" / "experiments" / "p1db_waveform"))
import stage05 as s


def one(mod: str, seed: int, method: str, nd: int, dt: float, reference_air: float, timing):
    tx, labels, alphabet = s.waveform(mod, seed)
    db = 3
    env = s.E1 * 10 ** (db / 20) * tx
    if dt != 1e-9:
        env = resample_poly(env, 2, 1)
    out = s.atomic_intensity(env, dt=dt, nd=nd)
    if dt != 1e-9:
        out = out[::2]
    noise = np.random.default_rng(seed + 9991).normal(0, s.SIGMA, len(tx))
    rx = s.extract(out + noise, mod)
    air, ser, ber, evm, var, _ = s.split_eval(rx, labels, alphabet, mod, timing)
    return {"method": method, "modulation": mod, "seed": seed, "power_db": db, "Nd": nd, "dt_ns": dt*1e9, "AIR_reference_1ns_Nd1501": reference_air, "AIR_validation": air, "AIR_abs_difference": abs(air-reference_air), "SER_validation": ser, "BER_validation": ber, "EVM_validation": evm, "cal_variance": var, "status": "complete" if np.isfinite(air) else "nonfinite"}


def main():
    outdir = ROOT / "results" / "p1db_waveform_stage05"
    results = pd.read_csv(outdir / "02_all_results.csv")
    version = {"CE": "fair_v2", "QPSK": "fair_v2", "16QAM": "fair_v2", "OFDM": "fair_v4_randomized_balanced_ofdm_training"}
    cases = [("CE", 20262002), ("QPSK", 20262003), ("16QAM", 20262004), ("OFDM", 20262005)]
    rows = []
    path = outdir / "05_numerical_validation.csv"
    if path.exists():
        path.unlink()
    for mod, seed in cases:
        refrow = results[(results.config_version == version[mod]) & (results.modulation == mod) & (results.seed == seed) & (results.model == "atomic") & (results.Pavg_over_P1dB_dB == 3)].iloc[0]
        ref = float(refrow.AIR_native)
        timing = None if pd.isna(refrow.ce_timing_offset) else int(refrow.ce_timing_offset)
        for nd, dt, method in ((1501, .5e-9, "timestep_0p5ns"), (501, 1e-9, "doppler_Nd501")):
            t0 = time.perf_counter()
            row = one(mod, seed, method, nd, dt, ref, timing)
            row["runtime_s"] = time.perf_counter()-t0
            rows.append(row)
            with path.open("a", newline="", encoding="utf-8") as f:
                w = csv.DictWriter(f, fieldnames=list(row))
                if f.tell() == 0:
                    w.writeheader()
                w.writerow(row)
                f.flush()
            print(method, mod, "AIR", row["AIR_validation"], "diff", row["AIR_abs_difference"], "s", row["runtime_s"], flush=True)


if __name__ == "__main__":
    main()
