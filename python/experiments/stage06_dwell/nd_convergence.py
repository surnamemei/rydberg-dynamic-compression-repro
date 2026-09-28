"""Step 1: thermal-quadrature (Nd) convergence on representative Stage-0.5 cases.

Tolerances are fixed in results/stage06_dwell_physics/00_preregistered_criteria.json.
"""
from __future__ import annotations

import time

import numpy as np
import pandas as pd

from common import OUT, s
import controlled_dwell as k  # archived Stage-0.5 Stage-K pair
from utils.transient_quantum import get_normal_quadrature

import sys

NDS = tuple(int(x) for x in sys.argv[1].split(",")) if len(sys.argv) > 1 else (3001, 2001, 1501, 1001, 751)  # highest first: it is the reference
CASES = [
    ("A_lowpower_QPSK_m10", "QPSK", 20262003, -10),
    ("B_P1dB_OFDM_0", "OFDM", 20262005, 0),
    ("C_plus3_16QAM", "16QAM", 20262004, 3),
    ("D_strongest_QPSK_p6", "QPSK", 20262003, 6),
]
REF_DB = -6


def baseband(intensity):
    t = np.arange(len(intensity)) / s.FS
    return s.old.sosfilt(s.old.sos, (intensity - s.old.dc) * np.exp(-2j * np.pi * s.old.IF * t))[::s.UPSAMPLE]


def nmse(a, ref):
    return float(np.sum(np.abs(a - ref) ** 2) / np.sum(np.abs(ref) ** 2))


def comm_case(mod, seed, db, nd, cache):
    key = (mod, seed, db, nd)
    if key in cache:
        return cache[key]
    tx, labels, alphabet = s.waveform(mod, seed)
    env = s.E1 * 10 ** (db / 20) * tx
    noise = np.random.default_rng(seed + 9991).normal(0, s.SIGMA, len(tx))
    t0 = time.perf_counter()
    y = s.atomic_intensity(env, nd=nd)
    raw = s.extract(y + noise, mod)
    air, ser, ber, evm, var, _ = s.split_eval(raw, labels, alphabet, mod, None)
    ac = y - np.mean(y)
    cache[key] = {"AIR": air, "SER": ser, "EVM": evm, "bb": baseband(y), "probe_ac": ac[::10],
                  "out_ac_rms": float(np.sqrt(np.mean(ac ** 2))), "runtime_s": time.perf_counter() - t0}
    return cache[key]


def main():
    rows, cache, traces = [], {}, {}
    for name, mod, seed, db in CASES:
        for nd in NDS:
            r = comm_case(mod, seed, db, nd, cache)
            ref = comm_case(mod, seed, REF_DB, nd, cache)
            top = cache[(mod, seed, db, NDS[0])]
            rel = 100 * (ref["AIR"] - r["AIR"]) / ref["AIR"]
            rows.append({"case": name, "modulation": mod, "seed": seed, "Pavg_over_P1dB_dB": db, "Nd": nd,
                         "Nd_retained_within_3sigma": int(np.count_nonzero(np.abs(get_normal_quadrature(nd)[0]) <= 3)),
                         "AIR": r["AIR"], "AIR_ref_m6dB": ref["AIR"], "AIR_loss_vs_m6dB": ref["AIR"] - r["AIR"],
                         "relative_AIR_loss_pct": rel, "EVM": r["EVM"], "SER": r["SER"], "output_AC_RMS": r["out_ac_rms"],
                         "probe_AC_NMSE_vs_top": nmse(r["probe_ac"], top["probe_ac"]),
                         "baseband_NMSE_vs_top": nmse(r["bb"], top["bb"]), "runtime_s": r["runtime_s"]})
            traces[f"{name}_Nd{nd}"] = r["bb"]
            print({k2: v for k2, v in rows[-1].items()}, flush=True)
    df = pd.DataFrame(rows)
    for name, g in df.groupby("case"):
        top = g[g.Nd == NDS[0]].iloc[0]
        df.loc[g.index, "AIR_diff_vs_top"] = g.AIR - top.AIR
        df.loc[g.index, "rel_loss_diff_vs_top_pp"] = g.relative_AIR_loss_pct - top.relative_AIR_loss_pct
        df.loc[g.index, "out_rms_rel_diff_vs_top"] = g.output_AC_RMS / top.output_AC_RMS - 1

    # Case E: the archived Stage-0.5 controlled dwell pair (physics-level metric only).
    a, b = k.make_pair()
    ktop = {}
    for nd in NDS:
        vals = {}
        for wname, env in (("short", a), ("long", b)):
            y = s.atomic_intensity(env, nd=nd)
            bb = baseband(y)
            high = np.abs(env)[::s.UPSAMPLE] > s.E1
            vals[wname] = (bb, float((np.mean(np.abs(bb[high])) - np.mean(np.abs(bb[~high]))) / s.E1))
        if nd == NDS[0]:
            ktop = vals
        dc = vals["long"][1] - vals["short"][1]
        dct = ktop["long"][1] - ktop["short"][1]
        rows_e = {"case": "E_stage05_dwell_pair", "modulation": "two-level", "seed": 20260926, "Pavg_over_P1dB_dB": -0.90, "Nd": nd,
                  "contrast_short": vals["short"][1], "contrast_long": vals["long"][1], "pair_contrast_difference": dc,
                  "pair_contrast_difference_diff_vs_top": dc - dct,
                  "baseband_NMSE_vs_top": max(nmse(vals[w][0], ktop[w][0]) for w in ("short", "long"))}
        df = pd.concat([df, pd.DataFrame([rows_e])], ignore_index=True)
        print(rows_e, flush=True)
    df.to_csv(OUT / "02_nd_convergence.csv", index=False)
    np.savez_compressed(OUT / "nd_convergence_baseband_traces.npz", **traces)


if __name__ == "__main__":
    main()
