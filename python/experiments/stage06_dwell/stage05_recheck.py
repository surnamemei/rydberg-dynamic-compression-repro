"""Re-check the Stage-0.5 +3 dB headline with converged quadrature and all-zone controls.

Same waveforms, seeds, noise, receiver and absolute field as Stage-0.5 (+3 dB and the
-6 dB reference in the Stage-0.5 coordinate, i.e. relative to E1dB(Nd=1501)).
Full model at Nd=4001; controls from the Nd=4001 tables (F1 = Stage-0.5 fundamental-only
convention; MH = all-zone).
"""
from __future__ import annotations

import json

import numpy as np
import pandas as pd

from common import OUT, STAGE05, s
import models06 as M

MODS = ("CE", "QPSK", "16QAM", "OFDM")
SEEDS = [20262001 + 100 * r for r in range(8)]
ND = 4001


def main():
    f1 = M.Controls.from_npz(OUT / "artifacts" / "controls_Nd4001.npz")
    mh = M.ControlsMH(OUT / "artifacts" / "controls_MH_Nd4001.npz", f1)
    rows = []
    path = OUT / "09_stage05_recheck.csv"
    for mod in MODS:
        for base in SEEDS:
            seed = base + {"CE": 1, "QPSK": 2, "16QAM": 3, "OFDM": 4}[mod]
            tx, labels, alphabet = s.waveform(mod, seed)
            noise = np.random.default_rng(seed + 9991).normal(0, s.SIGMA, len(tx))
            for db in (-6, 3):
                saved = json.loads((STAGE05 / "jobs" / f"{s.job_id(mod, seed, db, 1501)}.json").read_text())
                old_rows = {r["model"]: r for r in saved["rows"]}
                env = s.E1 * 10 ** (db / 20) * tx
                outs = {"atomic": M.atomic(env, ND), "static": f1.static(env), "static_lti": f1.static_lti(env),
                        "static_mh": mh.static_mh(env), "static_lti_mh": mh.static_lti_mh(env)}
                for m, y in outs.items():
                    raw = s.extract(y + noise, mod)
                    air, ser, ber, evm, var, _ = s.split_eval(raw, labels, alphabet, mod, old_rows["atomic"]["ce_timing_offset"])
                    rows.append({"modulation": mod, "seed": seed, "Pavg_stage05_coord_dB": db,
                                 "Pavg_over_P1dB_converged_dB": db + 20 * np.log10(s.E1 / f1.e1db), "model": m, "Nd_full": ND,
                                 "AIR": air, "EVM": evm, "AIR_stage05_archived": old_rows.get(m, {}).get("AIR_native", np.nan)})
                print(mod, seed, db, {m: round(r["AIR"], 4) for m, r in zip(outs, rows[-len(outs):])}, flush=True)
            pd.DataFrame(rows).to_csv(path, index=False)


if __name__ == "__main__":
    main()
