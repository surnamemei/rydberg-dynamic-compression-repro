"""Converged regeneration of the Stage-0.5 fair-waveform communication grid
(protocol: results/final_validation/regen_stage05/00_protocol.json).

Unchanged from Stage-0.5: waveforms, seeds, noise draw, receiver chain, split_eval, CE timing rule.
Changed: Nd=4001; power axis relative to the converged E1dB; converged controls incl. all-zone.
Archived results/p1db_waveform_stage05 is not touched.
"""
from __future__ import annotations

import csv
import json
import time

import numpy as np

from fv_common import FV, OUT06, s
import models06 as M

OUT = FV / "regen_stage05"
ND = 4001
MODS = ("CE", "QPSK", "16QAM", "OFDM")
POWERS_8 = (-6, -3, 0, 3, 6)
POWERS_4 = (-15, -10)
MOD_OFFSET = {"CE": 1, "QPSK": 2, "16QAM": 3, "OFDM": 4}


def plan():
    jobs = []
    for mod in MODS:
        for r in range(8):
            seed = 20262001 + 100 * r + MOD_OFFSET[mod]
            powers = sorted((POWERS_4 if r < 4 else ()) + POWERS_8)
            for p in powers:
                jobs.append((mod, r, seed, p))
    return jobs


def main():
    (OUT / "jobs").mkdir(parents=True, exist_ok=True)
    f1 = M.Controls.from_npz(OUT06 / "artifacts" / "controls_Nd4001.npz")
    mh = M.ControlsMH(OUT06 / "artifacts" / "controls_MH_Nd4001.npz", f1)
    e1 = f1.e1db
    jobs = plan()
    timing, cache, t_start, durs = {}, {}, time.perf_counter(), []
    rows_path = OUT / "rows.csv"
    for i, (mod, r, seed, p) in enumerate(jobs, 1):
        jid = f"regen_{mod}_s{seed}_p{p:+03d}_Nd{ND}"
        jp = OUT / "jobs" / f"{jid}.json"
        t0 = time.perf_counter()
        if jp.exists():
            saved = json.loads(jp.read_text())
            timing[(mod, seed)] = saved["ce_timing_offset"]
        else:
            if (mod, seed) not in cache:
                cache.clear()
                cache[(mod, seed)] = s.waveform(mod, seed)
            tx, labels, alphabet = cache[(mod, seed)]
            env = e1 * 10 ** (p / 20) * tx
            noise = np.random.default_rng(seed + 9991).normal(0, s.SIGMA, len(tx))
            outs = {"atomic": M.atomic(env, ND), "static": f1.static(env), "static_lti": f1.static_lti(env),
                    "static_mh": mh.static_mh(env), "static_lti_mh": mh.static_lti_mh(env), "linear": f1.linear(env)}
            if not np.isfinite(outs["atomic"]).all():
                raise RuntimeError(f"nonfinite atomic output {jid}")
            tm = timing.get((mod, seed))
            if mod == "CE" and tm is None:
                *_, tm = s.split_eval(s.extract(outs["atomic"] + noise, mod), labels, alphabet, mod, None)
            timing[(mod, seed)] = tm
            mag = np.abs(tx)
            papr_db = float(10 * np.log10(np.max(mag ** 2) / np.mean(mag ** 2)))
            rows = []
            for m, y in outs.items():
                air, ser, ber, evm, var, _ = s.split_eval(s.extract(y + noise, mod), labels, alphabet, mod, tm)
                rows.append({"job_id": jid, "modulation": mod, "seed": seed, "realization": r, "model": m, "Pavg_over_P1dB_dB": p,
                             "Ppeak_over_P1dB_dB": p + papr_db, "E1dB_Vpm": e1, "E_rms_Vpm": e1 * 10 ** (p / 20), "PAPR_dB": papr_db,
                             "AIR_native": air, "SER": ser, "BER": ber, "EVM": evm, "cal_var": var, "ce_timing_offset": tm,
                             "information_symbol_rate_Hz": s.rate_hz(mod), "Nd": ND, "dt_s": 1e-9})
            tmp = jp.with_suffix(".tmp")
            tmp.write_text(json.dumps({"job_id": jid, "status": "complete", "ce_timing_offset": tm, "rows": rows}, default=float))
            tmp.replace(jp)
            with rows_path.open("a", newline="") as f:
                w = csv.DictWriter(f, fieldnames=list(rows[0]))
                if f.tell() == 0:
                    w.writeheader()
                w.writerows(rows)
        durs.append(time.perf_counter() - t0)
        el = time.perf_counter() - t_start
        line = f"[{i}/{len(jobs)}] {100 * i / len(jobs):5.1f}% | {mod} | r{r} | {p:+d} dB | elapsed {el:.0f} s | avg {np.mean(durs):.1f} s/job | ETA {np.mean(durs) * (len(jobs) - i):.0f} s"
        print(line, flush=True)
        (OUT / "progress.txt").write_text(line + "\n")


if __name__ == "__main__":
    main()
