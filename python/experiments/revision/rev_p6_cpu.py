"""P6 (i)-(ii): double-precision (upstream CPU C++ RK4, unmodified) closure runs.

(i)  OP1 CW and REF (zero-drift, seed 20260701) at a_H, 110 us, identical to the GPU P2 runs except device.
(ii) Long-dwell family, realization 0, unmodulated carrier, xi = 32, +3 dB (2.4 ms), through run06.run_job with only
     the full-model call switched to the CPU double-precision path; compared with the archived GPU job.
"""
from __future__ import annotations

import json
import time

import numpy as np

from rev_common import REV
import rev_phase as RP
import run06 as R
import models06 as M
from common import s


def atomic_cpu(env, nd, dt=1e-9):
    n = len(env)
    t = np.arange(n) * dt
    cfg = s.old.SimConfig(300, n, dt, n * dt, t)
    res = s.old.sim.run(s.old.raqr, cfg, s.old.raqr.A_LO + env * np.exp(2j * np.pi * s.old.IF * t), Nd=nd,
                        device="cpu", use_cpp_rk=True, n_jobs=30)
    return 10 ** (res.probeResponse / 10)


def main():
    log = {}
    for name, spec in (("OP1_CW_H_fp64", {"kind": "cw", "op": "OP1", "level": RP.LEVEL_H}),
                       ("OP1_REF_s20260701_H_fp64", {"kind": "qpsk", "op": "OP1", "level": RP.LEVEL_H, "seed": 20260701, "Ts": 1e-6, "Tp": 300e-9})):
        t0 = time.perf_counter()
        RP.run(name, spec, "p6", device="cpu")
        log[name] = time.perf_counter() - t0
        print(name, f"{log[name]:.1f} s", flush=True)
    R.N, R.XIS, R.CW_CARRIER = 2_400_000, [1.0, 4.0, 8.0, 16.0, 32.0], True
    R.controls()
    R.OUT = REV / "p6"
    M.atomic = atomic_cpu
    t0 = time.perf_counter()
    rows = R.run_job("dwell", 0, "xi_32", 3, 4001, 1e-9, "p6_fp64_long", {}, False)
    log["long_dwell_xi32_fp64"] = time.perf_counter() - t0
    (REV / "p6" / "long_dwell_xi32_fp64_rows.json").write_text(json.dumps(rows, default=float))
    (REV / "p6" / "cpu_wall_times_s.json").write_text(json.dumps(log, indent=1))
    print("done", log, flush=True)


if __name__ == "__main__":
    main()
