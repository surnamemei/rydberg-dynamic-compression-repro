"""GPU run plan for P2-P6 (preregistered in results/revision/00_revision_preregistration.json).

Order: P2 sequences (OP2, OP1) -> P3 amplitude sweep -> OP3 controls and runs -> P4/P5 offset grid -> P6 extended segment.
Every run is skipped if its output exists, so the plan is resumable.
"""
from __future__ import annotations

import json
import time

import numpy as np

from rev_common import REV
import rev_phase as RP


def build_op3_controls():
    out = REV / "artifacts" / "controls_LO0.425_Nd4001.npz"
    if out.exists():
        return
    from fv_controls import cw_table, lti_kernel
    t0 = time.perf_counter()
    d = {**cw_table(4001, 0.425), **lti_kernel(4001, 0.425), "Nd": 4001, "dt_s": 1e-9, "LO_Vpm": 0.425}
    out.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(out, **d)
    print("OP3 controls: E1dB''", d["e1db"], "settle", d["cw_settle_rel"], "lin_check", d["lin_check"], f"{time.perf_counter() - t0:.0f} s", flush=True)


def plan():
    H = RP.LEVEL_H
    q = lambda op, lvl, seed, case, **kw: {"kind": "qpsk", "op": op, "level": lvl, "seed": seed, "Ts": RP.CASES[case][0], "Tp": RP.CASES[case][1], **kw}  # noqa: E731
    p = []
    for op in ("OP2", "OP1"):                                               # P2
        p.append((f"{op}_CW_H", {"kind": "cw", "op": op, "level": H}, "p2"))
        for seed in RP.SEEDS:
            for case in RP.CASES:
                p.append((f"{op}_{case}_s{seed}_H", q(op, H, seed, case), "p2"))
    for op in ("OP1", "OP2"):                                               # P3
        for db in RP.LEVELS_DB:
            lvl = 10 ** (db / 20)
            p.append((f"{op}_CW_p{db:+d}", {"kind": "cw", "op": op, "level": lvl}, "p3"))
            for seed in RP.SEEDS[:3]:
                p.append((f"{op}_REF_s{seed}_p{db:+d}", q(op, lvl, seed, "REF"), "p3"))
    p.append(("OP3_CONTROLS", None, None))                                  # P4 third operating point
    p.append(("OP3_CW_H", {"kind": "cw", "op": "OP3", "level": H}, "p4"))
    for seed in RP.SEEDS[:3]:
        p.append((f"OP3_REF_s{seed}_H", q("OP3", H, seed, "REF"), "p4"))
    for mhz in (0.125, RP.OFFSETS_MHZ[-5], -RP.OFFSETS_MHZ[-5], RP.OFFSETS_MHZ[-2], -RP.OFFSETS_MHZ[-2]):
        p.append((f"OP3_OFF{mhz:+.3f}_H", {"kind": "offset", "op": "OP3", "level": H, "delta_Hz": mhz * 1e6}, "p4"))
    levels = [("H", H)] + [(f"p{db:+d}", 10 ** (db / 20)) for db in RP.LEVELS_DB]
    for op in ("OP1", "OP2"):                                               # P4 sweep (H slice) and P5 grid
        for tag, lvl in levels:
            for mhz in RP.OFFSETS_MHZ:
                p.append((f"{op}_OFF{mhz:+.3f}_{tag}", {"kind": "offset", "op": op, "level": lvl, "delta_Hz": mhz * 1e6}, "p45"))
    p.append(("OP1_REF_s20260701_H_tmod200", q("OP1", H, 20260701, "REF", t_mod=200e-6), "p6"))   # P6 (iii)
    return p


def main():
    pl = plan()
    (REV / "gpu_plan.json").write_text(json.dumps([[n, s, d] for n, s, d in pl], indent=0, default=float))
    t_start, done = time.perf_counter(), 0
    for i, (name, spec, sub) in enumerate(pl, 1):
        t0 = time.perf_counter()
        if name == "OP3_CONTROLS":
            build_op3_controls()
        else:
            RP.run(name, spec, sub)
        done += 1
        el = time.perf_counter() - t_start
        line = f"[{i}/{len(pl)}] {name} {time.perf_counter() - t0:.1f} s | elapsed {el / 60:.1f} min | ETA {el / done * (len(pl) - i) / 60:.1f} min"
        print(line, flush=True)
        (REV / "gpu_plan_progress.txt").write_text(line + "\n")


if __name__ == "__main__":
    main()
