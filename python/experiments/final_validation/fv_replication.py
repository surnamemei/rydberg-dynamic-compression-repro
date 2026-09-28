"""Second-operating-point replication (preregistered: results/final_validation/replication/00_preregistration.json).

Constant-amplitude phase/frequency cases at LO' = 0.35 V/m (a' = 2.4468 E1dB'), zero-drift rate controls
at the original point, and two Nd=8001 closure checks.  Nothing beyond the preregistered run matrix.
"""
from __future__ import annotations

import json
import time

import numpy as np

from fv_common import DT, FV, IF, OUT06, run_field
import models06 as M
import pm_run as P

LEVEL = 3 * 0.5774 * 10 ** (3 / 20)  # 2.4468: a_H / E1dB at the original point
RUNS = FV / "replication" / "runs"
T0, T_END = P.T0, P.T_END
SEED_A, SEED_B = P.SEED_A, P.SEED_B


def zd_transitions(case):
    d = P.transitions(case).copy()
    sign = 1.0
    for i in np.flatnonzero(np.isclose(np.abs(d), np.pi)):
        d[i] = sign * np.pi
        sign = -sign
    return d


def phase(case, t):
    if not case.get("zd"):
        return P.phase(case, t)
    tt = t - T0
    phi = np.zeros_like(t)
    Ts, Tp = case["Ts"], case["Tp"]
    for k, dk in enumerate(zd_transitions(case), start=1):
        if dk == 0:
            continue
        tb = k * Ts
        if Tp == 0:
            phi[tt >= tb] += dk
        else:
            u = np.clip((tt - (tb - Tp / 2)) / Tp, 0, 1)
            phi += dk * (0.5 - 0.5 * np.cos(np.pi * u))
    return phi


def cases(e1_prime):
    a2 = LEVEL * e1_prime
    q = lambda seed, Ts, Tp: {"kind": "qpsk", "seed": seed, "Ts": Ts, "Tp": Tp, "zd": True}  # noqa: E731
    c = {
        "LO2_CW": {"kind": "cw"},
        "LO2_OFF_+0.125": {"kind": "offset", "delta_Hz": 0.125e6},
        "LO2_OFF_+1.31": {"kind": "offset", "delta_Hz": P.F_PI2},
        "LO2_OFF_-1.31": {"kind": "offset", "delta_Hz": -P.F_PI2},
        "LO2_QPSK_REF_zd": q(SEED_A, 1e-6, 300e-9),
        "LO2_QPSK_REF_zd_s2": q(SEED_B, 1e-6, 300e-9),
        "LO2_QPSK_JUMP_zd": q(SEED_A, 1e-6, 0.0),
        "LO2_QPSK_LONGRAMP_zd": q(SEED_A, 1e-6, 900e-9),
        "LO2_QPSK_SLOW_zd": q(SEED_A, 4e-6, 300e-9),
        "LO2_QPSK_FAST_zd": q(SEED_A, 0.5e-6, 300e-9),
        "LO2_POWERSTEP": {"kind": "cw", "powerstep": True},
        "LO1_QPSK_SLOW_zd": q(SEED_A, 4e-6, 300e-9),
        "LO1_QPSK_REF_zd": q(SEED_A, 1e-6, 300e-9),
        "LO1_QPSK_FAST_zd": q(SEED_A, 0.5e-6, 300e-9),
        "LO1_OFF_+0.125": {"kind": "offset", "delta_Hz": 0.125e6},
    }
    for name, v in c.items():
        v["lo"] = 0.35 if name.startswith("LO2") else 0.5
        v["amp"] = a2 if name.startswith("LO2") else P.A_H
    return c


def stats(case, dt=DT):
    t = np.arange(int(round(T_END / dt))) * dt
    phi = phase(case, t)
    m = t >= T0
    f = (np.diff(phi) / (2 * np.pi * dt))[m[1:]]
    ntr = int(np.count_nonzero((zd_transitions(case) if case.get("zd") else P.transitions(case)))) if case["kind"] == "qpsk" else 0
    return {"n_transitions": ntr, "transition_rate_per_us": ntr / 60.0, "mean_abs_freq_offset_Hz": float(np.mean(np.abs(f))),
            "mean_signed_freq_offset_Hz": float(np.mean(f)), "peak_abs_freq_offset_Hz": float(np.max(np.abs(f)))}


def run_case(name, case, nd, ctrl):
    n = int(round(T_END / DT))
    t = np.arange(n) * DT
    phi = phase(case, t)
    amp = np.full(n, case["amp"])
    if case.get("powerstep"):
        amp[t < T0] = case["amp"] / 3
    env = amp * np.exp(1j * phi)
    y, probe_db, rhovec = run_field(case["lo"] + env * np.exp(2j * np.pi * IF * t), nd, case["lo"])
    ylin = ctrl.linear(env)
    z = {m: P.demod(v, phi, t) for m, v in (("full", y), ("linear", ylin))}
    vx, pv, pops = P.populations(rhovec, nd)
    dec = 20
    np.savez_compressed(RUNS / f"{name}_Nd{nd}.npz", t=t[::dec], phi=phi[::dec], amp=amp[::dec], probe_db=probe_db[::dec],
                        z_full=z["full"][::dec], z_linear=z["linear"][::dec], vx=vx, pv=pv, pops=pops, Nd=nd, lo=case["lo"], a=case["amp"])


def long_dwell_check():
    """xi=32 long-dwell member (r0, CW carrier, +3 dB) at Nd=8001, written under final_validation only."""
    import run06 as R
    R.N, R.XIS, R.CW_CARRIER = 2_400_000, [1.0, 4.0, 8.0, 16.0, 32.0], True
    R.controls()          # load control tables from the Stage-0.6 artifacts before redirecting outputs
    R.OUT = FV / "replication"
    rows = R.run_job("dwell", 0, "xi_32", 3, 8001, 1e-9, "fv_nd8001_long", {}, False)
    (FV / "replication" / "long_dwell_xi32_Nd8001_rows.json").write_text(json.dumps(rows, default=float))


def main():
    RUNS.mkdir(parents=True, exist_ok=True)
    ctl2 = np.load(FV / "artifacts" / "controls_LO0.35_Nd4001.npz")
    e1p = float(ctl2["e1db"])
    ctrl2 = M.Controls.from_npz(FV / "artifacts" / "controls_LO0.35_Nd4001.npz")
    ctrl1 = M.Controls.from_npz(OUT06 / "artifacts" / "controls_Nd4001.npz")
    cs = cases(e1p)
    cfg = {k: {**{kk: vv for kk, vv in v.items()}, **stats(v)} for k, v in cs.items()}
    (FV / "replication" / "02_run_configs.json").write_text(json.dumps({"E1dB_prime_Vpm": e1p, "a_prime_Vpm": LEVEL * e1p, "cases": cfg}, indent=1, default=float))
    plan = [(k, 4001) for k in cs if k != "LO1_OFF_+0.125"] + [("LO1_OFF_+0.125", 8001), ("LONG_XI32", 8001)]
    t_start, durs = time.perf_counter(), []
    for i, (name, nd) in enumerate(plan, 1):
        t0 = time.perf_counter()
        if name == "LONG_XI32":
            if not (FV / "replication" / "long_dwell_xi32_Nd8001_rows.json").exists():
                long_dwell_check()
        elif not (RUNS / f"{name}_Nd{nd}.npz").exists():
            run_case(name, cs[name], nd, ctrl2 if name.startswith("LO2") else ctrl1)
        durs.append(time.perf_counter() - t0)
        el = time.perf_counter() - t_start
        line = f"[{i}/{len(plan)}] {100 * i / len(plan):5.1f}% | {name} | Nd={nd} | elapsed {el:.0f} s | avg {np.mean(durs):.1f} s/run | ETA {np.mean(durs) * (len(plan) - i):.0f} s"
        print(line, flush=True)
        (FV / "replication" / "progress.txt").write_text(line + "\n")


if __name__ == "__main__":
    main()
