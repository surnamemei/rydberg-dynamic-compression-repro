"""Constant-amplitude phase/frequency runs for the revision (P2-P6).

Generalizes final_validation/fv_replication.run_case: identical timeline (0-50 us CW, then the case-specific
modulation), zero-drift transition rule, demodulator (phase_mechanism/pm_run.demod) and linear reference, with a
configurable modulated duration and device (the upstream CPU path is double precision; the GPU kernel is single).
"""
from __future__ import annotations

import dataclasses
import json
import time

import numpy as np

from rev_common import REV
from common import OUT as OUT06, ROOT, s
import models06 as M
import pm_run as P

FV = ROOT / "results" / "final_validation"
IF = 5e6
DT = 1e-9
T0 = 50e-6
LEVEL_H = 3 * 0.5774 * 10 ** (3 / 20)          # a_H / E1dB = 2.4468 (7.77 dB)
SEEDS = [20260701, 20260702, 20260703, 20260704, 20260705, 20260706]
OFFSETS_MHZ = [0.125, 0.25, 0.5, 0.75, 1.0, P.F_PI2 / 1e6, 1.75, 2.0, P.F_PI / 1e6, 3.0]
OFFSETS_MHZ = sorted([-v for v in OFFSETS_MHZ] + OFFSETS_MHZ)
LEVELS_DB = [-6, -3, 0, 3, 6, 8]
OPS = {"OP1": {"lo": 0.5, "controls": OUT06 / "artifacts" / "controls_Nd4001.npz"},
       "OP2": {"lo": 0.35, "controls": FV / "artifacts" / "controls_LO0.35_Nd4001.npz"},
       "OP3": {"lo": 0.425, "controls": REV / "artifacts" / "controls_LO0.425_Nd4001.npz"}}
CASES = {"REF": (1e-6, 300e-9), "FAST": (0.5e-6, 300e-9), "SLOW": (4e-6, 300e-9), "JUMP": (1e-6, 0.0), "LONGRAMP": (1e-6, 900e-9)}
_CTRL = {}


def ctrl(op):
    if op not in _CTRL:
        _CTRL[op] = M.Controls.from_npz(OPS[op]["controls"])
    return _CTRL[op]


def e1db(op):
    return float(np.load(OPS[op]["controls"])["e1db"])


def zd_transitions(seed, Ts, t_mod):
    """Zero-drift QPSK transitions: fv_replication.zd_transitions applied to pm_run.transitions, any duration."""
    nsym = int(round(t_mod / Ts))
    labels = np.random.default_rng(seed).integers(0, 4, nsym + 1)
    d = np.angle(np.exp(1j * np.pi / 2 * np.diff(labels)))
    d[np.isclose(np.abs(d), np.pi)] = np.pi
    d = d[: nsym - 1]
    sign = 1.0
    for i in np.flatnonzero(np.isclose(np.abs(d), np.pi)):
        d[i] = sign * np.pi
        sign = -sign
    return d


def phase(spec, t):
    tt = t - T0
    phi = np.zeros_like(t)
    if spec["kind"] == "cw":
        return phi
    if spec["kind"] == "offset":
        m = tt >= 0
        phi[m] = 2 * np.pi * spec["delta_Hz"] * tt[m]
        return phi
    Ts, Tp = spec["Ts"], spec["Tp"]
    for k, dk in enumerate(zd_transitions(spec["seed"], Ts, spec["t_mod"]), start=1):
        if dk == 0:
            continue
        tb = k * Ts
        if Tp == 0:
            phi[tt >= tb] += dk
        else:
            u = np.clip((tt - (tb - Tp / 2)) / Tp, 0, 1)
            phi += dk * (0.5 - 0.5 * np.cos(np.pi * u))
    return phi


def run(name, spec, subdir, device="cuda", n_jobs=30, force=False):
    """spec: kind in {cw, offset, qpsk}; op; level (a/E1dB); seed/Ts/Tp or delta_Hz; t_mod (default 60 us)."""
    out = REV / subdir / f"{name}.npz"
    if out.exists() and not force:
        return out
    out.parent.mkdir(parents=True, exist_ok=True)
    op = spec["op"]
    lo, amp = OPS[op]["lo"], spec["level"] * e1db(op)
    spec = {**spec, "t_mod": spec.get("t_mod", 60e-6)}
    n = int(round((T0 + spec["t_mod"]) / DT))
    t = np.arange(n) * DT
    phi = phase(spec, t)
    env = amp * np.exp(1j * phi)
    field = lo + env * np.exp(2j * np.pi * IF * t)
    raqr = dataclasses.replace(s.old.raqr, A_LO=lo)
    cfg = s.old.SimConfig(300, n, DT, n * DT, t)
    t0 = time.perf_counter()
    if device == "cuda":
        res = s.old.sim.run(raqr, cfg, field, Nd=spec.get("Nd", 4001), device="cuda")
    else:
        res = s.old.sim.run(raqr, cfg, field, Nd=spec.get("Nd", 4001), device="cpu", use_cpp_rk=True, n_jobs=n_jobs)
    wall = time.perf_counter() - t0
    y = 10 ** (res.probeResponse / 10)
    ylin = ctrl(op).linear(env)
    z = {m: P.demod(v, phi, t) for m, v in (("full", y), ("linear", ylin))}
    dec = 20
    np.savez_compressed(out, t=t[::dec], phi=phi[::dec], z_full=z["full"][::dec], z_linear=z["linear"][::dec],
                        probe_db=res.probeResponse[::dec], lo=lo, amp=amp, wall_s=wall, device=device,
                        spec=json.dumps({k: (float(v) if isinstance(v, (np.floating, float)) else v) for k, v in spec.items()}))
    return out


def g_final(path, t_mod=60e-6, window=(20.5e-6, 0.5e-6)):
    """Median |z_full|/|z_linear| over [T_end-20.5 us, T_end-0.5 us) (pm_analyze / fv_analyze_replication windows)."""
    z = np.load(path)
    t = z["t"]
    t_end = T0 + t_mod
    m = (t >= t_end - window[0]) & (t < t_end - window[1])
    return float(np.median(np.abs(z["z_full"][m]) / np.abs(z["z_linear"][m])))


def g_window(path, a_us, b_us):
    z = np.load(path)
    t = z["t"] * 1e6
    m = (t >= a_us) & (t < b_us)
    return float(np.median(np.abs(z["z_full"][m]) / np.abs(z["z_linear"][m])))
