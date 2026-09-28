"""Phase-mechanism dissection runner.

Pre-registered in results/stage06_dwell_physics/phase_mechanism/00_preregistration.json.
Constant |E| = a_H throughout; 0-50 us CW (common settled state), 50-110 us case-specific
phase/frequency modulation; the run ends at the plateau end so rhovec_final is the
end-of-plateau atomic state.
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))

import numpy as np
import pandas as pd

from common import OUT as OUT06, s
import models06 as M
from utils.transient_quantum import get_normal_quadrature, get_physical_constant

PM = OUT06 / "phase_mechanism"
IF = 5e6
E1C = float(np.load(OUT06 / "artifacts" / "controls_Nd4001.npz")["e1db"])
A_H = 3 * 0.5774 * 10 ** (3 / 20) * E1C
T0, T_END = 50e-6, 110e-6
T_MOD = T_END - T0
SEED_A, SEED_B = 20260701, 20260702
F_PI2 = (np.pi / 2) / (4 * 300e-9)  # peak inst. freq of a pi/2 step with a 300 ns raised-cosine ramp: dphi/(4 Tp) [Hz]
F_PI = 2 * F_PI2

CASES = {
    "CW": {"kind": "cw"},
    "OFF_+0.125": {"kind": "offset", "delta_Hz": 0.125e6},
    "OFF_+1.31": {"kind": "offset", "delta_Hz": F_PI2},
    "OFF_-1.31": {"kind": "offset", "delta_Hz": -F_PI2},
    "OFF_+2.62": {"kind": "offset", "delta_Hz": F_PI},
    "OFF_-2.62": {"kind": "offset", "delta_Hz": -F_PI},
    "QPSK_REF": {"kind": "qpsk", "seed": SEED_A, "Ts": 1e-6, "Tp": 300e-9},
    "QPSK_REF_s2": {"kind": "qpsk", "seed": SEED_B, "Ts": 1e-6, "Tp": 300e-9},
    "QPSK_SLOW": {"kind": "qpsk", "seed": SEED_A, "Ts": 4e-6, "Tp": 300e-9},
    "QPSK_FAST": {"kind": "qpsk", "seed": SEED_A, "Ts": 0.5e-6, "Tp": 300e-9},
    "QPSK_JUMP": {"kind": "qpsk", "seed": SEED_A, "Ts": 1e-6, "Tp": 0.0},
    "QPSK_LONGRAMP": {"kind": "qpsk", "seed": SEED_A, "Ts": 1e-6, "Tp": 900e-9},
    "TOGGLE": {"kind": "toggle", "Ts": 1e-6, "Tp": 300e-9},
}
CHECKS = ["CW", "QPSK_REF", "QPSK_JUMP", "OFF_+1.31"]


def transitions(case):
    nsym = int(round(T_MOD / case["Ts"]))
    if case["kind"] == "qpsk":
        labels = np.random.default_rng(case["seed"]).integers(0, 4, nsym + 1)
        d = np.angle(np.exp(1j * np.pi / 2 * np.diff(labels)))
        d[np.isclose(np.abs(d), np.pi)] = np.pi  # campaign convention: pi steps taken as +pi
    else:
        d = np.resize(np.array([0.0, np.pi / 2, -np.pi / 2, np.pi]), nsym)
    return d[: nsym - 1]  # one transition per interior symbol boundary


def phase(case, t):
    """phi(t) relative to the CW carrier (zero before T0), defined in continuous time."""
    tt = t - T0
    phi = np.zeros_like(t)
    if case["kind"] == "cw":
        return phi
    if case["kind"] == "offset":
        m = tt >= 0
        phi[m] = 2 * np.pi * case["delta_Hz"] * tt[m]
        return phi
    Ts, Tp = case["Ts"], case["Tp"]
    for k, dk in enumerate(transitions(case), start=1):
        if dk == 0:
            continue
        tb = k * Ts
        if Tp == 0:
            phi[tt >= tb] += dk
        else:
            u = np.clip((tt - (tb - Tp / 2)) / Tp, 0, 1)
            phi += dk * (0.5 - 0.5 * np.cos(np.pi * u))
    return phi


def case_stats(name, case, dt=1e-9):
    t = np.arange(int(round(T_END / dt))) * dt
    phi = phase(case, t)
    m = t >= T0
    f = np.diff(phi) / (2 * np.pi * dt)
    f = f[m[1:]]
    if case["kind"] in ("qpsk", "toggle"):
        d = transitions(case)
        n_tr, exc = int(np.count_nonzero(d)), float(np.sum(np.abs(d)))
        ramp = case["Tp"]
        sym = case["Ts"]
    else:
        n_tr, exc, ramp, sym = 0, float(abs(2 * np.pi * case.get("delta_Hz", 0) * T_MOD)), np.nan, np.nan
    return {"case": name, "kind": case["kind"], "symbol_period_s": sym, "ramp_duration_s": ramp, "n_transitions_nonzero": n_tr,
            "transition_rate_per_us": n_tr / (T_MOD * 1e6), "total_phase_excursion_rad": exc,
            "mean_abs_freq_offset_Hz": float(np.mean(np.abs(f))), "rms_freq_offset_Hz": float(np.sqrt(np.mean(f ** 2))),
            "peak_abs_freq_offset_Hz": float(np.max(np.abs(f))), "mean_signed_freq_offset_Hz": float(np.mean(f)),
            "amplitude_Vpm": A_H, "amplitude_max_dev": 0.0, "modulated_duration_s": T_MOD, "carrier_center_Hz": IF}


def run_full(case, nd, dt):
    n = int(round(T_END / dt))
    t = np.arange(n) * dt
    env = A_H * np.exp(1j * phase(case, t))
    field = s.old.raqr.A_LO + env * np.exp(2j * np.pi * IF * t)
    res = s.old.sim.run(s.old.raqr, s.old.SimConfig(300, n, dt, n * dt, t), field, Nd=nd, device="cuda")
    step = int(round(1e-9 / dt))
    return 10 ** (res.probeResponse[::step] / 10), res.probeResponse[::step], res.rhovec_final


def populations(rhovec, nd):
    x, p = get_normal_quadrature(nd)
    keep = np.abs(x) <= 3
    x, p = x[keep], p[keep] / p[keep].sum()
    m_cs = (132.9 / 1e3) / 6.02e23
    vx = np.sqrt(get_physical_constant("Boltzmann") * 300 / m_cs) * x
    pops = np.real(rhovec[[0, 5, 10, 15], :])
    return vx, p, pops


def demod(y, phi, t):
    k = np.ones(200) / 200
    base = np.convolve(y, k, "same")
    z = 2 * (y - base) * np.exp(-1j * (2 * np.pi * IF * t + phi))
    return np.convolve(z.real, k, "same") + 1j * np.convolve(z.imag, k, "same")


def main():
    PM.mkdir(parents=True, exist_ok=True)
    (PM / "runs").mkdir(exist_ok=True)
    f1 = M.Controls.from_npz(OUT06 / "artifacts" / "controls_Nd4001.npz")
    mh = M.ControlsMH(OUT06 / "artifacts" / "controls_MH_Nd4001.npz", f1)
    pd.DataFrame([case_stats(k, v) for k, v in CASES.items()]).to_csv(PM / "02_phase_mechanism_configs.csv", index=False)
    plan = [(k, 4001, 1e-9) for k in CASES] + [(k, 8001, 1e-9) for k in CHECKS] + [(k, 4001, 0.5e-9) for k in CHECKS]
    t1 = np.arange(int(round(T_END / 1e-9))) * 1e-9
    t_start, durs = time.perf_counter(), []
    for i, (name, nd, dt) in enumerate(plan, 1):
        out = PM / "runs" / f"{name}_Nd{nd}_dt{dt * 1e9:g}.npz"
        if not out.exists():
            t0 = time.perf_counter()
            case = CASES[name]
            y, probe_db, rhovec = run_full(case, nd, dt)
            phi1 = phase(case, t1)
            env1 = A_H * np.exp(1j * phi1)
            ylin = f1.linear(env1)
            models = {"full": y, "linear": ylin}
            if (nd, dt) == (4001, 1e-9):
                models.update({"static_mh": mh.static_mh(env1), "static_lti_mh": mh.static_lti_mh(env1), "dsh": M.dsh(env1, 2.535e-6, mh)})
            z = {m: demod(v, phi1, t1) for m, v in models.items()}
            vx, pv, pops = populations(rhovec, nd)
            dec = 20
            np.savez_compressed(out, t=t1[::dec], phi=phi1[::dec], probe_db=probe_db[::dec], y_full=y[::dec],
                                **{f"z_{m}": v[::dec] for m, v in z.items()}, vx=vx, pv=pv, pops=pops, Nd=nd, dt=dt)
            durs.append(time.perf_counter() - t0)
        el = time.perf_counter() - t_start
        avg = np.mean(durs) if durs else 0.0
        line = f"[{i}/{len(plan)}] {100 * i / len(plan):5.1f}% | {name} | Nd={nd} dt={dt * 1e9:g} ns | elapsed {el:.0f} s | avg {avg:.1f} s/run | ETA {avg * (len(plan) - i):.0f} s"
        print(line, flush=True)
        (PM / "progress.txt").write_text(line + "\n")


if __name__ == "__main__":
    main()
