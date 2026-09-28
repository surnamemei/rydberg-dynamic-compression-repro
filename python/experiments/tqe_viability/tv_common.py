"""Shared code of the TQE-viability pass (decision rules: results/tqe_viability/00_decision_rules.json)."""
from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "revision2"))

import numpy as np  # noqa: E402

from rev2_common import (DT, KERNEL_OP1, LEVEL_H, OMEGA_P_WEAK, REV2, ROOT, SEEDS3, T0, Kernel, R, RP, raqr, s)  # noqa: E402,F401
from utils.transient_quantum import get_normal_quadrature  # noqa: E402

TV = ROOT / "results" / "tqe_viability"
FV = ROOT / "results" / "final_validation"
REV = ROOT / "results" / "revision"
KERNELS = {0.35: FV / "artifacts" / "controls_LO0.35_Nd4001.npz", 0.425: REV / "artifacts" / "controls_LO0.425_Nd4001.npz", 0.5: KERNEL_OP1}
HBAR = 6.626e-34 / (2.0 * np.pi)
MU_MW = 1443.45 * 1.6e-19 * 5.2918e-11
P_GRID_NS = (400, 360, 320, 300, 280, 260, 250, 240, 230, 220, 210, 200, 190, 180, 170, 160, 150, 140, 130, 125, 120, 115, 110, 105,
             100, 95, 90, 85, 80)


def lo_rabi_hz(a_lo):
    return MU_MW * a_lo / HBAR / (2 * np.pi)


def n_box(f_if):
    """Smallest boxcar of at least 200 samples spanning an integer number of IF periods (200 samples at 5, 10 and 15 MHz;
    P ceil(200/P) for the integer periods P of the P2 grid)."""
    for n in range(200, 2001):
        cyc = n * f_if * 1e-9
        if abs(cyc - round(cyc)) < 1e-9:
            return n
    raise ValueError(f_if)


def demod(y, phi, t, f_if):
    n = n_box(f_if)
    k = np.ones(n) / n
    base = np.convolve(y, k, "same")
    z = 2 * (y - base) * np.exp(-1j * (2 * np.pi * f_if * t + phi))
    return np.convolve(z.real, k, "same") + 1j * np.convolve(z.imag, k, "same")


_VQ = {}


def velocity_weights(temperature=300):
    if temperature not in _VQ:
        m_cs = (132.9 / 1e3) / 6.02e23
        from utils.transient_quantum import get_physical_constant
        sig = np.sqrt(get_physical_constant("Boltzmann") * temperature / m_cs)
        x, p = get_normal_quadrature(4001)
        keep = np.abs(x) <= 3.0
        _VQ[temperature] = (sig * x[keep], p[keep] / p[keep].sum())
    return _VQ[temperature]


def state_diagnostics(rhovec):
    """P5: trace deviation and positivity of the final density matrices (16 x Nv, column-major 4x4 per class)."""
    vx, p = velocity_weights()
    rho = np.asarray(rhovec)
    mats = rho.T.reshape(-1, 4, 4).transpose(0, 2, 1)            # k = i + 4 j  ->  rho[i, j]
    tr = np.real(np.trace(mats, axis1=1, axis2=2))
    herm = 0.5 * (mats + np.conj(np.transpose(mats, (0, 2, 1))))
    ev = np.linalg.eigvalsh(herm)
    mbar = np.tensordot(p, mats, axes=1)
    evbar = np.linalg.eigvalsh(0.5 * (mbar + mbar.conj().T))
    return {"trace_dev_max": float(np.max(np.abs(tr - 1))), "min_eig_min": float(ev.min()),
            "thermal_trace_dev": float(abs(np.real(np.trace(mbar)) - 1)), "thermal_min_eig": float(evbar.min()),
            "hermiticity_dev_max": float(np.max(np.abs(mats - np.conj(np.transpose(mats, (0, 2, 1))))))}


def full_diag(field, cfg_raqr, nd=4001):
    """Probe intensity of the unmodified simulator plus P5 diagnostics of the final density matrices."""
    n = len(field)
    t = np.arange(n) * DT
    res = s.old.sim.run(cfg_raqr, s.old.SimConfig(300, n, DT, n * DT, t), field, Nd=nd, device="cuda")
    return 10 ** (res.probeResponse / 10), state_diagnostics(res.rhovec_final)


def phase_run(out, spec, cfg_raqr, f_if, amp_vpm, kern):
    """Constant-amplitude run (0-50 us CW, 50-110 us modulation) at amplitude amp_vpm; saves z_full/z_linear and diagnostics."""
    import json
    import time
    if out.exists():
        return out
    out.parent.mkdir(parents=True, exist_ok=True)
    n = int(round((T0 + 60e-6) / DT))
    t = np.arange(n) * DT
    phi = RP.phase({**spec, "t_mod": 60e-6}, t)
    env = amp_vpm * np.exp(1j * phi)
    t0 = time.perf_counter()
    y, diag = full_diag(cfg_raqr.A_LO + env * np.exp(2j * np.pi * f_if * t), cfg_raqr)
    wall = time.perf_counter() - t0
    ylin = kern.linear(env, f_if)
    z, zl = demod(y, phi, t, f_if), demod(ylin, phi, t, f_if)
    dec = 20
    np.savez_compressed(out, t=t[::dec], phi=phi[::dec], z_full=z[::dec], z_linear=zl[::dec], IF=f_if, amp=amp_vpm, A_LO=cfg_raqr.A_LO,
                        Omega_p=cfg_raqr.Omega_p, wall_s=wall, spec=json.dumps(spec), diag=json.dumps(diag))
    return out
