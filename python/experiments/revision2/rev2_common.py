"""Shared code for the final receiver-physics validation pass (results/revision2/00_decision_rules.json).

IF- and probe-parametrized versions of the archived pipelines. At IF = 5 MHz with the standard probe they perform the
same operations as the archived code (validated in rev2_run.py): the CW-table rule of final_validation/fv_controls,
the field -> probe conversion of stage05.atomic_intensity, the demodulator of phase_mechanism/pm_run, the receiver
chain and gain fit of stage06_dwell/run06, and the zero-drift phase trajectories of revision/rev_phase.
"""
from __future__ import annotations

import dataclasses
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
EXP = HERE.parent
for p in (EXP / "stage06_dwell", EXP / "stage06_dwell" / "phase_mechanism", EXP / "final_validation", EXP / "revision"):
    sys.path.insert(0, str(p))

import numpy as np  # noqa: E402
from scipy.signal import fftconvolve, sosfilt  # noqa: E402

from common import OUT as OUT06, ROOT, s  # noqa: E402
import run06 as R  # noqa: E402
import rev_phase as RP  # noqa: E402
from build_controls06 import AMPS  # noqa: E402

REV2 = ROOT / "results" / "revision2"
DT = 1e-9
T0 = RP.T0
LEVEL_H = RP.LEVEL_H
SEEDS3 = RP.SEEDS[:3]
FIVE_OFFSETS_MHZ = (0.125, RP.OFFSETS_MHZ[-5], -RP.OFFSETS_MHZ[-5], RP.OFFSETS_MHZ[-2], -RP.OFFSETS_MHZ[-2])
OMEGA_P_WEAK = 2.0 * np.pi * 2.02e6
KERNEL_OP1 = OUT06 / "artifacts" / "controls_Nd4001.npz"


def raqr(lo=0.5, omega_p=None):
    kw = {"A_LO": lo}
    if omega_p is not None:
        kw["Omega_p"] = omega_p
    return dataclasses.replace(s.old.raqr, **kw)


def full(field, cfg_raqr, nd=4001):
    """Probe intensity of the unmodified simulator (same call and conversion as stage05.atomic_intensity)."""
    n = len(field)
    t = np.arange(n) * DT
    res = s.old.sim.run(cfg_raqr, s.old.SimConfig(300, n, DT, n * DT, t), field, Nd=nd, device="cuda")
    return 10 ** (res.probeResponse / 10)


def demod(y, phi, t, f_if):
    """pm_run.demod with the IF as a parameter (200 ns boxcar: an integer number of IF periods at 5, 10, 15 MHz)."""
    k = np.ones(200) / 200
    base = np.convolve(y, k, "same")
    z = 2 * (y - base) * np.exp(-1j * (2 * np.pi * f_if * t + phi))
    return np.convolve(z.real, k, "same") + 1j * np.convolve(z.imag, k, "same")


class Kernel:
    """Small-signal field kernel (hR, hI) and zero-field probe level; the linear model at any IF (models06.Controls.linear)."""

    def __init__(self, hR, hI, dc):
        self.hR, self.hI, self.dc = hR, hI, float(dc)

    @classmethod
    def from_npz(cls, path):
        d = np.load(path)
        return cls(d["hR"], d["hI"], d["means"][0])

    def linear(self, env, f_if):
        t = np.arange(len(env)) / 1e9                      # models06 time base (FS = 1e9)
        field = env * np.exp(2j * np.pi * f_if * t)
        return self.dc + fftconvolve(field.real, self.hR)[:len(env)] + fftconvolve(field.imag, self.hI)[:len(env)]


def cw_table(cfg_raqr, f_if, nd=4001, n=80_000, tail=10_000, kmax=3):
    """fv_controls.cw_table with the IF as a parameter, also returning harmonic zones k = 0..kmax (k = 0: mean level)."""
    t = np.arange(n) * DT
    C, early = [], []
    for a in AMPS:
        y = full(cfg_raqr.A_LO + a * np.exp(2j * np.pi * f_if * t), cfg_raqr, nd)
        row = [np.mean(y[-tail:])]
        for k in range(1, kmax + 1):
            ph = np.exp(-2j * np.pi * k * f_if * t)
            row.append(2 * np.mean(y[-tail:] * ph[-tail:]))
        C.append(row)
        ph1 = np.exp(-2j * np.pi * f_if * t)
        early.append(2 * np.mean(y[-2 * tail:-tail] * ph1[-2 * tail:-tail]))
    C = np.array(C)
    harm, means = C[:, 1], C[:, 0].real
    slope = np.mean(np.abs(harm[1:3]) / AMPS[1:3])
    gain = np.full(len(AMPS), np.nan)
    gain[1:] = 20 * np.log10(np.abs(harm[1:]) / (slope * AMPS[1:]))
    e1 = e1db_rule(AMPS, gain)
    settle = float(np.max(np.abs(np.array(early) - harm)[1:] / np.abs(harm[1:])))
    return {"amps": AMPS, "C": C, "harm": harm, "means": means, "gain_db": gain, "e1db": e1, "slope": slope,
            "cw_settle_rel": settle, "IF_Hz": f_if}


def e1db_rule(amps, gain):
    """Archived rule: first -1 dB crossing of the fundamental gain re the weak-tone slope (NaN if none on the grid)."""
    hit = np.flatnonzero(gain[1:] <= -1)
    if not len(hit):
        return float("nan")
    j = hit[0] + 1
    return float(np.interp(-1, gain[j - 1:j + 1][::-1], amps[j - 1:j + 1][::-1]))


def extend_table(cfg_raqr, f_if, tab, nd=4001, n=80_000, tail=10_000, kmax=3, a_max=20.0):
    """Deviation D1 (02_deviations.json): amplitudes 0.80 x 1.15^j V/m above the archived grid, until >= 3 E1dB(IF)."""
    t = np.arange(n) * DT
    amps, C = list(tab["amps"]), [list(r) for r in tab["C"]]
    e1 = float(tab["e1db"])
    j = 0
    while True:
        a = 0.80 * 1.15 ** j
        if (np.isfinite(e1) and a >= 3 * e1) or a > a_max:
            break
        y = full(cfg_raqr.A_LO + a * np.exp(2j * np.pi * f_if * t), cfg_raqr, nd)
        row = [np.mean(y[-tail:])] + [2 * np.mean(y[-tail:] * np.exp(-2j * np.pi * k * f_if * t[-tail:])) for k in range(1, kmax + 1)]
        amps.append(a)
        C.append(row)
        j += 1
        if not np.isfinite(e1):
            g = 20 * np.log10(np.abs(np.array(C)[1:, 1]) / (float(tab["slope"]) * np.array(amps)[1:]))
            e1 = e1db_rule(np.array(amps), np.r_[np.nan, g])
    amps, C = np.array(amps), np.array(C)
    gain = np.full(len(amps), np.nan)
    gain[1:] = 20 * np.log10(np.abs(C[1:, 1]) / (float(tab["slope"]) * amps[1:]))
    return {**{k: v for k, v in tab.items()}, "amps": amps, "C": C, "harm": C[:, 1], "means": C[:, 0].real, "gain_db": gain,
            "e1db": e1db_rule(amps, gain), "n_archived_grid": len(tab["amps"]), "a_max_needed": 3 * e1db_rule(amps, gain)}


def lti_kernel(cfg_raqr, nd=4001, n=62_000, pulse=1000, amp=2.0):
    """fv_controls.lti_kernel with the RAQR configuration as a parameter."""
    lo = cfg_raqr.A_LO

    def resp(spike):
        f = np.full(n, lo, complex)
        f[pulse] += spike
        return full(f, cfg_raqr, nd)
    hp, hm, ip, im = resp(amp), resp(-amp), resp(1j * amp), resp(-1j * amp)
    return {"hR": (hp[pulse:] - hm[pulse:]) / (2 * amp), "hI": (ip[pulse:] - im[pulse:]) / (2 * amp),
            "lin_check": float(np.max(np.abs(hp[pulse:] + hm[pulse:] - 2 * np.mean(hp[:pulse]))) / np.max(np.abs(hp[pulse:] - hm[pulse:])))}


def zone_interp(table, a, k):
    C = table["C"][:, k]
    return np.interp(a, table["amps"], C.real) + 1j * np.interp(a, table["amps"], C.imag)


def static_zone(env, table, k, f_if):
    """Memoryless CW-matched output of harmonic zone k (k = 0: baseline shift relative to zero field)."""
    a = np.abs(env)
    if k == 0:
        return (zone_interp(table, a, 0) - table["C"][0, 0]).real
    u = np.ones(len(env), complex)
    np.divide(env, a, out=u, where=a > 0)
    t = np.arange(len(env)) * DT
    return np.real(zone_interp(table, a, k) * u ** k * np.exp(2j * np.pi * k * f_if * t))


def bb(y, dc, f_if):
    """run06.bb with the IF as a parameter."""
    t = np.arange(len(y)) / 1e9
    return sosfilt(R.SOS, (y - dc) * np.exp(-2j * np.pi * f_if * t))[::R.DEC]


def fitted_gain(yb, xs, w0, w1):
    X = np.column_stack([xs, np.ones_like(xs)])
    coef, *_ = np.linalg.lstsq(X, yb[w0:w1], rcond=None)
    return complex(coef[0])
