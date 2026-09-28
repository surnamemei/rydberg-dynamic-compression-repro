"""Rebuild the Stage-0.5 control ingredients at a converged thermal quadrature.

Same definitions as p1db_reference.py / p1db_linear_kernel.py (CW 5 MHz fundamental of
probe intensity; P1dB where the fundamental gain drops 1 dB below the weak-tone slope;
LTI kernel from +-2 V/m 1 ns field impulses on the real and imaginary LO quadratures),
with a longer CW settle (80 us, fundamental over the last 10 us), a denser amplitude grid
and a 60 us kernel.  Stage-0.5 artifacts are not modified.
"""
from __future__ import annotations

import sys
import time

import numpy as np

from common import OUT, s

DT = 1e-9
IF = 5e6
AMPS = np.r_[0, .005, .01, .02, .03, .04, .05, .06, .07, .08, .085, .09, .095, .10, .11, .12, .13, .14, .15, .17, .20, .23, .26, .30, .35, .40, .50, .60, .70]


def run(field, nd):
    n = len(field)
    t = np.arange(n) * DT
    res = s.old.sim.run(s.old.raqr, s.old.SimConfig(300, n, DT, n * DT, t), field, Nd=nd, device="cuda")
    return 10 ** (res.probeResponse / 10)


def cw_table(nd, n=80_000, tail=10_000):
    t = np.arange(n) * DT
    harm, means, harm_early = [], [], []
    for a in AMPS:
        y = run(s.old.raqr.A_LO + a * np.exp(2j * np.pi * IF * t), nd)
        ph = np.exp(-2j * np.pi * IF * t)
        harm.append(2 * np.mean(y[-tail:] * ph[-tail:]))
        harm_early.append(2 * np.mean(y[-2 * tail:-tail] * ph[-2 * tail:-tail]))
        means.append(np.mean(y[-tail:]))
        print("CW", nd, a, abs(harm[-1]), flush=True)
    harm, means = np.array(harm), np.array(means)
    slope = np.mean(np.abs(harm[1:3]) / AMPS[1:3])
    gain = np.full(len(AMPS), np.nan)
    gain[1:] = 20 * np.log10(np.abs(harm[1:]) / (slope * AMPS[1:]))
    j = np.flatnonzero(gain[1:] <= -1)[0] + 1
    e1 = float(np.interp(-1, gain[j - 1:j + 1][::-1], AMPS[j - 1:j + 1][::-1]))
    settle = float(np.max(np.abs(np.array(harm_early) - harm)[1:] / np.abs(harm[1:])))
    return {"amps": AMPS, "harm": harm, "means": means, "gain_db": gain, "e1db": e1, "slope": slope, "cw_settle_rel": settle}


def lti_kernel(nd, n=62_000, pulse=1000, amp=2.0):
    def resp(spike):
        f = np.full(n, s.old.raqr.A_LO, complex)
        f[pulse] += spike
        return run(f, nd)
    hp, hm, ip, im = resp(amp), resp(-amp), resp(1j * amp), resp(-1j * amp)
    return {"hR": (hp[pulse:] - hm[pulse:]) / (2 * amp), "hI": (ip[pulse:] - im[pulse:]) / (2 * amp),
            "lin_check": float(np.max(np.abs(hp[pulse:] + hm[pulse:] - 2 * np.mean(hp[:pulse]))) / np.max(np.abs(hp[pulse:] - hm[pulse:])))}


def main():
    nd = int(sys.argv[1]) if len(sys.argv) > 1 else 4001
    kernel = "--no-kernel" not in sys.argv
    t0 = time.perf_counter()
    cw = cw_table(nd)
    out = {**cw, "Nd": nd, "dt_s": DT, "LO_Vpm": s.old.raqr.A_LO}
    if kernel:
        out.update(lti_kernel(nd))
    (OUT / "artifacts").mkdir(exist_ok=True)
    np.savez_compressed(OUT / "artifacts" / f"controls_Nd{nd}.npz", **out)
    print("E1dB", nd, cw["e1db"], "slope", cw["slope"], "settle", cw["cw_settle_rel"], "seconds", time.perf_counter() - t0, flush=True)


if __name__ == "__main__":
    main()
