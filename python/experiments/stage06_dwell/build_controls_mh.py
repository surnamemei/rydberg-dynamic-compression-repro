"""Multi-harmonic (all-zone) CW tables and the small-signal baseline kernel at Nd=4001.

For each CW amplitude a:  c_k(a) = 2 <I(t) e^{-j k w t}> over the last 10 us of an 80 us
tone (k >= 1), c_0(a) = <I> - dc.  A memoryless envelope model built from these tables
reproduces every CW steady-state intensity waveform exactly (all harmonics, incl. the
rectified baseline).  The zone-0 LTI kernel is the derivative of the normalized probe
response to a +1% LO step (unit DC gain), i.e. the small-signal baseline dynamics at the
same LO operating point used for the Stage-0.5 field kernel.
"""
from __future__ import annotations

import sys

import numpy as np

from common import OUT
from build_controls06 import AMPS as A0, DT, IF, run
from common import s

K = 8
AMPS = np.r_[A0, .8, .9, 1.0]


def main():
    nd = int(sys.argv[1]) if len(sys.argv) > 1 else 4001
    n, tail = 80_000, 10_000
    t = np.arange(n) * DT
    dc = None
    C = np.zeros((len(AMPS), K + 1), complex)
    for i, a in enumerate(AMPS):
        y = run(s.old.raqr.A_LO + a * np.exp(2j * np.pi * IF * t), nd)
        yt, tt = y[-tail:], t[-tail:]
        if i == 0:
            dc = float(np.mean(yt))
        C[i, 0] = np.mean(yt) - dc
        for k in range(1, K + 1):
            C[i, k] = 2 * np.mean(yt * np.exp(-2j * np.pi * k * IF * tt))
        # reconstruction check of the periodic steady state from K harmonics
        rec = dc + C[i, 0].real + sum(np.real(C[i, k] * np.exp(2j * np.pi * k * IF * tt)) for k in range(1, K + 1))
        print("MH", nd, a, "zone mags", np.round(np.abs(C[i]), 7).tolist(), "recon_rel_err", float(np.std(yt - rec) / max(np.std(yt), 1e-15)), flush=True)
    m, n0 = 90_000, 10_000
    f = np.full(m, s.old.raqr.A_LO, complex)
    f[n0:] *= 1.01
    y = run(f, nd)
    step = (y[n0 - 1:] - np.mean(y[n0 - 2000:n0])) / (np.mean(y[-5000:]) - np.mean(y[n0 - 2000:n0]))
    h0 = np.diff(step)
    h0 = h0[:60_000]
    h0 /= h0.sum()
    np.savez_compressed(OUT / "artifacts" / f"controls_MH_Nd{nd}.npz", amps=AMPS, C=C, dc=dc, h0=h0, K=K, Nd=nd)
    print("done", flush=True)


if __name__ == "__main__":
    main()
