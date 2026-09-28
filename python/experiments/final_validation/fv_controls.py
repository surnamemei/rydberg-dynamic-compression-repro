"""CW table, E1dB and small-signal kernel at a second LO (preregistered protocol, identical to
stage06_dwell/build_controls06.py except for the LO amplitude)."""
from __future__ import annotations

import sys
import time

import numpy as np

from fv_common import DT, FV, IF, run_field
from build_controls06 import AMPS


def cw_table(nd, lo, n=80_000, tail=10_000):
    t = np.arange(n) * DT
    harm, means, early = [], [], []
    for a in AMPS:
        y, _, _ = run_field(lo + a * np.exp(2j * np.pi * IF * t), nd, lo)
        ph = np.exp(-2j * np.pi * IF * t)
        harm.append(2 * np.mean(y[-tail:] * ph[-tail:]))
        early.append(2 * np.mean(y[-2 * tail:-tail] * ph[-2 * tail:-tail]))
        means.append(np.mean(y[-tail:]))
        print("CW", lo, nd, a, abs(harm[-1]), flush=True)
    harm, means = np.array(harm), np.array(means)
    slope = np.mean(np.abs(harm[1:3]) / AMPS[1:3])
    gain = np.full(len(AMPS), np.nan)
    gain[1:] = 20 * np.log10(np.abs(harm[1:]) / (slope * AMPS[1:]))
    j = np.flatnonzero(gain[1:] <= -1)[0] + 1
    e1 = float(np.interp(-1, gain[j - 1:j + 1][::-1], AMPS[j - 1:j + 1][::-1]))
    settle = float(np.max(np.abs(np.array(early) - harm)[1:] / np.abs(harm[1:])))
    return {"amps": AMPS, "harm": harm, "means": means, "gain_db": gain, "e1db": e1, "slope": slope, "cw_settle_rel": settle}


def lti_kernel(nd, lo, n=62_000, pulse=1000, amp=2.0):
    def resp(spike):
        f = np.full(n, lo, complex)
        f[pulse] += spike
        return run_field(f, nd, lo)[0]
    hp, hm, ip, im = resp(amp), resp(-amp), resp(1j * amp), resp(-1j * amp)
    return {"hR": (hp[pulse:] - hm[pulse:]) / (2 * amp), "hI": (ip[pulse:] - im[pulse:]) / (2 * amp),
            "lin_check": float(np.max(np.abs(hp[pulse:] + hm[pulse:] - 2 * np.mean(hp[:pulse]))) / np.max(np.abs(hp[pulse:] - hm[pulse:])))}


def main():
    lo = float(sys.argv[1]) if len(sys.argv) > 1 else 0.35
    nd = 4001
    t0 = time.perf_counter()
    out = {**cw_table(nd, lo), **lti_kernel(nd, lo), "Nd": nd, "dt_s": DT, "LO_Vpm": lo}
    (FV / "artifacts").mkdir(parents=True, exist_ok=True)
    np.savez_compressed(FV / "artifacts" / f"controls_LO{lo:g}_Nd{nd}.npz", **out)
    print("E1dB", lo, out["e1db"], "slope", out["slope"], "settle", out["cw_settle_rel"], "lin_check", out["lin_check"],
          "seconds", time.perf_counter() - t0, flush=True)


if __name__ == "__main__":
    main()
