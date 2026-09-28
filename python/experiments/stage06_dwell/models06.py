"""Receiver-model outputs (probe intensity) for Stage-0.6.

Conventions are those validated in Stage-0.5 (stage05.static_intensity /
static_lti_intensity / p1db_comm_stage0.linear_intensity):
  M_STATIC     memoryless complex CW-fundamental map  c(|E|) * E/|E|  on the IF envelope;
  M_LTI_STATIC Hammerstein: static complex gain G(|E|) = c(|E|) / (|E| * H_IF) applied to
               the envelope, then the small-signal LTI kernel (static -> LTI ordering),
               so that a CW tone reproduces the CW curve exactly;
  M_LINEAR     small-signal LTI kernel only (no compression), used as the AIR reference;
  M_FULL       full thermal transient simulator.
Only the source of the tables differs: Stage-0.6 uses tables rebuilt at a converged Nd.
"""
from __future__ import annotations

import numpy as np
from scipy.signal import fftconvolve

from common import s

FS = 1e9
IF = 5e6


class Controls:
    def __init__(self, amps, harm, dc, hR, hI, e1db, hif_n=80_000, hif_tail=10_000, label=""):
        self.amps, self.harm, self.dc, self.hR, self.hI, self.e1db, self.label = amps, harm, float(dc), hR, hI, float(e1db), label
        n = hif_n
        t = np.arange(n) / FS
        probe = self.linear(np.full(n, .01, complex))
        self.hif = 2 * np.mean((probe[-hif_tail:] - self.dc) * np.exp(-2j * np.pi * IF * t[-hif_tail:])) / .01
        g = np.zeros(len(amps), complex)
        g[1:] = harm[1:] / (amps[1:] * self.hif)
        g[0] = g[1]
        self.gain = g

    @classmethod
    def from_npz(cls, path, **kw):
        d = np.load(path)
        return cls(d["amps"], d["harm"], d["means"][0], d["hR"], d["hI"], d["e1db"], label=str(path), **kw)

    @classmethod
    def stage05(cls):
        o = s.old
        return cls(o.AMPS, o.HARM, o.dc, o.HR, o.HI, o.E1, hif_n=40_000, hif_tail=4_000, label="stage05_Nd1501")

    def _cw(self, a, table):
        return np.interp(a, self.amps, table.real) + 1j * np.interp(a, self.amps, table.imag)

    def static(self, env):
        a = np.abs(env)
        out = np.zeros(len(env), complex)
        np.divide(self._cw(a, self.harm) * env, a, out=out, where=a > 0)
        t = np.arange(len(env)) / FS
        return self.dc + np.real(out * np.exp(2j * np.pi * IF * t))

    def linear(self, env):
        t = np.arange(len(env)) / FS
        field = env * np.exp(2j * np.pi * IF * t)
        return self.dc + fftconvolve(field.real, self.hR)[:len(env)] + fftconvolve(field.imag, self.hI)[:len(env)]

    def static_lti(self, env):
        return self.linear(self._cw(np.abs(env), self.gain) * env)


class ControlsMH:
    """All-zone CW-matched controls.

    M_STATIC_MH      I = dc + c0(a) + sum_k Re(c_k(a) u^k e^{jkwt}),  a=|E|, u=E/a  (memoryless;
                     reproduces every CW periodic steady state, incl. the rectified baseline).
    M_LTI_STATIC_MH  zone 1: Stage-0.5 Hammerstein (static gain -> field LTI kernel);
                     zone 0: static c0(a) -> small-signal baseline kernel h0 (unit DC gain);
                     zones k>=2: static.  A CW tone again reproduces the CW steady state.
    """

    def __init__(self, path, f1: Controls):
        d = np.load(path)
        self.amps, self.C, self.dc, self.h0, self.K = d["amps"], d["C"], float(d["dc"]), d["h0"], int(d["K"])
        self.f1 = f1

    def _zone(self, a, k):
        return np.interp(a, self.amps, self.C[:, k].real) + 1j * np.interp(a, self.amps, self.C[:, k].imag)

    def _zones(self, env, kmin):
        a = np.abs(env)
        u = np.ones(len(env), complex)
        np.divide(env, a, out=u, where=a > 0)
        t = np.arange(len(env)) / FS
        out = np.zeros(len(env))
        for k in range(kmin, self.K + 1):
            out += np.real(self._zone(a, k) * u ** k * np.exp(2j * np.pi * k * IF * t))
        return out, a

    def static_mh(self, env):
        out, a = self._zones(env, 1)
        return self.dc + self._zone(a, 0).real + out

    def static_lti_mh(self, env):
        high, a = self._zones(env, 2)
        base = fftconvolve(self._zone(a, 0).real, self.h0)[:len(env)]
        return self.f1.static_lti(env) - self.f1.dc + self.dc + base + high


def dsh(env, tau_s, mh: ControlsMH):
    """Single-slow-state dynamic Hammerstein (pre-registered FU-E control, not fitted).

    P_s' = (|E|^2 - P_s)/tau_s, a_eff = sqrt(P_s); the fundamental zone uses the CW-matched
    Hammerstein gain evaluated at a_eff instead of |E|; zones 0 and k>=2 are identical to
    M_LTI_STATIC_MH.  tau_s -> 0 recovers M_LTI_STATIC_MH exactly; CW tones are reproduced.
    """
    from scipy.signal import lfilter
    f1 = mh.f1
    p = np.abs(env) ** 2
    alpha = (1 / FS) / tau_s
    ps = lfilter([alpha], [1, -(1 - alpha)], p, zi=[(1 - alpha) * np.mean(p)])[0]
    g = f1._cw(np.sqrt(np.maximum(ps, 0)), f1.gain)
    high, a = mh._zones(env, 2)
    base = fftconvolve(mh._zone(a, 0).real, mh.h0)[:len(env)]
    return f1.linear(g * env) - f1.dc + mh.dc + base + high


def atomic(env, nd, dt=1e-9):
    return s.atomic_intensity(env, dt=dt, nd=nd)
