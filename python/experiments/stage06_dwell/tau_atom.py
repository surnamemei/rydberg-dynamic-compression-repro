"""Step 2: physical atomic relaxation timescale from full-thermal step responses.

Definitions (all on the normalized excess e(t) = (y(t) - y_final) / (y(t0-) - y_final),
t measured from the step, envelope E(t) = max_{t' >= t} |e(t')|):
  tau_1e : first t with E(t) <= 1/e
  tau_90 : first t with E(t) <= 0.1   (90% of the relaxation complete)
  tau_dom: single-exponential slope of log|e| over the tail window 0.3 >= E >= 0.03,
           reported with its RMS log-residual; a two-exponential fit is also reported.
The IF-envelope responses are demodulated with a one-IF-period (200 ns) boxcar,
which nulls DC and 2*IF without adding a filter pole; this limits resolution to ~0.2 us.
"""
from __future__ import annotations

import json

import numpy as np
import pandas as pd
from scipy.optimize import curve_fit

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from common import OUT, ROOT, s

DT = 1e-9
IF = s.old.IF
NDS = (1501, 4001, 8001)
PERIOD = int(round(1 / (IF * DT)))  # 200 samples


def run_field(field, nd):
    n = len(field)
    t = np.arange(n) * DT
    cfg = s.old.SimConfig(300, n, DT, n * DT, t)
    res = s.old.sim.run(s.old.raqr, cfg, field, Nd=nd, device="cuda")
    return 10 ** (res.probeResponse / 10)


def if_envelope(y):
    t = np.arange(len(y)) * DT
    z = 2 * (y - np.mean(y)) * np.exp(-2j * np.pi * IF * t)
    k = np.ones(PERIOD) / PERIOD
    zr = np.convolve(z.real, k, mode="same")
    zi = np.convolve(z.imag, k, mode="same")
    return zr + 1j * zi


def relax_metrics(e, t):
    """e: normalized excess starting at the step (e[0] ~ 1)."""
    env = np.maximum.accumulate(np.abs(e)[::-1])[::-1]
    def first(th):
        idx = np.flatnonzero(env <= th)
        return float(t[idx[0]]) if len(idx) else np.nan
    out = {"tau_1e_s": first(np.exp(-1)), "tau_90_s": first(0.1)}
    tail = (env <= 0.3) & (env >= 0.03) & (np.abs(e) > 1e-6)
    if np.count_nonzero(tail) > 50:
        slope, icpt = np.polyfit(t[tail], np.log(np.abs(e[tail])), 1)
        resid = np.log(np.abs(e[tail])) - (slope * t[tail] + icpt)
        out["tau_dom_s"] = float(-1 / slope) if slope < 0 else np.nan
        out["tau_dom_fit_rms_log_resid"] = float(np.sqrt(np.mean(resid ** 2)))
        out["tau_dom_window_s"] = f"{t[tail][0]:.3e}-{t[tail][-1]:.3e}"
    # Two-exponential fit on the signed excess after the first 100 ns (fast optical transients excluded).
    m = t >= 100e-9
    def f2(tt, a1, t1, a2, t2):
        return a1 * np.exp(-tt / t1) + a2 * np.exp(-tt / t2)
    try:
        p, _ = curve_fit(f2, t[m], e[m], p0=(0.5, 0.3e-6, 0.5, 3e-6), bounds=([-5, 1e-9, -5, 1e-8], [5, 1e-4, 5, 1e-3]), maxfev=20000)
        if p[1] > p[3]:
            p = p[[2, 3, 0, 1]]
        r = e[m] - f2(t[m], *p)
        out.update({"fit2_a_fast": p[0], "fit2_tau_fast_s": p[1], "fit2_a_slow": p[2], "fit2_tau_slow_s": p[3],
                    "fit2_rms_resid": float(np.sqrt(np.mean(r ** 2)))})
    except Exception as exc:  # noqa: BLE001
        out["fit2_error"] = repr(exc)
    return out


def experiments(nd):
    A = s.old.raqr.A_LO
    res, traces = [], {}
    # (1) small LO amplitude step: +1% and -1% (linearity check).
    n0, n = 10_000, 90_000
    for sign in (+1, -1):
        field = np.full(n, A, complex)
        field[n0:] *= 1 + sign * 0.01
        y = run_field(field, nd)
        yf = np.mean(y[-5000:])
        e = (y[n0:] - yf) / (np.mean(y[n0 - 500:n0]) - yf)
        t = np.arange(len(e)) * DT
        res.append({"experiment": f"LO_step_{'+' if sign > 0 else '-'}1pct", "Nd": nd, **relax_metrics(e, t)})
        traces[res[-1]["experiment"]] = (t, e)
    # (2) IF-envelope steps: linear point and around P1dB; (3) large-signal plateau.
    warm = 60_000
    for name, ea, eb, hold in (("IF_step_small_0.05to0.10E1", .05, .10, None),
                               ("IF_step_P1dB_1.00to1.05E1", 1.00, 1.05, None),
                               ("IF_step_P1dB_1.05to1.00E1", 1.05, 1.00, None),
                               ("IF_plateau_0.5to1.5E1_on", .5, 1.5, 30_000)):
        total = warm + (hold or 0) + 60_000
        tt = np.arange(total) * DT
        amp = np.full(total, ea * s.E1)
        amp[warm:] = eb * s.E1
        if hold:
            amp[warm + hold:] = ea * s.E1
        # Start CW from t=0; the steady-state init only sees the t=0 field, so warm up for 60 us.
        field = A + amp * np.exp(2j * np.pi * IF * tt)
        y = run_field(field, nd)
        z = if_envelope(y)
        segs = [(name, warm, warm + hold if hold else total)]
        if hold:
            segs = [(name, warm, warm + hold), (name.replace("_on", "_off_recovery"), warm + hold, total)]
        for sname, i0, i1 in segs:
            pre = np.mean(z[i0 - 5 * PERIOD:i0 - PERIOD])
            fin = np.mean(z[i1 - 5 * PERIOD - PERIOD:i1 - PERIOD])
            u = (pre - fin) / abs(pre - fin)
            e = np.real((z[i0:i1 - PERIOD // 2] - fin) * np.conj(u)) / abs(pre - fin)
            # The boxcar spreads the step over +-100 ns; start at the step centre.
            t = np.arange(len(e)) * DT
            gain_pre, gain_fin = abs(pre) / amp[i0 - 1], abs(fin) / amp[i1 - 1]
            res.append({"experiment": sname, "Nd": nd, "envelope_gain_before": gain_pre, "envelope_gain_after": gain_fin,
                        "pre_step_drift_rel": float(abs(z[i0 - PERIOD] - z[i0 - 20 * PERIOD]) / abs(pre)),
                        **relax_metrics(e, t)})
            traces[sname] = (t, e)
    return res, traces


def main():
    rows, all_tr = [], {}
    for nd in NDS:
        r, tr = experiments(nd)
        rows += r
        all_tr[nd] = tr
        for x in r:
            print(json.dumps({k: (round(v, 9) if isinstance(v, float) else v) for k, v in x.items()}), flush=True)
    # (4) archived weak-impulse LTI kernel (Stage-0.5, Nd=1501): envelope of the IF-demodulated impulse response.
    kern = np.load(ROOT / "artifacts" / "p1db_lti_kernel.npz")
    h = kern["hR"]
    t = np.arange(len(h)) * DT
    hb = if_envelope(np.r_[np.zeros(PERIOD), h])[PERIOD:]
    mag = np.abs(hb)
    i0 = int(np.argmax(mag))
    e = mag[i0:] / mag[i0]
    rows.append({"experiment": "weak_impulse_kernel_IF_envelope(Stage-0.5 artifact)", "Nd": 1501, **relax_metrics(e, t[:len(e)])})
    all_tr.setdefault(1501, {})["weak_impulse"] = (t[:len(e)], e)
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "03_tau_atom_results.csv", index=False)
    np.savez_compressed(OUT / "tau_atom_traces.npz", **{f"{nd}__{k}__t": v[0][::10] for nd, d in all_tr.items() for k, v in d.items()},
                        **{f"{nd}__{k}__e": v[1][::10] for nd, d in all_tr.items() for k, v in d.items()})
    print(df.to_string(), flush=True)


if __name__ == "__main__":
    main()
