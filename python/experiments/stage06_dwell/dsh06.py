"""FU-E: evaluate the pre-registered single-slow-state control (DSH) on Stage-0.6 waveforms.

CPU only.  For every (design, realization, variant, power) of a given tag the waveform is
regenerated exactly as run06.py did, DSH outputs are computed for tau in TAUS, and the same
run06.analyze metrics are stored (rows_<tag>_<design>_dsh.csv).  For realizations with
saved traces the waveform-level NMSE against the full atomic output is also stored.
"""
from __future__ import annotations

import argparse
import os

os.environ.setdefault("OMP_NUM_THREADS", "2")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "2")

import numpy as np
import pandas as pd
from joblib import Parallel, delayed

TAUS = (2.535e-6, 1.94e-6, 4.0e-6)


def task(design, r, powers, tag, cw, n_samples, xis):
    import run06 as R
    import models06 as M
    R.CW_CARRIER = cw
    if n_samples:
        R.N = n_samples
    if xis:
        R.XIS = xis
    ctrl = R.controls()
    mh = R._CTRL["mh"]
    amps, meta, car, labels, noise = R.realization(design, r)
    rows = []
    for v, a in amps.items():
        for p in powers:
            env = ctrl.e1db * 10 ** (p / 20) * a * car
            env_cp = np.r_[env[-R.CP:], env]
            outs = {f"dsh_tau{t * 1e6:.3f}us": M.dsh(env_cp, t, mh) for t in TAUS}
            res, ys, x = R.analyze(outs, env_cp, noise, labels, ctrl, a)
            jid = f"{design}_r{r}_{v}_p{p:+d}_Nd4001_dt1"
            tr = R.OUT / "traces" / tag / f"{jid}.npz"
            nm = {}
            if tr.exists():
                z = np.load(tr)
                w0, w1 = (R.CP + R.SKIP) // R.DEC, len(x) - 2000 // R.DEC
                full = z["y_atomic"][w0:w1]
                cands = {**{m: ys[m][w0:w1] for m in ys}, **{m: z["y_" + m][w0:w1] for m in ("static_lti_mh", "static_mh", "static_lti", "static")}}
                for m, c in cands.items():
                    g = np.vdot(c, full) / np.vdot(c, c)
                    nm[m] = float(np.mean(np.abs(full - g * c) ** 2) / np.mean(np.abs(full) ** 2))
            for m, val in res.items():
                rows.append({"job_id": jid, "tag": tag, "design": design, "realization": r, "variant": v, "Pavg_over_P1dB_dB": p,
                             "model": m, **val, "NMSE_vs_full_trace": nm.get(m, np.nan),
                             **({f"NMSE_vs_full_trace_{k}": vv for k, vv in nm.items() if not k.startswith("dsh")} if m == next(iter(res)) else {})})
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", default="main")
    ap.add_argument("--design", choices=("dwell", "shuffle"), required=True)
    ap.add_argument("--reals", default="0-7")
    ap.add_argument("--powers", default="-6,0,3,6")
    ap.add_argument("--cw-carrier", action="store_true")
    ap.add_argument("--n-samples", type=int, default=None)
    ap.add_argument("--xis", default=None)
    ap.add_argument("--jobs", type=int, default=8)
    a = ap.parse_args()
    lo, hi = (int(v) for v in a.reals.split("-"))
    powers = [int(p) for p in a.powers.split(",")]
    xis = [float(v) for v in a.xis.split(",")] if a.xis else None
    out = Parallel(n_jobs=a.jobs)(delayed(task)(a.design, r, powers, a.tag, a.cw_carrier, a.n_samples, xis) for r in range(lo, hi + 1))
    from common import OUT
    df = pd.DataFrame([row for rows in out for row in rows])
    df.to_csv(OUT / f"rows_{a.tag}_{a.design}_dsh.csv", index=False)
    print(df.groupby(["Pavg_over_P1dB_dB", "variant", "model"]).D_BLA.mean().unstack().round(3).to_string(), flush=True)


if __name__ == "__main__":
    main()
