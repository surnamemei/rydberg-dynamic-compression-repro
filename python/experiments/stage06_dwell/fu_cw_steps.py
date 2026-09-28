"""Follow-up A/B (pre-registered in 11_followup_preregistration.json).

FU-A: tau sensitivity - +5% IF-envelope step starting at 0.8/1.0/1.2 x converged E1dB.
FU-B: one 60 us plateau a_L -> a_H -> a_L at the dwell-family levels (0 and +3 dB) with an
      unmodulated carrier and with the 1 Msym/s QPSK phase carrier.
Gains are |IF envelope| / |E| (one-IF-period boxcar demodulation), normalized by the
CW small-signal gain |c1(a)|/a at a -> 0 from the all-zone table (Nd=4001).
"""
from __future__ import annotations

import json

import numpy as np
import pandas as pd

import tau_atom as T
import waveforms06 as W
from common import OUT, s

ND = 4001
MH = np.load(OUT / "artifacts" / "controls_MH_Nd4001.npz")
E1C = float(np.load(OUT / "artifacts" / "controls_Nd4001.npz")["e1db"])
G0 = abs(MH["C"][1, 1]) / MH["amps"][1]
JOBS = OUT / "jobs" / "fu_cw_steps"
JOBS.mkdir(parents=True, exist_ok=True)


def cw_gain(a):
    c1 = np.interp(a, MH["amps"], np.abs(MH["C"][:, 1]))
    return c1 / a / G0


def mark(name, payload):
    tmp = JOBS / f"{name}.tmp"
    tmp.write_text(json.dumps(payload, default=float))
    tmp.replace(JOBS / f"{name}.json")


def tau_steps():
    rows = []
    for f in (0.8, 1.0, 1.2):
        name = f"tau_{f:.1f}xE1"
        if (JOBS / f"{name}.json").exists():
            rows.append(json.loads((JOBS / f"{name}.json").read_text()))
            continue
        warm, total = 60_000, 120_000
        tt = np.arange(total) * T.DT
        amp = np.full(total, f * E1C)
        amp[warm:] *= 1.05
        z = T.if_envelope(T.run_field(s.old.raqr.A_LO + amp * np.exp(2j * np.pi * T.IF * tt), ND))
        pre = np.mean(z[warm - 5 * T.PERIOD:warm - T.PERIOD])
        fin = np.mean(z[total - 6 * T.PERIOD:total - T.PERIOD])
        u = (pre - fin) / abs(pre - fin)
        e = np.real((z[warm:total - T.PERIOD // 2] - fin) * np.conj(u)) / abs(pre - fin)
        row = {"experiment": name, "E_start_Vpm": f * E1C, "E_start_over_E1dB_converged": f, "Nd": ND, **T.relax_metrics(e, np.arange(len(e)) * T.DT)}
        mark(name, row)
        rows.append(row)
        print(row, flush=True)
    pd.DataFrame(rows).to_csv(OUT / "03b_tau_atom_sensitivity.csv", index=False)


def plateaus():
    rows, traces = [], {}
    warm, hold, tail = 60_000, 60_000, 60_000
    n = warm + hold + tail
    tt = np.arange(n) * T.DT
    for pdb in (0, 3):
        a_l = 0.5774 * 10 ** (pdb / 20) * E1C
        a_h = 3 * a_l
        amp = np.full(n, a_l)
        amp[warm:warm + hold] = a_h
        for carrier in ("cw", "qpsk"):
            name = f"plateau_p{pdb:+d}_{carrier}"
            if carrier == "cw":
                car = np.ones(n, complex)
            else:
                car, _ = W.phase_carrier(np.random.default_rng(20260607), n, 1000, 300)
            y = T.run_field(s.old.raqr.A_LO + amp * car * np.exp(2j * np.pi * T.IF * tt), ND)
            z = T.if_envelope(y)
            g = np.abs(z) / amp / G0
            us = lambda t_us: int(t_us * 1000)  # noqa: E731
            seg_h = g[warm:warm + hold]
            seg_r = g[warm + hold:]
            g_h_end = float(np.median(seg_h[-us(10):-T.PERIOD]))
            g_l_pre = float(np.median(g[warm - us(10):warm - T.PERIOD]))
            g_l_end = float(np.median(seg_r[-us(10):-T.PERIOD]))

            def settle(seg, final, tol=0.02):
                env = np.maximum.accumulate(np.abs(seg - final)[::-1])[::-1]
                i = np.flatnonzero(env <= tol * final)
                return float(i[0] * T.DT) if len(i) else np.nan
            row = {"experiment": name, "Pavg_over_P1dB_dB": pdb, "carrier": carrier, "a_L_Vpm": a_l, "a_H_Vpm": a_h, "Nd": ND,
                   "gain_low_before": g_l_pre, "gain_plateau_first_us": float(np.median(seg_h[T.PERIOD:us(1)])),
                   "gain_plateau_last10us": g_h_end, "cw_static_gain_aH": float(cw_gain(a_h)), "cw_static_gain_aL": float(cw_gain(a_l)),
                   "plateau_vs_cw_static_rel": g_h_end / float(cw_gain(a_h)) - 1,
                   "plateau_settle_2pct_s": settle(seg_h[T.PERIOD:], g_h_end),
                   "gain_low_first_us_after": float(np.median(seg_r[T.PERIOD:us(1)])), "gain_low_1_to_5us_after": float(np.median(seg_r[us(1):us(5)])),
                   "gain_low_end": g_l_end, "recovery_settle_2pct_s": settle(seg_r[T.PERIOD:], g_l_end)}
            rows.append(row)
            traces[name] = g[::20]
            mark(name, row)
            print(row, flush=True)
    pd.DataFrame(rows).to_csv(OUT / "12_plateau_steady_state.csv", index=False)
    np.savez_compressed(OUT / "plateau_gain_traces.npz", dt_s=20e-9, warm_s=warm * 1e-9, hold_s=hold * 1e-9, **traces)


if __name__ == "__main__":
    tau_steps()
    plateaus()
