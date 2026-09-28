"""POST-HOC diagnostic on already-saved phase-mechanism outputs (no new simulation).

Question (for claims C6/C9): after a phase transition at constant amplitude, does the full
model's IF (beat) response LAG the new drive phase and relax over microseconds (an
operational "loss and recovery of phase lock"), or does only its magnitude change?

theta(t) = arg(z_full * conj(z_linear)): phase of the full response relative to the
small-signal model driven by the identical input (removes the small-signal lag).
If the full response kept its old phase, theta would jump by -dphi at the transition.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
RUNS = ROOT / "results" / "stage06_dwell_physics" / "phase_mechanism" / "runs"
OUT = Path(__file__).resolve().parent
T0 = 50.0  # us


def load(case):
    z = np.load(RUNS / f"{case}_Nd4001_dt1.npz")
    t = z["t"] * 1e6
    zf, zl = z["z_full"], z["z_linear"]
    g = np.abs(zf) / np.abs(zl)
    u = zf * np.conj(zl) / np.maximum(np.abs(zf * np.conj(zl)), 1e-30)  # unit phasor of the relative phase
    return t, z["phi"], g, u, np.abs(zf), np.abs(zl)


def cang(u):
    """circular mean phase of unit phasors"""
    return float(np.angle(np.mean(u)))


def cstd(u):
    return float(np.sqrt(max(-2 * np.log(max(abs(np.mean(u)), 1e-12)), 0)))


def smooth(x, n):
    return np.convolve(x, np.ones(n) / n, "same")


def transition_table(case, ts_us):
    t, phi, g, u, af, al = load(case)
    gs = smooth(g, 10)  # 200 ns
    us = smooth(u.real, 10) + 1j * smooth(u.imag, 10)
    rows, traces = [], []
    nb = int(round(60 / ts_us))
    for k in range(1, nb):
        tb = T0 + k * ts_us
        i0 = np.searchsorted(t, tb)
        dphi = phi[min(i0 + 50, len(phi) - 1)] - phi[max(i0 - 50, 0)]
        if abs(dphi) < 1e-6:
            continue
        pre = (t >= tb - 0.6) & (t < tb - 0.2)
        post = (t >= tb) & (t < min(tb + ts_us - 0.3, 109.4))
        g_pre, th_pre = np.median(g[pre]), cang(u[pre])
        af_pre, al_pre = np.median(af[pre]), np.median(al[pre])
        dg = gs[post] - g_pre
        dth = np.angle(us[post] * np.exp(-1j * th_pre))  # wrapped phase deviation
        af_rel, al_rel = smooth(af, 10)[post] / af_pre, smooth(al, 10)[post] / al_pre
        tt = t[post] - tb
        imin = int(np.argmin(dg))
        imax_th = int(np.argmax(np.abs(dth)))
        # recovery of the gain dip: time after the minimum until |dg| < |dg_min|/e (NaN if not within window)
        rec = np.flatnonzero(np.abs(dg[imin:]) <= abs(dg[imin]) / np.e)
        rows.append({"case": case, "t_boundary_us": tb, "dphi_rad": dphi, "g_pre": g_pre, "dip_depth": float(dg[imin]), "dip_time_us": float(tt[imin]),
                     "dip_recovery_1e_us": float(tt[imin + rec[0]] - tt[imin]) if len(rec) else np.nan,
                     "max_abs_dtheta_rad": float(abs(dth[imax_th])), "dtheta_at_max_rad": float(dth[imax_th]), "time_max_dtheta_us": float(tt[imax_th]),
                     "window_us": float(tt[-1]),
                     "min_absfull_rel_first_0p5us": float(np.min(af_rel[tt < .5])), "min_abslin_rel_first_0p5us": float(np.min(al_rel[tt < .5])),
                     "max_absfull_rel_first_0p5us": float(np.max(af_rel[tt < .5])), "max_abslin_rel_first_0p5us": float(np.max(al_rel[tt < .5]))})
        traces.append((dphi, tt, dg, dth))
    return pd.DataFrame(rows), traces


def main():
    tab, tr = transition_table("QPSK_SLOW", 4.0)
    tab.to_csv(OUT / "posthoc_isolated_transitions.csv", index=False)
    print(tab.round(3).to_string())
    # steady-state relative phase per case (last 20 us) vs CW
    rows = []
    for c in ("CW", "OFF_+0.125", "OFF_+1.31", "OFF_-1.31", "OFF_+2.62", "OFF_-2.62", "QPSK_SLOW", "QPSK_REF", "QPSK_REF_s2", "QPSK_FAST", "QPSK_JUMP", "QPSK_LONGRAMP", "TOGGLE"):
        t, phi, g, u, af, al = load(c)
        m = (t >= 89.5) & (t < 109.5)
        pre = (t >= 40) & (t < 50)
        rows.append({"case": c, "dtheta_final_vs_pre_rad": float(np.angle(np.mean(u[m]) * np.exp(-1j * cang(u[pre])))),
                     "circ_std_final_rad": cstd(u[m]), "g_final": float(np.median(g[m])), "g_std_last20": float(np.std(g[m]))})
    ss = pd.DataFrame(rows)
    ss.to_csv(OUT / "posthoc_steady_relative_phase.csv", index=False)
    print(ss.round(3).to_string())

    fig, ax = plt.subplots(1, 2, figsize=(11, 3.8), sharex=True)
    for dphi, tt, dg, dth in tr:
        c = "#eb6834" if abs(abs(dphi) - np.pi) < .1 else "#2a78d6"
        ax[0].plot(tt, dg, color=c, lw=.9)
        ax[1].plot(tt, dth / np.pi, color=c, lw=.9)
    for a in ax:
        a.axhline(0, color="k", lw=.6)
        a.set_xlabel("time after phase transition (µs)")
    ax[0].plot([], [], color="#2a78d6", label="±π/2 steps")
    ax[0].plot([], [], color="#eb6834", label="π steps")
    ax[0].legend(fontsize=7, frameon=False)
    ax[0].set_ylabel("gain change vs pre-transition")
    ax[1].set_ylabel("response phase change vs pre-transition (units of π)")
    fig.suptitle("POST-HOC: isolated transitions (QPSK_SLOW, 300 ns ramps), full response relative to the small-signal response")
    fig.tight_layout()
    fig.savefig(OUT / "posthoc_isolated_transition_phase.png", dpi=150)


if __name__ == "__main__":
    main()
