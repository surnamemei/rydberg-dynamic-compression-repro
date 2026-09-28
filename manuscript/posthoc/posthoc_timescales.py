"""POST-HOC (saved data only; no new simulation) - timescale evidence for claim C9.

1. Compare the transient after a POWER step and after PHASE steps at the same amplitude
   a_H = 0.179 V/m (+3 dB dwell-family level):
     power step : FU-B plateau_p+3_cw / p+0_cw (a_L -> a_H, unmodulated carrier; gain traces)
     phase step : phase-mechanism QPSK_SLOW isolated transitions at constant a_H
   Metrics: dominant oscillation period in 0.3-3.7 us (zero-padded FFT of the detrended excess)
   and envelope decay (1/e of the running-max envelope).
2. Nd-independence of the oscillation period in the Step-2 step responses (Nd 1501/4001/8001):
   a velocity-quadrature recurrence would scale with Nd; a physical mode would not.
"""
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
R06 = ROOT / "results" / "stage06_dwell_physics"
HERE = Path(__file__).resolve().parent


def dominant_period(x, dt_us, lo=0.8, hi=4.0):
    x = x - np.polyval(np.polyfit(np.arange(len(x)), x, 2), np.arange(len(x)))
    n = 1 << 14
    S = np.abs(np.fft.rfft(x * np.hanning(len(x)), n)) ** 2
    f = np.fft.rfftfreq(n, dt_us)  # cycles per us
    m = (f >= 1 / hi) & (f <= 1 / lo)
    return float(1 / f[m][np.argmax(S[m])])


def env_decay(x, t):
    env = np.maximum.accumulate(np.abs(x)[::-1])[::-1]
    i = np.flatnonzero(env <= env[0] / np.e)
    return float(t[i[0]]) if len(i) else np.nan


def power_vs_phase():
    rows = []
    z = np.load(R06 / "plateau_gain_traces.npz")
    dt = float(z["dt_s"]) * 1e6
    for p in ("+0", "+3"):
        g = z[f"plateau_p{p}_cw"]
        i0 = int(round(float(z["warm_s"]) * 1e6 / dt))
        seg = g[i0: i0 + int(12 / dt)]
        tt = np.arange(len(seg)) * dt
        ex = seg - np.median(g[i0 + int(40 / dt): i0 + int(58 / dt)])
        m = (tt >= 0.3) & (tt <= 8)
        rows.append({"transient": f"power step a_L->a_H ({p} dB level), CW carrier", "osc_period_us": dominant_period(ex[m], dt),
                     "osc_period_gain_us": np.nan, "envelope_1e_us": env_decay(ex[tt >= 0.3], tt[tt >= 0.3] - 0.3)})
    tab = pd.read_csv(HERE / "posthoc_isolated_transitions.csv")
    pm = np.load(R06 / "phase_mechanism" / "runs" / "QPSK_SLOW_Nd4001_dt1.npz")
    t = pm["t"] * 1e6
    g = np.abs(pm["z_full"]) / np.abs(pm["z_linear"])
    u = pm["z_full"] * np.conj(pm["z_linear"])
    dtp = t[1] - t[0]
    for r in tab.itertuples():
        post = (t >= r.t_boundary_us + 0.3) & (t < r.t_boundary_us + 3.7)
        pre = (t >= r.t_boundary_us - 0.6) & (t < r.t_boundary_us - 0.2)
        th = np.angle(u[post] * np.conj(np.mean(u[pre])))
        gg = g[post] - np.median(g[pre])
        kind = "pi/2" if abs(abs(r.dphi_rad) - np.pi / 2) < .1 else "pi"
        rows.append({"transient": f"phase step {kind} at t={r.t_boundary_us:.0f} us (constant a_H)", "osc_period_us": dominant_period(th, dtp),
                     "osc_period_gain_us": dominant_period(gg, dtp), "envelope_1e_us": env_decay(gg, t[post] - t[post][0])})
    d = pd.DataFrame(rows)
    d.to_csv(HERE / "posthoc_timescales.csv", index=False)
    return d


def period_vs_nd():
    z = np.load(R06 / "tau_atom_traces.npz")
    rows = []
    for exp in ("IF_step_P1dB_1.00to1.05E1", "IF_step_P1dB_1.05to1.00E1", "IF_plateau_0.5to1.5E1_on", "IF_plateau_0.5to1.5E1_off_recovery", "LO_step_+1pct"):
        for nd in (1501, 4001, 8001):
            t = z[f"{nd}__{exp}__t"] * 1e6
            e = z[f"{nd}__{exp}__e"]
            m = (t >= .3) & (t <= 8)
            rows.append({"experiment": exp, "Nd": nd, "osc_period_us": dominant_period(e[m], t[1] - t[0]),
                         "envelope_1e_us": env_decay(e[t >= .3], t[t >= .3] - .3)})
    d = pd.DataFrame(rows)
    d.to_csv(HERE / "posthoc_oscillation_vs_Nd.csv", index=False)
    rec = {nd: {"coupling_comb_us": 1 / (1 / 510e-9 * 7 * 137 / (nd - 2)) * 1e6, "two_photon_comb_us": 1 / ((1 / 510e-9 - 1 / 852e-9) * 7 * 137 / (nd - 2)) * 1e6}
           for nd in (1501, 4001, 8001)}
    return d, rec


if __name__ == "__main__":
    pd.set_option("display.width", 200)
    print(power_vs_phase().round(2).to_string())
    d, rec = period_vs_nd()
    print(d.pivot_table(index="experiment", columns="Nd", values="osc_period_us").round(3).to_string())
    print("velocity-grid recurrence periods an artifact would follow (us):", {k: {kk: round(vv, 2) for kk, vv in v.items()} for k, v in rec.items()})
