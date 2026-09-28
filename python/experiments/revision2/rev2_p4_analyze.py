"""P4 analysis: internal density-matrix traces around an isolated +pi/2 phase step (rules: 00_decision_rules.json, P4).

rho vector layout (upstream, column-major): k = i + 4 j for rho[i, j] (0-based), e.g. rho21 -> 1, rho43 -> 11.
IF components of each element x(t) on the 10 ns record, with the demodulator used for the output (pm_run.demod: the 200 ns
boxcar base is removed before mixing): c+-(t) = 200 ns boxcar of (x - base)(t) exp(-+j(2 pi f_IF t + phi(t))); for the probe
observable Im(rho21) also its real-signal IF component 2 boxcar((Im x - base) exp(-j(...))). Slow component: 200 ns boxcar of x.
Output: demodulated IF envelope of the continuous run (same demodulator as P1).
"""
from __future__ import annotations

import json

import numpy as np
import pandas as pd

from rev2_common import REV2, demod

F_IF = 5e6
ELEMENTS = {"rho21 (probe)": 1, "rho32 (coupling)": 6, "rho43 (RF)": 11, "rho31 (two-photon)": 2, "rho41 (three-photon)": 3,
            "rho42": 7, "rho22": 5, "rho33": 10, "rho44": 15}
PRE = (49.2e-6, 49.7e-6)
# Windows end at 54.75 us: within 0.2 us of the record end (55 us) the 200 ns base removal and demodulation boxcars run
# off the record (the demodulated magnitude jumps by >10x there), the same reason the X window excludes its final 0.5 us.
SEARCH = (49.8e-6, 54.75e-6)
CORR = (50.0e-6, 54.75e-6)


def boxcar(x, n=20):
    k = np.ones(n) / n
    if np.iscomplexobj(x):
        return np.convolve(x.real, k, "same") + 1j * np.convolve(x.imag, k, "same")
    return np.convolve(x, k, "same")


def track(t, m, m_out):
    pre = (t >= PRE[0]) & (t < PRE[1])
    srch = (t >= SEARCH[0]) & (t < SEARCH[1])
    cor = (t >= CORR[0]) & (t < CORR[1])
    mn = m / np.mean(m[pre])
    mo = m_out / np.mean(m_out[pre])
    i, io = np.argmin(np.where(srch, mn, np.inf)), np.argmin(np.where(srch, mo, np.inf))
    r = float(np.corrcoef(mn[cor], mo[cor])[0, 1])
    dip = float(1 - mn[i])
    return {"dip": dip, "t_min_us": t[i] * 1e6, "dt_min_vs_output_us": (t[i] - t[io]) * 1e6, "r_vs_output": r,
            "tracks": bool(dip >= .2 and abs(t[i] - t[io]) <= .3e-6 and r >= .8)}, mn, mo


def phase_extreme(t, dph):
    """Largest |deviation| of the output phase from the drive phase after the step (demodulated with the input phase)."""
    srch = (t >= SEARCH[0]) & (t < SEARCH[1])
    i = np.argmax(np.where(srch, np.abs(dph), -np.inf))
    return {"output_phase_max_abs_dev_rad": float(abs(dph[i])), "output_phase_dev_signed_rad": float(dph[i]), "output_phase_t_max_us": float(t[i] * 1e6)}


def analyze(name):
    d = np.load(REV2 / "p4" / f"{name}.npz")
    t_y, y, phi_y = d["t_y"], d["y_cont"], d["phi"]
    z_out = demod(y, phi_y, t_y, F_IF)
    t = d["t_rec"][:-1]                                   # uniform 10 ns grid (the last segment is shorter)
    idx = np.searchsorted(t_y, t - 1e-13)
    phi = phi_y[idx]
    zo = z_out[idx]
    m_out = np.abs(zo)
    ph_out = np.unwrap(np.angle(zo))
    rows, traces = [], {"t": t, "phi_in": phi, "out_mag": m_out, "out_phase": ph_out}
    pre = (t >= PRE[0]) & (t < PRE[1])
    for cls, rho in (("thermal", d["rho_bar"][:-1]), ("v0", d["rho_v0"][:-1])):
        for lab, k in ELEMENTS.items():
            x = rho[:, k]
            slow = boxcar(x)
            comps = [(+1, "+IF"), (-1, "-IF")] + ([(0, "Im part, IF")] if k == 1 else [])
            for sgn, cname in comps:
                if sgn == 0:
                    c = 2 * boxcar((x.imag - boxcar(x.imag)) * np.exp(-1j * (2 * np.pi * F_IF * t + phi)))
                else:
                    c = boxcar((x - slow) * np.exp(-sgn * 1j * (2 * np.pi * F_IF * t + phi)))
                res, mn, mo = track(t, np.abs(c), m_out)
                ph = np.unwrap(np.angle(c))
                rows.append({"case": name, "class": cls, "element": lab, "component": cname,
                             "pre_mag": float(np.mean(np.abs(c[pre]))), "phase_excursion_rad": float(np.max(np.abs(ph - np.mean(ph[pre])))),
                             **res})
                traces[f"{cls}|{lab}|{cname}"] = mn
                traces[f"{cls}|{lab}|{cname}|phase"] = ph - np.mean(ph[pre])
            sv = slow.real if k in (5, 10, 15) else np.abs(slow)
            rel = sv / np.mean(sv[pre])
            srch = (t >= SEARCH[0]) & (t < SEARCH[1])
            j = np.argmax(np.where(srch, np.abs(rel - 1), -np.inf))
            rows.append({"case": name, "class": cls, "element": lab, "component": "slow", "pre_mag": float(np.mean(sv[pre])),
                         "dip": float(1 - rel[j]), "t_min_us": t[j] * 1e6, "dt_min_vs_output_us": np.nan, "r_vs_output": np.nan, "tracks": False})
            traces[f"{cls}|{lab}|slow"] = rel
    traces["out_mag_norm"] = m_out / np.mean(m_out[pre])
    traces["out_phase_rel"] = ph_out - np.mean(ph_out[pre])
    np.savez_compressed(REV2 / "p4" / f"{name}_traces.npz", **traces)
    out_dip = float(1 - np.min(traces["out_mag_norm"][(t >= SEARCH[0]) & (t < SEARCH[1])]))
    return pd.DataFrame(rows), {"case": name, "chain_dev": float(d["chain_dev"]), "output_dip": out_dip,
                                "output_t_min_us": float(t[np.argmin(np.where((t >= SEARCH[0]) & (t < SEARCH[1]), m_out, np.inf))] * 1e6),
                                **phase_extreme(t, ph_out - np.mean(ph_out[pre]))}


def main():
    all_rows, meta = [], []
    for name in ("step_aH", "step_smallsignal"):
        df, m = analyze(name)
        all_rows.append(df)
        meta.append(m)
    df = pd.concat(all_rows)
    df.to_csv(REV2 / "40_p4_internal_tracking.csv", index=False)
    th = df[(df["class"] == "thermal") & (df.case == "step_aH")]
    internal = th[th.component != "Im part, IF"]
    verdict = {"meta": meta,
               "probe_observable_consistency": th[th.component == "Im part, IF"][["dip", "dt_min_vs_output_us", "r_vs_output"]].to_dict("records"),
               "tracking_components_aH_thermal": internal[internal.tracks][["element", "component", "dip", "dt_min_vs_output_us", "r_vs_output"]].to_dict("records"),
               "any_internal_IF_component_tracks": bool(internal.tracks.any())}
    (REV2 / "41_p4_verdict.json").write_text(json.dumps(verdict, indent=1, default=float))
    pd.set_option("display.width", 250)
    print(df[df["class"] == "thermal"].round(3).to_string())
    print(json.dumps(verdict, indent=1, default=float))


if __name__ == "__main__":
    main()
