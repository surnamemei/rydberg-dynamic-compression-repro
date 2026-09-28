"""P1 analysis: IF dependence (rules: results/revision2/00_decision_rules.json, P1; 01_rule_clarifications.json)."""
from __future__ import annotations

import json

import numpy as np
import pandas as pd
from scipy import stats

from rev2_common import FIVE_OFFSETS_MHZ, REV2, ROOT, SEEDS3, R, bb, fitted_gain, static_zone
import rev_phase as RP

IFS = (5, 10, 15)
REV = ROOT / "results" / "revision"
JOBS = ROOT / "results" / "stage06_dwell_physics" / "jobs" / "main"


def ci(v):
    v = np.asarray(v, float)
    n = len(v)
    h = stats.t.ppf(.975, n - 1) * v.std(ddof=1) / np.sqrt(n) if n > 1 else np.nan
    return float(v.mean()), float(v.mean() - h), float(v.mean() + h), n


def g(path):
    return RP.g_final(path)


def validation():
    """5 MHz: IF-parametrized code vs archived runs (same GPU kernel, same inputs)."""
    rows = []
    d5 = REV2 / "p1" / "IF5"
    pairs = [("CW_H", REV / "p2" / "OP1_CW_H.npz")] + [(f"REF_s{s}_H", REV / "p2" / f"OP1_REF_s{s}_H.npz") for s in SEEDS3]
    pairs += [(f"OFF{m:+.3f}_H", REV / "p45" / f"OP1_OFF{m:+.3f}_H.npz") for m in FIVE_OFFSETS_MHZ]
    for name, arch in pairs:
        a, b = np.load(arch), np.load(d5 / "runs" / f"{name}.npz")
        rows.append({"item": f"phase run {name}", "z_full_identical": bool(np.array_equal(a["z_full"], b["z_full"])),
                     "max_rel_diff_z_full": float(np.max(np.abs(a["z_full"] - b["z_full"])) / np.max(np.abs(a["z_full"]))),
                     "max_rel_diff_z_linear": float(np.max(np.abs(a["z_linear"] - b["z_linear"])) / np.max(np.abs(a["z_linear"]))),
                     "archived": g(arch), "new": g(d5 / "runs" / f"{name}.npz")})
    for r in range(4):
        for v in ("xi_0.05", "xi_4"):
            arch = {x["model"]: x for x in json.loads((JOBS / f"dwell_r{r}_{v}_p+3_Nd4001_dt1.json").read_text())["rows"]}
            new = np.load(d5 / "dwell" / f"dwell_r{r}_{v}_p+3.npz")
            for m_new, m_arch in (("full", "atomic"), ("linear", "linear"), ("static", "static")):
                a, b = arch[m_arch]["gain_abs"], abs(complex(new[f"g_{m_new}"]))
                rows.append({"item": f"dwell r{r} {v} {m_new} gain_abs", "archived": a, "new": b, "rel_diff": abs(a - b) / a})
    df = pd.DataFrame(rows)
    df.to_csv(REV2 / "10_p1_validation_5MHz.csv", index=False)
    return df


def load_table(f_mhz):
    """5 MHz: archived grid; 10 and 15 MHz: grid extended above 0.70 V/m (deviation D1, 02_deviations.json)."""
    p = REV2 / "tables" / (f"OP1_IF{f_mhz}.npz" if f_mhz == 5 else f"OP1_IF{f_mhz}_ext.npz")
    return dict(np.load(p))


_REAL = {}


def realization(r):
    if r not in _REAL:
        _REAL.clear()
        _REAL[r] = R.realization("dwell", r)
    return _REAL[r]


def static_gain(f_mhz, tab, r, v):
    """Fitted gain of the CW-matched static surrogate (fundamental zone; models06.Controls.static) on the dwell job."""
    amps, meta, car, labels, noise = realization(r)
    env = tab["e1db"] * 10 ** (3 / 20) * amps[v] * car
    env_cp = np.r_[env[-R.CP:], env]
    a = np.abs(env_cp)
    u = np.ones(len(env_cp), complex)
    np.divide(env_cp, a, out=u, where=a > 0)
    c1 = np.interp(a, tab["amps"], tab["harm"].real) + 1j * np.interp(a, tab["amps"], tab["harm"].imag)
    y = tab["means"][0] + np.real(c1 * u * np.exp(2j * np.pi * f_mhz * 1e6 * np.arange(len(env_cp)) / 1e9))
    x = R.sosfilt(R.SOS, env_cp)[::R.DEC]
    w0, w1 = (R.CP + R.SKIP) // R.DEC, len(x) - 2000 // R.DEC
    return fitted_gain(bb(y, tab["means"][0], f_mhz * 1e6), x[w0:w1], w0, w1), float(a.max()), float(tab["amps"][-1])


def leakage(f_mhz, tab):
    """In-band power (receiver chain at the IF) of the memoryless CW-matched zones 0, 2, 3 relative to zone 1."""
    rows = []
    for r in range(4):
        amps, meta, car, labels, noise = R.realization("dwell", r)
        for v in ("xi_0.05", "xi_4"):
            env = tab["e1db"] * 10 ** (3 / 20) * amps[v] * car
            env_cp = np.r_[env[-R.CP:], env]
            n = len(env_cp) // R.DEC
            w0, w1 = (R.CP + R.SKIP) // R.DEC, n - 2000 // R.DEC
            p = {k: float(np.mean(np.abs(bb(static_zone(env_cp, tab, k, f_mhz * 1e6), 0.0, f_mhz * 1e6)[w0:w1]) ** 2)) for k in range(4)}
            rows.append({"IF_MHz": f_mhz, "realization": r, "variant": v, **{f"zone{k}_rel_zone1_dB": 10 * np.log10(p[k] / p[1]) for k in (0, 2, 3)},
                         "all_other_zones_rel_zone1": (p[0] + p[2] + p[3]) / p[1]})
    return rows


def main():
    val = validation()
    rows, dwell_rows, leak_rows, p2x = [], [], [], []
    for f in IFS:
        d = REV2 / "p1" / f"IF{f}"
        tab = load_table(f)
        gcw = g(d / "runs" / "CW_H.npz")
        xs = [1 - g(d / "runs" / f"REF_s{s}_H.npz") / gcw for s in SEEDS3]
        mx, lo, hi, n = ci(xs)
        R_off = {m: g(d / "runs" / f"OFF{m:+.3f}_H.npz") / gcw for m in FIVE_OFFSETS_MHZ}
        # dwell: fitted gain relative to the linear model, paired change xi 4 - xi 0.05
        gr = {}
        for r in range(4):
            for v in ("xi_0.05", "xi_4"):
                z = np.load(d / "dwell" / f"dwell_r{r}_{v}_p+3.npz")
                assert abs(float(z["e1db"]) - float(tab["e1db"])) <= 1e-12 * float(tab["e1db"])
                g_st, a_peak, a_grid = static_gain(f, tab, r, v)
                if f == 5:                                   # same table and code path as the run itself
                    assert abs(g_st - complex(z["g_static"])) <= 1e-12 * abs(g_st)
                gr[("full", r, v)] = abs(complex(z["g_full"])) / abs(complex(z["g_linear"]))
                gr[("static", r, v)] = abs(g_st) / abs(complex(z["g_linear"]))
                dwell_rows.append({"IF_MHz": f, "realization": r, "variant": v, "gain_rel_full": gr[("full", r, v)], "gain_rel_static": gr[("static", r, v)],
                                   "peak_field_Vpm": a_peak, "cw_grid_max_Vpm": a_grid})
        dfull = [gr[("full", r, "xi_4")] - gr[("full", r, "xi_0.05")] for r in range(4)]
        dst = [gr[("static", r, "xi_4")] - gr[("static", r, "xi_0.05")] for r in range(4)]
        df_m, df_lo, df_hi, _ = ci(dfull)
        ds_m, ds_lo, ds_hi, _ = ci(dst)
        cw_gain_db = tab["gain_db"]
        a_h = RP.LEVEL_H * tab["e1db"]
        rows.append({"IF_MHz": f, "E1dB_Vpm": float(tab["e1db"]), "small_signal_slope": float(tab["slope"]), "cw_settle_rel": float(tab["cw_settle_rel"]),
                     "CW_gain_at_aH_rel_small_signal": gcw,
                     "CW_table_gain_at_aH_dB": float(np.interp(a_h, tab["amps"][1:], cw_gain_db[1:])),
                     "X_REF_mean": mx, "X_REF_lo": lo, "X_REF_hi": hi, "X_REF_n": n, **{f"X_REF_s{s}": x for s, x in zip(SEEDS3, xs)},
                     "phase_effect_present": bool(mx >= .05 and lo > 0),
                     "R_+0.125": R_off[0.125], "S5": max(R_off.values()) - 1, **{f"R_{m:+.3f}": v for m, v in R_off.items()},
                     "gain_rel_full_xi0.05": float(np.mean([gr[("full", r, "xi_0.05")] for r in range(4)])),
                     "gain_rel_full_xi4": float(np.mean([gr[("full", r, "xi_4")] for r in range(4)])),
                     "d_gain_full": df_m, "d_gain_full_lo": df_lo, "d_gain_full_hi": df_hi,
                     "d_gain_static": ds_m, "d_gain_static_lo": ds_lo, "d_gain_static_hi": ds_hi,
                     "dwell_effect_present": bool(df_m < 0 and df_hi < 0 and abs(df_m) >= 3 * abs(ds_m))})
        leak_rows += leakage(f, tab)
        # coherent / magnitude decomposition of the new reference runs (P2 definitions)
        import rev2_p2 as P2
        c = P2.metrics(d / "runs" / "CW_H.npz")
        for s in SEEDS3:
            q = P2.metrics(d / "runs" / f"REF_s{s}_H.npz")
            X = {k: 1 - q[f"g_{k}"] / c[f"g_{k}"] for k in ("coh", "mag", "rms", "med")}
            p2x.append({"IF_MHz": f, "seed": s, **{f"X_{k}": v for k, v in X.items()},
                        "share_phase": np.log(q["kappa"] / c["kappa"]) / np.log(1 - X["coh"]), "phase_circ_sd_rad": q["phase_circ_sd_rad"]})
    summ = pd.DataFrame(rows)
    b5 = summ.set_index("IF_MHz")
    ratio_x = b5.loc[15, "X_REF_mean"] / b5.loc[5, "X_REF_mean"]
    ratio_d = b5.loc[15, "d_gain_full"] / b5.loc[5, "d_gain_full"]
    both = bool(summ.phase_effect_present.all() and summ.dwell_effect_present.all())
    absent_high = bool(not (b5.loc[[10, 15], "phase_effect_present"].all() and b5.loc[[10, 15], "dwell_effect_present"].all()))
    if absent_high:
        cls = "DOES_NOT_SURVIVE"
    elif both and ratio_x >= .5 and ratio_d >= .5:
        cls = "SURVIVES_IF_CHANGE"
    elif both:
        cls = "WEAKENS"
    else:
        cls = "DOES_NOT_SURVIVE"
    summ.to_csv(REV2 / "11_p1_if_summary.csv", index=False)
    pd.DataFrame(dwell_rows).to_csv(REV2 / "12_p1_dwell_per_realization.csv", index=False)
    leak = pd.DataFrame(leak_rows)
    leak.to_csv(REV2 / "13_p1_zone_leakage.csv", index=False)
    pd.DataFrame(p2x).to_csv(REV2 / "14_p1_ref_decomposition.csv", index=False)
    verdict = {"classification": cls, "X_REF_ratio_15_over_5": ratio_x, "d_gain_ratio_15_over_5": ratio_d,
               "validation_5MHz_all_phase_runs_identical": bool(val[val.item.str.startswith("phase")].z_full_identical.all()),
               "validation_5MHz_max_dwell_gain_rel_diff": float(val[val.item.str.startswith("dwell")].rel_diff.max())}
    (REV2 / "15_p1_verdict.json").write_text(json.dumps(verdict, indent=1, default=float))
    pd.set_option("display.width", 250)
    print(val.to_string())
    print(summ.T.to_string())
    print(leak.groupby(["IF_MHz", "variant"]).mean(numeric_only=True).drop(columns="realization").round(2).to_string())
    print(pd.DataFrame(p2x).groupby("IF_MHz").mean().round(4).to_string())
    print(json.dumps(verdict, indent=1, default=float))


if __name__ == "__main__":
    main()
