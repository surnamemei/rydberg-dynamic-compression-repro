"""P2-P6 analysis (preregistered rules in results/revision/00_revision_preregistration.json)."""
from __future__ import annotations

import json

import numpy as np
import pandas as pd
from scipy import stats

from rev_common import REV, ROOT
import rev_phase as RP
import pm_run as P

FV = ROOT / "results" / "final_validation"
PM = ROOT / "results" / "stage06_dwell_physics" / "phase_mechanism"
H = RP.LEVEL_H
SEEDS = RP.SEEDS
CASES = list(RP.CASES)


def ci(v):
    v = np.asarray(v, float)
    n = len(v)
    h = stats.t.ppf(.975, n - 1) * v.std(ddof=1) / np.sqrt(n) if n > 1 else np.nan
    return float(v.mean()), float(v.mean() - h), float(v.mean() + h), n


def resolved(lo, hi):
    return bool(lo > 0 or hi < 0)


def path(sub, name):
    return REV / sub / f"{name}.npz"


# ---------------------------------------------------------------------------------------------------------- P2
def p2():
    rows, con = [], []
    X = {}
    for op in ("OP1", "OP2"):
        gcw = RP.g_final(path("p2", f"{op}_CW_H"))
        for case in CASES:
            for seed in SEEDS:
                x = 1 - RP.g_final(path("p2", f"{op}_{case}_s{seed}_H")) / gcw
                X[(op, case, seed)] = x
                rows.append({"op": op, "case": case, "seed": seed, "X": x, "g_CW": gcw})
    per = pd.DataFrame(rows)
    summ = []
    for (op, case), g in per.groupby(["op", "case"]):
        m, lo, hi, n = ci(g.X)
        summ.append({"op": op, "case": case, "X_mean": m, "ci_lo": lo, "ci_hi": hi, "n": n, "sd": float(g.X.std(ddof=1)),
                     "min": float(g.X.min()), "max": float(g.X.max())})
    summ = pd.DataFrame(summ)
    for op in ("OP1", "OP2"):
        for case in ("JUMP", "LONGRAMP", "SLOW", "FAST"):
            d = [X[(op, case, s)] - X[(op, "REF", s)] for s in SEEDS]
            m, lo, hi, n = ci(d)
            con.append({"contrast": f"{case}-REF", "op": op, "mean": m, "ci_lo": lo, "ci_hi": hi, "n": n, "resolved": resolved(lo, hi)})
    for case in CASES:
        d = [X[("OP1", case, s)] - X[("OP2", case, s)] for s in SEEDS]
        m, lo, hi, n = ci(d)
        con.append({"contrast": f"OP1-OP2 {case}", "op": "both", "mean": m, "ci_lo": lo, "ci_hi": hi, "n": n, "resolved": resolved(lo, hi)})
    con = pd.DataFrame(con)
    # determinism check against archived zero-drift runs (identical cases)
    arch = pd.read_csv(FV / "replication" / "03_replication_results.csv")
    chk = []
    pairs = {("OP2", "REF", 20260701): "LO2_QPSK_REF_zd", ("OP2", "REF", 20260702): "LO2_QPSK_REF_zd_s2", ("OP2", "JUMP", 20260701): "LO2_QPSK_JUMP_zd",
             ("OP2", "LONGRAMP", 20260701): "LO2_QPSK_LONGRAMP_zd", ("OP2", "SLOW", 20260701): "LO2_QPSK_SLOW_zd", ("OP2", "FAST", 20260701): "LO2_QPSK_FAST_zd",
             ("OP1", "REF", 20260701): "LO1_QPSK_REF_zd", ("OP1", "SLOW", 20260701): "LO1_QPSK_SLOW_zd", ("OP1", "FAST", 20260701): "LO1_QPSK_FAST_zd"}
    for (op, case, seed), name in pairs.items():
        a = np.load(FV / "replication" / "runs" / f"{name}_Nd4001.npz")
        b = np.load(path("p2", f"{op}_{case}_s{seed}_H"))
        r = arch[arch.case == name]
        chk.append({"case": name, "z_full_identical": bool(np.array_equal(a["z_full"], b["z_full"])),
                    "max_abs_diff_z_full": float(np.max(np.abs(a["z_full"] - b["z_full"]))),
                    "X_archived": float(r.X.iloc[0]) if "X" in r and len(r) else np.nan, "X_new": X[(op, case, seed)]})
    return per, summ, con, pd.DataFrame(chk), X


# ---------------------------------------------------------------------------------------------------------- P3
def p3():
    rows, summ = [], []
    for op in ("OP1", "OP2"):
        g_small = None
        for db in RP.LEVELS_DB:
            gcw = RP.g_final(path("p3", f"{op}_CW_p{db:+d}"))
            xs = [1 - RP.g_final(path("p3", f"{op}_REF_s{s}_p{db:+d}")) / gcw for s in SEEDS[:3]]
            for s, x in zip(SEEDS[:3], xs):
                rows.append({"op": op, "level_dB": db, "seed": s, "X": x, "g_CW": gcw})
            m, lo, hi, n = ci(xs)
            summ.append({"op": op, "level_dB": db, "g_CW_rel_small_signal": gcw, "X_mean": m, "ci_lo": lo, "ci_hi": hi, "n": n,
                         "materially_nonzero": bool(m >= .05 and lo > 0)})
    summ = pd.DataFrame(summ)
    onset = {}
    for op, g in summ.groupby("op"):
        g = g.sort_values("level_dB")
        flags = g.materially_nonzero.values
        lv = None
        for i in range(len(flags)):
            if flags[i:].all():
                lv = int(g.level_dB.values[i])
                break
        onset[op] = lv
    return pd.DataFrame(rows), summ, onset


# ---------------------------------------------------------------------------------------------------------- P4
def offset_metrics(op, name, cw_name):
    z = np.load(path(*name))
    zc = np.load(path(*cw_name))
    t = z["t"]
    m = (t >= RP.T0 + 60e-6 - 20.5e-6) & (t < RP.T0 + 60e-6 - .5e-6)
    g = np.median(np.abs(z["z_full"][m]) / np.abs(z["z_linear"][m]))
    gcw = np.median(np.abs(zc["z_full"][m]) / np.abs(zc["z_linear"][m]))
    return {"R": float(g / gcw), "X_offset": float(1 - g / gcw),
            "large_signal_out_vs_CW": float(np.median(np.abs(z["z_full"][m])) / np.median(np.abs(zc["z_full"][m]))),
            "small_signal_out_vs_CW": float(np.median(np.abs(z["z_linear"][m])) / np.median(np.abs(zc["z_linear"][m]))),
            "C": complex(np.median(z["z_full"][m].real) + 1j * np.median(z["z_full"][m].imag)),
            "C_cw": complex(np.median(zc["z_full"][m].real) + 1j * np.median(zc["z_full"][m].imag))}


def p4(p2X):
    sweep = []
    for op in ("OP1", "OP2"):
        sweep.append({"op": op, "offset_MHz": 0.0, "R": 1.0, "X_offset": 0.0, "large_signal_out_vs_CW": 1.0, "small_signal_out_vs_CW": 1.0})
        for mhz in RP.OFFSETS_MHZ:
            r = offset_metrics(op, ("p45", f"{op}_OFF{mhz:+.3f}_H"), ("p2", f"{op}_CW_H"))
            sweep.append({"op": op, "offset_MHz": mhz, **{k: v for k, v in r.items() if k not in ("C", "C_cw")}})
    sweep = pd.DataFrame(sweep).sort_values(["op", "offset_MHz"])
    five = [0.125, RP.OFFSETS_MHZ[-5], -RP.OFFSETS_MHZ[-5], RP.OFFSETS_MHZ[-2], -RP.OFFSETS_MHZ[-2]]
    op3 = []
    for mhz in five:
        r = offset_metrics("OP3", ("p4", f"OP3_OFF{mhz:+.3f}_H"), ("p4", "OP3_CW_H"))
        op3.append({"op": "OP3", "offset_MHz": mhz, **{k: v for k, v in r.items() if k not in ("C", "C_cw")}})
    op3 = pd.DataFrame(op3)
    S5 = {}
    for op in ("OP1", "OP2"):
        g = sweep[(sweep.op == op) & (np.isclose(sweep.offset_MHz.values[:, None], np.array(five)[None, :]).any(axis=1))]
        S5[op] = float(g.R.max() - 1)
    S5["OP3"] = float(op3.R.max() - 1)
    gcw3 = RP.g_final(path("p4", "OP3_CW_H"))
    x3 = [1 - RP.g_final(path("p4", f"OP3_REF_s{s}_H")) / gcw3 for s in SEEDS[:3]]
    xref = {"OP3": ci(x3)}
    for op in ("OP1", "OP2"):
        xref[op] = ci([p2X[(op, "REF", s)] for s in SEEDS[:3]])
    e1 = {op: RP.e1db(op) for op in ("OP1", "OP2", "OP3")}
    lo_ = {op: RP.OPS[op]["lo"] for op in e1}
    between = lambda v3, v1, v2: min(v1, v2) - .05 <= v3 <= max(v1, v2) + .05  # noqa: E731
    mono = between(S5["OP3"], S5["OP1"], S5["OP2"]) and between(xref["OP3"][0], xref["OP1"][0], xref["OP2"][0])
    iso = S5["OP3"] < .5 * S5["OP1"] and abs(S5["OP3"] - S5["OP2"]) <= .10
    cls = "MONOTONIC_IN_LO" if mono else ("OP1_ISOLATED" if iso else "OTHER")
    ops = pd.DataFrame([{"op": op, "A_LO_Vpm": lo_[op], "E1dB_Vpm": e1[op], "S5": S5[op], "X_REF_mean_3seeds": xref[op][0],
                         "X_REF_ci_lo": xref[op][1], "X_REF_ci_hi": xref[op][2]} for op in ("OP2", "OP3", "OP1")])
    return sweep, op3, ops, cls


# ---------------------------------------------------------------------------------------------------------- P5
def calibration(op, level_tag, cw_name):
    """C(df) at one level: complex steady-state IF output on the offset grid (df = 0 from the CW run)."""
    f = [0.0]
    C = [None]
    for mhz in RP.OFFSETS_MHZ:
        r = offset_metrics(op, ("p45", f"{op}_OFF{mhz:+.3f}_{level_tag}"), cw_name)
        f.append(mhz * 1e6)
        C.append(r["C"])
        C[0] = r["C_cw"]
    o = np.argsort(f)
    return np.array(f)[o], np.array(C)[o]


def surrogate_g(run_path, fgrid, Cgrid):
    """Memoryless g(a, f_i) surrogate on the exact phase trajectory of a stored run; returns g_sur(t) (20 ns grid)."""
    z = np.load(run_path)
    spec = json.loads(str(z["spec"]))
    t_mod = spec.get("t_mod", 60e-6)
    n = int(round((RP.T0 + t_mod) / RP.DT))
    t = np.arange(n) * RP.DT
    phi = RP.phase(spec, t)
    fi = np.r_[0.0, np.diff(phi) / (2 * np.pi * RP.DT)]
    fi = np.clip(fi, fgrid[0], fgrid[-1])
    Cs = np.interp(fi, fgrid, Cgrid.real) + 1j * np.interp(fi, fgrid, Cgrid.imag)
    y = np.real(Cs * np.exp(1j * (2 * np.pi * RP.IF * t + phi)))
    zs = P.demod(y, phi, t)[::20]
    return np.abs(zs) / np.abs(z["z_linear"]), z["t"]


def g_sur_final(run_path, fgrid, Cgrid, t_mod=60e-6):
    g, t = surrogate_g(run_path, fgrid, Cgrid)
    m = (t >= RP.T0 + t_mod - 20.5e-6) & (t < RP.T0 + t_mod - .5e-6)
    return float(np.median(g[m]))


def p5(p2X, p3_rows):
    rows = []
    Xs = {}
    for op in ("OP1", "OP2"):
        fgrid, Cgrid = calibration(op, "H", ("p2", f"{op}_CW_H"))
        gcw = g_sur_final(path("p2", f"{op}_CW_H"), fgrid, Cgrid)
        for case in CASES:
            for seed in SEEDS:
                xs = 1 - g_sur_final(path("p2", f"{op}_{case}_s{seed}_H"), fgrid, Cgrid) / gcw
                Xs[(op, case, seed)] = xs
                rows.append({"op": op, "case": case, "seed": seed, "X_full": p2X[(op, case, seed)], "X_sur": xs})
    per = pd.DataFrame(rows)
    t1 = []
    for op in ("OP1", "OP2"):
        xf = np.mean([p2X[(op, "REF", s)] for s in SEEDS])
        xs = np.mean([Xs[(op, "REF", s)] for s in SEEDS])
        t1.append({"op": op, "X_full_REF": xf, "X_sur_REF": xs, "pass": bool(abs(xs - xf) <= .25 * xf)})
    t2 = []
    for op in ("OP1", "OP2"):
        for case in ("SLOW", "FAST", "JUMP", "LONGRAMP"):
            df_ = [p2X[(op, case, s)] - p2X[(op, "REF", s)] for s in SEEDS]
            ds_ = [Xs[(op, case, s)] - Xs[(op, "REF", s)] for s in SEEDS]
            m, lo, hi, _ = ci(df_)
            eligible = resolved(lo, hi)
            ds = float(np.mean(ds_))
            t2.append({"op": op, "contrast": f"{case}-REF", "d_full": m, "d_full_ci_lo": lo, "d_full_ci_hi": hi, "d_sur": ds, "eligible": eligible,
                       "pass": bool(eligible and np.sign(ds) == np.sign(m) and abs(ds - m) <= .25 * abs(m)) if eligible else None})
    t3 = []
    for op in ("OP1", "OP2"):
        for db in RP.LEVELS_DB:
            fgrid, Cgrid = calibration(op, f"p{db:+d}", ("p3", f"{op}_CW_p{db:+d}"))
            gcw = g_sur_final(path("p3", f"{op}_CW_p{db:+d}"), fgrid, Cgrid)
            xs = [1 - g_sur_final(path("p3", f"{op}_REF_s{s}_p{db:+d}"), fgrid, Cgrid) / gcw for s in SEEDS[:3]]
            xf = p3_rows[(p3_rows.op == op) & (p3_rows.level_dB == db)].X.mean()
            t3.append({"op": op, "level_dB": db, "X_full_mean": float(xf), "X_sur_mean": float(np.mean(xs))})
    t1, t2, t3 = pd.DataFrame(t1), pd.DataFrame(t2), pd.DataFrame(t3)
    el = t2[t2.eligible]
    frac = float(el["pass"].mean()) if len(el) else np.nan
    verdict = "REPRODUCES" if (t1["pass"].all() and len(el) and frac >= .75) else "FAILS"
    return per, t1, t2, t3, verdict


# ---------------------------------------------------------------------------------------------------------- P6
def p6(p2X):
    out = {}
    g_cw64 = RP.g_final(path("p6", "OP1_CW_H_fp64"))
    x64 = 1 - RP.g_final(path("p6", "OP1_REF_s20260701_H_fp64")) / g_cw64
    x32 = p2X[("OP1", "REF", 20260701)]
    out["phase_fp64"] = {"X_fp64": x64, "X_fp32_gpu": x32, "abs_dX": abs(x64 - x32), "pass_tol_0.01": bool(abs(x64 - x32) <= .01),
                         "g_CW_fp64": g_cw64, "g_CW_fp32": RP.g_final(path("p2", "OP1_CW_H"))}
    rows = pd.DataFrame(json.loads((REV / "p6" / "long_dwell_xi32_fp64_rows.json").read_text())).set_index("model")
    arch = pd.read_csv(ROOT / "results" / "stage06_dwell_physics" / "rows_fu_long_cw_dwell.csv")
    arch = arch[arch.job_id == "dwell_r0_xi_32_p+3_Nd4001_dt1"].set_index("model")
    gap = lambda r: (r.loc["atomic", "D_ref_linear"] - r.loc["static_mh", "D_ref_linear"]) / r.loc["static_mh", "D_ref_linear"]  # noqa: E731
    out["long_dwell_fp64"] = {"D_ref_full_fp64": float(rows.loc["atomic", "D_ref_linear"]), "D_ref_full_fp32": float(arch.loc["atomic", "D_ref_linear"]),
                              "D_ref_static": float(rows.loc["static_mh", "D_ref_linear"]), "rel_gap_fp64": float(gap(rows)), "rel_gap_fp32": float(gap(arch)),
                              "gap_change_pp": float(100 * abs(gap(rows) - gap(arch))), "pass_tol_1pp": bool(100 * abs(gap(rows) - gap(arch)) <= 1.0),
                              "D_BLA_full_fp64": float(rows.loc["atomic", "D_BLA"]), "D_BLA_full_fp32": float(arch.loc["atomic", "D_BLA"])}
    gcw = RP.g_final(path("p2", "OP1_CW_H"))
    ext = path("p6", "OP1_REF_s20260701_H_tmod200")
    win = []
    for k in range(7):
        a = 50 + 40 + 20 * k
        b = a + 20 - (0.5 if k == 6 else 0)
        win.append({"window_after_onset_us": f"{40 + 20 * k}-{60 + 20 * k}", "X": 1 - RP.g_window(ext, a, b) / gcw})
    w = pd.DataFrame(win)
    out["extended_segment"] = {"windows": win, "max_abs_dev_from_first": float(np.max(np.abs(w.X - w.X.iloc[0]))),
                               "pass_tol_0.03": bool(np.max(np.abs(w.X - w.X.iloc[0])) <= .03)}
    timing = {"gpu_phase_run_110us_s": float(np.load(path("p2", "OP1_REF_s20260701_H"))["wall_s"]),
              "gpu_phase_run_250us_s": float(np.load(ext)["wall_s"]),
              "cpu_fp64_phase_run_110us_s": float(np.load(path("p6", "OP1_REF_s20260701_H_fp64"))["wall_s"])}
    cpu = json.loads((REV / "p6" / "cpu_wall_times_s.json").read_text())
    timing["cpu_fp64_long_dwell_2.44ms_s"] = cpu.get("long_dwell_xi32_fp64")
    timing["gpu_long_dwell_2.44ms_archived_runtime_s"] = float(arch.loc["atomic", "runtime_atomic_s"])
    p1rows = pd.read_csv(REV / "p1" / "rows_p1.csv")
    timing["gpu_dwell_job_840us_median_s"] = float(p1rows[p1rows.model == "atomic"].runtime_atomic_s.median())
    out["timing"] = timing
    return out


def p5_posthoc_mean_based():
    """POST HOC sensitivity check: X from the window MEAN (instead of the preregistered median) of g, full vs g(a,f)."""
    out = []
    for op in ("OP1", "OP2"):
        fgrid, Cgrid = calibration(op, "H", ("p2", f"{op}_CW_H"))

        def xmean(g, t):
            m = (t >= RP.T0 + 60e-6 - 20.5e-6) & (t < RP.T0 + 60e-6 - .5e-6)
            return float(np.mean(g[m]))
        gs, t = surrogate_g(path("p2", f"{op}_CW_H"), fgrid, Cgrid)
        gs_cw = xmean(gs, t)
        z = np.load(path("p2", f"{op}_CW_H"))
        gf_cw = xmean(np.abs(z["z_full"]) / np.abs(z["z_linear"]), z["t"])
        for case in CASES:
            xs, xf = [], []
            for seed in SEEDS:
                pth = path("p2", f"{op}_{case}_s{seed}_H")
                g, t = surrogate_g(pth, fgrid, Cgrid)
                xs.append(1 - xmean(g, t) / gs_cw)
                z = np.load(pth)
                xf.append(1 - xmean(np.abs(z["z_full"]) / np.abs(z["z_linear"]), z["t"]) / gf_cw)
            out.append({"op": op, "case": case, "X_full_meanbased": float(np.mean(xf)), "X_sur_meanbased": float(np.mean(xs))})
    return pd.DataFrame(out)


def calibration_table():
    rows = []
    for op in ("OP1", "OP2"):
        for tag, cw in [("H", ("p2", f"{op}_CW_H"))] + [(f"p{db:+d}", ("p3", f"{op}_CW_p{db:+d}")) for db in RP.LEVELS_DB]:
            f, C = calibration(op, tag, cw)
            for fi, ci_ in zip(f, C):
                rows.append({"op": op, "level": tag, "offset_Hz": fi, "C_re": ci_.real, "C_im": ci_.imag, "C_abs": abs(ci_)})
    return pd.DataFrame(rows)


def main():
    per2, summ2, con2, chk2, p2X = p2()
    per3, summ3, onset = p3()
    sweep, op3, ops, cls = p4(p2X)
    per5, t1, t2, t3, v5 = p5(p2X, per3)
    o6 = p6(p2X)
    for name, df in (("10_p2_per_sequence", per2), ("11_p2_summary", summ2), ("12_p2_contrasts", con2), ("13_p2_determinism_check", chk2),
                     ("20_p3_per_sequence", per3), ("21_p3_summary", summ3), ("30_p4_offset_sweep", sweep), ("31_p4_op3_offsets", op3),
                     ("32_p4_operating_points", ops), ("40_p5_per_sequence", per5), ("41_p5_T1", t1), ("42_p5_T2", t2), ("43_p5_T3", t3)):
        df.to_csv(REV / f"{name}.csv", index=False)
    p5_posthoc_mean_based().to_csv(REV / "45_p5_posthoc_mean_based.csv", index=False)
    calibration_table().to_csv(REV / "44_p5_calibration_C.csv", index=False)
    pd.DataFrame(o6["extended_segment"]["windows"]).to_csv(REV / "51_p6_extended_segment_windows.csv", index=False)
    summary = {"P3_onset_dB_re_P1dB": onset, "P4_classification": cls, "P5_verdict": v5, "P6": o6}
    (REV / "50_p2_p6_summary.json").write_text(json.dumps(summary, indent=1, default=float))
    pd.set_option("display.width", 200)
    print(summ2.to_string(), con2.to_string(), chk2.to_string(), summ3.to_string(), sep="\n\n")
    print(sweep.to_string(), op3.to_string(), ops.to_string(), cls, sep="\n\n")
    print(t1.to_string(), t2.to_string(), t3.to_string(), v5, sep="\n\n")
    print(json.dumps(o6, indent=1, default=float))


if __name__ == "__main__":
    main()
