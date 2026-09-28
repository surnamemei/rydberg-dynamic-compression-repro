"""Audit of the final figure set (no plotting, no simulation).

1. files: PDF/SVG/PNG exist; PDF page size (pdfinfo), all fonts embedded (pdffonts), PNG pixel size/dpi (PIL);
   SVG contains no <text> elements (text is converted to paths) and no file paths.
2. wording: pdftotext of every PDF is scanned for forbidden or internal strings (e.g. "out of lock", "phase lock",
   internal case/model/stage labels, filenames) and for ASCII hyphen-minus before digits.
3. provenance: the figure scripts import no simulation module.
4. numbers: plotted values recorded by each script (figures/audit/<fig>_plotted_values.csv) are compared with
   manuscript/numbers_ledger.csv (built independently from the same result files by build_numbers_ledger.py).
Writes figures/audit/figure_audit_results.json and prints a summary; exit code 1 if any check fails.
"""
from __future__ import annotations

import ast
import csv
import json
import re
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from PIL import Image

HERE = Path(__file__).resolve().parent
FIG = HERE.parent
MS = FIG.parent
AUD = FIG / "audit"
MAIN = [f"fig{i}" for i in range(1, 11)]
SUPP = [f"figS{i}" for i in range(1, 11)]

L = {r["id"]: r for r in csv.DictReader(open(MS / "numbers_ledger.csv"))}


def lv(key, idx=None):
    s = L[key]["value"]
    s = re.sub(r"np\.float64\(([^)]*)\)", r"\1", s)
    try:
        v = float(s)
    except ValueError:
        v = ast.literal_eval(s)
    return v if idx is None else v[idx]


def pv(fig):
    return pd.read_csv(AUD / f"{fig}_plotted_values.csv", dtype={"x": str, "panel": str}, keep_default_na=False)


def sel(df, panel, series, x=None):
    d = df[(df.panel == panel) & (df.series == series)]
    if x is not None:
        try:
            xf = float(x)
            d = d[[_isnum(v) and abs(float(v) - xf) < 1e-9 for v in d.x]]
        except (TypeError, ValueError):
            d = d[d.x == str(x)]
    return d


def _isnum(v):
    try:
        float(v)
        return True
    except ValueError:
        return False


RESULTS = {}


def check(fig, name, got, want, tol):
    ok = bool(np.all(np.abs(np.asarray(got, float) - np.asarray(want, float)) <= tol))
    RESULTS.setdefault(fig, {"numeric": [], "files": {}, "text": {}})["numeric"].append(
        {"check": name, "plotted": np.asarray(got, float).round(10).tolist(), "ledger_or_source": np.asarray(want, float).round(10).tolist(), "tol": tol, "pass": ok})
    return ok


def y1(df, panel, series, x=None):
    d = sel(df, panel, series, x)
    assert len(d) == 1, (panel, series, x, len(d))
    return float(d.y.iloc[0])


def ys(df, panel, series):
    return sel(df, panel, series).y.astype(float).values


EX = 1e-9
FMT = {"CE": "CE", "QPSK": "QPSK", "16-QAM": "16QAM", "OFDM": "OFDM"}
PS = {-6: "-6", 0: "+0", 3: "+3", 6: "+6"}

# ------------------------------------------------------------------ numbers
d = pv("fig1")
for n in (1501, 4001, 8001):
    check("fig1", f"(b) E1dB Nd={n}", y1(d, "b", f"E1dB Nd={n}"), lv(f"N_E1_{n}"), EX)
check("fig1", "(c) |c0/c1| nearest E1dB", y1(d, "c", "|c0/c1| nearest E1dB"), lv("N_c0c1_nearE1"), EX)
check("fig1", "(c) max |c0/c1| 0.75-3.2 E1dB", y1(d, "c", "max |c0/c1| 0.75-3.2 E1dB"), lv("N_c0c1_max"), EX)

d = pv("fig2")
for n in (1501, 4001, 8001):
    check("fig2", f"(b) E1dB Nd={n}", y1(d, "b", f"E1dB Nd={n}", n), lv(f"N_E1_{n}"), EX)
nm = d[(d.panel == "a") & d.series.str.startswith("NMSE")]
for n in (1501, 3001, 4001):
    check("fig2", f"(a) max NMSE at Nd={n}", nm[[float(v) == n for v in nm.x]].y.max(), lv(f"N_nd{n}_NMSE"), EX)

d = pv("fig3")
for f, lf in FMT.items():
    for p in (0, 3, 6):
        check("fig3", f"(a) loss {f} {p:+d} dB", y1(d, "a", f"loss {f}", p), lv(f"R_loss_{lf}_{PS[p]}"), EX)
        for c in ("static_mh", "static_lti_mh"):
            check("fig3", f"(b) excess {c} {f} {p:+d} dB", y1(d, "b", f"excess {c} {f}", p), lv(f"R_exc_{lf}_{PS[p]}_{c}"), EX)

d = pv("fig4")
for p in (-6, 0, 3, 6):
    check("fig4", f"(d) cycle-shuffle dD full {p:+d} dB", y1(d, "d", "dD cycle-shuffle atomic", p), lv(f"S_full_{PS[p]}"), EX)
    check("fig4", f"(d) cycle-shuffle dD static {p:+d} dB", y1(d, "d", "dD cycle-shuffle static_mh", p), lv(f"S_static_{PS[p]}"), EX)
for des in ("dwell", "shuffle"):
    check("fig4", f"(b) histograms identical ({des})", y1(d, "b", f"histogram identical {des}"), 0.0, 0.0)
for _, r in d[d.series.str.startswith("time-weighted")].iterrows():
    check("fig4", f"(b) {r.series} vs stored T_dwell_over_tau", r.y, float(r.note.split("=")[-1]), 5e-7)

d = pv("figS5")
for p in (-6, 0, 3):
    for xi in (0.1, 1, 4):
        check("figS5", f"(a) dD full {p:+d} dB xi={xi:g}", y1(d, "a", f"dD atomic {p:+d} dB", xi), lv(f"D_full_{PS[p]}_xi_{xi:g}"), EX)
    check("figS5", f"(a) dD LTI+static {p:+d} dB xi=4", y1(d, "a", f"dD static_lti_mh {p:+d} dB", 4), lv(f"D_lti_{PS[p]}_xi4"), EX)
for p in (0, 3):
    check("figS5", f"(a) dD static {p:+d} dB xi=4", y1(d, "a", f"dD static_mh {p:+d} dB", 4), lv(f"D_smh_{PS[p]}_xi4"), EX)
check("figS5", "(b) dwell contrasts within 25%", y1(d, "b", "dwell contrasts within 25%"), lv("E_dsh_dwell", 0), 0)
cl = sel(d, "b", "DSH vs full: shuffle cycles")
for k, p in enumerate((3, 6)):
    check("figS5", f"(b) clustering {p:+d} dB (DSH, full)", [float(cl.y.iloc[k]), float(cl.x.iloc[k])], [lv(f"E_dsh_clu_{PS[p]}", 0), lv(f"E_dsh_clu_{PS[p]}", 1)], EX)

d = pv("fig6")
check("fig6", "(a) tau_1e LO step (us)", y1(d, "a", "e(t) LO_step_+1pct"), lv("T_LO_1e") * 1e6, EX)
check("fig6", "(a) tau_1e linear-regime step (us)", y1(d, "a", "e(t) IF_step_small_0.05to0.10E1"), lv("T_lin_1e") * 1e6, EX)
check("fig6", "(a) tau_atom (+5% step, us)", y1(d, "a", "e(t) IF_step_P1dB_1.00to1.05E1"), lv("T_atom") * 1e6, EX)
check("fig6", "(b) plateau first us / last 10 us", [y1(d, "b", "window median 0.56 (first µs)"), y1(d, "b", "window median 0.43 (last 10 µs)")], list(lv("P_onset_+3")), EX)
check("fig6", "(b) after: first us / 1-5 us", [y1(d, "b", "window median 0.36 (first µs after)"), y1(d, "b", "window median 0.60 (1–5 µs after)")], list(lv("P_rec_+3")[:2]), EX)
tr = ys(d, "b", "trace")
check("fig6", "(b) plotted trace inside axis range 0-1.3 (no clipping)", [float(tr.min() >= 0), float(tr.max() <= 1.3)], [1, 1], 0)

d = pv("fig7")
for p in (0, 3):
    ref = y1(d, "a", f"CW-static gain {p:+d} dB")
    check("fig7", f"(a) plateau vs CW-static, unmodulated {p:+d} dB", y1(d, "a", f"late-plateau gain cw {p:+d} dB") / ref - 1, lv(f"P_QS_{PS[p]}"), 1e-6)
    check("fig7", f"(a) plateau vs CW-static, QPSK carrier {p:+d} dB", y1(d, "a", f"late-plateau gain qpsk {p:+d} dB") / ref - 1, lv(f"P_carrier_{PS[p]}"), 1e-6)
for car in ("cw", "qpsk"):
    for xi, key in ((0.05, "xi_0.05"), (4, "xi_4"), (32, "xi_32")):
        check("fig7", f"(b) D_ref full/static {car} xi={xi:g}", [y1(d, "b", f"D_ref atomic {car}", xi), y1(d, "b", f"D_ref static_mh {car}", xi)], list(lv(f"Q_ref_{car}_{key}")), 1e-9)
    for xi, key in ((0.05, "xi_0.05"), (32, "xi_32")):
        check("fig7", f"(c) gain full/static {car} xi={xi:g}", [y1(d, "c", f"gain atomic {car}", xi), y1(d, "c", f"gain static_mh {car}", xi)], list(lv(f"Q_gain_{car}_{key}")), 1e-9)
check("fig7", "(b) D_ref LTI+static cw xi=32", y1(d, "b", "D_ref static_lti_mh cw", 32), lv("Q_ref_lti_cw_xi_32"), 1e-9)
check("fig7", "(d) prespecified full-static gap xi=4 / 32", list(ys(d, "d", "full-static gap")), list(lv("Q_DBLA_cw")), EX)
check("fig7", "(d) static D at xi=0.05 (ill-conditioning)", y1(d, "d", "D (prespecified) static_mh cw", 0.05), lv("Q_DBLA_patho"), EX)

d = pv("fig8")
cases = ["CW", "OFF_+0.125", "OFF_+1.31", "OFF_-1.31", "OFF_+2.62", "OFF_-2.62", "QPSK_REF", "QPSK_REF_s2", "QPSK_SLOW", "QPSK_FAST", "QPSK_JUMP", "QPSK_LONGRAMP", "TOGGLE"]
for c in cases:
    check("fig8", f"(a) X {c}", y1(d, "a", f"X full {c}", c), lv(f"F_X_{c}"), EX)
chk = d[(d.panel == "a") & d.series.str.startswith("X check")]
dx = [abs(float(r.y) - lv(f"F_X_{r.x}")) for _, r in chk.iterrows()]
check("fig8", "(a) max |dX| numerical checks", max(dx), lv("F_numerics"), EX)
xl = d[(d.panel == "a") & d.series.str.startswith("X LTI+static")].y.astype(float)
g0, g1 = lv("F_gmodels")
check("fig8", "(a) surrogate |X| <= spread of pipeline gains 0.434-0.436", float(xl.abs().max() <= (g1 - g0) / g0 + 1e-12), 1, 0)
check("fig8", "(b) X original slow/ref/fast", list(ys(d, "b", "X original")), [lv("F_X_QPSK_SLOW"), lv("F_X_QPSK_REF"), lv("F_X_QPSK_FAST")], EX)
check("fig8", "(b) X six zero-drift sequences slow/ref/fast", list(ys(d, "b", "X zero-drift six sequences")), [lv(f"RV_X_OP1_{k}") for k in ("SLOW", "REF", "FAST")], EX)
c3 = sel(d, "c", "X vs peak |f_inst - f_IF| (phase)")
check("fig8", "(c) X longramp/ref/jump", list(c3.y.astype(float)), [lv("F_X_QPSK_LONGRAMP"), lv("F_X_QPSK_REF"), lv("F_X_QPSK_JUMP")], EX)
check("fig8", "(c) peak excursion jump/ref/longramp (MHz)", sorted(c3.x.astype(float), reverse=True), list(lv("F_peak")), 1e-9)
mf = ys(d, "d", "min |z| rel. pre, first 0.5 us, full")
ml = ys(d, "d", "min |z| rel. pre, first 0.5 us, small-signal")
lag = ys(d, "d", "max response-phase lag (rad)")
check("fig8", "(d) pi/2 min magnitude range, full", [mf.min(), mf.max()], list(lv("H_pi2_mag")), EX)
check("fig8", "(d) pi/2 min magnitude range, small-signal", [ml.min(), ml.max()], list(lv("H_pi2_lin")), EX)
check("fig8", "(d) pi/2 phase-lag range (magnitude > 0.4)", [lag[mf > 0.4].min(), lag[mf > 0.4].max()], list(lv("H_pi2_phase")), EX)
check("fig8", "(d) OP1 power-step periods", list(ys(d, "d", "osc period OP1 power (us)")), list(lv("H_per_power")), EX)
php = ys(d, "d", "osc period OP1 phase (us)")
check("fig8", "(d) OP1 phase-step period range", [php.min(), php.max()], list(lv("H_per_phase")), EX)
check("fig8", "(d) OP2 power / phase median period", [y1(d, "d", "osc period OP2 power (us)"), float(np.median(ys(d, "d", "osc period OP2 phase (us)")))], list(lv("V_c9a")), EX)
e1 = {"reference": "RV_X_OP1_REF", "slow": "RV_X_OP1_SLOW", "fast": "RV_X_OP1_FAST", "steps": "RV_X_OP1_JUMP",
      "900 ns ramps": "RV_X_OP1_LONGRAMP", "offset +0.125": "F_X_OFF_+0.125", "offset +1.31": "F_X_OFF_+1.31", "offset −1.31": "F_X_OFF_-1.31"}
e2 = {"reference": "RV_X_OP2_REF", "slow": "RV_X_OP2_SLOW", "fast": "RV_X_OP2_FAST", "steps": "RV_X_OP2_JUMP",
      "900 ns ramps": "RV_X_OP2_LONGRAMP", "offset +0.125": "V_X_LO2_OFF_+0.125", "offset +1.31": "V_X_LO2_OFF_+1.31", "offset −1.31": "V_X_LO2_OFF_-1.31"}
for k in e1:
    check("fig8", f"(e) OP1 {k}", y1(d, "e", f"X OP1 {k}"), lv(e1[k]), EX)
for k in e2:
    check("fig8", f"(e) OP2 {k}", y1(d, "e", f"X OP2 {k}"), lv(e2[k]), EX)
ev = d[d.panel == "e"].y.astype(float)
check("fig8", "(e) all values inside axis range -0.92..0.72 (no clipping)", [float(ev.min() > -0.92), float(ev.max() < 0.72)], [1, 1], 0)
av = d[(d.panel == "a") & d.series.str.startswith(("X full", "X check"))].y.astype(float)
check("fig8", "(a) all values inside axis range -2.2..0.85 (no clipping)", [float(av.min() > -2.2), float(av.max() < 0.85)], [1, 1], 0)

d = pv("figS1")
ad = d[(d.panel == "a") & d.series.str.startswith("AIR diff")]
for n in (1501, 3001, 4001):
    check("figS1", f"(a) max |AIR diff| at Nd={n}", ad[[float(v) == n for v in ad.x]].y.abs().max(), lv(f"N_nd{n}_AIR"), EX)
check("figS1", "(b) plateau-onset period at Nd 1501/4001/8001", list(ys(d, "b", "osc period plateau onset")), list(lv("N_osc_on")), EX)

d = pv("figS2")
W = pd.read_csv(MS.parent / "results" / "p1db_waveform_stage05" / "03_waveform_statistics_extended.csv")
W = W[W.config_version.isin(["fair_v2", "fair_v4_randomized_balanced_ofdm_training"]) & ~((W.modulation == "OFDM") & (W.config_version == "fair_v2"))]
per_seed = W.drop_duplicates(["modulation", "seed"]).groupby("modulation").B_occ_99_Hz.mean() / 1e6
lo, hi = [float(v) for v in L["R_B99"]["text"].replace("MHz", "").replace("–", " ").split()]
for f, lf in FMT.items():
    got = y1(d, "", f"B99 mean {f} (MHz)")
    check("figS2", f"B99 mean {f} (MHz) = equal-weight mean of stored per-seed B99", got, per_seed[lf], 1e-9)
    check("figS2", f"B99 mean {f} rounds into the quoted range {lo:.2f}-{hi:.2f} MHz", float(lo - 0.005 <= got < hi + 0.005), 1, 0)

d = pv("figS3")
for xi in (0.1, 1, 4):
    check("figS3", f"(a) dD full -6 dB xi={xi:g}", y1(d, "a", "dD atomic -6 dB", xi), lv(f"D_full_-6_xi_{xi:g}"), EX)
    check("figS3", f"(b) dD full +6 dB xi={xi:g}", y1(d, "b", "dD atomic +6 dB", xi), lv(f"D_full_+6_xi_{xi:g}"), EX)
check("figS3", "(b) dD LTI+static +6 dB xi=4", y1(d, "b", "dD static_lti_mh +6 dB", 4), lv("D_lti_+6_xi4"), EX)

d = pv("figS4")
check("figS4", "(a) Rydberg population change range, phase cases (%)", list(ys(d, "a", "range phase cases (%)")), [100 * v for v in lv("F_ryd")], 1e-7)
check("figS4", "(b) max |probe shift| phase cases (dB)", y1(d, "b", "max |shift| phase cases (dB)"), lv("F_probe"), EX)
check("figS4", "(b) max |probe shift| offsets (dB)", y1(d, "b", "max |shift| offsets (dB)"), lv("F_probe_off"), EX)

d = pv("fig5")
check("fig5", "(a) fitted gain QPSK +3 dB, full, xi 0.05 / 4", [y1(d, "a", "gain_rel atomic qpsk +3 dB", 0.05), y1(d, "a", "gain_rel atomic qpsk +3 dB", 4)],
      list(lv("RV_g_qpsk_ref_xi4")), EX)
check("fig5", "(b) fitted gain unmodulated +3 dB, full, xi 0.05 / 0.25 / 4",
      [y1(d, "b", "gain_rel atomic cw +3 dB", x) for x in (0.05, 0.25, 4)], list(lv("RV_g_cw_ref_xi4")), EX)
check("fig5", "(c) D_ref full, xi=4 minus reference = paired mean change",
      y1(d, "c", "D_ref atomic +3 dB", 4) - y1(d, "c", "D_ref atomic +3 dB", 0.05), lv("RV_p1_Dref_C1_full"), 1e-12)
check("fig5", "(c) D_ref LTI+static, xi=4 minus reference = paired mean change",
      y1(d, "c", "D_ref static_lti_mh +3 dB", 4) - y1(d, "c", "D_ref static_lti_mh +3 dB", 0.05), lv("RV_p1_Dref_C1_lti"), 1e-12)
check("fig5", "(d) declustering D_ref, full, +3 dB, 2-us span", y1(d, "d", "dD_ref shuffle atomic +3 dB", 2), lv("RV_p1_Dref_C2_full"), 1e-12)
check("fig5", "(d) declustering D_ref, full, +3 dB, span range 2-16 us",
      [min(ys(d, "d", "dD_ref shuffle atomic +3 dB")), max(ys(d, "d", "dD_ref shuffle atomic +3 dB"))], list(lv("R2_C2_Dref_span_range")), 1e-12)
check("fig5", "(d) declustering D_ref, full, +6 dB, 2-us span", y1(d, "d", "dD_ref shuffle atomic +6 dB", 2), lv("R2_C2_Dref_full_6dB"), 1e-12)

d = pv("figS6")
for op in ("OP1", "OP2"):
    for db in (-6, -3, 0, 3, 6, 8):
        check("figS6", f"(a) X REF {op} {db:+d} dB", y1(d, "a", f"X REF {op}", db), lv(f"RV_P3_{op}_{db:+d}"), EX)
check("figS6", "(a) max |X| of the g(a,f) surrogate", max(np.abs(ys(d, "a", "X REF g(a,f) surrogate OP1")).max(), np.abs(ys(d, "a", "X REF g(a,f) surrogate OP2")).max()),
      lv("RV_P5_T3max"), EX)
check("figS6", "(b) OP1 R at -0.25 / +0.125 MHz", [y1(d, "b", "R OP1", -0.25), y1(d, "b", "R OP1", 0.125)], list(lv("RV_R_OP1_near")), EX)
check("figS6", "(b) OP1 R at +3 MHz", y1(d, "b", "R OP1", 3.0), lv("RV_R_OP1_max"), EX)
check("figS6", "(b) OP3 R at +0.125 MHz", y1(d, "b", "R OP3", 0.125), lv("RV_R_OP3_near"), EX)
check("figS6", "(c) X REF (3 sequences) OP2 / OP3 / OP1", [y1(d, "c", f"X REF 3 sequences {op}") for op in ("OP2", "OP3", "OP1")], list(lv("RV_XREF3")), EX)

# ------------------------------------------------------------------ clipping (generic, recorded by figstyle.save)
for fig in MAIN + SUPP:
    a = pv(fig)
    row = a[a.panel == "_audit"]
    check(fig, "no data point outside its axis limits (lines, error bars, bars, fills, histograms)", float(row.y.iloc[0]), 0.0, 0)

# ------------------------------------------------------------------ files, fonts, wording
d = pv("figS8")
for f_ in (5, 10, 15):
    check("figS8", f"(b) X REF {f_} MHz", y1(d, "b", f"X REF {f_} MHz"), lv(f"R2_XREF_IF{f_}"), 1e-12)
    check("figS8", f"(c) dwell gain change full / static {f_} MHz", [y1(d, "c", f"d gain full {f_} MHz"), y1(d, "c", f"d gain static {f_} MHz")],
          [lv(f"R2_dg_IF{f_}", 0), lv(f"R2_dg_IF{f_}", 3)], 1e-12)
check("figS8", "(b) X REF weak probe", y1(d, "b", "X REF weak probe"), lv("R2_wp_X"), 1e-12)
check("figS8", "(c) dwell gain change full / static weak probe", [y1(d, "c", "d gain full weak probe"), y1(d, "c", "d gain static weak probe")],
      [lv("R2_wp_dg", 0), lv("R2_wp_dg", 3)], 1e-12)
check("figS8", "(d) S5 = max R - 1 at 5/10/15 MHz", [max(ys(d, "d", f"R {f_} MHz")) - 1 for f_ in (5, 10, 15)], list(lv("R2_S5_IF")), 1e-12)
check("figS8", "(a) matched level 10 MHz and weak probe (V/m)", [float(sel(d, "a", "matched level 10 MHz").x.iloc[0]), float(sel(d, "a", "matched level weak probe").x.iloc[0])],
      [lv("R2_aH_IF10"), lv("R2_wp_aH")], 1e-12)
check("figS8", "(e) slow time constants weak / standard probe (us)", [y1(d, "e", "slow time constant weak (us)"), y1(d, "e", "slow time constant standard (us)")],
      [lv("R2_wp_tau", 0) * 1e6, lv("R2_wp_tau", 2) * 1e6], 1e-12)

d = pv("fig9")
for cfg in ("C1", "C2", "C3", "C4"):
    check("fig9", f"(b) X at signal/LO 0.2/0.358/0.4 {cfg}", [y1(d, "b", f"X {cfg}", r) for r in (0.2, 0.358, 0.4)], list(lv(f"TV_Xsig_{cfg}")), 1e-12)
check("fig9", "(c) X at CW gain 0.9/0.85, 5 MHz", [y1(d, "c", "X C1", g) for g in (0.9, 0.85)], list(lv("TV_Xcw_C1")), 1e-12)
check("fig9", "(c) X at CW gain 0.9/0.85, 15 MHz", [y1(d, "c", "X C3", g) for g in (0.9, 0.85)], list(lv("TV_Xcw_C3")), 1e-12)
for cfg, op in (("C1", "OP1"), ("C5", "OP3"), ("C6", "OP2")):
    check("fig9", f"(c) X at CW gain 0.8/0.6/0.4 {op}", [y1(d, "c", f"X {cfg}", g) for g in (0.8, 0.6, 0.4)], list(lv(f"TV_op_{op}")), 1e-12)
check("fig9", "(e) X_coh maximum per LO field (0.35/0.425/0.5)", [max(ys(d, "e", f"X_coh vs r2, LO {lo}")) for lo in ("0.35", "0.425", "0.5")], list(lv("TV_xcoh_max")), 1e-12)
gmin = [sel(d, "d", f"g_CW vs r2, LO {lo}") for lo in ("0.35", "0.425", "0.5")]
check("fig9", "(d) r2 of the CW-gain minimum within the ledger range", [float(g.x.astype(float).values[g.y.astype(float).values.argmin()]) for g in gmin],
      [min(max(float(g.x.astype(float).values[g.y.astype(float).values.argmin()]), lv("TV_r2_gmin", 0)), lv("TV_r2_gmin", 1)) for g in gmin], 1e-12)

d = pv("fig10")
nm = lv("TV_gmp_nmse")
for fam in ("dwell", "shuffle", "phase", "fair"):
    best = min(y1(d, "a", f"median NMSE {m}", fam) for m in ("static_mh", "static_lti_mh", "dsh"))
    check("fig10", f"(a) median NMSE GMP / best surrogate, {fam}", [y1(d, "a", "median NMSE gmp", fam), best], list(nm[fam]), 1e-12)
check("fig10", "(b) GMP |AIR error| QPSK/16QAM/OFDM/CE", [y1(d, "b", "mean |AIR error| gmp", m) for m in ("QPSK", "16QAM", "OFDM", "CE")], list(lv("TV_gmp_air")), 1e-12)

d = pv("figS9")
check("figS9", "(f) X_coh maximum per LO field (0.35/0.425/0.5)", [max(ys(d, "f", f"X_coh vs r2, LO {lo}")) for lo in ("0.35", "0.425", "0.5")], list(lv("TV_xcoh_max")), 1e-12)

d = pv("figS10")
ofdm = [y1(d, "b", "static_mh minus full, +3 dB, OFDM", c) for c in (0.25, 0.5, 1.0, 2.0, 4.0)]
check("figS10", "(b) static minus full, OFDM +3 dB, range over c", [min(ofdm), max(ofdm)], list(lv("TV_noise_ofdm3")), 1e-12)
check("figS10", "(a) full-model loss +3 dB at c = 4, CE / OFDM", [y1(d, "a", "loss +3 dB CE", 4.0), y1(d, "a", "loss +3 dB OFDM", 4.0)], list(lv("TV_noise_loss4")), 1e-12)

d = pv("figS7")
check("figS7", "(a) a_H output minimum (rel. pre-step)", y1(d, "a", "a_H output minimum (rel. pre-step), 49.8-54.75 us"), lv("R2_p4_out", 0), 1e-12)
check("figS7", "(b) rho43 +IF phase at the window end (rad)", y1(d, "b", "a_H rho43 +IF phase at the window end (rad)"), lv("R2_p4_drift", 0), 1e-12)

FORBID = [r"out of lock", r"phase[- ]lock", r"(?<!symbol )phase transition", r"\bM_FULL\b", r"\bM_STATIC", r"\bM_LTI", r"\bDSH\b", r"\bTOGGLE\b",
          r"QPSK_", r"OFF_", r"_zd\b", r"\bzd\b", r"\bLO1\b", r"\bLO2\b", r"\bLO'", r"FU-[A-Z]", r"\bStage", r"\bstage\b", r"Step \d", r"fig_?\w*\d",
          r"\.csv", r"\.npz", r"\.json", r"\.png", r"candidate", r"universal", r"mechanism", r"\bC9a\b", r"D_BLA", r"results/", r"manuscript/",
          r"pre-?regist", r"extra compression", r"τatom", r"tau_atom"]
for fig in MAIN + SUPP:
    out = FIG / ("final" if fig in MAIN else "supplement")
    res = RESULTS.setdefault(fig, {"numeric": [], "files": {}, "text": {}})
    files = {ext: (out / f"{fig}.{ext}").exists() for ext in ("pdf", "svg", "png")}
    info = subprocess.run(["pdfinfo", str(out / f"{fig}.pdf")], capture_output=True, text=True).stdout
    m = re.search(r"Page size:\s+([\d.]+) x ([\d.]+) pts", info)
    w_in, h_in = float(m.group(1)) / 72, float(m.group(2)) / 72
    fonts = subprocess.run(["pdffonts", str(out / f"{fig}.pdf")], capture_output=True, text=True).stdout.splitlines()[2:]
    font_rows = [ln.split() for ln in fonts if ln.strip()]
    emb = all(("yes" in r[-5:-2][:1]) for r in font_rows) if font_rows else False
    fnames = sorted({r[0].split("+")[-1] for r in font_rows})
    im = Image.open(out / f"{fig}.png")
    dpi = im.info.get("dpi", (0, 0))[0]
    svg = (out / f"{fig}.svg").read_text()
    svg_text = len(re.findall(r"<text\b", svg))
    svg_paths = bool(re.search(r"/home/|results/|\.csv|\.npz", svg))
    txt = subprocess.run(["pdftotext", "-layout", str(out / f"{fig}.pdf"), "-"], capture_output=True, text=True).stdout
    hits = sorted({h for pat in FORBID for h in re.findall(pat, txt, flags=re.I)})
    ascii_minus = re.findall(r"(?<![\w)])-\d", txt)
    res["files"] = {"pdf_svg_png": files, "size_in": [round(w_in, 2), round(h_in, 2)], "size_mm": [round(w_in * 25.4, 1), round(h_in * 25.4, 1)],
                    "fonts": fnames, "fonts_embedded": emb, "png_px": list(im.size), "png_dpi": round(float(dpi)), "svg_text_elements": svg_text,
                    "svg_contains_paths_or_files": svg_paths}
    res["text"] = {"forbidden_hits": hits, "ascii_hyphen_minus_before_digit": len(ascii_minus)}

src_bad = {}
for f in HERE.glob("*.py"):
    s = f.read_text()
    bad = re.findall(r"^\s*(?:import|from)\s+(run06|common|stage05|p1db_comm_stage0|utils\.transient_quantum|dsh06|models06|tau_atom|nd_convergence)\b", s, flags=re.M)
    if bad:
        src_bad[f.name] = bad

# ------------------------------------------------------------------ verdicts
summary = {}
all_ok = not src_bad
for fig in MAIN + SUPP:
    r = RESULTS[fig]
    n_ok = sum(c["pass"] for c in r["numeric"])
    f = r["files"]
    ok = (n_ok == len(r["numeric"]) and all(f["pdf_svg_png"].values()) and f["fonts_embedded"] and f["svg_text_elements"] == 0
          and not f["svg_contains_paths_or_files"] and not r["text"]["forbidden_hits"] and r["text"]["ascii_hyphen_minus_before_digit"] == 0
          and f["png_dpi"] >= 600)
    all_ok &= ok
    summary[fig] = {"numeric_checks": f"{n_ok}/{len(r['numeric'])}", "size_in": f["size_in"], "size_mm": f["size_mm"], "fonts": f["fonts"],
                    "fonts_embedded": f["fonts_embedded"], "png": f"{f['png_px'][0]}x{f['png_px'][1]} px @ {f['png_dpi']} dpi",
                    "forbidden_text": r["text"]["forbidden_hits"], "ascii_minus": r["text"]["ascii_hyphen_minus_before_digit"],
                    "verdict": "PASS" if ok else "FAIL"}
    for c in r["numeric"]:
        if not c["pass"]:
            print("FAIL", fig, c)
out = {"summary": summary, "simulation_imports_in_figure_code": src_bad, "details": RESULTS, "all_pass": bool(all_ok)}
(AUD / "figure_audit_results.json").write_text(json.dumps(out, indent=1, ensure_ascii=False))
for fig, s in summary.items():
    print(f"{fig:6s} {s['verdict']}  numeric {s['numeric_checks']:>6s}  {s['size_mm'][0]}x{s['size_mm'][1]} mm  fonts {','.join(s['fonts'])} emb={s['fonts_embedded']}  "
          f"{s['png']}  forbidden={s['forbidden_text'] or 0} ascii-minus={s['ascii_minus']}")
print("simulation imports in figure code:", src_bad or "none")
print("ALL PASS" if all_ok else "FAILURES PRESENT")
sys.exit(0 if all_ok else 1)
