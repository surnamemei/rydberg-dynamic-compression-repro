"""P2: coherent vs magnitude vs RMS gain of all saved constant-amplitude phase runs (no new runs).

Rules: results/revision2/00_decision_rules.json, P2. Window: [T_end - 20.5 us, T_end - 0.5 us), as for the published X.
  g_coh = |mean(z_full conj(z_lin)/|z_lin|)| / mean|z_lin|      g_mag = mean|z_full| / mean|z_lin|
  g_rms = sqrt(mean|z_full|^2 / mean|z_lin|^2)                   X_q = 1 - g_q / g_q(CW);  X = published median metric
  kappa = |mean(|z_full| e^{j theta})| / mean|z_full|, theta = arg(z_full conj(z_lin)); g_coh = g_mag kappa
  shares of ln(1 - X_coh): magnitude ln(1 - X_mag), phase misalignment ln(kappa / kappa_CW)
"""
from __future__ import annotations

import json
import re

import numpy as np
import pandas as pd
from scipy import stats

from rev2_common import REV2, ROOT, T0

REV = ROOT / "results" / "revision"
PM = ROOT / "results" / "stage06_dwell_physics" / "phase_mechanism" / "runs"


def window(t):
    t_end = T0 + 60e-6
    return (t >= t_end - 20.5e-6) & (t < t_end - 0.5e-6)


def metrics(path):
    z = np.load(path)
    m = window(z["t"])
    zf, zl = z["z_full"][m], z["z_linear"][m]
    rot = zf * np.conj(zl) / np.abs(zl)
    theta = np.angle(rot)
    unit = np.mean(np.exp(1j * theta))
    return {"g_coh": abs(np.mean(rot)) / np.mean(np.abs(zl)), "g_mag": np.mean(np.abs(zf)) / np.mean(np.abs(zl)),
            "g_rms": np.sqrt(np.mean(np.abs(zf) ** 2) / np.mean(np.abs(zl) ** 2)),
            "g_med": float(np.median(np.abs(zf) / np.abs(zl))), "kappa": abs(np.mean(rot)) / np.mean(np.abs(zf)),
            "phase_mean_rad": float(np.angle(unit)), "phase_circ_sd_rad": float(np.sqrt(-2 * np.log(min(abs(unit), 1.0)))),
            "phase_p05_p95_rad": float(np.diff(np.percentile(np.unwrap(theta), [5, 95]))[0]),
            "mag_cv": float(np.std(np.abs(zf)) / np.mean(np.abs(zf))), "lin_mag_cv": float(np.std(np.abs(zl)) / np.mean(np.abs(zl)))}


def inventory():
    """(group, op, case, seed, level, path, cw_path) for every saved constant-amplitude phase run."""
    rows = []
    for p in sorted((REV / "p2").glob("OP*_*_s*_H.npz")):
        op, case, seed = re.match(r"(OP\d)_(\w+?)_s(\d+)_H", p.stem).groups()
        rows.append(("revision_P2", op, case, int(seed), "H", p, REV / "p2" / f"{op}_CW_H.npz"))
    for p in sorted((REV / "p4").glob("OP3_*_H.npz")):
        if p.stem == "OP3_CW_H":
            continue
        m = re.match(r"OP3_REF_s(\d+)_H", p.stem)
        case, seed = ("REF", int(m.group(1))) if m else (p.stem[4:-2], 0)
        rows.append(("revision_P4_OP3", "OP3", case, seed, "H", p, REV / "p4" / "OP3_CW_H.npz"))
    for p in sorted((REV / "p45").glob("OP*_OFF*.npz")):
        op, off, lev = re.match(r"(OP\d)_(OFF[+-][\d.]+)_(H|p[+-]\d+)", p.stem).groups()
        cw = REV / "p2" / f"{op}_CW_H.npz" if lev == "H" else REV / "p3" / f"{op}_CW_{lev}.npz"
        rows.append(("revision_P4P5_offsets", op, off, 0, lev, p, cw))
    for p in sorted((REV / "p3").glob("OP*_REF_s*_p*.npz")):
        op, seed, lev = re.match(r"(OP\d)_REF_s(\d+)_(p[+-]\d+)", p.stem).groups()
        rows.append(("revision_P3_levels", op, "REF", int(seed), lev, p, REV / "p3" / f"{op}_CW_{lev}.npz"))
    for p in sorted(PM.glob("*_Nd4001_dt1.npz")):
        if p.stem.startswith("CW_"):
            continue
        rows.append(("stage06_phase_mechanism(non-zero-drift)", "OP1", p.stem.replace("_Nd4001_dt1", ""), 0, "H", p, PM / "CW_Nd4001_dt1.npz"))
    return rows


def ci(v):
    v = np.asarray(v, float)
    n = len(v)
    h = stats.t.ppf(.975, n - 1) * v.std(ddof=1) / np.sqrt(n) if n > 1 else np.nan
    return float(v.mean()), float(v.mean() - h), float(v.mean() + h), n


def main():
    cache, rows = {}, []
    for group, op, case, seed, lev, p, cw in inventory():
        if not cw.exists():
            continue
        if cw not in cache:
            cache[cw] = metrics(cw)
        c, r = cache[cw], metrics(p)
        row = {"group": group, "op": op, "case": case, "seed": seed, "level": lev, "run": str(p.relative_to(ROOT)), "cw": str(cw.relative_to(ROOT)), **r}
        for q in ("coh", "mag", "rms", "med"):
            row[f"X_{q}"] = 1 - r[f"g_{q}"] / c[f"g_{q}"]
        row["X"] = row.pop("X_med")
        row["kappa_rel_CW"] = r["kappa"] / c["kappa"]
        row["phase_lag_vs_CW_rad"] = float(np.angle(np.exp(1j * (r["phase_mean_rad"] - c["phase_mean_rad"]))))
        L = np.log(1 - row["X_coh"])
        row["share_magnitude"] = np.log(1 - row["X_mag"]) / L if abs(L) > 1e-9 else np.nan
        row["share_phase"] = np.log(row["kappa_rel_CW"]) / L if abs(L) > 1e-9 else np.nan
        row["kappa_CW"] = c["kappa"]
        rows.append(row)
    per = pd.DataFrame(rows)
    per.to_csv(REV2 / "20_p2_per_run.csv", index=False)
    # Summary per operating point and case (a_H, revision runs used by the manuscript).
    main_ = per[per.group.isin(["revision_P2", "revision_P4_OP3"]) & (per.case.isin(["REF", "FAST", "SLOW", "JUMP", "LONGRAMP"]))]
    summ = []
    for (op, case), g in main_.groupby(["op", "case"]):
        d = {"op": op, "case": case, "n": len(g)}
        for q in ("X", "X_coh", "X_mag", "X_rms", "kappa_rel_CW", "phase_lag_vs_CW_rad", "phase_circ_sd_rad", "share_magnitude", "share_phase", "mag_cv"):
            m, lo, hi, _ = ci(g[q])
            d.update({q: m, f"{q}_lo": lo, f"{q}_hi": hi})
        d["max_abs_X_minus_X_mag"] = float(np.max(np.abs(g.X - g.X_mag)))
        summ.append(d)
    summ = pd.DataFrame(summ)
    summ.to_csv(REV2 / "21_p2_summary.csv", index=False)
    ref = summ[summ.case == "REF"].set_index("op")
    verdict = {
        "reference_phase_share_by_op": ref.share_phase.to_dict(),
        "reference_magnitude_share_by_op": ref.share_magnitude.to_dict(),
        "phase_share_exceeds_0.25_at_any_op": bool((ref.share_phase > .25).any()),
        "X_coh_called_compression": bool(not (ref.share_phase > .25).any()),
        "X_vs_X_mag_max_abs_diff_by_op_all_cases": {op: float(np.max(np.abs(g.X - g.X_mag))) for op, g in main_.groupby("op")},
        "X_vs_X_mag_max_abs_diff_by_op_REF": {op: float(np.max(np.abs(g.X - g.X_mag))) for op, g in main_[main_.case == "REF"].groupby("op")},
    }
    verdict["extra_compression_name_kept"] = bool(all(v <= .05 for v in verdict["X_vs_X_mag_max_abs_diff_by_op_all_cases"].values()))
    (REV2 / "22_p2_verdict.json").write_text(json.dumps(verdict, indent=1, default=float))
    print(summ[["op", "case", "n", "X", "X_mag", "X_rms", "X_coh", "kappa_rel_CW", "phase_lag_vs_CW_rad", "phase_circ_sd_rad",
                "share_magnitude", "share_phase", "max_abs_X_minus_X_mag"]].round(4).to_string())
    print(json.dumps(verdict, indent=1, default=float))


if __name__ == "__main__":
    main()
