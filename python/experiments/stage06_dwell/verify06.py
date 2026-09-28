"""Independent audit of the Stage-0.6 campaign outputs (does not import analyze06)."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

OUT = Path(__file__).resolve().parents[3] / "results" / "stage06_dwell_physics"
MODELS = {"atomic", "static", "static_lti", "static_mh", "static_lti_mh", "linear"}
lines = []


def say(s=""):
    print(s, flush=True)
    lines.append(str(s))


def ci(x):
    x = np.asarray(x, float)
    h = stats.t.ppf(.975, len(x) - 1) * x.std(ddof=1) / np.sqrt(len(x))
    return x.mean(), x.mean() - h, x.mean() + h, len(x)


# 1. completeness against the campaign plan
plan = {("main", "dwell"): (8, 7, 4), ("main", "shuffle"): (8, 4, 4), ("num_Nd8001", "dwell"): (2, 4, 1),
        ("num_dt0p5", "dwell"): (2, 4, 1), ("num_Nd8001", "shuffle"): (2, 2, 1), ("num_dt0p5", "shuffle"): (2, 2, 1)}
say("== 1. completeness")
rows = {}
for (tag, design), (nr, nv, npow) in plan.items():
    raw = pd.read_csv(OUT / f"rows_{tag}_{design}.csv")
    d = raw.drop_duplicates(["job_id", "model"], keep="last")
    rows[(tag, design)] = d
    jobs = d.job_id.nunique()
    per_job_models = d.groupby("job_id").model.apply(set)
    ok_models = all(s == MODELS for s in per_job_models)
    nfiles = len(list((OUT / "jobs" / tag).glob(f"{design}_*.json")))
    dup = len(raw) - len(d)
    conflicting = raw.groupby(["job_id", "model"]).AIR.nunique().max()
    say(f"{tag:11s} {design:7s} jobs={jobs} expected={nr * nv * npow} checkpoint_files={nfiles} all_6_models={ok_models} "
        f"duplicate_rows_dropped={dup} max_distinct_AIR_per_(job,model)={conflicting} finite={np.isfinite(d[['AIR', 'D_BLA']]).all().all()}")

# 2. checkpoint JSON vs CSV for a sample of jobs
say("\n== 2. checkpoint JSON vs CSV (5 random jobs)")
d = rows[("main", "dwell")]
rng = np.random.default_rng(0)
for jid in rng.choice(d.job_id.unique(), 5, replace=False):
    js = {r["model"]: r for r in json.loads((OUT / "jobs" / "main" / f"{jid}.json").read_text())["rows"]}
    diffs = [abs(js[m]["D_BLA"] - d[(d.job_id == jid) & (d.model == m)].D_BLA.iloc[0]) for m in MODELS]
    say(f"{jid}: max |D_BLA json - csv| = {max(diffs):.2e}")

# 3. independent recomputation of the headline pair statistics
say("\n== 3. independent recomputation: dwell pair Delta D_BLA (B - xi_0.05), paired over realizations")
w = d.pivot_table(index=["Pavg_over_P1dB_dB", "variant", "realization"], columns="model", values="D_BLA")
for p in (-6, 0, 3, 6):
    for v in ("xi_1", "xi_4"):
        a, b = w.loc[(p, "xi_0.05")], w.loc[(p, v)]
        full = ci(b.atomic - a.atomic)
        exc = ci((b.atomic - a.atomic) - (b.static_lti_mh - a.static_lti_mh))
        say(f"{p:+d} dB {v}: dFull={full[0]:+.3f} [{full[1]:+.3f},{full[2]:+.3f}] n={full[3]}   full-LTI_MH={exc[0]:+.3f} [{exc[1]:+.3f},{exc[2]:+.3f}]")
s = rows[("main", "shuffle")]
ws = s.pivot_table(index=["Pavg_over_P1dB_dB", "variant", "realization"], columns="model", values="D_BLA")
say("shuffle BLOCK_SHUFFLED_CYCLES - ORIGINAL:")
for p in (-6, 0, 3, 6):
    a, b = ws.loc[(p, "ORIGINAL")], ws.loc[(p, "BLOCK_SHUFFLED_CYCLES")]
    f, e, st = ci(b.atomic - a.atomic), ci((b.atomic - a.atomic) - (b.static_lti_mh - a.static_lti_mh)), ci(b.static_mh - a.static_mh)
    say(f"{p:+d} dB: dFull={f[0]:+.3f} [{f[1]:+.3f},{f[2]:+.3f}]  full-LTI_MH={e[0]:+.3f} [{e[1]:+.3f},{e[2]:+.3f}]  dStaticMH={st[0]:+.3f} [{st[1]:+.3f},{st[2]:+.3f}]  "
        f"mean D full orig/cyc = {a.atomic.mean():.3f}/{b.atomic.mean():.3f} ({100 * (b.atomic.mean() / a.atomic.mean() - 1):+.1f}%)")

# 4. sanity: linear model floor, spike-count confound, controls saturation
say("\n== 4. sanity checks")
say(f"linear-model D_BLA range: {d[d.model == 'linear'].D_BLA.min():.4f} .. {d[d.model == 'linear'].D_BLA.max():.4f}")
cfg = pd.read_csv(OUT / "04_controlled_pair_configs.csv")
sp = cfg[cfg.design == "dwell"].groupby("variant")[["meta_n_spikes", "meta_spike_share_of_high_time", "T_dwell_over_tau"]].mean()
say("spikes per member (mean over realizations):\n" + sp.round(3).to_string())
for p in (0, 3):
    m = w.loc[p].groupby("variant").mean()
    say(f"{p:+d} dB mean D_BLA by variant:\n" + m[["atomic", "static_lti_mh", "static_mh", "linear"]].round(3).to_string())

# 5. plateau-interior gain from the realization-0 traces (summary wording check)
say("\n== 5. plateau-interior gain, dwell xi_4, +3 dB, realization 0 (relative to the linear small-signal fit)")
z = np.load(OUT / "traces" / "main" / "dwell_r0_xi_4_p+3_Nd4001_dt1.npz")
cp = int(z["cp"])
x = z["x"][cp:]
a = z["a"]
g = np.vdot(x, z["y_linear"][cp:]) / np.vdot(x, x)
hi = a > (a.min() + a.max()) / 2
edge = np.diff(np.r_[False, hi, False].astype(int))
st, en = np.flatnonzero(edge == 1), np.flatnonzero(edge == -1)
long_runs = [(i, j) for i, j in zip(st, en) if (j - i) > 100]  # > 5 us at 20 MHz
for m in ("atomic", "static_mh", "static_lti_mh"):
    gi = np.abs(z["y_" + m][cp:]) / np.maximum(np.abs(g * x), 1e-15)
    first = np.median(np.concatenate([gi[i:i + 20] for i, j in long_runs]))
    last = np.median(np.concatenate([gi[(i + j) // 2:j - 5] for i, j in long_runs]))
    low_after = np.median(np.concatenate([gi[j + 5:j + 40] for i, j in long_runs if j + 40 < len(gi)]))
    say(f"{m:14s} plateau first 1 us: {first:.3f}   plateau 2nd half: {last:.3f}   low level 0.25-2 us after plateau: {low_after:.3f}")

(OUT / "10_verification_report.txt").write_text("\n".join(lines) + "\n")
