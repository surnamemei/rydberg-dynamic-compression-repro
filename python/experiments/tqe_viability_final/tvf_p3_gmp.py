"""Part III: exactly one standard memory model (GMP) vs the CW-derived surrogates (plan: results/tqe_viability_final/00_plan.md, C).

Pipeline, data loaders and basis are those of the prior pass (python/experiments/tqe_viability/tv_p3_gmp.py), with the campaign
grid: nonlinear order K in {3, 5, 7}; memory depth D in {0.25, 0.5, 1, 2, 4} us (M = 20 D taps at 20 MHz for aligned and
harmonic-zone terms, L_b = M lagging lags, K_b = 2); relative ridge in {1e-8, 1e-5}. Split fixed in the plan before fitting.
"""
from __future__ import annotations

import itertools
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

from tvf_common import FIN, ROOT, git_state, now

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "tqe_viability"))
import tv_p3_gmp as G  # noqa: E402

G.KA_MAX, G.MA_MAX, G.KB, G.LB_MAX = 7, 80, 2, 80
OUT = FIN / "03_memory_model"
GRID = [dict(K=k, D=d, M=int(round(20 * d)), lam=l) for k, d, l in itertools.product((3, 5, 7), (0.25, 0.5, 1.0, 2.0, 4.0), (1e-8, 1e-5))]
COMPARATORS = ("static_mh", "static_lti_mh", "dsh")
MODELS = ("gmp",) + COMPARATORS
LABEL = {"gmp": "GMP", "static_mh": "all-zone static", "static_lti_mh": "all-zone LTI+static", "dsh": "single-slow-state", "full": "full model"}
REVP3 = ROOT / "results" / "revision" / "p3"
LEVELS = (-6, -3, 0, 3, 6, 8)


def n_coef(K, M):
    return (K + 9) * M + 2


def items_all():
    items = G.list_waveforms()
    for lev in LEVELS:
        tag = f"p{lev:+d}"
        items.append({"family": "phase_level", "path": REVP3 / f"OP1_CW_{tag}.npz", "seed": 0, "case": "CW", "level": lev, "split": "heldout_level"})
        for seed in (20260701, 20260702, 20260703):
            items.append({"family": "phase_level", "path": REVP3 / f"OP1_REF_s{seed}_{tag}.npz", "seed": seed, "case": "REF", "level": lev,
                          "split": "heldout_level"})
    return items


def load(it, comps=False):
    return G.load({**it, "family": "phase" if it["family"] == "phase_level" else it["family"]}, comps=comps)


def xmetrics(yb, lin, xwin):
    a, b = xwin
    y, l_ = yb[a:b], lin[a:b]
    return {"g_med": float(np.median(np.abs(y) / np.abs(l_))), "g_mag": float(np.mean(np.abs(y)) / np.mean(np.abs(l_))),
            "g_coh": float(abs(np.mean(y * np.conj(l_) / np.abs(l_))) / np.mean(np.abs(l_)))}


def fit():
    items = items_all()
    fam_n = {}
    for it in items:
        if it["split"] == "train":
            fam_n[it["family"]] = fam_n.get(it["family"], 0) + 1
    A = b = names = None
    t0 = time.perf_counter()
    for it in items:
        if it["split"] != "train":
            continue
        d = load(it)
        Phi, names = G.basis_max(d["x"])
        a0, a1 = d["w"]
        P_, y = Phi[a0:a1], d["y"][a0:a1]
        w = 1.0 / (fam_n[it["family"]] * (a1 - a0))
        A = (0 if A is None else A) + w * (P_.conj().T @ P_)
        b = (0 if b is None else b) + w * (P_.conj().T @ y)
    print(f"{now()} Gram matrix from {sum(fam_n.values())} training waveforms {fam_n} in {time.perf_counter() - t0:.0f} s", flush=True)
    sols = {}
    for c in GRID:
        idx = G.subset(names, c["K"], c["M"], c["M"])
        assert len(idx) == n_coef(c["K"], c["M"]), (c, len(idx))
        As, bs = A[np.ix_(idx, idx)], b[idx]
        sc = 1 / np.sqrt(np.real(np.diag(As)))
        Asc = As * sc[:, None] * sc[None, :]
        th = np.linalg.solve(Asc + c["lam"] * np.trace(Asc).real / len(idx) * np.eye(len(idx)), bs * sc) * sc
        sols[json.dumps(c)] = (idx, th)
    val = {k: {} for k in sols}
    for it in items:
        if it["split"] != "val":
            continue
        d = load(it)
        Phi, _ = G.basis_max(d["x"])
        for k, (idx, th) in sols.items():
            val[k].setdefault(it["family"], []).append(G.nmse(Phi[:, idx] @ th, d["y"], d["w"]))
    rows = []
    for k, fam in val.items():
        c = json.loads(k)
        rows.append({**c, "memory_taps": c["M"], "lagging_lags": c["M"], "n_coef_complex": n_coef(c["K"], c["M"]), "n_real_params": 2 * n_coef(c["K"], c["M"]),
                     **{f"val_nmse_{f}": float(np.median(v)) for f, v in fam.items()},
                     "val_nmse_family_mean": float(np.mean([np.median(v) for v in fam.values()]))})
    cfgdf = pd.DataFrame(rows).sort_values("val_nmse_family_mean").reset_index(drop=True)
    cfgdf["selected"] = False
    cfgdf.loc[0, "selected"] = True
    cfgdf.to_csv(OUT / "gmp_config.csv", index=False)
    top = cfgdf.iloc[0]
    key = [k for k in sols if json.loads(k) == {"K": int(top.K), "D": float(top.D), "M": int(top.M), "lam": float(top.lam)}][0]
    idx, th = sols[key]
    np.savez_compressed(OUT / "gmp_coefficients.npz", idx=idx, theta=th, selected=key,
                        names=np.array([f"{a}|{k_}|{m}" for a, k_, m in names], dtype=object))
    print(f"{now()} selected {key}: {len(idx)} complex coefficients; validation family-mean NMSE {top.val_nmse_family_mean:.4g}", flush=True)
    return items, idx, th, key


def test(items, idx, th):
    rows = []
    for it in items:
        is_ref = it["family"] == "phase" and it.get("case") == "CW"
        if it["split"] not in ("test", "heldout_level") and not is_ref:
            continue
        d = load(it, comps=True)
        Phi, _ = G.basis_max(d["x"])
        models = {"full": d["y"], "gmp": Phi[:, idx] @ th, **d["comps"]}
        g_lin, _, f_lin = G.gain_and_dref(d["lin"], d["xs"], d["w"])
        base = {"family": it["family"], "id": d["id"], "split": "reference" if is_ref else it["split"]}
        for k in ("variant", "p", "mod", "seed", "case"):
            if k in d:
                base[k] = d[k]
        if "level" in it:
            base["level"] = it["level"]
            base["case"] = it["case"]
        for m, yb in models.items():
            g, r_pow, _ = G.gain_and_dref(yb, d["xs"], d["w"])
            row = {**base, "model": m, "nmse_vs_full": np.nan if m == "full" else G.nmse(yb, d["y"], d["w"]), "gain_rel": abs(g) / abs(g_lin),
                   "D_ref": r_pow / f_lin}
            if it["family"] in ("phase", "phase_level"):
                row.update(xmetrics(yb, d["lin"], d["xwin"]))
            if it["family"] == "fair":
                tm = json.loads((ROOT / "results" / "final_validation" / "regen_stage05" / "jobs" / f"{d['id']}.json").read_text())["ce_timing_offset"]
                row["AIR"] = G.air(yb, d, tm)
            rows.append(row)
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "gmp_results.csv", index=False)
    return df


def realization_of(i):
    return int(i.split("_r")[1].split("_")[0])


def evaluate(df, key):
    t = df[df.split == "test"]
    out = {"nmse": {}, "contrasts": {}, "air_err": {}, "levels": {}}
    nm = t[t.model != "full"].groupby(["family", "model"]).nmse_vs_full.median().unstack()
    out["nmse"] = {fam: {m: float(r[m]) for m in MODELS} for fam, r in nm.iterrows()}
    dw = t[(t.family == "dwell") & (t.p == 3)].copy()
    dw["r"] = dw.id.map(realization_of)
    sh = t[(t.family == "shuffle") & (t.p == 3)].copy()
    sh["r"] = sh.id.map(realization_of)
    ph = df[df.family == "phase"]
    cw = ph[ph.case == "CW"].set_index("model")
    for m in ("full",) + MODELS:
        piv = dw[dw.model == m].pivot_table(index="r", columns="variant", values="gain_rel")
        out["contrasts"].setdefault("dwell_gain", {})[m] = float((piv["xi_4"] - piv["xi_0.05"]).mean())
        piv = sh[sh.model == m].pivot_table(index="r", columns="variant", values="D_ref")
        out["contrasts"].setdefault("decluster_Dref", {})[m] = float((piv["BLOCK_SHUFFLED_CYCLES"] - piv["ORIGINAL"]).mean())
        for case in ("REF", "FAST", "SLOW", "JUMP", "LONGRAMP"):
            sub = ph[(ph.case == case) & (ph.split == "test") & (ph.model == m)]
            for g in ("med", "coh", "mag"):
                nmx = {"med": "X", "coh": "X_coh", "mag": "X_mag"}[g]
                out["contrasts"].setdefault(f"{nmx}_{case}", {})[m] = float((1 - sub[f"g_{g}"] / cw.loc[m, f"g_{g}"]).mean())
    fa = t[(t.family == "fair") & (t.p == 3)]
    for mod, g in fa.groupby("mod"):
        full = g[g.model == "full"].set_index("id").AIR
        out["air_err"][mod] = {m: float(np.mean(np.abs(g[g.model == m].set_index("id").AIR - full))) for m in MODELS}
    lv = df[df.family == "phase_level"]
    for lev, g in lv.groupby("level"):
        cwl = g[g.case == "CW"].set_index("model")
        ref = g[g.case == "REF"]
        out["levels"][int(lev)] = {m: {nmx: float((1 - ref[ref.model == m][f"g_{gk}"] / cwl.loc[m, f"g_{gk}"]).mean())
                                       for gk, nmx in (("med", "X"), ("coh", "X_coh"), ("mag", "X_mag"))} for m in ("full",) + MODELS}

    def rep(k, m):
        full, mod = out["contrasts"][k]["full"], out["contrasts"][k][m]
        return bool(np.sign(mod) == np.sign(full) and abs(mod - full) <= .25 * abs(full))

    def keys_of(m):
        return {"dwell gain contrast": rep("dwell_gain", m), "declustering D_ref contrast": rep("decluster_Dref", m),
                "phase REF (X and X_coh)": rep("X_REF", m) and rep("X_coh_REF", m),
                "fair AIR at +3 dB": bool(all(e[m] <= .05 for e in out["air_err"].values()))}
    keys = {m: keys_of(m) for m in MODELS}
    best = {fam: min(v[m] for m in COMPARATORS) for fam, v in out["nmse"].items()}
    half = {fam: bool(out["nmse"][fam]["gmp"] <= .5 * best[fam]) for fam in best}
    p8 = {fam: bool(out["nmse"][fam]["gmp"] <= .8 * best[fam]) for fam in best}
    unique = [k for k, v in keys["gmp"].items() if v and not any(keys[c][k] for c in COMPARATORS)]
    if all(keys["gmp"].values()) and all(half.values()):
        cls = "SUCCEEDS"
    elif sum(p8.values()) >= 2 or unique:
        cls = "PARTIAL SUCCESS"
    else:
        cls = "FAILS"
    # held-out drive levels (S3 condition)
    lev_ok = {}
    for lev, v in out["levels"].items():
        if lev == 8 or v["full"]["X"] < .05:
            continue
        lev_ok[lev] = bool(abs(v["gmp"]["X"] - v["full"]["X"]) <= .25 * abs(v["full"]["X"])
                           and abs(v["gmp"]["X_coh"] - v["full"]["X_coh"]) <= .25 * abs(v["full"]["X_coh"]))
    s3 = bool(cls == "SUCCEEDS" and lev_ok and all(lev_ok.values()))
    # prior-pass P3 rule on the same fit (X only for the phase contrast)
    pk = {"dwell": rep("dwell_gain", "gmp"), "decluster": rep("decluster_Dref", "gmp"), "X_REF": rep("X_REF", "gmp"),
          "AIR": bool(all(e["gmp"] <= .05 for e in out["air_err"].values()))}
    prior = ("GMP_SUCCEEDS" if all(pk.values()) and all(half.values()) else
             ("PARTIAL_SUCCESS" if sum(pk.values()) >= 2 or sum(half.values()) >= 2 else "GMP_FAILS"))
    return {"classification": cls, "selected": json.loads(key), "key_contrasts_reproduced": keys, "gmp_nmse_le_half_best": half,
            "gmp_nmse_le_0.8_best": p8, "gmp_unique_contrasts": unique, "heldout_levels_reproduced": lev_ok, "stop_rule_S3": s3,
            "prior_pass_P3_rule": {"classification": prior, "keys": pk}, **out}


def figure(v):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    col = {"full": "#2a78d6", "gmp": "#4a3aa7", "static_mh": "#eb6834", "static_lti_mh": "#1baf7a", "dsh": "#eda100"}
    fig, axs = plt.subplots(1, 4, figsize=(13, 3.4))
    fams = list(v["nmse"])
    x = np.arange(len(fams))
    for j, m in enumerate(MODELS):
        axs[0].bar(x + (j - 1.5) * 0.2, [v["nmse"][f][m] for f in fams], 0.19, color=col[m], label=LABEL[m])
    axs[0].set_yscale("log")
    axs[0].set_xticks(x, fams)
    axs[0].set_ylabel("median test NMSE vs full model")
    axs[0].legend(fontsize=6.5, frameon=False)
    keys = ["dwell_gain", "decluster_Dref", "X_REF", "X_coh_REF", "X_FAST", "X_SLOW", "X_JUMP", "X_LONGRAMP"]
    xs = np.arange(len(keys))
    for j, m in enumerate(("full",) + MODELS):
        axs[1].plot(xs + (j - 2) * 0.08, [v["contrasts"][k][m] for k in keys], marker="o", lw=0, ms=4, color=col[m], label=LABEL[m])
    axs[1].set_xticks(xs, ["dwell Δgain", "decl. ΔD_ref", "X REF", "X_coh REF", "X fast", "X slow", "X steps", "X 900ns"], rotation=45, ha="right", fontsize=7)
    axs[1].axhline(0, color="#898781", lw=0.6)
    axs[1].set_title("held-out contrasts (test set means)", fontsize=8)
    mods = list(v["air_err"])
    xm = np.arange(len(mods))
    for j, m in enumerate(MODELS):
        axs[2].bar(xm + (j - 1.5) * 0.2, [v["air_err"][md][m] for md in mods], 0.19, color=col[m])
    axs[2].axhline(0.05, color="#898781", lw=0.6, ls="--")
    axs[2].set_xticks(xm, mods)
    axs[2].set_ylabel("mean |AIR − AIR_full| at +3 dB (bit)")
    levs = sorted(v["levels"])
    for m in ("full",) + MODELS:
        axs[3].plot(levs, [v["levels"][l_][m]["X_coh"] for l_ in levs], marker="o", ms=3.5, color=col[m], label=LABEL[m])
    axs[3].axhline(0, color="#898781", lw=0.6)
    axs[3].set_xlabel("drive level re P1dB (dB); +8 dB ≈ training level")
    axs[3].set_ylabel(r"$X_\mathrm{coh}$, REF case (held-out levels)")
    fig.suptitle(f"Part III: GMP ({v['selected']}) vs CW-derived surrogates — {v['classification']}", fontsize=9)
    fig.tight_layout()
    fig.savefig(OUT / "fig_model_comparison.png", dpi=160)
    plt.close(fig)


def summary(v, cfgdf):
    s = v["selected"]
    L = ["# Part III — one standard memory model (GMP) vs the CW-derived surrogates", "",
         f"**Classification: {v['classification']}** (rule: `00_plan.md`, section C; generated {v['generated']}, code {v['code_version'][:7]}).", "",
         f"**Selected model** (minimum family-averaged validation NMSE over {len(cfgdf)} grid points): nonlinear order K = {s['K']}, memory depth "
         f"{s['D']} µs (M = {s['M']} taps at 20 MHz; {s['M']} lagging lags, K_b = 2), ridge {s['lam']:g}; {n_coef(s['K'], s['M'])} complex coefficients "
         f"({2 * n_coef(s['K'], s['M'])} real parameters). Split (fixed before fitting): train fair/dwell/shuffle r0–r1 + phase seeds 20260705–06 + CW; "
         "validation r2–r3 + seed 20260704; test r4–r7 + seeds 20260701–03; held-out drive levels: revision P3 OP1 REF seeds 1–3 at −6…+6 dB re P1dB "
         "(+8 dB ≈ training level, reported separately).", "",
         "## Median test NMSE vs the full model", "", "| Family | GMP | all-zone static | all-zone LTI+static | single-slow-state | GMP ≤ 0.5 best | GMP ≤ 0.8 best |",
         "|---|---:|---:|---:|---:|---|---|"]
    for fam, r in v["nmse"].items():
        L.append(f"| {fam} | {r['gmp']:.4f} | {r['static_mh']:.4f} | {r['static_lti_mh']:.4f} | {r['dsh']:.4f} | {v['gmp_nmse_le_half_best'][fam]} | {v['gmp_nmse_le_0.8_best'][fam]} |")
    L += ["", "## Held-out contrasts (test-set means)", "", "| Quantity | Full model | GMP | all-zone static | all-zone LTI+static | single-slow-state |", "|---|---:|---:|---:|---:|---:|"]
    for k, d in v["contrasts"].items():
        L.append(f"| {k} | {d['full']:+.4f} | " + " | ".join(f"{d[m]:+.4f}" for m in MODELS) + " |")
    L += ["", "## Mean |AIR − AIR_full| at +3 dB (bit; archived noise)", "", "| Format | GMP | all-zone static | all-zone LTI+static | single-slow-state |", "|---|---:|---:|---:|---:|"]
    for mod, e in v["air_err"].items():
        L.append(f"| {mod} | {e['gmp']:.3f} | {e['static_mh']:.3f} | {e['static_lti_mh']:.3f} | {e['dsh']:.3f} |")
    L += ["", "## Held-out drive levels (REF case, seeds 1–3; X / X_coh)", "", "| level (dB re P1dB) | full | GMP | all-zone static | all-zone LTI+static | single-slow-state |", "|---|---|---|---|---|---|"]
    for lev, d in sorted(v["levels"].items()):
        L.append(f"| {lev:+d} | " + " | ".join(f"{d[m]['X']:+.3f} / {d[m]['X_coh']:+.3f}" for m in ("full",) + MODELS) + " |")
    L += ["", "## Decision", "", f"- Key contrasts reproduced (GMP): {v['key_contrasts_reproduced']['gmp']}.",
          f"- Key contrasts reproduced by comparators: { {m: v['key_contrasts_reproduced'][m] for m in COMPARATORS} }.",
          f"- Contrasts reproduced only by the GMP: {v['gmp_unique_contrasts']}.",
          f"- Held-out levels reproduced (X and X_coh within 25%, levels with X_full ≥ 0.05): {v['heldout_levels_reproduced']}.",
          f"- **Class: {v['classification']}.** Stop rule S3: {v['stop_rule_S3']}. Prior-pass P3 rule on this fit: {v['prior_pass_P3_rule']['classification']}.",
          "", "Grid (validation, family-mean NMSE; lower is better): see `gmp_config.csv`."]
    (OUT / "gmp_summary.md").write_text("\n".join(L) + "\n")


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    t0 = time.perf_counter()
    items, idx, th, key = fit()
    df = test(items, idx, th)
    v = evaluate(df, key)
    rev, dirty = git_state()
    v.update({"generated": now(), "code_version": rev, "code_dirty": dirty, "runtime_s": time.perf_counter() - t0})
    (OUT / "gmp_verdict.json").write_text(json.dumps(v, indent=1, default=float))
    cfgdf = pd.read_csv(OUT / "gmp_config.csv")
    figure(v)
    summary(v, cfgdf)
    print(json.dumps({k: v[k] for k in ("classification", "selected", "key_contrasts_reproduced", "gmp_nmse_le_half_best", "gmp_unique_contrasts",
                                        "heldout_levels_reproduced", "stop_rule_S3", "prior_pass_P3_rule", "runtime_s")}, indent=1, default=float))


if __name__ == "__main__":
    main()
