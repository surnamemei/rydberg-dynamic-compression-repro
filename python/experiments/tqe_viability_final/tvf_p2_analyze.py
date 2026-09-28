"""Part II analysis: LO / dressed-state resonance hypothesis (plan: results/tqe_viability_final/00_plan.md, section B).

First call (no 02_resonance/refinement.json): applies the plan's one-time refinement rule to the coarse scan and writes
refinement.json. If refinement is triggered, run `tvf_run.py IIr` and call this script again; the class is final only when
refinement.json exists and every listed refinement run is present.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from tvf_common import FIN, P_GRID_NS, SEEDS6, TV, git_state, lo_rabi_hz, now, raqr

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "revision2"))
import rev2_p2 as P2  # noqa: E402
from tvf_p1_analyze import paths as p1_paths, point as p1_point, wrap  # noqa: E402

OUT = FIN / "02_resonance"
LOS = (0.35, 0.425, 0.5)
P_NEW = (500, 450)
LO_COL = {0.35: "#8fb9ec", 0.425: "#4f8fdc", 0.5: "#1c4f9a"}
LO_MK = {0.35: "o", 0.425: "s", 0.5: "D"}
FC_MHZ = abs(raqr(0.5).Omega_c) / (2 * np.pi) / 1e6
FP_MHZ = abs(raqr(0.5).Omega_p) / (2 * np.pi) / 1e6
NORMS = {"f": ("IF (MHz)", "abs"), "r1": (r"$r_1=f/(\Omega_\mathrm{LO}/2\pi)$", "ratio"), "r2": (r"$r_2=f/(\Omega_\mathrm{LO}/4\pi)$", "ratio"),
         "d1": (r"$d_1=f-\Omega_\mathrm{LO}/2\pi$ (MHz)", "diff"), "d2": (r"$d_2=f-\Omega_\mathrm{LO}/4\pi$ (MHz)", "diff"),
         "s3": (r"$s_3=f/(\Omega_3/4\pi)$", "ratio"), "d3": (r"$d_3=f-\Omega_3/4\pi$ (MHz)", "diff")}
FEATURES = (("g_min", "g_CW", "min"), ("S_max", "S", "max"), ("Xcoh_max", "X_coh", "max"))


def point_files(lo, P):
    prior = TV / "p2" / f"LO{lo:g}"
    new = OUT / "runs" / f"LO{lo:g}"
    base = prior if (prior / f"P{P}_CW.npz").exists() else new
    return base / f"P{P}_CW.npz", base / f"P{P}_REF_s{SEEDS6[0]}.npz", "prior-pass P2" if base == prior else "campaign"


def selectivity(f, g, delta_max=1.31, step=0.01):
    lf = np.log(g)
    out = np.full(len(f), np.nan)
    for i, fi in enumerate(f):
        fd = fi + np.arange(-delta_max, delta_max + step / 2, step)
        m = (fd >= f.min()) & (fd <= f.max())
        out[i] = np.max(np.abs(np.interp(fd[m], f, lf) - lf[i]))
    return out


def scan(refinement):
    rows = []
    for lo in LOS:
        f_lo = lo_rabi_hz(lo) / 1e6
        f3 = np.hypot(FC_MHZ, f_lo)
        periods = sorted(set(P_GRID_NS) | set(P_NEW) | set(refinement.get(f"{lo:g}", [])), reverse=True)
        for P in periods:
            cw, ref, src = point_files(lo, P)
            if not (cw.exists() and ref.exists()):
                rows.append({"A_LO_Vpm": lo, "P_ns": P, "IF_MHz": 1e3 / P, "status": "missing"})
                continue
            c, q = P2.metrics(cw), P2.metrics(ref)
            zc, zq = np.load(cw), np.load(ref)
            mc, mq = P2.window(zc["t"]), P2.window(zq["t"])
            f = 1e3 / P
            dq, dc = json.loads(str(zq["diag"])), json.loads(str(zc["diag"]))
            rows.append({"A_LO_Vpm": lo, "P_ns": P, "IF_MHz": f, "f": f, "status": "ok", "source": src, "LO_Rabi_MHz": f_lo, "Omega3_MHz": f3,
                         "amp_Vpm": float(zc["amp"]), "g_CW": c["g_med"],
                         "X": 1 - q["g_med"] / c["g_med"], "X_coh": 1 - q["g_coh"] / c["g_coh"], "X_mag": 1 - q["g_mag"] / c["g_mag"],
                         "X_rms": 1 - q["g_rms"] / c["g_rms"], "lag_rad": wrap(q["phase_mean_rad"] - c["phase_mean_rad"]),
                         "jitter_rad": q["phase_circ_sd_rad"],
                         "env_red": 1 - np.mean(np.abs(zq["z_full"][mq])) / np.mean(np.abs(zc["z_full"][mc])),
                         "r1": f / f_lo, "r2": f / (f_lo / 2), "d1": f - f_lo, "d2": f - f_lo / 2, "s3": f / (f3 / 2), "d3": f - f3 / 2,
                         "thermal_trace_dev": max(dq["thermal_trace_dev"], dc["thermal_trace_dev"]),
                         "min_class_eig": min(dq["min_eig_min"], dc["min_eig_min"])})
    df = pd.DataFrame(rows)
    ok = df[df.status == "ok"].copy()
    ok["S"] = np.nan
    for lo, g in ok.groupby("A_LO_Vpm"):
        g = g.sort_values("IF_MHz")
        ok.loc[g.index, "S"] = selectivity(g.IF_MHz.values, g.g_CW.values)
    return df, ok


def features(ok):
    feats = {}
    for lo, g in ok.groupby("A_LO_Vpm"):
        g = g.sort_values("IF_MHz").reset_index(drop=True)
        feats[lo] = {}
        for name, col, kind in FEATURES:
            i = int(g[col].idxmin() if kind == "min" else g[col].idxmax())
            edge = i in (0, len(g) - 1)
            feats[lo][name] = None if edge else {n: float(g.loc[i, n]) for n in NORMS} | {"value": float(g.loc[i, col])}
            feats[lo][name + "_edge"] = edge
    return feats


def aligned(feats, name, norm):
    locs = [feats[lo][name] for lo in LOS]
    if any(x is None for x in locs):
        return False, None
    v = np.array([x[norm] for x in locs])
    f = np.array([x["f"] for x in locs])
    if NORMS[norm][1] == "ratio":
        sp, spf = float(v.max() / v.min()) if v.min() > 0 else np.inf, float(f.max() / f.min())
        return bool(sp <= 1.08 and spf >= 1.25), {"spread": sp, "spread_f": spf}
    rg, rgf = float(v.max() - v.min()), float(f.max() - f.min())
    return bool(rg <= 0.4 and rgf >= 1.0), {"range_MHz": rg, "range_f_MHz": rgf}


def collapse(ok, col, norm, log=False):
    cur = {}
    for lo, g in ok.groupby("A_LO_Vpm"):
        g = g.sort_values(norm)
        y = np.log(g[col].values) if log else g[col].values
        cur[lo] = (g[norm].values, y)
    lo_, hi_ = max(a.min() for a, _ in cur.values()), min(a.max() for a, _ in cur.values())
    if hi_ <= lo_:
        return np.nan
    grid = np.linspace(lo_, hi_, 200)
    stack = np.array([np.interp(grid, a, v) for a, v in cur.values()])
    return float(np.sqrt(np.mean(np.var(stack, axis=0))))


def decide_refinement(ok):
    out = {"decided": now(), "rule": "00_plan.md section B: X_coh max >= 0.2, not at an edge, adjacent spacing > 0.2 MHz -> up to 4 IFs "
           "(integer-ns periods nearest to +-0.1 and +-0.2 MHz)", "periods": {}, "per_LO": {}}
    for lo, g in ok.groupby("A_LO_Vpm"):
        g = g.sort_values("IF_MHz").reset_index(drop=True)
        i = int(g.X_coh.idxmax())
        f0, x0 = float(g.IF_MHz[i]), float(g.X_coh[i])
        edge = i in (0, len(g) - 1)
        sp = max(float(g.IF_MHz[i] - g.IF_MHz[i - 1]) if i > 0 else 0, float(g.IF_MHz[i + 1] - g.IF_MHz[i]) if i < len(g) - 1 else 0)
        trig = bool(x0 >= .2 and not edge and sp > .2)
        have = set(g.P_ns.astype(int))
        new = []
        if trig:
            for df_ in (-0.2, -0.1, 0.1, 0.2):
                P = int(round(1e3 / (f0 + df_)))
                if P not in have and P not in new:
                    new.append(P)
            new = new[:4]
            out["periods"][f"{lo:g}"] = new
        out["per_LO"][f"{lo:g}"] = {"Xcoh_max": x0, "at_IF_MHz": f0, "edge": edge, "max_adjacent_spacing_MHz": sp, "triggered": trig, "new_periods": new}
    out["triggered"] = bool(out["periods"])
    return out


def seq_sd():
    """Sequence-to-sequence SD at 5 MHz, signal/LO 0.358 (Part I, six sequences) for C1 (0.5), C5 (0.425), C6 (0.35)."""
    plan = json.loads((FIN / "01_matching" / "plan.json").read_text())
    out = {}
    for cfg, lo in (("C1", 0.5), ("C5", 0.425), ("C6", 0.35)):
        p = [x for x in plan["points"] if x["config"] == cfg and x["rule"] == "matched_signal_to_LO" and x["target"] == 0.358][0]
        _, per = p1_point(p)
        out[lo] = {k: float(np.std([per[s_][k] for s_ in per], ddof=1)) for k in ("X", "X_coh", "X_mag")}
    return out


def oscillations():
    obs = {0.5: {"power steps": (1 / 1.62, 1 / 1.60), "phase steps": (1 / 1.74, 1 / 1.63)}, 0.35: {"power steps": (1 / 0.560,), "phase steps": (1 / 0.558,)}}
    rows = []
    for lo, kinds in obs.items():
        fl = lo_rabi_hz(lo) / 1e6
        f3 = np.hypot(FC_MHZ, fl)
        cand = {"Omega_LO/2pi": fl, "Omega_LO/4pi": fl / 2, "Omega_c/2pi": FC_MHZ, "Omega_p/2pi": FP_MHZ, "Omega_3/2pi": f3, "Omega_3/4pi": f3 / 2,
                "|f_IF - Omega_LO/2pi|": abs(5 - fl), "|f_IF - Omega_LO/4pi|": abs(5 - fl / 2), "|2 f_IF - Omega_LO/2pi|": abs(10 - fl)}
        for kind, fs in kinds.items():
            for name, val in cand.items():
                rel = [abs(val - f_) / f_ for f_ in fs]
                rows.append({"A_LO_Vpm": lo, "observation": kind, "observed_MHz": "/".join(f"{f_:.3f}" for f_ in fs), "candidate": name,
                             "candidate_MHz": val, "min_rel_diff": min(rel), "match_10pct": bool(min(rel) <= .10)})
    df = pd.DataFrame(rows)
    cons = []
    for (name, kind), g in df.groupby(["candidate", "observation"]):
        if g.match_10pct.all() and len(g) == 2:
            cons.append(f"{name} ({kind})")
    return df, cons


def main():
    refp = OUT / "refinement.json"
    ref = json.loads(refp.read_text()) if refp.exists() else None
    df, ok = scan(ref["periods"] if ref else {})
    missing = df[df.status == "missing"]
    if ref is None:
        ref = decide_refinement(ok)
        refp.write_text(json.dumps(ref, indent=1, default=float))
        if ref["triggered"]:                       # provisional until the refinement runs exist (tvf_run.py IIr)
            df, ok = scan(ref["periods"])
            missing = df[df.status == "missing"]
    final = bool(len(missing) == 0)
    feats = features(ok)
    align, coll = {}, {}
    for norm in NORMS:
        if norm == "f":
            continue
        align[norm] = {name: aligned(feats, name, norm) for name, _, _ in FEATURES}
        coll[norm] = {"ln_g_CW": collapse(ok, "g_CW", norm, log=True) / collapse(ok, "g_CW", "f", log=True),
                      "X_coh": collapse(ok, "X_coh", norm) / collapse(ok, "X_coh", "f")}
    clear = [n for n in align if align[n]["Xcoh_max"][0] and (align[n]["g_min"][0] or align[n]["S_max"][0])
             and coll[n]["ln_g_CW"] <= .5 and coll[n]["X_coh"] <= .5]
    sugg = [n for n in align if any(a[0] for a in align[n].values()) or coll[n]["ln_g_CW"] <= .75 or coll[n]["X_coh"] <= .75]
    cls = "CLEAR SCALING" if clear else ("SUGGESTIVE ONLY" if sugg else "NO SIMPLE SCALING")
    # prior-pass P2 rule (r1/r2; features g_min, X_max (median X), S_max on raw curves incl. edges)
    pf = {}
    for lo, g in ok.groupby("A_LO_Vpm"):
        pf[lo] = {"g_min": g.loc[g.g_CW.idxmin()], "X_max": g.loc[g.X.idxmax()], "S_max": g.loc[g.S.idxmax()]}
    pv = {}
    for ax in ("r1", "r2"):
        okf = [ft for ft in ("g_min", "X_max", "S_max") if max(pf[lo][ft][ax] for lo in LOS) / min(pf[lo][ft][ax] for lo in LOS) <= 1.08
               and max(pf[lo][ft]["IF_MHz"] for lo in LOS) / min(pf[lo][ft]["IF_MHz"] for lo in LOS) >= 1.25]
        red = {"ln_g": collapse(ok, "g_CW", ax, True) / collapse(ok, "g_CW", "f", True), "X": collapse(ok, "X", ax) / collapse(ok, "X", "f")}
        pv[ax] = {"features_aligned": okf, "rms_ratio": red}
    pclear = any(len(v["features_aligned"]) >= 2 and all(r <= .5 for r in v["rms_ratio"].values()) for v in pv.values())
    psugg = any(len(v["features_aligned"]) >= 1 or any(r <= .75 for r in v["rms_ratio"].values()) for v in pv.values())
    pcls = "CLEAR_RESONANCE_SCALING" if pclear else ("SUGGESTIVE_ONLY" if psugg else "NO_SIMPLE_SCALING")
    osc, consistent = oscillations()
    osc.to_csv(OUT / "oscillation_comparison.csv", index=False)
    sd = seq_sd()
    df_out = df.merge(ok[["A_LO_Vpm", "P_ns", "S"]], on=["A_LO_Vpm", "P_ns"], how="left")
    df_out.to_csv(OUT / "resonance_scan.csv", index=False)
    rev, dirty = git_state()
    verdict = {"generated": now(), "code_version": rev, "code_dirty": dirty, "final": final, "classification": cls,
               "normalizations_meeting_CLEAR": clear, "normalizations_meeting_SUGGESTIVE": sugg,
               "features": {f"{lo:g}": feats[lo] for lo in LOS},
               "alignment": {n: {k: {"aligned": v[0], **(v[1] or {})} for k, v in a.items()} for n, a in align.items()},
               "collapse_ratio_vs_absolute_f": coll, "refinement": ref, "sequence_sd_at_5MHz_ratio0.358": sd,
               "prior_pass_P2_rule": {"classification": pcls, "per_normalization": pv},
               "oscillation_consistent_candidates": consistent, "Omega_c_MHz": FC_MHZ, "Omega_p_MHz": FP_MHZ}
    (OUT / "resonance_verdict.json").write_text(json.dumps(verdict, indent=1, default=float))
    figures(ok, sd)
    summary(ok, verdict, osc)
    print(json.dumps({k: v for k, v in verdict.items() if k not in ("features",)}, indent=1, default=float))


def figures(ok, sd):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    def draw(ax, norm, col, err=None):
        for lo, g in ok.groupby("A_LO_Vpm"):
            g = g.sort_values(norm)
            ax.plot(g[norm], g[col], "-", marker=LO_MK[lo], ms=3.2, lw=1.0, color=LO_COL[lo], label=f"$A_\\mathrm{{LO}}$ = {lo:g} V/m")
            if err is not None:
                e = sd[lo][err]
                ax.fill_between(g[norm], g[col] - e, g[col] + e, color=LO_COL[lo], alpha=0.15, lw=0)
    # IF scan
    fig, axs = plt.subplots(4, 1, figsize=(7.5, 9), sharex=True)
    for ax, (col, lab, err) in zip(axs, (("g_CW", "CW gain / small-signal", None), ("S", "selectivity S (±1.31 MHz)", None),
                                        ("X_coh", r"$X_\mathrm{coh}$ (±1 seq. SD)", "X_coh"), ("X_mag", r"$X_\mathrm{mag}$ (±1 seq. SD)", "X_mag"))):
        draw(ax, "IF_MHz", col, err)
        for lo in LOS:
            fl = lo_rabi_hz(lo) / 1e6
            ax.axvline(fl / 2, color=LO_COL[lo], lw=0.7, ls=":")
            ax.axvline(fl, color=LO_COL[lo], lw=0.7, ls="--")
        if col.startswith("X"):
            ax.axhline(0, color="#898781", lw=0.6)
        ax.set_ylabel(lab)
    axs[-1].set_xlabel("IF (MHz)   (dotted: $\\Omega_\\mathrm{LO}/4\\pi$; dashed: $\\Omega_\\mathrm{LO}/2\\pi$)")
    axs[0].legend(fontsize=7, frameon=False)
    fig.suptitle("Part II IF scan at matched signal/LO = 0.358 (one sequence per IF)", fontsize=9)
    fig.tight_layout()
    fig.savefig(OUT / "fig_if_scan.png", dpi=160)
    plt.close(fig)
    norms = [n for n in NORMS if n != "f"]
    for fname, cols in (("fig_scaled_if.png", (("g_CW", "CW gain / small-signal"), ("S", "selectivity S"))),
                        ("fig_phase_effect_vs_scaled_if.png", (("X_coh", r"$X_\mathrm{coh}$"), ("X_mag", r"$X_\mathrm{mag}$")))):
        fig, axs = plt.subplots(len(cols), len(norms), figsize=(2.6 * len(norms), 2.4 * len(cols)), sharey="row")
        for i, (col, lab) in enumerate(cols):
            for j, n in enumerate(norms):
                ax = axs[i, j]
                draw(ax, n, col, col if col.startswith("X") else None)
                ax.axvline(1.0 if NORMS[n][1] == "ratio" else 0.0, color="#898781", lw=0.6, ls=":")
                if col.startswith("X"):
                    ax.axhline(0, color="#898781", lw=0.6)
                if i == len(cols) - 1:
                    ax.set_xlabel(NORMS[n][0], fontsize=8)
                if j == 0:
                    ax.set_ylabel(lab)
        axs[0, 0].legend(fontsize=6.5, frameon=False)
        fig.tight_layout()
        fig.savefig(OUT / fname, dpi=150)
        plt.close(fig)


def summary(ok, v, osc):
    L = ["# Part II — LO / dressed-state resonance hypothesis", "",
         f"**Classification: {v['classification']}**{'' if v['final'] else ' (PROVISIONAL: refinement runs missing)'} "
         f"(rule: `00_plan.md`, section B; generated {v['generated']}, code {v['code_version'][:7]}).", "",
         "Matched state: signal/LO = 0.358 at every IF and LO field (standard probe); CW + one reference sequence (20260701) per IF; "
         f"{ok.groupby('A_LO_Vpm').size().to_dict()} IFs per LO field. Sequence-to-sequence SD at 5 MHz (Part I, six sequences): "
         + "; ".join(f"LO {lo:g}: X {s['X']:.3f}, X_coh {s['X_coh']:.3f}, X_mag {s['X_mag']:.3f}" for lo, s in v["sequence_sd_at_5MHz_ratio0.358"].items()) + ".", "",
         "## Feature locations (grid points; edge features ineligible)", "", "| LO (V/m) | Ω_LO/2π (MHz) | feature | IF (MHz) | r1 | r2 | d1 | d2 | s3 | d3 | value |",
         "|---|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|"]
    for lo, fts in v["features"].items():
        for name, _, _ in FEATURES:
            x = fts[name]
            if x is None:
                L.append(f"| {lo} | {lo_rabi_hz(float(lo)) / 1e6:.3f} | {name} | edge | | | | | | | |")
            else:
                L.append(f"| {lo} | {lo_rabi_hz(float(lo)) / 1e6:.3f} | {name} | {x['f']:.3f} | {x['r1']:.3f} | {x['r2']:.3f} | {x['d1']:+.3f} | "
                         f"{x['d2']:+.3f} | {x['s3']:.3f} | {x['d3']:+.3f} | {x['value']:.3f} |")
    L += ["", "## Alignment and collapse", "", "| normalization | g_min aligned | S_max aligned | X_coh-max aligned | collapse ratio ln g_CW | collapse ratio X_coh |",
          "|---|---|---|---|---:|---:|"]
    for n, a in v["alignment"].items():
        c = v["collapse_ratio_vs_absolute_f"][n]
        L.append(f"| {n} | {a['g_min']['aligned']} | {a['S_max']['aligned']} | {a['Xcoh_max']['aligned']} | {c['ln_g_CW']:.2f} | {c['X_coh']:.2f} |")
    L += ["", f"Normalizations meeting CLEAR: {v['normalizations_meeting_CLEAR']}; meeting SUGGESTIVE: {v['normalizations_meeting_SUGGESTIVE']}.",
          f"Prior-pass P2 rule: {v['prior_pass_P2_rule']['classification']}.", "",
          "## II-C transient oscillation frequencies vs candidate scales (descriptive; match = within 10%)", "",
          "| LO (V/m) | observation | observed (MHz) | candidate | candidate (MHz) | rel. diff | match |", "|---|---|---|---|---:|---:|---|"]
    for r in osc.itertuples():
        L.append(f"| {r.A_LO_Vpm:g} | {r.observation} | {r.observed_MHz} | {r.candidate} | {r.candidate_MHz:.3f} | {r.min_rel_diff:.2f} | {'yes' if r.match_10pct else ''} |")
    L += ["", f"Candidates matching at both LO fields (consistent): {v['oscillation_consistent_candidates'] or 'none'}. No match is called causal.", "",
          f"Refinement: {json.dumps(v['refinement']['per_LO'], default=float)}."]
    (OUT / "resonance_summary.md").write_text("\n".join(L) + "\n")


if __name__ == "__main__":
    main()
