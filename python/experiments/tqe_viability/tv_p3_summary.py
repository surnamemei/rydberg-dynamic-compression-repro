"""P3 summary and classification (rules: results/tqe_viability/00_decision_rules.json, P3)."""
from __future__ import annotations

import json

import numpy as np
import pandas as pd

from tv_common import TV

MODELS = ("gmp", "static_mh", "static_lti_mh", "dsh")
LABEL = {"gmp": "GMP", "static_mh": "all-zone static", "static_lti_mh": "all-zone LTI+static", "dsh": "single-slow-state", "full": "full model"}


def realization_of(i):
    return int(i.split("_r")[1].split("_")[0])


def main():
    df = pd.read_csv(TV / "03_memory_model_per_waveform.csv")
    sel = json.loads(str(np.load(TV / "03_gmp_coefficients.npz", allow_pickle=True)["selected"]))
    test = df[df.split == "test"]
    rows = []
    # (1) NMSE by family
    nm = test[test.model != "full"].groupby(["family", "model"]).nmse_vs_full.median().unstack()
    for fam, r in nm.iterrows():
        best = min(r[m] for m in MODELS if m != "gmp")
        rows.append({"metric": "median NMSE vs full", "family": fam, **{LABEL[m]: r[m] for m in MODELS}, "gmp_over_best_comparator": r["gmp"] / best})
    # (2) dwell gain contrast (xi 4 minus 0.05, +3 dB)
    dw = test[(test.family == "dwell") & (test.p == 3)].copy()
    dw["r"] = dw.id.map(realization_of)
    contrasts = {}
    for m in ("full",) + MODELS:
        piv = dw[dw.model == m].pivot_table(index="r", columns="variant", values="gain_rel")
        contrasts.setdefault("dwell_gain", {})[m] = (piv["xi_4"] - piv["xi_0.05"])
    # (3) declustering D_ref contrast (+3 dB)
    sh = test[(test.family == "shuffle") & (test.p == 3)].copy()
    sh["r"] = sh.id.map(realization_of)
    for m in ("full",) + MODELS:
        piv = sh[sh.model == m].pivot_table(index="r", columns="variant", values="D_ref")
        contrasts.setdefault("decluster_Dref", {})[m] = (piv["BLOCK_SHUFFLED_CYCLES"] - piv["ORIGINAL"])
    # (4) phase X of the reference case (baseband definition, identical for every model)
    ph = df[df.family == "phase"]
    cw = ph[ph.case == "CW"].set_index("model").g_bb
    for case in ("REF", "FAST", "SLOW", "JUMP", "LONGRAMP"):
        sub = ph[(ph.case == case) & (ph.split == "test")]
        for m in ("full",) + MODELS:
            v = sub[sub.model == m].set_index("seed").g_bb
            contrasts.setdefault(f"X_{case}", {})[m] = 1 - v / cw[m]
    # (5) AIR at +3 dB per format (archived noise level)
    fa = test[(test.family == "fair") & (test.p == 3)]
    air_err = {}
    for mod, g in fa.groupby("mod"):
        full = g[g.model == "full"].set_index("id").AIR
        for m in MODELS:
            air_err.setdefault(mod, {})[m] = float(np.mean(np.abs(g[g.model == m].set_index("id").AIR - full)))
    for k, d in contrasts.items():
        full = d["full"].mean()
        rows.append({"metric": f"{k} (mean over test set)", "family": k.split("_")[0] if not k.startswith("X_") else "phase",
                     "full model": full, **{LABEL[m]: d[m].mean() for m in MODELS},
                     "gmp_rel_error": abs(d["gmp"].mean() - full) / abs(full) if full != 0 else np.nan})
    for mod, e in air_err.items():
        rows.append({"metric": f"mean |AIR - AIR_full| at +3 dB, {mod} (bit)", "family": "fair", **{LABEL[m]: e[m] for m in MODELS}})
    res = pd.DataFrame(rows)
    res.to_csv(TV / "03_memory_model_results.csv", index=False)
    # classification
    def reproduced(key, m="gmp"):
        full, mod = contrasts[key]["full"].mean(), contrasts[key][m].mean()
        return bool(np.sign(mod) == np.sign(full) and abs(mod - full) <= .25 * abs(full))
    keys = {"dwell gain contrast": reproduced("dwell_gain"), "declustering D_ref contrast": reproduced("decluster_Dref"),
            "phase X (REF)": reproduced("X_REF"), "fair AIR at +3 dB": bool(all(e["gmp"] <= .05 for e in air_err.values()))}
    nmse_ok = {fam: bool(r["gmp"] <= .5 * min(r[m] for m in MODELS if m != "gmp")) for fam, r in nm.iterrows()}
    comp_keys = {m: {"dwell": reproduced("dwell_gain", m), "decluster": reproduced("decluster_Dref", m), "X_REF": reproduced("X_REF", m)} for m in MODELS}
    if all(keys.values()) and all(nmse_ok.values()):
        cls = "GMP_SUCCEEDS"
    elif sum(keys.values()) >= 2 or sum(nmse_ok.values()) >= 2:
        cls = "PARTIAL_SUCCESS"
    else:
        cls = "GMP_FAILS"
    verdict = {"classification": cls, "selected_hyperparameters": sel, "key_contrasts_reproduced": keys, "nmse_half_of_best_comparator": nmse_ok,
               "comparators_key_contrasts": comp_keys}
    (TV / "03_memory_model_verdict.json").write_text(json.dumps(verdict, indent=1, default=float))
    val = pd.read_csv(TV / "03_memory_model_validation.csv")
    n_coef = int(val.iloc[0].n_coef)
    lines = ["# P3: one standard memory model (GMP) vs the CW-derived surrogates", "",
             f"**Classification: {cls}** (rules: `00_decision_rules.json`, P3).", "",
             f"**Model.** Generalized memory polynomial, selected on the validation set: nonlinear order K_a = {sel['Ka']}, aligned memory "
             f"M_a = {sel['Ma']} taps ({sel['Ma'] * 50} ns at 20 MHz), lagging cross terms K_b = 2 with L_b = {sel['Lb']} lags "
             f"({sel['Lb'] * 50 / 1000:.1f} us), ridge {sel['lam']:g}; harmonic-zone terms for zones 0 and 2; {n_coef} complex coefficients. "
             "Trained on fair r0-r1, dwell/shuffle r0-r1 and phase seeds 20260705-06 (+ CW); validated on r2-r3 and seed 20260704; "
             "tested on r4-r7 and seeds 20260701-03.", "",
             "**Key contrasts on the held-out test set** (full model vs models):", ""]
    lines.append("| Quantity | Full model | GMP | All-zone static | All-zone LTI+static | Single-slow-state |")
    lines.append("|---|---:|---:|---:|---:|---:|")
    for k, lab in (("dwell_gain", "Dwell gain contrast (xi 4 - 0.05, +3 dB)"), ("decluster_Dref", "Declustering D_ref contrast (+3 dB)"),
                   ("X_REF", "X, reference phase case (baseband)"), ("X_FAST", "X, fast rate"), ("X_SLOW", "X, slow rate"),
                   ("X_JUMP", "X, steps"), ("X_LONGRAMP", "X, 900 ns ramps")):
        d = contrasts[k]
        lines.append(f"| {lab} | {d['full'].mean():+.4f} | " + " | ".join(f"{d[m].mean():+.4f}" for m in MODELS) + " |")
    lines += ["", "**Median test NMSE vs the full model** (lower is better):", "", "| Family | GMP | All-zone static | All-zone LTI+static | Single-slow-state |",
              "|---|---:|---:|---:|---:|"]
    for fam, r in nm.iterrows():
        lines.append(f"| {fam} | {r['gmp']:.4f} | {r['static_mh']:.4f} | {r['static_lti_mh']:.4f} | {r['dsh']:.4f} |")
    lines += ["", "**Mean |AIR error| at +3 dB (bit), archived noise level:**", "", "| Format | GMP | All-zone static | All-zone LTI+static | Single-slow-state |",
              "|---|---:|---:|---:|---:|"]
    for mod, e in air_err.items():
        lines.append(f"| {mod} | {e['gmp']:.3f} | {e['static_mh']:.3f} | {e['static_lti_mh']:.3f} | {e['dsh']:.3f} |")
    lines += ["", f"Key contrasts reproduced by the GMP: {keys}.", f"GMP NMSE at most half of the best comparator's: {nmse_ok}.", "",
              "Floor: a 20-tap linear FIR in the same pipeline fits the linear small-signal model with NMSE ~0.01 (the 1% slow component of the "
              "small-signal response lies beyond the prespecified aligned memory)."]
    (TV / "03_memory_model_summary.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))
    print(json.dumps(verdict, indent=1, default=float))


if __name__ == "__main__":
    main()
