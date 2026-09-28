"""Runs of the TQE-viability pass (decision rules: results/tqe_viability/00_decision_rules.json). Resumable.

  p1  fair cross-configuration matching (matched CW gain; matched signal/LO ratio)
  p2  resonance scan: 3 LO fields x 29 IFs at signal/LO 0.358 (CW + one reference sequence)
  p4  re-run of the 192 fair-waveform jobs, storing each model's noiseless receiver baseband
Usage: python3 tv_run.py [p1|p2|p4 ...]
"""
from __future__ import annotations

import json
import sys
import time

import numpy as np

from tv_common import (DT, KERNELS, KERNEL_OP1, OMEGA_P_WEAK, P_GRID_NS, REV2, SEEDS3, TV, Kernel, full_diag, phase_run, raqr, s)

REF_SPEC = {"kind": "qpsk", "Ts": 1e-6, "Tp": 300e-9}
CONFIGS = {"IF5": (5e6, None), "IF10": (10e6, None), "IF15": (15e6, None), "IF5_weak": (5e6, OMEGA_P_WEAK)}
TABLES = {"IF5": "OP1_IF5.npz", "IF10": "OP1_IF10_ext.npz", "IF15": "OP1_IF15.npz", "IF5_weak": "OP1_IF5_weakprobe_ext.npz"}
A_LO = 0.5


def log(msg):
    print(msg, flush=True)
    with open(TV / "progress.txt", "a") as fh:
        fh.write(msg + "\n")


def kernel_for(cfg_name):
    if cfg_name == "IF5_weak":
        k = np.load(REV2 / "tables" / "kernel_weakprobe.npz")
        return Kernel(k["hR"], k["hI"], np.load(REV2 / "tables" / TABLES[cfg_name])["means"][0])
    return Kernel.from_npz(KERNEL_OP1)


def target_amplitude(table, target):
    """First crossing of g_CW/g_small = target on the CW table (linear interpolation in dB); None if none."""
    a, g = table["amps"][1:], table["gain_db"][1:]
    tdb = 20 * np.log10(target)
    hit = np.flatnonzero(g <= tdb)
    if not len(hit) or hit[0] == 0:
        return None
    j = hit[0]
    return float(np.interp(tdb, [g[j], g[j - 1]], [a[j], a[j - 1]]))


def p1_plan():
    rules = json.loads((TV / "00_decision_rules.json").read_text())["P1_fair_matching"]
    plan, infeasible = [], []
    for cname in CONFIGS:
        tab = np.load(REV2 / "tables" / TABLES[cname])
        for tgt in rules["A_matched_CW_gain"]["targets"]:
            a = target_amplitude(tab, tgt)
            ok = a is not None and a / A_LO <= 1.0 and a <= float(tab["amps"][-1])
            (plan if ok else infeasible).append({"config": cname, "rule": "matched_CW_gain", "target": tgt, "amp_Vpm": a})
        for ratio in rules["B_matched_signal_to_LO"]["ratios"]:
            plan.append({"config": cname, "rule": "matched_signal_to_LO", "target": ratio, "amp_Vpm": ratio * A_LO})
    return plan, infeasible


def stage_p1():
    plan, infeasible = p1_plan()
    (TV / "p1").mkdir(parents=True, exist_ok=True)
    (TV / "p1" / "plan.json").write_text(json.dumps({"run": plan, "infeasible": infeasible}, indent=1, default=float))
    log(f"p1: {len(plan)} matched points to run; infeasible: {[(x['config'], x['target']) for x in infeasible]}")
    for pt in plan:
        f_if, op = CONFIGS[pt["config"]]
        cfg = raqr(A_LO, op)
        kern = kernel_for(pt["config"])
        tag = f"{pt['config']}_{pt['rule']}_{pt['target']:g}"
        d = TV / "p1" / "runs"
        phase_run(d / f"{tag}_CW.npz", {"kind": "cw"}, cfg, f_if, pt["amp_Vpm"], kern)
        for seed in SEEDS3:
            phase_run(d / f"{tag}_REF_s{seed}.npz", {**REF_SPEC, "seed": seed}, cfg, f_if, pt["amp_Vpm"], kern)
        log(f"p1 {tag}: a = {pt['amp_Vpm']:.4f} V/m done")


def stage_p2():
    for lo in (0.35, 0.425, 0.5):
        cfg = raqr(lo)
        kern = Kernel.from_npz(KERNELS[lo])
        amp = 0.358 * lo
        for P in P_GRID_NS:
            f_if = 1e9 / P
            d = TV / "p2" / f"LO{lo:g}"
            phase_run(d / f"P{P}_CW.npz", {"kind": "cw"}, cfg, f_if, amp, kern)
            phase_run(d / f"P{P}_REF_s{SEEDS3[0]}.npz", {**REF_SPEC, "seed": SEEDS3[0]}, cfg, f_if, amp, kern)
        log(f"p2 LO {lo:g} V/m: {len(P_GRID_NS)} IFs done")


def stage_p4():
    """Re-run of results/final_validation/regen_stage05 (same waveforms, seeds, noise draw, receiver) storing noiseless baseband."""
    sys.path.insert(0, str(TV.parents[1] / "python" / "experiments" / "final_validation"))
    import regen_stage05 as G
    import models06 as M
    from common import OUT as OUT06
    import run06
    f1 = M.Controls.from_npz(OUT06 / "artifacts" / "controls_Nd4001.npz")
    mh = M.ControlsMH(OUT06 / "artifacts" / "controls_MH_Nd4001.npz", f1)
    e1 = f1.e1db
    old = s.old
    out_dir = TV / "p4" / "jobs"
    out_dir.mkdir(parents=True, exist_ok=True)
    cache = {}
    t_start = time.perf_counter()
    jobs = G.plan()
    for i, (mod, r, seed, p) in enumerate(jobs, 1):
        jid = f"regen_{mod}_s{seed}_p{p:+03d}_Nd4001"
        out = out_dir / f"{jid}.npz"
        if out.exists():
            continue
        if (mod, seed) not in cache:
            cache.clear()
            cache[(mod, seed)] = s.waveform(mod, seed)
        tx, labels, alphabet = cache[(mod, seed)]
        env = e1 * 10 ** (p / 20) * tx
        n = len(env)
        field = old.raqr.A_LO + env * np.exp(2j * np.pi * old.IF * (np.arange(n) * DT))
        y_atomic, diag = full_diag(field, old.raqr)
        outs = {"atomic": y_atomic, "static": f1.static(env), "static_lti": f1.static_lti(env), "static_mh": mh.static_mh(env),
                "static_lti_mh": mh.static_lti_mh(env), "linear": f1.linear(env), "dsh": M.dsh(env, run06.TAU_S, mh)}
        tt = np.arange(n) / s.FS
        mix = np.exp(-2j * np.pi * old.IF * tt)
        bbs = {m: old.sosfilt(old.sos, (y - old.dc) * mix)[::s.UPSAMPLE] for m, y in outs.items()}
        noise = np.random.default_rng(seed + 9991).normal(0, s.SIGMA, len(tx))
        bb_noise = old.sosfilt(old.sos, noise * mix)[::s.UPSAMPLE]
        # validation of the decomposition at c = 1: original receiver path for the atomic output
        air_orig = s.split_eval(s.extract(y_atomic + noise, mod), labels, alphabet, mod,
                                json.loads((G.OUT / "jobs" / f"{jid}.json").read_text())["ce_timing_offset"])[0]
        np.savez_compressed(out, **{f"bb_{m}": v.astype(np.complex128) for m, v in bbs.items()}, bb_noise=bb_noise,
                            air_atomic_orig_c1=air_orig, diag=json.dumps(diag), mod=mod, r=r, seed=seed, p=p)
        el = time.perf_counter() - t_start
        if i % 8 == 0 or i == len(jobs):
            log(f"p4 [{i}/{len(jobs)}] {mod} r{r} {p:+d} dB, elapsed {el:.0f} s")
    log("p4 done")


def main():
    TV.mkdir(parents=True, exist_ok=True)
    for st in sys.argv[1:] or ["p1", "p2", "p4"]:
        {"p1": stage_p1, "p2": stage_p2, "p4": stage_p4}[st]()


if __name__ == "__main__":
    main()
