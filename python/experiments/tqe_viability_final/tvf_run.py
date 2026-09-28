"""Campaign runner (plan: results/tqe_viability_final/00_plan.md). Resumable; logs every job (Step 2) and the run manifest.

  I   Part I matched states (01_matching/plan.json): CW + REF seeds 20260701-06; prior 3-seed points get seeds 04-06
  II  Part II new IFs (P = 500, 450 ns) at signal/LO 0.358, 3 LO fields: CW + REF seed 20260701
  IIr Part II refinement IFs listed in 02_resonance/refinement.json (written by the Part II analysis if the plan's rule fires)
  V   Part V-A prefix runs: phase case (55/110/165/220 us) and long-dwell case (105/210/420/840 us)
Usage: python3 tvf_run.py STAGE [STAGE ...]
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

from tvf_common import (CONFIGS, DT, FIN, KERNEL_OP1, PRIOR_PID_FILE, R, REF_SPEC, ROOT, RP, SEEDS6, T0, TV, Kernel, KERNELS, Progress, full_diag,
                        kernel, manifest_row, now, phase_run, raqr, raqr_of)

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "revision2"))
import rev2_p2 as P2  # noqa: E402

COMMAND = "python3 python/experiments/tqe_viability_final/tvf_run.py " + " ".join(sys.argv[1:])
P_NEW = (500, 450)
LOS = (0.35, 0.425, 0.5)


def wait_for_prior_runner(prog=None):
    """Never two GPU jobs at once: wait until the prior-pass runner (tv_run.py) has exited."""
    if PRIOR_PID_FILE is None:
        return
    try:
        pid = int(PRIOR_PID_FILE.read_text().strip())
    except Exception:
        return
    first = True
    while Path(f"/proc/{pid}").exists():
        if first:
            print(f"{now()} waiting for prior-pass runner PID {pid} to exit", flush=True)
            (FIN / "00_progress.txt").write_text(f"[0/?] waiting for the prior-pass runner (PID {pid}, P4 reruns) to exit before using the GPU\nupdated: {now()}\n")
            first = False
        time.sleep(15)


# ------------------------------------------------------------------ job builders
def jobs_part1():
    plan = json.loads((FIN / "01_matching" / "plan.json").read_text())
    out = []
    for p in plan["points"]:
        cfg = p["config"]
        d = FIN / "01_matching" / "runs"
        seeds = [s_ for s_ in SEEDS6 if s_ not in p["seeds_prior"]]
        base = {"stage": "I", "config": cfg, "matching_rule": p["rule"], "parameter_value": p["target"], "amp": p["amp_Vpm"],
                "source": "python/experiments/tqe_viability_final/tvf_run.py (phase_run from tqe_viability/tv_common.py)"}
        if not p["prior_tag"]:
            out.append({**base, "job_id": f"{p['tag']}_CW", "spec": {"kind": "cw"}, "out": d / f"{p['tag']}_CW.npz", "seed": ""})
        for s_ in seeds:
            out.append({**base, "job_id": f"{p['tag']}_REF_s{s_}", "spec": {**REF_SPEC, "seed": s_}, "out": d / f"{p['tag']}_REF_s{s_}.npz",
                        "seed": s_, "cw_path": (TV / "p1" / "runs" / f"{p['prior_tag']}_CW.npz") if p["prior_tag"] else d / f"{p['tag']}_CW.npz"})
    return out


def jobs_part2(periods):
    out = []
    for lo in LOS:
        d = FIN / "02_resonance" / "runs" / f"LO{lo:g}"
        for P in periods:
            base = {"stage": "II", "config": f"LO{lo:g}", "matching_rule": "matched_signal_to_LO_0.358", "parameter_value": f"IF {1e3 / P:.4f} MHz",
                    "lo": lo, "P": P, "amp": 0.358 * lo, "source": "python/experiments/tqe_viability_final/tvf_run.py (phase_run)"}
            out.append({**base, "job_id": f"LO{lo:g}_P{P}_CW", "spec": {"kind": "cw"}, "out": d / f"P{P}_CW.npz", "seed": ""})
            out.append({**base, "job_id": f"LO{lo:g}_P{P}_REF_s{SEEDS6[0]}", "spec": {**REF_SPEC, "seed": SEEDS6[0]},
                        "out": d / f"P{P}_REF_s{SEEDS6[0]}.npz", "seed": SEEDS6[0], "cw_path": d / f"P{P}_CW.npz"})
    return out


def jobs_part5():
    out = []
    for T in (55e-6, 110e-6, 165e-6, 220e-6):
        out.append({"stage": "V", "job_id": f"VA_phase_T{T * 1e6:.0f}us", "config": "C1", "matching_rule": "a_H", "parameter_value": f"T={T * 1e6:.0f}us",
                    "kind": "va_phase", "T": T, "out": FIN / "05_numerical" / "runs" / f"VA_phase_T{T * 1e6:.0f}us.npz", "seed": SEEDS6[0],
                    "source": "python/experiments/tqe_viability_final/tvf_run.py (run_va_phase)"})
    for T in (105e-6, 210e-6, 420e-6, 840e-6):
        out.append({"stage": "V", "job_id": f"VA_dwell_T{T * 1e6:.0f}us", "config": "C1", "matching_rule": "dwell xi=4 r0 +3dB",
                    "parameter_value": f"T={T * 1e6:.0f}us", "kind": "va_dwell", "T": T,
                    "out": FIN / "05_numerical" / "runs" / f"VA_dwell_T{T * 1e6:.0f}us.npz", "seed": "r0",
                    "source": "python/experiments/tqe_viability_final/tvf_run.py (run_va_dwell)"})
    return out


# ------------------------------------------------------------------ executors
def run_phase_job(job):
    if job["stage"] == "I":
        cfg = job["config"]
        f_if = CONFIGS[cfg][0]
        phase_run(job["out"], job["spec"], raqr_of(cfg), f_if, job["amp"], kernel(cfg))
        job["cfg_json"] = {"config": cfg, "IF_Hz": f_if, "A_LO": CONFIGS[cfg][1], "Omega_p": CONFIGS[cfg][2], "amp_Vpm": job["amp"],
                           "spec": job["spec"], "rule": job["matching_rule"], "target": job["parameter_value"]}
    else:
        lo, P = job["lo"], job["P"]
        f_if = 1e9 / P
        phase_run(job["out"], job["spec"], raqr(lo), f_if, job["amp"], Kernel.from_npz(KERNELS[lo]))
        job["cfg_json"] = {"A_LO": lo, "IF_Hz": f_if, "amp_Vpm": job["amp"], "spec": job["spec"], "rule": "matched signal/LO 0.358"}


def latest_phase(job):
    z = np.load(job["out"])
    diag = json.loads(str(z["diag"]))
    if job["spec"]["kind"] == "cw":
        c = P2.metrics(job["out"])
        return f"{job['job_id']}: g_CW = {c['g_med']:.3f}; thermal trace dev {diag['thermal_trace_dev']:.1e}"
    cw = job.get("cw_path")
    if cw is not None and Path(cw).exists():
        c, q = P2.metrics(cw), P2.metrics(job["out"])
        return (f"{job['job_id']}: X = {1 - q['g_med'] / c['g_med']:+.3f}, X_coh = {1 - q['g_coh'] / c['g_coh']:+.3f}; "
                f"thermal trace dev {diag['thermal_trace_dev']:.1e}")
    return job["job_id"]


def a_h():
    return float(np.load(ROOT / "results" / "revision" / "p2" / "OP1_CW_H.npz")["amp"])


def run_va_phase(job):
    """Prefix runs of one field: CW 0-50 us, then REF modulation generated for the longest duration (220 us) and truncated."""
    job["out"].parent.mkdir(parents=True, exist_ok=True)
    t_long = np.arange(int(round(220e-6 / DT))) * DT
    phi = RP.phase({**REF_SPEC, "seed": SEEDS6[0], "t_mod": 220e-6 - T0}, t_long)
    amp, cfg = a_h(), raqr_of("C1")
    n = int(round(job["T"] / DT))
    t = t_long[:n]
    field = cfg.A_LO + amp * np.exp(1j * phi[:n]) * np.exp(2j * np.pi * 5e6 * t)
    t0 = time.perf_counter()
    y, diag = full_diag(field, cfg)
    np.savez_compressed(job["out"], diag=json.dumps(diag), T=job["T"], amp=amp, wall_s=time.perf_counter() - t0, y_tail=y[-20000:][::20])
    job["cfg_json"] = {"case": "C1 REF seed 20260701 at a_H, modulation generated for 220 us, truncated", "T_s": job["T"], "amp_Vpm": amp}
    return f"{job['job_id']}: trace dev max {diag['trace_dev_max']:.2e}, min eig {diag['min_eig_min']:.2e}, thermal {diag['thermal_trace_dev']:.1e}"


def run_va_dwell(job):
    e1 = float(np.load(KERNEL_OP1)["e1db"])            # the E1dB used by the archived dwell runs (stage-06 controls)
    job["out"].parent.mkdir(parents=True, exist_ok=True)
    amps, meta, car, labels, noise = R.realization("dwell", 0)
    env = e1 * 10 ** (3 / 20) * amps["xi_4"] * car
    env_cp = np.r_[env[-R.CP:], env]
    n = int(round(job["T"] / DT))
    cfg = raqr_of("C1")
    t = np.arange(n) * DT
    field = cfg.A_LO + env_cp[:n] * np.exp(2j * np.pi * 5e6 * t)
    t0 = time.perf_counter()
    y, diag = full_diag(field, cfg)
    np.savez_compressed(job["out"], diag=json.dumps(diag), T=job["T"], wall_s=time.perf_counter() - t0, n_total=len(env_cp), y_tail=y[-20000:][::20])
    job["cfg_json"] = {"case": "dwell xi=4, r0, +3 dB (stage06 realization), cyclic prefix, truncated", "T_s": job["T"], "n_full": len(env_cp)}
    return f"{job['job_id']}: trace dev max {diag['trace_dev_max']:.2e}, min eig {diag['min_eig_min']:.2e}, thermal {diag['thermal_trace_dev']:.1e}"


def main():
    stages = sys.argv[1:]
    jobs = []
    for st in stages:
        if st == "I":
            jobs += jobs_part1()
        elif st == "II":
            jobs += jobs_part2(P_NEW)
        elif st == "IIr":
            ref = json.loads((FIN / "02_resonance" / "refinement.json").read_text())
            for lo, periods in ref.get("periods", {}).items():
                jobs += [j for j in jobs_part2(periods) if j["lo"] == float(lo)]
        elif st == "V":
            jobs += jobs_part5()
        else:
            raise SystemExit(f"unknown stage {st}")
    todo = [j for j in jobs if not j["out"].exists()]
    wait_for_prior_runner()
    prog = Progress(total=len(jobs), done0=len(jobs) - len(todo))
    print(f"{now()} {len(jobs)} jobs, {len(todo)} to run", flush=True)
    for job in todo:
        prog.start(job)
        t0 = time.perf_counter()
        try:
            if job.get("kind") == "va_phase":
                latest = run_va_phase(job)
            elif job.get("kind") == "va_dwell":
                latest = run_va_dwell(job)
            else:
                run_phase_job(job)
                latest = latest_phase(job)
            status = "done"
        except Exception as exc:  # recorded, not hidden
            latest, status = f"{job['job_id']}: FAILED {type(exc).__name__}: {exc}", "failed"
        runtime = time.perf_counter() - t0
        if status == "done":
            manifest_row(job, job["out"], runtime, COMMAND)
        prog.finish(job, status, latest)
        print(f"{now()} {latest} ({runtime:.1f} s)", flush=True)
    print(f"{now()} stages {stages} finished", flush=True)


if __name__ == "__main__":
    main()
