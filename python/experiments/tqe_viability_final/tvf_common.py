"""Shared code of the final adversarial-validation campaign (plan: results/tqe_viability_final/00_plan.md).

Reuses the TQE-viability pass code (python/experiments/tqe_viability) unchanged; adds the campaign configurations C1-C6, the
table-derived state descriptors (E/E1dB, CW gain, local CW slope, pathology flag), progress logging and the run manifest.
"""
from __future__ import annotations

import csv
import hashlib
import json
import os
import subprocess
import sys
import threading
import time
from datetime import datetime
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "tqe_viability"))

from tv_common import (DT, KERNELS, KERNEL_OP1, OMEGA_P_WEAK, P_GRID_NS, REV2, ROOT, SEEDS3, T0, TV, FV, REV, Kernel, R, RP,  # noqa: E402,F401
                       demod, full_diag, lo_rabi_hz, n_box, phase_run, raqr, s, state_diagnostics)

FIN = ROOT / "results" / "tqe_viability_final"
SEEDS6 = tuple(int(x) for x in RP.SEEDS[:6])
REF_SPEC = {"kind": "qpsk", "Ts": 1e-6, "Tp": 300e-9}
# PID file of an earlier GPU runner to wait for (never two GPU jobs at once); optional, set via the environment.
PRIOR_PID_FILE = Path(os.environ["TVF_PRIOR_PID_FILE"]) if os.environ.get("TVF_PRIOR_PID_FILE") else None
F_RF_CARRIER_HZ = 6.946e9   # Cs 47D5/2 -> 48P3/2 transition of the model's source configuration (rotating frame; not a simulator input)

# name: (IF Hz, A_LO V/m, Omega_p or None, CW table, kernel spec, description)
CONFIGS = {
    "C1": (5e6, 0.5, None, REV2 / "tables" / "OP1_IF5.npz", ("op1",), "5 MHz, standard probe, A_LO 0.5 V/m (reference)"),
    "C2": (10e6, 0.5, None, REV2 / "tables" / "OP1_IF10_ext.npz", ("op1",), "10 MHz, standard probe, A_LO 0.5 V/m"),
    "C3": (15e6, 0.5, None, REV2 / "tables" / "OP1_IF15_ext.npz", ("op1",), "15 MHz, standard probe, A_LO 0.5 V/m"),
    "C4": (5e6, 0.5, OMEGA_P_WEAK, REV2 / "tables" / "OP1_IF5_weakprobe_ext.npz", ("weak",), "5 MHz, weaker probe (2.02 MHz), A_LO 0.5 V/m"),
    "C5": (5e6, 0.425, None, KERNELS[0.425], ("lo", 0.425), "5 MHz, standard probe, A_LO 0.425 V/m (OP3)"),
    "C6": (5e6, 0.35, None, KERNELS[0.35], ("lo", 0.35), "5 MHz, standard probe, A_LO 0.35 V/m (OP2)"),
}
# prior-pass directory names of C1-C4 (results/tqe_viability/p1/runs)
PRIOR_NAME = {"C1": "IF5", "C2": "IF10", "C3": "IF15", "C4": "IF5_weak"}


def table(cfg):
    z = np.load(CONFIGS[cfg][3])
    return {"amps": np.asarray(z["amps"], float), "gain_db": np.asarray(z["gain_db"], float), "e1db": float(z["e1db"])}


def kernel(cfg):
    spec = CONFIGS[cfg][4]
    if spec[0] == "op1":
        return Kernel.from_npz(KERNEL_OP1)
    if spec[0] == "weak":
        k = np.load(REV2 / "tables" / "kernel_weakprobe.npz")
        return Kernel(k["hR"], k["hI"], np.load(CONFIGS[cfg][3])["means"][0])
    return Kernel.from_npz(KERNELS[spec[1]])


def raqr_of(cfg):
    _, lo, op, *_ = CONFIGS[cfg]
    return raqr(lo, op)


# ------------------------------------------------------------------ table-derived state descriptors
def target_amplitude(tab, target):
    """First crossing of g_CW/g_small = target (linear interpolation in dB vs amplitude); None if no crossing on the grid."""
    a, g = tab["amps"][1:], tab["gain_db"][1:]
    tdb = 20 * np.log10(target)
    hit = np.flatnonzero(g <= tdb)
    if not len(hit) or hit[0] == 0:
        return None
    j = hit[0]
    return float(np.interp(tdb, [g[j], g[j - 1]], [a[j], a[j - 1]]))


def cw_gain_table(tab, a):
    """CW fundamental gain re the weak-tone slope at amplitude a (linear interpolation in dB)."""
    return float(10 ** (np.interp(a, tab["amps"][1:], tab["gain_db"][1:]) / 20))


def local_slope(tab, a):
    """Local CW slope d ln|output| / d ln a at amplitude a (1 = linear, 0 = saturated, < 0 = output falls with amplitude);
    central differences on the table grid, interpolated linearly in amplitude."""
    aa, g = tab["amps"][1:], tab["gain_db"][1:]
    ln_out = np.log(aa) + g / 20 * np.log(10)
    sl = np.gradient(ln_out, np.log(aa))
    return float(np.interp(a, aa, sl))


def pathology(tab, a, tol_db=0.02):
    """Plan rule I-B: a matched amplitude is pathological if the CW gain curve has a local extremum within [0.85 a, 1.15 a]
    (grid points inside the window plus the nearest grid point outside on each side; steps below tol_db count as flat)."""
    aa, g = tab["amps"][1:], tab["gain_db"][1:]
    inside = np.flatnonzero((aa >= 0.85 * a) & (aa <= 1.15 * a))
    lo_i = np.flatnonzero(aa < 0.85 * a)
    hi_i = np.flatnonzero(aa > 1.15 * a)
    idx = sorted(set(inside.tolist() + ([lo_i[-1]] if len(lo_i) else []) + ([hi_i[0]] if len(hi_i) else [])))
    d = np.diff(g[idx])
    d = d[np.abs(d) >= tol_db]
    return bool(len(d) and (d.min() < 0 < d.max()))


def describe(cfg, a):
    f_if, lo, *_ = CONFIGS[cfg]
    tab = table(cfg)
    inside = a <= tab["amps"][-1]
    return {"amp_Vpm": a, "signal_to_LO": a / lo, "E_over_E1dB": a / tab["e1db"], "signal_gt_LO": bool(a > lo),
            "table_CW_gain": cw_gain_table(tab, a) if inside else None, "local_CW_slope": local_slope(tab, a) if inside else None,
            "pathological": pathology(tab, a) if inside else None, "inside_table": bool(inside)}


# ------------------------------------------------------------------ provenance
def git_state():
    try:
        rev = subprocess.run(["git", "-C", str(ROOT), "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
        dirty = subprocess.run(["git", "-C", str(ROOT), "status", "--porcelain", "--", "python/experiments/tqe_viability",
                                "python/experiments/tqe_viability_final"], capture_output=True, text=True).stdout.strip()
        return rev, bool(dirty)
    except Exception:  # pragma: no cover
        return "unknown", True


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def now():
    return datetime.now().astimezone().isoformat(timespec="seconds")


class Progress:
    """Step 2 of the plan: after every job rewrite 00_progress.txt and append a row to 00_progress.csv. A watchdog thread logs a
    slow-job check (process and GPU state) when a job exceeds 5x the median runtime of its stage; it never kills anything."""
    COLS = ["timestamp", "stage", "job_id", "config", "matching_rule", "parameter_value", "elapsed_s", "avg_job_s", "eta_s", "status"]

    def __init__(self, total, done0=0):
        self.total, self.done = total, done0
        self.t_start = time.perf_counter()
        self.durations, self.by_stage = [], {}
        self.cur = None
        self.csv = FIN / "00_progress.csv"
        if not self.csv.exists():
            with open(self.csv, "w", newline="") as fh:
                csv.writer(fh).writerow(self.COLS)
        self._stop = threading.Event()
        threading.Thread(target=self._watch, daemon=True).start()

    def _row(self, job, elapsed, avg, eta, status):
        with open(self.csv, "a", newline="") as fh:
            csv.writer(fh).writerow([now(), job["stage"], job["job_id"], job.get("config", ""), job.get("matching_rule", ""),
                                     job.get("parameter_value", ""), f"{elapsed:.1f}", f"{avg:.1f}", f"{eta:.0f}", status])

    def start(self, job):
        self.cur = (job, time.perf_counter(), False)

    def finish(self, job, status, latest):
        elapsed = time.perf_counter() - self.cur[1]
        self.cur = None
        self.done += 1
        self.durations.append(elapsed)
        self.by_stage.setdefault(job["stage"], []).append(elapsed)
        avg = float(np.mean(self.durations[-10:]))
        eta = avg * (self.total - self.done)
        self._row(job, elapsed, avg, eta, status)
        txt = (f"[{self.done}/{self.total}] {100 * self.done / self.total:.1f}%\n"
               f"stage: {job['stage']}\ncase: {job['job_id']}\n"
               f"elapsed: {time.perf_counter() - self.t_start:.0f} s (campaign runner)\n"
               f"rolling s/job (last 10): {avg:.1f}\nETA: {eta / 60:.1f} min\nlatest result: {latest}\nupdated: {now()}\n")
        (FIN / "00_progress.txt").write_text(txt)

    def note(self, stage, job_id, status):
        self._row({"stage": stage, "job_id": job_id}, 0.0, 0.0, 0.0, status)

    def _watch(self):
        while not self._stop.is_set():
            time.sleep(30)
            cur = self.cur
            if cur is None or cur[2]:
                continue
            job, t0, _ = cur
            hist = self.by_stage.get(job["stage"], [])
            if len(hist) >= 5 and time.perf_counter() - t0 > 5 * float(np.median(hist)):
                gpu = subprocess.run(["nvidia-smi", "--query-gpu=utilization.gpu,memory.used", "--format=csv,noheader"],
                                     capture_output=True, text=True).stdout.strip()
                self.cur = (job, t0, True)
                self._row(job, time.perf_counter() - t0, float(np.median(hist)), 0.0, f"SLOW (>5x median); process alive; GPU {gpu}")
                with open(FIN / "00_progress.txt", "a") as fh:
                    fh.write(f"SLOW job {job['job_id']}: {time.perf_counter() - t0:.0f} s > 5x median {np.median(hist):.1f} s; GPU {gpu}\n")


def manifest_row(job, out, runtime, command):
    rev, dirty = git_state()
    path = FIN / "00_run_manifest.csv"
    new = not path.exists()
    with open(path, "a", newline="") as fh:
        w = csv.writer(fh)
        if new:
            w.writerow(["timestamp", "job_id", "stage", "command", "git_commit", "code_dirty", "config_json", "seed", "runtime_s",
                        "source", "output", "output_sha256"])
        w.writerow([now(), job["job_id"], job["stage"], command, rev, dirty, json.dumps(job.get("cfg_json", {}), default=float),
                    job.get("seed", ""), f"{runtime:.1f}", job.get("source", ""), str(Path(out).relative_to(ROOT)), sha256(out)])


def md_table(df, floatfmt=".3g"):
    """Markdown table without optional dependencies."""
    cols = list(df.columns)

    def fmt(v):
        if isinstance(v, (float, np.floating)):
            return "" if np.isnan(v) else format(v, floatfmt)
        return str(v)
    lines = ["| " + " | ".join(str(c) for c in cols) + " |", "|" + "---|" * len(cols)]
    for r in df.itertuples(index=False):
        lines.append("| " + " | ".join(fmt(v) for v in r) + " |")
    return "\n".join(lines)
