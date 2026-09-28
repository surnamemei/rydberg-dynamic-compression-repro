"""Read-only progress monitor for the Stage-0.6 campaign (run_campaign.sh).

Never touches the running jobs: completion is detected from the per-job checkpoint JSONs
that run06.py writes atomically at the end of each job (jobs/<tag>/<job_id>.json).
The Stage-0.5 re-check (first campaign step) has no per-job files; its 64 full-atomic
runs are counted as completed from 09_stage05_recheck.csv.

ETA weights: Nd=8001 jobs ~3.5x and dt=0.5 ns jobs ~2x a main Nd=4001 job (Step-1 timings).
"""
from __future__ import annotations

import csv
import json
import re
import statistics
import sys
import time
from datetime import datetime
from pathlib import Path

import pandas as pd

OUT = Path(__file__).resolve().parents[3] / "results" / "stage06_dwell_physics"
CFG = json.loads((OUT / "stage06_config.json").read_text())
START = datetime.fromisoformat(sys.argv[1]).timestamp() if len(sys.argv) > 1 else None
# Optional plan file: list of {tag, job_id, weight, family, realization, xi, power_db}.
PLAN_FILE = OUT / sys.argv[2] if len(sys.argv) > 2 else None
TXT, CSV = OUT / "progress_stage06.txt", OUT / "progress_stage06.csv"
DWELL = ["xi_0.05"] + [f"xi_{x:g}" for x in CFG["xis"]]
SHUF = ["ORIGINAL", "BLOCK_SHUFFLED_CYCLES", "BLOCK_SHUFFLED_FIXED", "SHUFFLED"]
R8 = range(8)


def plan():
    """(tag, job_id, weight) for every run06 job in run_campaign.sh, in execution order."""
    jobs = []

    def add(design, reals, powers, variants, nd=4001, dt=1, tag="main", w=1.0):
        for r in reals:
            for p in powers:
                for v in variants:
                    jobs.append((tag, f"{design}_r{r}_{v}_p{p:+d}_Nd{nd}_dt{dt:g}", w))
    add("dwell", R8, (3, 0, 6), DWELL)
    add("shuffle", R8, (3, 0, 6), SHUF)
    add("dwell", R8, (-6,), DWELL)
    add("shuffle", R8, (-6,), SHUF)
    num_d, num_s = ["xi_0.05", "xi_0.25", "xi_1", "xi_4"], ["ORIGINAL", "BLOCK_SHUFFLED_CYCLES"]
    add("dwell", range(2), (3,), num_d, nd=8001, tag="num_Nd8001", w=3.5)
    add("dwell", range(2), (3,), num_d, dt=0.5, tag="num_dt0p5", w=2.0)
    add("shuffle", range(2), (3,), num_s, nd=8001, tag="num_Nd8001", w=3.5)
    add("shuffle", range(2), (3,), num_s, dt=0.5, tag="num_dt0p5", w=2.0)
    return jobs


def parse(jid):
    m = re.match(r"(dwell|shuffle)_r(\d+)_(.+)_p([+-]\d+)_Nd(\d+)_dt([\d.]+)$", jid)
    design, r, v, p, nd, dt = m.groups()
    xi = v.split("_", 1)[1] if v.startswith("xi_") else v
    return design, int(r), xi, int(p), int(nd), dt


def fmt(sec):
    sec = int(max(sec, 0))
    return f"{sec // 3600:d}h{sec % 3600 // 60:02d}m{sec % 60:02d}s"


def main():
    if PLAN_FILE:
        spec = json.loads(PLAN_FILE.read_text())
        jobs = [(j["tag"], j["job_id"], j["weight"]) for j in spec]
        info = {(j["tag"], j["job_id"]): j for j in spec}
        recheck_n, prev_t = 0, START
    else:
        jobs, info = plan(), {}
        recheck_n = len(pd.read_csv(OUT / "09_stage05_recheck.csv")) // 5 if (OUT / "09_stage05_recheck.csv").exists() else 0
        prev_t = (OUT / "09_stage05_recheck.csv").stat().st_mtime
    total = recheck_n + len(jobs)
    seen, durs = set(), []  # durs: (elapsed_s, weight)
    if not CSV.exists():
        with CSV.open("w", newline="") as f:
            csv.writer(f).writerow(["timestamp", "job_id", "family", "realization", "xi", "power_db", "elapsed_s", "avg_job_s", "eta_s"])
    while True:
        done_now = []
        for tag, jid, w in jobs:
            p = OUT / "jobs" / tag / f"{jid}.json"
            if (tag, jid) not in seen and p.exists():
                done_now.append((p.stat().st_mtime, tag, jid, w))
        for mt, tag, jid, w in sorted(done_now):
            seen.add((tag, jid))
            el = mt - prev_t
            prev_t = mt
            durs.append((el, w))
            recent = durs[-10:]
            avg_norm = sum(e / ww for e, ww in recent) / len(recent)
            remaining_w = sum(ww for t2, j2, ww in jobs if (t2, j2) not in seen)
            eta = remaining_w * avg_norm
            done = recheck_n + len(seen)
            if (tag, jid) in info:
                j = info[(tag, jid)]
                fam, r, xi, pw = j["family"], j["realization"], j["xi"], j["power_db"]
            else:
                design, r, xi, pw, nd, dt = parse(jid)
                fam = f"{design}/{tag}" + (f" Nd{nd}" if nd != 4001 else "") + (f" dt{dt}ns" if dt != "1" else "")
            since = (mt - START) if START else 0
            line = (f"[{done}/{total}] {100 * done / total:5.1f}% | {fam} | r{r} | xi={xi} | {pw:+d} dB | "
                    f"elapsed {fmt(since)} | rolling avg {sum(e for e, _ in recent) / len(recent):.1f} s/job | ETA {fmt(eta)}")
            print(line, flush=True)
            with CSV.open("a", newline="") as f:
                csv.writer(f).writerow([datetime.fromtimestamp(mt).isoformat(timespec="seconds"), jid, fam, r, xi, pw,
                                        round(el, 1), round(sum(e for e, _ in recent) / len(recent), 1), round(eta)])
            all_el = [e for e, _ in durs if e > 0]
            TXT.write_text(line + "\n"
                           f"completed {done} / total {total}" + (f" (incl. {recheck_n} Stage-0.5 re-check runs)" if recheck_n else "") + f", remaining {total - done}\n"
                           f"mean {statistics.mean(all_el):.1f} s/job, median {statistics.median(all_el):.1f} s/job over {len(all_el)} run06 jobs\n"
                           f"updated {datetime.now().isoformat(timespec='seconds')}\n")
        if len(seen) == len(jobs):
            break
        time.sleep(10)


if __name__ == "__main__":
    main()
