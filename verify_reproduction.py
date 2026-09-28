"""One-command reproduction check (no simulation, no GPU).

Run after `python3 fetch_upstream.py`. It
  1. checks every shipped file against SHA256SUMS;
  2. snapshots the shipped outputs;
  3. rebuilds all manuscript and supplementary figures (with the figure audit), the numbers ledger, the revision
     tables P2-P6, the P1 tables that derive from the shipped per-job metrics, and the final-check (revision2)
     tables from the shipped run files;
  4. compares every regenerated file with the shipped one (byte-identical; otherwise PDFs must be pixel-identical
     at 150 dpi with identical text, and CSV/JSON tables and NPZ arrays numerically identical to 1e-9 relative with
     identical text values);
  5. restores the shipped files, so the checkout stays identical to the release.
Exit code 0 means everything reproduced.
"""
from __future__ import annotations

import hashlib
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
PY = sys.executable


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def run(cmd, cwd):
    r = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True)
    return r.returncode, r.stdout + r.stderr


def same_pdf(a, b, tmp):
    if shutil.which("pdftoppm") is None:
        return False, "bytes differ (pdftoppm unavailable for a pixel comparison)"
    from PIL import Image, ImageChops
    ia, ib = tmp / "a", tmp / "b"
    subprocess.run(["pdftoppm", "-r", "150", "-png", "-singlefile", str(a), str(ia)], check=True)
    subprocess.run(["pdftoppm", "-r", "150", "-png", "-singlefile", str(b), str(ib)], check=True)
    pa, pb = Image.open(f"{ia}.png").convert("RGB"), Image.open(f"{ib}.png").convert("RGB")
    ta = subprocess.run(["pdftotext", str(a), "-"], capture_output=True, text=True).stdout
    tb = subprocess.run(["pdftotext", str(b), "-"], capture_output=True, text=True).stdout
    ok = pa.size == pb.size and ImageChops.difference(pa, pb).getbbox() is None and ta == tb
    return ok, "pixel- and text-identical (bytes differ)" if ok else "DIFFERENT"


def same_csv(a, b, rtol=1e-9):
    """Tables re-derived from shipped CSV inputs can differ from the shipped tables (computed in memory) in the last
    floating-point digits. Accept identical shape, columns and non-numeric cells with numbers equal to rtol."""
    import numpy as np
    import pandas as pd
    x, y = pd.read_csv(a), pd.read_csv(b)
    if list(x.columns) != list(y.columns) or x.shape != y.shape:
        return False, "DIFFERENT (shape or columns)"
    num = x.select_dtypes("number").columns
    other = [c for c in x.columns if c not in num]
    if not (x[other].fillna("<NA>").astype(str).values == y[other].fillna("<NA>").astype(str).values).all():
        return False, "DIFFERENT (non-numeric cells)"
    xv, yv = x[num].to_numpy(float), y[num].to_numpy(float)
    ok = np.allclose(xv, yv, rtol=rtol, atol=0, equal_nan=True)
    rel = np.nanmax(np.abs(xv - yv) / np.maximum(np.abs(yv), 1e-300)) if xv.size else 0.0
    return ok, f"numerically identical (max rel. diff {rel:.1e}; floating-point round trip)" if ok else f"DIFFERENT (max rel. diff {rel:.1e})"


def same_json(a, b, rtol=1e-9):
    import json
    import math

    volatile = {"generated", "code_version", "code_dirty", "runtime_s"}   # generation metadata of the final-campaign verdicts

    def eq(x, y):
        if isinstance(x, dict):
            return isinstance(y, dict) and x.keys() == y.keys() and all(eq(x[k], y[k]) for k in x if k not in volatile)
        if isinstance(x, list):
            return isinstance(y, list) and len(x) == len(y) and all(eq(u, v) for u, v in zip(x, y))
        if isinstance(x, bool) or isinstance(y, bool):
            return x == y
        if isinstance(x, (int, float)) and isinstance(y, (int, float)):
            if math.isnan(x) or math.isnan(y):
                return math.isnan(x) and math.isnan(y)
            return abs(x - y) <= rtol * max(abs(x), abs(y), 1e-300)
        return x == y
    ok = eq(json.loads(a.read_text()), json.loads(b.read_text()))
    return ok, "numerically identical to 1e-9 (floating-point round trip)" if ok else "DIFFERENT"


def same_npz(a, b, rtol=1e-9):
    import numpy as np
    x, y = np.load(a), np.load(b)
    if sorted(x.files) != sorted(y.files):
        return False, "DIFFERENT (array names)"
    for k in x.files:
        u, v = x[k], y[k]
        if u.shape != v.shape or (u.dtype.kind in "fc" and not np.allclose(u, v, rtol=rtol, atol=0, equal_nan=True)) or (u.dtype.kind not in "fc" and not np.array_equal(u, v)):
            return False, f"DIFFERENT ({k})"
    return True, "arrays identical to 1e-9 relative (bytes differ)"


def main():
    if not (HERE / "python" / "utils" / "transient_quantum.py").exists():
        sys.exit("Upstream RydbergComms not found: run `python3 fetch_upstream.py` first.")
    bad = []
    for line in (HERE / "SHA256SUMS").read_text().splitlines():
        h, rel = line.split("  ", 1)
        if sha(HERE / rel) != h:
            bad.append(rel)
    print(f"[1] SHA256SUMS: {'all files match' if not bad else f'{len(bad)} mismatches: {bad[:5]}'}")
    outputs = (sorted((HERE / "manuscript" / "figures" / "final").glob("fig*.pdf"))
               + sorted((HERE / "manuscript" / "figures" / "supplement").glob("figS*.pdf"))
               + sorted((HERE / "manuscript" / "figures" / "audit").glob("*_plotted_values.csv"))
               + [HERE / "manuscript" / "figures" / "audit" / "figure_audit_results.json", HERE / "manuscript" / "numbers_ledger.csv"]
               + sorted((HERE / "results" / "revision").glob("[1-5][0-9]_*.csv"))
               + [HERE / "results" / "revision" / "50_p2_p6_summary.json"]
               + [HERE / "results" / "revision" / "p1" / f for f in ("02_p1_contrasts.csv", "04_p1_decisive_table.csv",
                                                                      "06_p1_posthoc_gain_residual_decomposition.csv",
                                                                      "07_p1_posthoc_gain_vs_dwell.csv")]
               + sorted(p for p in (HERE / "results" / "revision2").glob("[1-4][0-9]_*") if p.suffix in (".csv", ".json"))
               + sorted((HERE / "results" / "revision2" / "p4").glob("*_traces.npz"))
               + [p for d in ("01_matching", "02_resonance", "03_memory_model", "04_noise", "05_numerical")
                  for p in sorted((HERE / "results" / "tqe_viability_final" / d).glob("*")) if p.suffix in (".csv", ".json")])
    restore_only = [p for d in ("01_matching", "02_resonance", "03_memory_model", "04_noise", "05_numerical")
                    for p in sorted((HERE / "results" / "tqe_viability_final" / d).glob("*")) if p.suffix in (".md", ".png", ".npz")]
    # The identified-memory-model fit (tvf_p3_gmp.py) also needs the 128 revision P1 dwell/shuffle traces (about 320 MB), which are
    # not shipped (regenerated bit-identically by python/experiments/revision/rev_p1_rerun.py on a GPU). Without them the step is skipped.
    have_traces = (HERE / "results" / "revision" / "p1" / "traces" / "p1").is_dir()
    if not have_traces:
        gmp = [p for p in outputs if p.parent.name == "03_memory_model"]
        outputs = [p for p in outputs if p not in gmp]
        restore_only += gmp
    before = {p for p in HERE.rglob("*") if p.is_file()}
    tmp = Path(tempfile.mkdtemp(prefix="verify_"))
    snap = {p: tmp / "snap" / p.relative_to(HERE) for p in outputs + restore_only}
    for p, q in snap.items():
        q.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(p, q)
    steps = {}
    rc, out = run([PY, "make_all.py"], HERE / "manuscript" / "figures" / "src")
    steps["figures + figure audit (make_all.py)"] = rc == 0 and "ALL PASS" in out
    rc, out = run([PY, "build_numbers_ledger.py"], HERE / "manuscript")
    steps["numbers ledger (build_numbers_ledger.py)"] = rc == 0
    rc, out = run([PY, "rev_p2_p6_analyze.py"], HERE / "python" / "experiments" / "revision")
    steps["revision tables P2-P6 (rev_p2_p6_analyze.py)"] = rc == 0
    code = ("import pandas as pd, rev_p1_analyze as A\n"
            "df = pd.read_csv(A.P1 / '01_p1_metric_rows.csv')\n"
            "con = A.contrasts(df); con.to_csv(A.P1 / '02_p1_contrasts.csv', index=False)\n"
            "tab, v = A.classify(con); tab.to_csv(A.P1 / '04_p1_decisive_table.csv', index=False)\n"
            "A.posthoc_decomposition(df).to_csv(A.P1 / '06_p1_posthoc_gain_residual_decomposition.csv', index=False)\n"
            "A.posthoc_gain_vs_dwell().to_csv(A.P1 / '07_p1_posthoc_gain_vs_dwell.csv', index=False)\n"
            "print('P1 verdict', v)\n")
    rc, out = run([PY, "-c", code], HERE / "python" / "experiments" / "revision")
    steps["revision tables P1 from shipped per-job metrics"] = rc == 0
    for script in ("rev2_p2.py", "rev2_p1_analyze.py", "rev2_p3_analyze.py", "rev2_p4_analyze.py"):
        rc, out = run([PY, script], HERE / "python" / "experiments" / "revision2")
        steps[f"final-check tables ({script})"] = rc == 0
    for script in ("tvf_p1_analyze.py", "tvf_p2_analyze.py", "tvf_p3_gmp.py", "tvf_p4_analyze.py", "tvf_p5_analyze.py"):
        if script == "tvf_p3_gmp.py" and not have_traces:
            steps[f"final-campaign tables ({script}; needs the revision P1 traces)"] = None
            continue
        rc, out = run([PY, script], HERE / "python" / "experiments" / "tqe_viability_final")
        steps[f"final-campaign tables ({script})"] = rc == 0
    for k, v in steps.items():
        print(f"[3] {k}: {'SKIPPED' if v is None else ('OK' if v else 'FAILED')}")
    diffs, notes = [], []
    for p, q in snap.items():
        if p in restore_only or sha(p) == sha(q):
            continue
        if p.suffix in (".pdf", ".csv", ".json", ".npz"):
            # the identified-memory-model fit solves ill-conditioned normal equations; multithreaded BLAS reduction order changes
            # its outputs at about 1e-8 relative, so its tables are compared at 1e-6 (all reported values have at most 4 digits)
            rt = 1e-6 if p.parent.name == "03_memory_model" else 1e-9
            cmp = {".pdf": lambda: same_pdf(p, q, tmp), ".csv": lambda: same_csv(p, q, rt), ".json": lambda: same_json(p, q, rt), ".npz": lambda: same_npz(p, q, rt)}
            ok, msg = cmp[p.suffix]()
            notes.append(f"{p.relative_to(HERE)}: {msg}")
            if ok:
                continue
        diffs.append(str(p.relative_to(HERE)))
    n_cmp = len(snap) - len(restore_only)
    print(f"[4] compared {n_cmp} regenerated files: {n_cmp - len(diffs)} reproduced, {len(diffs)} different ({len(restore_only)} summaries/plots restored without comparison)")
    for n in notes:
        print("    note:", n)
    for d in diffs:
        print("    DIFFERENT:", d)
    for p, q in snap.items():
        shutil.copy2(q, p)
    shutil.rmtree(tmp)
    extra = [p for p in HERE.rglob("*") if p.is_file() and p not in before]
    for p in extra:
        p.unlink()
    for d in sorted((d for d in HERE.rglob("__pycache__") if d.is_dir()), reverse=True):
        shutil.rmtree(d, ignore_errors=True)
    print(f"[5] shipped files restored; {len(extra)} by-products removed (SVG/PNG copies, caches)")
    ok = not bad and all(v is not False for v in steps.values()) and not diffs
    print("REPRODUCTION:", "PASS" if ok else "FAIL")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
