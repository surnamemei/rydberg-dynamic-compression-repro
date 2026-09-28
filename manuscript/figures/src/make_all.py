"""Regenerate all final figures from stored results, then run the audit.

Usage (from manuscript/figures/src):  python3 make_all.py
No simulation is run; every script reads stored CSV/JSON/NPZ results (see each script's docstring).
"""
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
SCRIPTS = [f"fig{i}.py" for i in range(1, 11)] + [f"figS{i}.py" for i in range(1, 11)]

if __name__ == "__main__":
    for s in SCRIPTS:
        print(f"== {s}", flush=True)
        subprocess.run([sys.executable, str(HERE / s)], cwd=HERE, check=True)
    sys.exit(subprocess.run([sys.executable, str(HERE / "audit_figures.py")], cwd=HERE).returncode)
