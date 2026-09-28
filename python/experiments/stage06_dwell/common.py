"""Shared Stage-0.6 setup: reuse the archived Stage-0.5 receiver chain unchanged.

Stage-0.5 modules call the Windows-only ``os.add_dll_directory`` at import time.
On Linux the CUDA runtime is found through the rpath baked into
``python/utils/bin/cu_ryd.so``, so a no-op shim keeps the archived code untouched.
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

if not hasattr(os, "add_dll_directory"):
    os.add_dll_directory = lambda path: None

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "results" / "stage06_dwell_physics"
STAGE05 = ROOT / "results" / "p1db_waveform_stage05"
sys.path.insert(0, str(ROOT / "python" / "experiments" / "p1db_waveform"))
sys.path.insert(0, str(ROOT / "python"))

import stage05 as s  # noqa: E402

old = s.old
E1 = s.E1
FS = s.FS
