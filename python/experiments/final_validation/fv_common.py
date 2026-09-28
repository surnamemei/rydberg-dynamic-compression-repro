"""Shared setup for the final-validation runs (reuses the Stage-0.6 code paths unchanged)."""
from __future__ import annotations

import dataclasses
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
S06 = HERE.parent / "stage06_dwell"
sys.path.insert(0, str(S06))
sys.path.insert(0, str(S06 / "phase_mechanism"))

import numpy as np  # noqa: E402

from common import OUT as OUT06, ROOT, s  # noqa: E402

FV = ROOT / "results" / "final_validation"
DT = 1e-9
IF = 5e6


def raqr_with_lo(lo):
    return dataclasses.replace(s.old.raqr, A_LO=lo)


def run_field(field, nd, lo, dt=DT):
    """Full thermal transient model; returns probe intensity, probe dB and final density matrices."""
    n = len(field)
    t = np.arange(n) * dt
    res = s.old.sim.run(raqr_with_lo(lo), s.old.SimConfig(300, n, dt, n * dt, t), field, Nd=nd, device="cuda")
    return 10 ** (res.probeResponse / 10), res.probeResponse, res.rhovec_final
