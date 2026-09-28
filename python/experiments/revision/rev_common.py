"""Shared paths for the pre-submission revision (preregistered in results/revision/00_revision_preregistration.json).

Archived code paths are imported unchanged; all outputs are redirected to results/revision/.
"""
from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
EXP = HERE.parent
for p in (EXP / "stage06_dwell", EXP / "stage06_dwell" / "phase_mechanism", EXP / "final_validation"):
    sys.path.insert(0, str(p))

from common import ROOT  # noqa: E402

REV = ROOT / "results" / "revision"
REV.mkdir(parents=True, exist_ok=True)
