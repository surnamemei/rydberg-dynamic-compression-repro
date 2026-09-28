"""Deterministic reconstruction of the stored-seed input envelopes (no atomic simulation).

Imports only waveforms06.py (numpy-only) and reproduces run06.realization() amplitude envelopes exactly
from results/stage06_dwell_physics/stage06_config.json. verify() checks bit-identity against the decimated
inputs stored with every saved trace (traces/main/*.npz, key 'a' = a[::50]).
"""
from __future__ import annotations

import json
import sys

import numpy as np

from figstyle import ROOT, RES06

sys.path.insert(0, str(ROOT / "python" / "experiments" / "stage06_dwell"))
import waveforms06 as W  # noqa: E402  (numpy only)

CFG = json.loads((RES06 / "stage06_config.json").read_text())
N = CFG["n_samples"]
TAU_SAMPLES = CFG["tau_atom_s"] * 1e9
DEC = 50


def envelopes(design: str, r: int = 0) -> dict:
    """Unit-RMS amplitude envelopes at 1 GHz for one realization (same RNG streams as run06.realization)."""
    rng = np.random.default_rng(CFG["seed_base"] + 1000 * r + (0 if design == "dwell" else 500))
    if design == "dwell":
        amps, _ = W.dwell_family(rng, N, TAU_SAMPLES, CFG["xis"], CFG["xi_ref"], CFG["edge_samples"], CFG["occupancy"], CFG["level_ratio"])
        return {f"xi_{k:g}": v for k, v in amps.items()}
    amps, _ = W.shuffle_set(rng, N, TAU_SAMPLES, **CFG["shuffle"])
    return amps


def verify(design: str, amps: dict, r: int = 0) -> dict:
    """Max |difference| between regenerated a[::50] and the stored decimated input of each saved trace."""
    out = {}
    for v, a in amps.items():
        f = RES06 / "traces" / "main" / f"{design}_r{r}_{v}_p+3_Nd4001_dt1.npz"
        if f.exists():
            out[v] = float(np.max(np.abs(np.load(f)["a"] - a[::DEC])))
    return out


def dwell_runs(a, thr):
    return W.dwell_runs(a, thr)
