"""Write stage06_config.json from the Step-2 tau definition (run after tau_atom.py, before any dwell job)."""
import json
import pandas as pd
from common import OUT

d = pd.read_csv(OUT / "03_tau_atom_results.csv")
tau = float(d[(d.experiment == "IF_step_P1dB_1.00to1.05E1") & (d.Nd == 4001)].tau_1e_s.iloc[0])
cfg = {"tau_atom_s": tau,
       "tau_definition": "tau_1e of the normalized IF-envelope excess after a +5% amplitude step from E1dB (full thermal, Nd=4001, dt=1 ns); fixed before any dwell run",
       "decisive_Nd": 4001, "controls_Nd": 4001, "highest_check_Nd": 8001, "n_samples": 800_000, "cyclic_prefix_samples": 40_000,
       "analysis_skip_samples": 5_000, "symbol_samples": 1000, "phase_ramp_samples": 300, "seed_base": 20260600,
       "xis": [0.1, 0.25, 0.5, 1, 2, 4], "xi_ref": 0.05, "edge_samples": 10, "occupancy": 0.25, "level_ratio": 3.0,
       "shuffle": {"fast_xi": 0.25, "slow_xi": 2.0, "sigma_m": 0.6, "block_xi": 0.25}}
(OUT / "stage06_config.json").write_text(json.dumps(cfg, indent=2))
print(tau)
