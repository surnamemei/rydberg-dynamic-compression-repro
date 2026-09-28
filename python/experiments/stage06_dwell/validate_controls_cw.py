"""Validate all-zone controls on off-grid CW tones against the full model (Nd=4001)."""
import time
import numpy as np
import pandas as pd
import models06 as M
from common import OUT

f1 = M.Controls.from_npz(OUT / "artifacts/controls_Nd4001.npz")
mh = M.ControlsMH(OUT / "artifacts/controls_MH_Nd4001.npz", f1)
rows = []
for a in (0.065, 0.115, 0.185, 0.28):
    env = np.full(80_000, a, complex)
    t0 = time.perf_counter()
    y = M.atomic(env, 4001)[-10_000:]
    rt = time.perf_counter() - t0
    for name, f in (("static_F1", f1.static), ("static_lti_F1", f1.static_lti), ("static_mh", mh.static_mh), ("static_lti_mh", mh.static_lti_mh)):
        z = f(env)[-10_000:]
        rows.append({"amp_Vpm": a, "model": name, "rel_rms_err_vs_full": float(np.std(y - z) / np.std(y)), "full_runtime_s": rt})
        print(rows[-1], flush=True)
pd.DataFrame(rows).to_csv(OUT / "07b_control_cw_validation.csv", index=False)
