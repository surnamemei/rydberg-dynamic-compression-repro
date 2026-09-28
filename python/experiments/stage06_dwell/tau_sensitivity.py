"""tau_atom sensitivity: the Step-2 'P1dB' steps used the Stage-0.5 E1dB (0.0877 V/m).
Repeat the +5% IF-envelope step at the converged E1dB (Nd=4001) and at 0.8x / 1.2x of it."""
import numpy as np
import pandas as pd
import tau_atom as T
from common import OUT, s

E1c = float(np.load(OUT / "artifacts/controls_Nd4001.npz")["e1db"])
rows = []
for f in (0.8, 1.0, 1.2):
    ea = f * E1c / s.E1
    warm, total = 60_000, 120_000
    tt = np.arange(total) * T.DT
    amp = np.full(total, ea * s.E1)
    amp[warm:] = 1.05 * ea * s.E1
    y = T.run_field(s.old.raqr.A_LO + amp * np.exp(2j * np.pi * T.IF * tt), 4001)
    z = T.if_envelope(y)
    pre = np.mean(z[warm - 5 * T.PERIOD:warm - T.PERIOD])
    fin = np.mean(z[total - 6 * T.PERIOD:total - T.PERIOD])
    u = (pre - fin) / abs(pre - fin)
    e = np.real((z[warm:total - T.PERIOD // 2] - fin) * np.conj(u)) / abs(pre - fin)
    rows.append({"experiment": f"IF_step_{f:.1f}xE1conv_+5pct", "E_start_Vpm": ea * s.E1, "E1dB_converged_Vpm": E1c, "Nd": 4001,
                 **T.relax_metrics(e, np.arange(len(e)) * T.DT)})
    print(rows[-1], flush=True)
pd.DataFrame(rows).to_csv(OUT / "03b_tau_atom_sensitivity.csv", index=False)
