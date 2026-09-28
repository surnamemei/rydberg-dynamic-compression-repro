"""POST-HOC (no simulation): spectral structure of the constant-amplitude test signals over the
modulated window, from the known phase sequences. Tests whether signal power near the IF
(where the nonlinear response is sharply selective) orders the extra compression X."""
from pathlib import Path
import sys
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "python" / "experiments" / "stage06_dwell" / "phase_mechanism"))
sys.path.insert(0, str(ROOT / "python" / "experiments" / "stage06_dwell"))
import pm_run as P  # noqa: E402

res = pd.read_csv(ROOT / "results/stage06_dwell_physics/phase_mechanism/03_phase_mechanism_results.csv")
res = res[(res.Nd == 4001) & (res.dt_ns == 1)].set_index("case")
t = np.arange(int(P.T0 * 1e9), int(P.T_END * 1e9)) * 1e-9
rows = []
for c, case in P.CASES.items():
    s = np.exp(1j * P.phase(case, t))
    n = 1 << 22
    S = np.abs(np.fft.fft(s, n)) ** 2
    f = np.fft.fftfreq(n, 1e-9)
    S /= S.sum()
    row = {"case": c, "X": float(res.loc[c, "X"])}
    for bw in (25e3, 125e3, 250e3, 500e3, 1e6):
        row[f"power_within_{bw / 1e3:g}kHz"] = float(S[np.abs(f) <= bw].sum())
    row["spectral_centroid_kHz"] = float(np.sum(S * f) / 1e3)
    row["max_bin_fraction"] = float(S.max())  # large for line spectra
    rows.append(row)
d = pd.DataFrame(rows)
d.to_csv(Path(__file__).resolve().parent / "posthoc_spectral_structure.csv", index=False)
pd.set_option("display.width", 200)
print(d.round(4).to_string())
ph = d[~d.case.str.startswith("OFF")]
for col in [c for c in d.columns if c.startswith("power_within")]:
    r = np.corrcoef(ph[col], ph.X)[0, 1]
    print(f"phase cases (incl. CW): corr(X, {col}) = {r:+.2f}")
