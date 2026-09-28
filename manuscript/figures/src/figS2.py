"""Supplementary Figure S2: spectra of the four fair-waveform formats.

Input waveforms regenerated from the stored seeds with fmtwave.py (verbatim copy of stage05.waveform; no atomic simulation);
every regenerated waveform reproduces its stored 99% bandwidth and PAPR (03_waveform_statistics_extended.csv).
PSD: Welch average (Hann, 16384-sample segments at 1 GHz) over the 8 realizations of each format, normalized to the peak.
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import welch

from figstyle import (fs, COL1, INK2, FMT, FMT_LABEL, ST05, panel, rec, save, relpath)
import fmtwave as F

st = F.stored_stats()
fig = plt.figure(figsize=(COL1, 2.4))
ax = fig.subplots()
ax.set_title("fair-waveform formats (input envelopes)", fontsize=fs(7.5), pad=3)
ax.axvspan(-3, 3, color="#f3f2ee", lw=0, zorder=0)
ax.text(0, -136, "receiver band ±3 MHz", ha="center", va="bottom", fontsize=fs(6.0), color=INK2)
for mod, key in (("CE", "CE"), ("QPSK", "QPSK"), ("16QAM", "16-QAM"), ("OFDM", "OFDM")):
    ps, b99 = [], []
    for s in F.seeds(mod):
        x = F.waveform(mod, s)
        b = F.occupied99(x)
        assert abs(b - st.loc[(mod, s)].B_occ_99_Hz) < 1e-3, (mod, s)
        b99.append(b)
        f, p = welch(x, fs=F.FS, window="hann", nperseg=16384, return_onesided=False, detrend=False)
        ps.append(p)
    f = np.fft.fftshift(f) / 1e6
    p = np.fft.fftshift(np.mean(ps, axis=0))
    pdb = 10 * np.log10(p / p.max())
    sel = np.abs(f) <= 12
    ax.plot(f[sel], pdb[sel], color=FMT[key], lw=0.9, label=f"{FMT_LABEL[key]}")
    rec("", f"B99 mean {key} (MHz)", "B99", float(np.mean(b99)) / 1e6, source=relpath(ST05 / "03_waveform_statistics_extended.csv"),
        note="regenerated per-seed B99 equal to stored values")
ax.set_xlim(-12, 12)
ax.set_ylim(-140, 5)
ax.set_yticks([0, -20, -40, -60, -80, -100, -120, -140])
ax.set_xlabel("frequency offset from the IF (MHz)")
ax.set_ylabel("normalized PSD (dB)")
ax.legend(loc="upper right", fontsize=fs(6.0), borderaxespad=0.2, handlelength=1.4)
ax.text(-11.5, 0, "99% bandwidths\n5.36–5.37 MHz", ha="left", va="top", fontsize=fs(6.0))

save(fig, "figS2", supplement=True)
