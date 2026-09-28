"""Seeded fair-waveform generator: a verbatim numpy/scipy copy of stage05.waveform() (python/experiments/p1db_waveform/stage05.py).

stage05.py imports p1db_comm_stage0, which constructs a TransientQuantumSimulator at import time; the generator is copied
here so that no simulator code is loaded. verify() compares the 99% bandwidth and PAPR of every regenerated waveform with
the stored per-seed statistics (results/p1db_waveform_stage05/03_waveform_statistics_extended.csv).
"""
from __future__ import annotations

import sys

import numpy as np
import pandas as pd
from scipy.signal import fftconvolve, resample_poly

from figstyle import ROOT, ST05

sys.path.insert(0, str(ROOT / "python"))
from utils.sim_SingleCarrier import SingleCarrier  # noqa: E402  (numpy/scipy only)

FS = 1e9
UPSAMPLE = 50
N_SYMBOLS = 1020
SPS = 4
N_SC = N_SYMBOLS * SPS
NFFT = 256
N_ACTIVE = 68
CP = 16
WOLA = 16
N_FRAMES = N_SYMBOLS // N_ACTIVE
BLOCK = NFFT + CP
FRAME = BLOCK + WOLA
OFDM_BINS = np.r_[np.arange(1, N_ACTIVE // 2 + 1), np.arange(NFFT - N_ACTIVE // 2, NFFT)]
RRC = SingleCarrier._rrcosdesign(None, beta=.2, span=10, sps=SPS)
ALPHA_CE = .36
QPSK = np.exp(1j * (np.pi / 4 + np.pi / 2 * np.arange(4)))
_grid = np.array([-3, -1, 1, 3])
QAM = (_grid[:, None] + 1j * _grid[None, :]).ravel() / np.sqrt(10)
FMT_INDEX = {"CE": 2, "QPSK": 3, "16QAM": 4, "OFDM": 5}


def seeds(mod):
    return [20262000 + 100 * r + FMT_INDEX[mod] for r in range(8)]


def _norm(x):
    return x / np.sqrt(np.mean(np.abs(x) ** 2))


def waveform(mod: str, seed: int):
    rng = np.random.default_rng(seed)
    if mod in ("QPSK", "16QAM"):
        alphabet = QPSK if mod == "QPSK" else QAM
        labels = rng.integers(0, len(alphabet), N_SYMBOLS)
        up = np.zeros(N_SC, complex)
        up[::SPS] = alphabet[labels]
        tx20 = fftconvolve(up, RRC, mode="same")
        return _norm(resample_poly(_norm(tx20), UPSAMPLE, 1))
    if mod == "CE":
        labels = rng.integers(0, 2, N_SYMBOLS)
        symbols = 2 * labels - 1
        impulses = np.zeros(N_SC)
        impulses[::SPS] = symbols
        t = np.arange(-8 * SPS, 8 * SPS + 1) / SPS
        alpha = np.sqrt(np.log(2)) / ALPHA_CE
        h = np.exp(-2 * (np.pi * t / alpha) ** 2)
        h /= h.sum()
        freq = fftconvolve(impulses, h, mode="same") * SPS
        freq_ns = np.repeat(freq, UPSAMPLE)
        phase = np.cumsum(freq_ns * (np.pi * 5e6 / (2 * FS)))
        return np.exp(1j * phase)
    if mod == "OFDM":
        labels = rng.integers(0, len(QPSK), (N_FRAMES, N_ACTIVE))
        for k in range(N_ACTIVE):
            labels[:4, k] = rng.permutation(4)
        freq = np.zeros((N_FRAMES, NFFT), complex)
        freq[:, OFDM_BINS] = QPSK[labels]
        time = np.fft.ifft(freq, axis=1) * np.sqrt(NFFT)
        window = np.ones(FRAME)
        ramp = np.sin(np.pi / 2 * (np.arange(WOLA) + .5) / WOLA) ** 2
        window[:WOLA] = ramp
        window[-WOLA:] = ramp[::-1]
        tx20 = np.zeros(N_FRAMES * BLOCK + WOLA, complex)
        for j in range(N_FRAMES):
            block = np.r_[time[j, -CP:], time[j], time[j, :WOLA]]
            tx20[j * BLOCK:j * BLOCK + FRAME] += block * window
        tx20 = tx20[:N_FRAMES * BLOCK]
        return _norm(resample_poly(_norm(tx20), UPSAMPLE, 1))
    raise ValueError(mod)


def occupied99(x):
    nfft = 1 << int(np.ceil(np.log2(len(x) * 8)))
    spec = np.fft.fftshift(np.fft.fft(x, nfft))
    p = np.abs(spec) ** 2
    c = np.cumsum(p) / np.sum(p)
    f = np.fft.fftshift(np.fft.fftfreq(nfft, 1 / FS))
    return float(f[np.searchsorted(c, .995)] - f[np.searchsorted(c, .005)])


def stored_stats():
    W = pd.read_csv(ST05 / "03_waveform_statistics_extended.csv")
    W = W[W.config_version.isin(["fair_v2", "fair_v4_randomized_balanced_ofdm_training"]) & ~((W.modulation == "OFDM") & (W.config_version == "fair_v2"))]
    return W.drop_duplicates(["modulation", "seed"]).set_index(["modulation", "seed"])
