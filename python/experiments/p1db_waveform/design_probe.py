from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
from scipy.signal import fftconvolve, resample_poly

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
import p1db_comm_stage0 as old

FS_OUT = 1e9
N_SYM = 1020
N_SPS = 4
N_SC_SAMPLES = N_SYM * N_SPS


def occupied99(x: np.ndarray, fs: float) -> float:
    nfft = 1 << int(np.ceil(np.log2(len(x) * 8)))
    spec = np.fft.fftshift(np.fft.fft(x, nfft))
    p = np.abs(spec) ** 2
    c = np.cumsum(p) / p.sum()
    f = np.fft.fftshift(np.fft.fftfreq(nfft, 1 / fs))
    return float(f[np.searchsorted(c, .995)] - f[np.searchsorted(c, .005)])


def sc(mod: str, seed: int) -> tuple[np.ndarray, np.ndarray]:
    rng = np.random.default_rng(seed)
    if mod == "QPSK":
        alphabet = old.QPSK
    else:
        alphabet = old.QAM
    labels = rng.integers(0, len(alphabet), N_SYM)
    up = np.zeros(N_SC_SAMPLES, complex)
    up[::N_SPS] = alphabet[labels]
    tx = fftconvolve(up, old.SingleCarrier._rrcosdesign(None, beta=.2, span=10, sps=4), mode="same")
    return tx / np.sqrt(np.mean(np.abs(tx) ** 2)), labels


def ofdm(seed: int, nactive: int = 68, cp: int = 16, w: int = 16):
    rng = np.random.default_rng(seed)
    nfft = 256
    bins = np.r_[np.arange(1, nactive // 2 + 1), np.arange(nfft - nactive // 2, nfft)]
    nframes = N_SYM // nactive
    labs = rng.integers(0, 4, (nframes, nactive))
    f = np.zeros((nframes, nfft), complex)
    f[:, bins] = old.QPSK[labs]
    time = np.fft.ifft(f, axis=1) * np.sqrt(nfft)
    step = nfft + cp
    frame = step + w
    win = np.ones(frame)
    ramp = np.sin(np.pi / 2 * (np.arange(w) + .5) / w) ** 2
    win[:w] = ramp
    win[-w:] = ramp[::-1]
    x = np.zeros(nframes * step + w, complex)
    for j in range(nframes):
        blk = np.r_[time[j, -cp:], time[j], time[j, :w]]
        x[j * step:j * step + frame] += blk * win
    x = x[:nframes * step]
    fs = 20e6
    x = resample_poly(x, 50, 1)
    return x / np.sqrt(np.mean(np.abs(x) ** 2)), labs, bins, FS_OUT


def ce(seed: int, bt: float = .3):
    rng = np.random.default_rng(seed)
    bits = rng.integers(0, 2, N_SYM)
    symbols = 2 * bits - 1
    impulses = np.zeros(N_SC_SAMPLES)
    impulses[::N_SPS] = symbols
    # Gaussian frequency pulse with 4-symbol span, normalized to unit DC gain.
    t = np.arange(-8 * N_SPS, 8 * N_SPS + 1) / N_SPS
    alpha = np.sqrt(np.log(2)) / bt
    h = np.exp(-2 * (np.pi * t / alpha) ** 2)
    h /= h.sum()
    freq = fftconvolve(impulses, h, mode="same") * N_SPS
    dphi = np.pi / (2 * N_SPS) * freq
    phase = np.cumsum(dphi)
    x = np.exp(1j * phase)
    return x, bits


def main():
    seed = 20261002
    sigs = {}
    for m in ("QPSK", "16QAM"):
        x, _ = sc(m, seed)
        sigs[m] = (x, 20e6)
    x, _ = ce(seed)
    sigs["CE"] = (x, 20e6)
    for nactive, cp, w in [(68, 16, 16)]:
        x, labels, bins, fs = ofdm(seed, nactive, cp, w)
        mod = f"OFDM_N{nactive}_CP{cp}_W{w}"
        sigs[mod] = (x, fs)
    for name, (x, fs) in sigs.items():
        bw = occupied99(x, fs)
        papr = 10 * np.log10(np.max(np.abs(x) ** 2) / np.mean(np.abs(x) ** 2))
        print(name, f"Fs={fs:.3f}", f"B99={bw/1e6:.6f} MHz", f"PAPR={papr:.3f} dB", f"peak={np.max(np.abs(x)):.4f}", f"n={len(x)}")


if __name__ == "__main__":
    main()




