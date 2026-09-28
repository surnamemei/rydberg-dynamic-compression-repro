"""Stage-0.6 controlled amplitude envelopes and the constant-envelope information carrier.

All amplitude envelopes are dimensionless with unit RMS; callers scale them to a field.

Dwell family (exact amplitude multiset)
---------------------------------------
Two levels a_L < a_H joined by raised-cosine edges of Te samples.  Every member of a
realization contains the same number n_A of rising (R) and falling (F) edges, the same
number of a_H plateau samples and the same number of a_L samples, so the sorted sample
vectors are bitwise identical.  Only the arrangement differs:
  * reference A: n_A plateaus of d0 samples each (shortest dwell);
  * member xi:   plateau samples regrouped into runs of ~xi*tau, the surplus R/F pairs
                 placed back-to-back as brief spikes (above threshold for Te samples).
Because mean dwell = high time / number of edges, an exact multiset with a fixed edge
shape cannot change the per-event mean dwell; the physically relevant statistic is the
time-weighted dwell  T_dw = sum(d_i^2) / sum(d_i)  (dwell of the event containing a
randomly chosen above-threshold sample), which is reported alongside plateau length.

Shuffle test (exact amplitude multiset)
---------------------------------------
ORIGINAL: clustered continuous envelope |z_fast| * exp(sigma_m * s_slow).
BLOCK-SHUFFLED: permutation of whole excursion cycles cut at upward median crossings
                (continuous at the joins; excursion shapes and dwell kept; clustering removed).
BLOCK-FIXED:    permutation of fixed-length blocks (cuts through excursions; adds jumps).
SHUFFLED:       random sample permutation (white envelope; strictest, least physical).
"""
from __future__ import annotations

import numpy as np

FS = 1e9


def rc_edge(te: int) -> np.ndarray:
    i = np.arange(te)
    return 0.5 - 0.5 * np.cos(np.pi * (i + 0.5) / te)  # symmetric: r[i] + r[te-1-i] = 1


def _split(total: int, n: int, rng, shape: float = 4.0) -> np.ndarray:
    """Split an integer total into n non-negative integer parts with Gamma(shape) weights."""
    if n <= 0:
        return np.zeros(0, int)
    w = rng.gamma(shape, 1.0, n)
    parts = np.floor(total * w / w.sum()).astype(int)
    rem = total - parts.sum()
    parts[rng.choice(n, rem, replace=False)] += 1
    return parts


def dwell_family(rng, n: int, tau_samples: float, xis, xi_ref: float, te: int, p: float = 0.25, ratio: float = 3.0):
    """Return {xi: unit-RMS amplitude}, metadata.  xi_ref sets the reference plateau d0."""
    d0 = int(round(xi_ref * tau_samples))
    high = int(round(p * n))
    n_a = high // (d0 + te)
    p_tot = high - n_a * te  # plateau samples; above-threshold = p_tot + n_a*te = high exactly
    low_tot = n - p_tot - 2 * n_a * te
    edge = rc_edge(te)
    a_l, a_h = 1.0, ratio
    rise = a_l + (a_h - a_l) * edge
    fall = rise[::-1]
    out, meta = {}, {}
    for xi in [xi_ref, *[x for x in xis if x != xi_ref]]:
        if xi == xi_ref:
            plateaus = np.full(n_a, d0)
            plateaus[: p_tot - n_a * d0] += 1
        else:
            n_b = max(1, int(round(p_tot / (xi * tau_samples))))
            plateaus = np.full(n_b, p_tot // n_b)
            plateaus[: p_tot - plateaus.sum()] += 1
            plateaus = np.r_[plateaus, np.zeros(n_a - n_b, int)]
        rng.shuffle(plateaus)
        gaps = _split(low_tot, n_a + 1, rng)
        pieces = [np.full(gaps[0], a_l)]
        for k, pl in enumerate(plateaus):
            pieces += [rise, np.full(pl, a_h), fall, np.full(gaps[k + 1], a_l)]
        a = np.concatenate(pieces)
        assert len(a) == n
        out[xi] = a
        d_above = plateaus + te
        meta[xi] = {"n_events": int(n_a), "n_plateaus": int(np.count_nonzero(plateaus)), "n_spikes": int(np.count_nonzero(plateaus == 0)),
                    "plateau_len_samples": float(np.mean(plateaus[plateaus > 0])),
                    "tw_dwell_samples": float(np.sum(d_above.astype(float) ** 2) / np.sum(d_above)),
                    "spike_share_of_high_time": float(np.count_nonzero(plateaus == 0) * te / high)}
    scale = np.sqrt(np.mean(out[xi_ref] ** 2))
    ref_sorted = np.sort(out[xi_ref])
    for xi in out:
        assert np.array_equal(np.sort(out[xi]), ref_sorted), "multiset not preserved"
        out[xi] = out[xi] / scale
    return out, meta


def _gauss_process(rng, n: int, corr_samples: float, complex_: bool = False) -> np.ndarray:
    """Stationary Gaussian process with Gaussian autocorrelation exp(-(lag/corr)^2), unit variance."""
    f = np.fft.rfftfreq(n) if not complex_ else np.fft.fftfreq(n)
    psd = np.exp(-(np.pi * f * corr_samples) ** 2)
    if complex_:
        w = rng.normal(size=n) + 1j * rng.normal(size=n)
        x = np.fft.ifft(np.fft.fft(w) * np.sqrt(psd))
        return x / np.sqrt(np.mean(np.abs(x) ** 2))
    w = rng.normal(size=n)
    x = np.fft.irfft(np.fft.rfft(w) * np.sqrt(psd), n)
    return x / np.std(x)


def shuffle_set(rng, n: int, tau_samples: float, fast_xi: float = 0.25, slow_xi: float = 2.0, sigma_m: float = 0.6, block_xi: float = 0.25):
    fast = np.abs(_gauss_process(rng, n, fast_xi * tau_samples, complex_=True))
    slow = _gauss_process(rng, n, slow_xi * tau_samples)
    orig = fast * np.exp(sigma_m * slow)
    orig /= np.sqrt(np.mean(orig ** 2))
    med = np.median(orig)
    up = np.flatnonzero((orig[:-1] < med) & (orig[1:] >= med)) + 1
    cycles = np.split(orig, up)
    head, body = cycles[0], cycles[1:]
    order = rng.permutation(len(body))
    block = np.concatenate([head, *[body[i] for i in order]])
    b = int(round(block_xi * tau_samples))
    nb = n // b
    blocks = orig[: nb * b].reshape(nb, b)[rng.permutation(nb)].ravel()
    blockfix = np.r_[blocks, orig[nb * b:]]
    shuffled = orig[rng.permutation(n)]
    out = {"ORIGINAL": orig, "BLOCK_SHUFFLED_CYCLES": block, "BLOCK_SHUFFLED_FIXED": blockfix, "SHUFFLED": shuffled}
    ref = np.sort(orig)
    for k, v in out.items():
        assert len(v) == n and np.array_equal(np.sort(v), ref), k
    return out, {"n_cycles": len(body), "block_len_samples": b}


def phase_carrier(rng, n: int, ts: int, tp: int):
    """Constant-envelope QPSK: phase steps smoothed by a raised cosine of tp samples at each boundary."""
    nsym = n // ts
    labels = rng.integers(0, 4, nsym)
    phases = np.pi / 4 + np.pi / 2 * labels
    phi = np.repeat(phases, ts).astype(float)
    ramp = rc_edge(tp)
    for k in range(1, nsym):
        d = np.angle(np.exp(1j * (phases[k] - phases[k - 1])))
        if abs(abs(d) - np.pi) < 1e-9:
            d = np.pi
        i0 = k * ts - tp // 2
        phi[i0:i0 + tp] = phases[k - 1] + d * ramp
    phi = np.r_[phi, np.full(n - len(phi), phi[-1])]
    return np.exp(1j * phi), labels


def dwell_runs(a: np.ndarray, thr: float) -> np.ndarray:
    mask = a > thr
    edge = np.diff(np.r_[False, mask, False].astype(np.int8))
    return np.flatnonzero(edge == -1) - np.flatnonzero(edge == 1)


def env_corr_1e(a: np.ndarray) -> float:
    ac = a - a.mean()
    nfft = 1 << int(np.ceil(np.log2(2 * len(ac) - 1)))
    sp = np.fft.rfft(ac, nfft)
    c = np.fft.irfft(sp * np.conj(sp), nfft)[: len(ac)]
    c /= c[0]
    cut = np.flatnonzero(c <= np.exp(-1))
    return float(cut[0]) if len(cut) else float(len(a))


def occupied99(x: np.ndarray) -> float:
    nfft = 1 << int(np.ceil(np.log2(len(x) * 2)))
    p = np.abs(np.fft.fftshift(np.fft.fft(x, nfft))) ** 2
    c = np.cumsum(p) / p.sum()
    f = np.fft.fftshift(np.fft.fftfreq(nfft, 1 / FS))
    return float(f[np.searchsorted(c, .995)] - f[np.searchsorted(c, .005)])


def ks_distance(a: np.ndarray, b: np.ndarray) -> float:
    grid = np.sort(np.r_[a, b])
    ca = np.searchsorted(np.sort(a), grid, side="right") / len(a)
    cb = np.searchsorted(np.sort(b), grid, side="right") / len(b)
    return float(np.max(np.abs(ca - cb)))
