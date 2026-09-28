"""Pilot Stage-1: thermal atomic QPSK, intensity readout, auxiliary-channel AIR.

Run from the project root. The Gaussian detectors are deliberately simple;
results are diagnostic and are not calibrated photodiode performance claims.
"""
import argparse
import json
import os
import sys
import time
from pathlib import Path

import numpy as np
from scipy.special import logsumexp

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "upstream-RydbergComms" / "python"))
os.add_dll_directory(str(ROOT / ".cuda-env" / "Library" / "bin"))
from utils.transient_quantum import SimConfig, TransientQuantumSimulator, configure_raqr


def simulate(symbols, *, ts_ns, amplitude, nd, if_cycles):
    raqr = configure_raqr("Transit")
    samples = np.repeat(np.exp(0.5j * np.pi * symbols), ts_ns)
    t = np.arange(len(samples)) * 1e-9
    field = raqr.A_LO + amplitude * samples * np.exp(2j * np.pi * if_cycles * t / (ts_ns * 1e-9))
    cfg = SimConfig(300.0, len(t), 1e-9, len(t) * 1e-9, t)
    start = time.perf_counter()
    result = TransientQuantumSimulator().run(raqr, cfg, field, init_method="SteadyState", Nd=nd, device="cuda")
    intensity = 10 ** (result.probeResponse / 10)
    return intensity.reshape(len(symbols), ts_ns), time.perf_counter() - start


def features(wave, nbin=8):
    # Eight averages preserve within-symbol transient shape in the measured intensity.
    edges = np.linspace(0, wave.shape[1], nbin + 1, dtype=int)
    return np.column_stack([wave[:, edges[i]:edges[i+1]].mean(axis=1) for i in range(nbin)])


def fit_gaussian(y, labels, nclass, variance_floor=0.05):
    global_mean = y.mean(axis=0)
    global_var = np.mean((y - global_mean) ** 2)
    mu = np.vstack([y[labels == k].mean(axis=0) if np.any(labels == k) else global_mean for k in range(nclass)])
    residual = y - mu[labels]
    variance = max(np.mean(residual ** 2), variance_floor * global_var, 1e-10)
    return mu, variance


def log_emissions(y, mu, variance):
    return -np.sum((y[:, None, :] - mu[None, :, :]) ** 2, axis=-1) / (2 * variance)


def gray_ber(pred, truth):
    code = np.array([0b00, 0b01, 0b11, 0b10], dtype=np.uint8)
    xor = np.bitwise_xor(code[pred], code[truth])
    return float(np.mean((xor & 1) + ((xor >> 1) & 1)) / 2)


def score_symbol(y, x, mu, variance):
    l = log_emissions(y, mu, variance)
    pred = np.argmax(l, axis=1)
    air = np.mean((l[np.arange(len(x)), x] - logsumexp(l, axis=1) + np.log(4)) / np.log(2))
    return float(np.mean(pred != x)), gray_ber(pred, x), float(air)


def score_sequence(y, x, mu, variance):
    # q(y_t | x_{t-1}, x_t); first symbol is excluded from both AIR and SER.
    emission = log_emissions(y[1:], mu, variance).reshape(-1, 4, 4)
    n = len(emission)
    prev, cur = x[:-1], x[1:]
    path_logq = float(np.sum(emission[np.arange(n), prev, cur]))
    log_alpha = np.full(4, -np.log(4))
    delta = np.full(4, -np.log(4))
    back = np.empty((n, 4), dtype=np.int8)
    for t in range(n):
        log_alpha = logsumexp(log_alpha[:, None] + emission[t] - np.log(4), axis=0)
        candidate = delta[:, None] + emission[t]
        back[t] = np.argmax(candidate, axis=0)
        delta = np.max(candidate, axis=0)
    log_den = float(logsumexp(log_alpha))
    decoded = np.empty(n, dtype=np.int8)
    state = int(np.argmax(delta))
    for t in range(n-1, -1, -1):
        decoded[t] = state
        state = int(back[t, state])
    air = (path_logq - log_den) / (n * np.log(2))
    return float(np.mean(decoded != cur)), gray_ber(decoded, cur), float(air)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--ts-ns", type=int, default=200)
    p.add_argument("--amplitude", type=float, default=0.2)
    p.add_argument("--nd", type=int, default=1501)
    p.add_argument("--train", type=int, default=600)
    p.add_argument("--test", type=int, default=600)
    p.add_argument("--snr-db", type=float, default=16)
    p.add_argument("--if-cycles", type=float, default=1.0)
    p.add_argument("--seed", type=int, default=20260925)
    p.add_argument("--save-waveforms", type=str, default="")
    p.add_argument("--load-waveforms", type=str, default="")
    args = p.parse_args()
    rng = np.random.default_rng(args.seed)
    x_train = rng.integers(0, 4, size=args.train)
    x_test = rng.integers(0, 4, size=args.test)
    if args.load_waveforms:
        data = np.load(args.load_waveforms)
        if (int(data["ts_ns"]) != args.ts_ns or float(data["amplitude"]) != args.amplitude
                or int(data["nd"]) != args.nd or float(data["if_cycles"]) != args.if_cycles):
            raise ValueError("Saved waveform settings differ from command-line settings")
        if not (np.array_equal(x_train, data["x_train"]) and np.array_equal(x_test, data["x_test"])):
            raise ValueError("Saved symbol sequences differ from requested seed and counts")
        train_wave, test_wave = data["train_intensity"], data["test_intensity"]
        train_seconds = test_seconds = 0.0
    else:
        train_wave, train_seconds = simulate(x_train, ts_ns=args.ts_ns, amplitude=args.amplitude, nd=args.nd, if_cycles=args.if_cycles)
        test_wave, test_seconds = simulate(x_test, ts_ns=args.ts_ns, amplitude=args.amplitude, nd=args.nd, if_cycles=args.if_cycles)
    if args.save_waveforms:
        np.savez_compressed(args.save_waveforms, x_train=x_train, x_test=x_test,
                            train_intensity=train_wave, test_intensity=test_wave,
                            ts_ns=args.ts_ns, amplitude=args.amplitude, if_cycles=args.if_cycles,
                            nd=args.nd)
    noise_sigma = np.std(train_wave) / 10 ** (args.snr_db / 20)
    y_train = features(train_wave + rng.normal(0, noise_sigma, train_wave.shape))
    y_test = features(test_wave + rng.normal(0, noise_sigma, test_wave.shape))
    # Normalization is learned only from training data.
    center = y_train.mean(axis=0)
    scale = max(float(np.std(y_train)), 1e-9)
    y_train = (y_train - center) / scale
    y_test = (y_test - center) / scale
    mu1, v1 = fit_gaussian(y_train[1:], x_train[1:], 4)
    labels4 = 4 * x_train[:-1] + x_train[1:]
    mu4, v4 = fit_gaussian(y_train[1:], labels4, 16)
    ser1, ber1, air1 = score_symbol(y_test[1:], x_test[1:], mu1, v1)
    ser4, ber4, air4 = score_sequence(y_test, x_test, mu4, v4)
    table = mu4.reshape(4, 4, -1)
    additive = table.mean(axis=1, keepdims=True) + table.mean(axis=0, keepdims=True) - table.mean(axis=(0, 1), keepdims=True)
    mu2 = additive.reshape(16, -1)
    v2 = max(float(np.mean((y_train[1:] - mu2[labels4]) ** 2)), 0.05 * float(np.var(y_train)), 1e-10)
    ser2, ber2, air2 = score_sequence(y_test, x_test, mu2, v2)
    interaction_fraction = float(np.sum((table-additive) ** 2) / max(np.sum((table-table.mean(axis=(0,1))) ** 2), 1e-10))
    out = {**vars(args), "train_seconds": train_seconds, "test_seconds": test_seconds,
           "if_hz": args.if_cycles / (args.ts_ns * 1e-9), "noise_sigma_intensity": float(noise_sigma),
           "intensity_mean": float(np.mean(train_wave)), "intensity_std": float(np.std(train_wave)),
           "M1_SER": ser1, "M1_BER": ber1, "M2_additive_SER": ser2, "M2_additive_BER": ber2,
           "M4_SER": ser4, "M4_BER": ber4, "M1_mismatched_AIR_bits_per_symbol": air1,
           "M2_additive_AIR_bits_per_symbol": air2, "transition_interaction_energy_fraction": interaction_fraction,
           "M4_mismatched_AIR_bits_per_symbol": air4, "Delta_AIR": air4-air1}
    print("STAGE1_RESULT " + json.dumps(out, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
