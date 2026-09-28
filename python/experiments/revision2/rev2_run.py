"""Runs of the final receiver-physics validation pass (decision rules: results/revision2/00_decision_rules.json).

Stages (each output is skipped if it exists, so the script is resumable):
  p1  IF dependence at 5, 10, 15 MHz: CW table, CW + reference QPSK (3 sequences), 5 constant offsets, dwell xi = 0.05/4
  p3  weak probe (Omega_p = 2pi x 2.02 MHz): CW table + kernel, CW + reference QPSK (3), dwell xi = 0.05/4, steps
  p4  internal density-matrix trace for an isolated +pi/2 phase step (and a small-signal control), chained segments
Usage: python3 rev2_run.py [p1|p3|p4 ...]
"""
from __future__ import annotations

import json
import sys
import time

import numpy as np

from rev2_common import (AMPS, DT, FIVE_OFFSETS_MHZ, KERNEL_OP1, LEVEL_H, OMEGA_P_WEAK, REV2, SEEDS3, T0, Kernel, R, RP,
                         bb, cw_table, demod, extend_table, fitted_gain, full, lti_kernel, raqr, s)

T_END = T0 + 60e-6


def log(msg):
    print(msg, flush=True)
    with open(REV2 / "progress.txt", "a") as fh:
        fh.write(msg + "\n")


def table(name, cfg_raqr, f_if):
    path = REV2 / "tables" / f"{name}.npz"
    if not path.exists():
        path.parent.mkdir(parents=True, exist_ok=True)
        t0 = time.perf_counter()
        d = cw_table(cfg_raqr, f_if)
        np.savez_compressed(path, **d)
        log(f"table {name}: E1dB = {d['e1db']:.6f} V/m, settle {d['cw_settle_rel']:.1e}, {time.perf_counter() - t0:.0f} s")
    return dict(np.load(path))


def kernel(name, cfg_raqr):
    path = REV2 / "tables" / f"kernel_{name}.npz"
    if not path.exists():
        d = lti_kernel(cfg_raqr)
        np.savez_compressed(path, **d)
        log(f"kernel {name}: lin_check {d['lin_check']:.3e}")
    return dict(np.load(path))


def phase_run(sub, name, spec, cfg_raqr, f_if, e1, kern):
    out = REV2 / sub / "runs" / f"{name}.npz"
    if out.exists():
        return out
    out.parent.mkdir(parents=True, exist_ok=True)
    n = int(round(T_END / DT))
    t = np.arange(n) * DT
    phi = RP.phase({**spec, "t_mod": 60e-6}, t)
    env = spec["level"] * e1 * np.exp(1j * phi)
    t0 = time.perf_counter()
    y = full(cfg_raqr.A_LO + env * np.exp(2j * np.pi * f_if * t), cfg_raqr)
    wall = time.perf_counter() - t0
    ylin = kern.linear(env, f_if)
    z, zl = demod(y, phi, t, f_if), demod(ylin, phi, t, f_if)
    dec = 20
    np.savez_compressed(out, t=t[::dec], phi=phi[::dec], z_full=z[::dec], z_linear=zl[::dec], IF=f_if, e1db=e1, wall_s=wall,
                        spec=json.dumps(spec))
    return out


def dwell_job(sub, r, variant, cfg_raqr, f_if, tab, kern):
    out = REV2 / sub / "dwell" / f"dwell_r{r}_{variant}_p+3.npz"
    if out.exists():
        return out
    out.parent.mkdir(parents=True, exist_ok=True)
    amps, meta, car, labels, noise = R.realization("dwell", r)
    env = tab["e1db"] * 10 ** (3 / 20) * amps[variant] * car
    env_cp = np.r_[env[-R.CP:], env]
    t = np.arange(len(env_cp)) * DT
    t0 = time.perf_counter()
    y_full = full(cfg_raqr.A_LO + env_cp * np.exp(2j * np.pi * f_if * t), cfg_raqr)
    wall = time.perf_counter() - t0
    y_lin = kern.linear(env_cp, f_if)
    a = np.abs(env_cp)
    u = np.ones(len(env_cp), complex)
    np.divide(env_cp, a, out=u, where=a > 0)
    c1 = np.interp(a, tab["amps"], tab["harm"].real) + 1j * np.interp(a, tab["amps"], tab["harm"].imag)
    y_st = tab["means"][0] + np.real(c1 * u * np.exp(2j * np.pi * f_if * np.arange(len(env_cp)) / 1e9))  # models06.static
    x = R.sosfilt(R.SOS, env_cp)[::R.DEC]
    w0, w1 = (R.CP + R.SKIP) // R.DEC, len(x) - 2000 // R.DEC
    yb = {m: bb(v, tab["means"][0], f_if) for m, v in (("full", y_full), ("linear", y_lin), ("static", y_st))}
    g = {m: fitted_gain(v, x[w0:w1], w0, w1) for m, v in yb.items()}
    np.savez_compressed(out, x=x, w0=w0, w1=w1, **{f"yb_{m}": v for m, v in yb.items()}, **{f"g_{m}": v for m, v in g.items()},
                        IF=f_if, e1db=tab["e1db"], wall_s=wall)
    return out


def stage_p1():
    kern = Kernel.from_npz(KERNEL_OP1)            # the field kernel does not depend on the IF
    cfg = raqr(0.5)
    for f_mhz in (5, 10, 15):
        f_if = f_mhz * 1e6
        sub = f"p1/IF{f_mhz}"
        tab = table(f"OP1_IF{f_mhz}", cfg, f_if)
        if f_mhz != 5:                             # deviation D1: grid extended above 0.70 V/m (02_deviations.json)
            ext = REV2 / "tables" / f"OP1_IF{f_mhz}_ext.npz"
            if not ext.exists():
                d = extend_table(cfg, f_if, tab)
                np.savez_compressed(ext, **d)
                log(f"extended table OP1_IF{f_mhz}: {len(d['amps'])} amplitudes up to {d['amps'][-1]:.3f} V/m, E1dB = {d['e1db']:.6f} V/m")
            tab = dict(np.load(ext))
        e1 = float(tab["e1db"])
        if f_mhz == 5:
            arch = float(np.load(KERNEL_OP1)["e1db"])
            log(f"VALIDATION E1dB at 5 MHz: new {e1!r} vs archived {arch!r} (rel. diff {abs(e1 - arch) / arch:.1e})")
        phase_run(sub, "CW_H", {"kind": "cw", "level": LEVEL_H}, cfg, f_if, e1, kern)
        for seed in SEEDS3:
            phase_run(sub, f"REF_s{seed}_H", {"kind": "qpsk", "level": LEVEL_H, "seed": seed, "Ts": 1e-6, "Tp": 300e-9}, cfg, f_if, e1, kern)
        for mhz in FIVE_OFFSETS_MHZ:
            phase_run(sub, f"OFF{mhz:+.3f}_H", {"kind": "offset", "level": LEVEL_H, "delta_Hz": mhz * 1e6}, cfg, f_if, e1, kern)
        for r in range(4):
            for v in ("xi_0.05", "xi_4"):
                dwell_job(sub, r, v, cfg, f_if, tab, kern)
        log(f"p1 IF {f_mhz} MHz done")


def step_run(sub, name, cfg_raqr, e1, start, factor, warm=60_000, total=120_000):
    """+(factor-1) IF-envelope step at start x E1dB (fu_cw_steps / tau_atom definitions), 5 MHz."""
    import tau_atom as T
    out = REV2 / sub / "steps" / f"{name}.json"
    if out.exists():
        return json.loads(out.read_text())
    out.parent.mkdir(parents=True, exist_ok=True)
    tt = np.arange(total) * DT
    amp = np.full(total, start * e1)
    amp[warm:] *= factor
    z = T.if_envelope(full(cfg_raqr.A_LO + amp * np.exp(2j * np.pi * 5e6 * tt), cfg_raqr))
    pre = np.mean(z[warm - 5 * T.PERIOD:warm - T.PERIOD])
    fin = np.mean(z[total - 6 * T.PERIOD:total - T.PERIOD])
    u = (pre - fin) / abs(pre - fin)
    e = np.real((z[warm:total - T.PERIOD // 2] - fin) * np.conj(u)) / abs(pre - fin)
    row = {"name": name, "start_over_E1dB": start, "factor": factor, "E1dB_Vpm": e1, **T.relax_metrics(e, np.arange(len(e)) * T.DT)}
    out.write_text(json.dumps(row, default=float, indent=1))
    np.savez_compressed(out.with_suffix(".npz"), e=e[::10], t=(np.arange(len(e)) * DT)[::10])
    return row


def stage_p3():
    cfg_w = raqr(0.5, OMEGA_P_WEAK)
    tab = table("OP1_IF5_weakprobe", cfg_w, 5e6)
    if not np.isfinite(tab["e1db"]) or 3 * tab["e1db"] > 0.70:      # deviation D3 (02_deviations.json)
        ext = REV2 / "tables" / "OP1_IF5_weakprobe_ext.npz"
        if not ext.exists():
            d = extend_table(cfg_w, 5e6, tab)
            np.savez_compressed(ext, **d)
            log(f"extended table weak probe: {len(d['amps'])} amplitudes up to {d['amps'][-1]:.3f} V/m, E1dB = {d['e1db']:.6f} V/m")
        tab = dict(np.load(ext))
    kd = kernel("weakprobe", cfg_w)
    kern = Kernel(kd["hR"], kd["hI"], tab["means"][0])
    e1 = float(tab["e1db"])
    sub = "p3"
    phase_run(sub, "CW_H", {"kind": "cw", "level": LEVEL_H}, cfg_w, 5e6, e1, kern)
    for seed in SEEDS3:
        phase_run(sub, f"REF_s{seed}_H", {"kind": "qpsk", "level": LEVEL_H, "seed": seed, "Ts": 1e-6, "Tp": 300e-9}, cfg_w, 5e6, e1, kern)
    for r in range(4):
        for v in ("xi_0.05", "xi_4"):
            dwell_job(sub, r, v, cfg_w, 5e6, tab, kern)
    e1_std = float(np.load(KERNEL_OP1)["e1db"])
    for probe, cfg, e in (("weak", cfg_w, e1), ("standard", raqr(0.5), e1_std)):
        step_run(sub, f"{probe}_linear_0.05to0.10", cfg, e, 0.05, 2.0)
        row = step_run(sub, f"{probe}_plus5pct_at_1.0", cfg, e, 1.0, 1.05)
        log(f"p3 step {probe}: tau_1e {row.get('tau_1e_s')}, slow {row.get('fit2_tau_slow_s')} (a {row.get('fit2_a_slow')})")
    log("p3 done")


def stage_p4():
    """Internal density matrices around an isolated +pi/2 step, via chained cudaThermalSim segments (unmodified upstream)."""
    from utils.cu_ryd import cudaThermalSim
    from utils.transient_quantum import get_normal_quadrature, get_physical_constant
    cfg = raqr(0.5)
    e1 = float(np.load(KERNEL_OP1)["e1db"])
    hbar = 6.626e-34 / (2.0 * np.pi)
    mu_MW = 1443.45 * 1.6e-19 * 5.2918e-11
    m_cs = (132.9 / 1e3) / 6.02e23
    sig = np.sqrt(get_physical_constant("Boltzmann") * 300 / m_cs)
    x, p = get_normal_quadrature(4001)
    keep = np.abs(x) <= 3.0
    vx, p = sig * x[keep], p[keep] / p[keep].sum()
    mid = len(vx) // 2
    for name, level in (("step_aH", LEVEL_H), ("step_smallsignal", 0.05)):
        out = REV2 / "p4" / f"{name}.npz"
        if out.exists():
            continue
        out.parent.mkdir(parents=True, exist_ok=True)
        n = int(round(55e-6 / DT))
        t = np.arange(n) * DT
        Tp = 300e-9
        uu = np.clip((t - (T0 - Tp / 2)) / Tp, 0, 1)
        phi = (np.pi / 2) * (0.5 - 0.5 * np.cos(np.pi * uu))
        field = cfg.A_LO + level * e1 * np.exp(1j * phi) * np.exp(2j * np.pi * 5e6 * t)
        tc = time.perf_counter()
        y_cont = full(field, cfg)
        wall_cont = time.perf_counter() - tc
        om = mu_MW * field / hbar
        rho = s.old.sim._init_steady_state(cfg, om, vx, len(vx))
        i0, seg = int(round(49e-6 / DT)), 10
        tc = time.perf_counter()

        def call(i, j, rho_in):
            m = j - i + 1
            sc = s.old.SimConfig(300, m, DT, m * DT, np.arange(m) * DT)
            r21, rf, _ = cudaThermalSim(om[i:j + 1], rho_in, vx, p, sc, cfg)
            return r21, rf
        r21, rho = call(0, i0, rho)
        rho21 = [r21]
        idx, rbar, rmid = [i0], [rho @ p], [rho[:, mid]]
        i = i0
        while i < n - 1:
            j = min(i + seg, n - 1)
            r21, rho = call(i, j, rho)
            rho21.append(r21[1:])
            idx.append(j)
            rbar.append(rho @ p)
            rmid.append(rho[:, mid])
            i = j
        wall_chain = time.perf_counter() - tc
        rho21 = np.concatenate(rho21)
        d, lam = cfg.d, cfg.lambda_p
        mu_12 = 4.5022 * 1.6e-19 * 5.2918e-11
        probe = (20.0 / np.log(10.0)) * (-np.pi * d / lam) * np.imag(-(2.0 * cfg.N0 * mu_12 ** 2) / (8.85e-12 * hbar * cfg.Omega_p) * rho21)
        y_chain = 10 ** (probe / 10)
        dev = float(np.max(np.abs(y_chain - y_cont)) / np.mean(np.abs(y_cont - y_cont.mean())))
        log(f"p4 {name}: chained vs continuous max |dy| / mean|y - <y>| = {dev:.2e}; walls {wall_cont:.1f} s / {wall_chain:.1f} s")
        np.savez_compressed(out, t_rec=np.array(idx) * DT, rho_bar=np.array(rbar), rho_v0=np.array(rmid), y_cont=y_cont[i0:], t_y=t[i0:],
                            phi=phi[i0:], level=level, e1db=e1, chain_dev=dev, v0=vx[mid])
    log("p4 done")


def main():
    REV2.mkdir(parents=True, exist_ok=True)
    stages = sys.argv[1:] or ["p1", "p4", "p3"]
    for st in stages:
        {"p1": stage_p1, "p3": stage_p3, "p4": stage_p4}[st]()


if __name__ == "__main__":
    main()
