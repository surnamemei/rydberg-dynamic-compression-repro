"""Supplementary Figure S1: additional convergence diagnostics (OP1).

(a) AIR difference vs Nd relative to Nd = 8001, prespecified +/-0.005 band (02_nd_convergence.csv);
(b) POST HOC: transient oscillation period vs Nd (manuscript/posthoc/posthoc_oscillation_vs_Nd.csv).
"""
import pandas as pd
import matplotlib.pyplot as plt

from figstyle import (fs, COL1, INK, INK2, MUTED, MODEL, FMT, ND, RES06, POSTHOC, panel, rec, save, relpath, zero_line)

d = pd.read_csv(RES06 / "02_nd_convergence.csv")
cases = [("A_lowpower_QPSK_m10", "QPSK −10 dB", FMT["QPSK"], "s", "white"), ("B_P1dB_OFDM_0", "OFDM 0 dB", FMT["OFDM"], "D", None),
         ("C_plus3_16QAM", "16-QAM +3 dB", FMT["16-QAM"], "^", None), ("D_strongest_QPSK_p6", "QPSK +6 dB", FMT["QPSK"], "s", None)]
fig = plt.figure(figsize=(COL1, 2.6))
axa, axb = fig.subplots(1, 2)

panel(axa, "a", "OP1")
axa.axhspan(-0.005, 0.005, color="#f3f2ee", lw=0, zorder=0)
axa.text(8500, -0.0085, "prespecified ±0.005", ha="right", va="top", fontsize=fs(6.0), color=INK2)
zero_line(axa)
for key, lab, col, mk, mfc in cases:
    g = d[d.case == key].sort_values("Nd")
    axa.plot(g.Nd, g.AIR_diff_vs_top, color=col, marker=mk, ms=3.2, mfc=mfc or col, mec=col if mfc else "white", mew=0.7 if mfc else 0.4, lw=0.9, label=lab)
    rec("a", f"AIR diff {lab}", g.Nd.values, g.AIR_diff_vs_top.values, source=relpath(RES06 / "02_nd_convergence.csv"))
axa.set_xscale("log")
axa.set_xticks([751, 1501, 4001, 8001])
axa.set_xticklabels(["751", "1501", "4001", "8001"])
axa.xaxis.set_minor_formatter(plt.NullFormatter())
axa.set_xlim(650, 9000)
axa.set_ylim(-0.15, 0.03)
axa.set_xlabel(ND)
axa.set_ylabel(f"AIR − AIR({ND} = 8001) (bit/symbol)")
axa.legend(loc="lower right", fontsize=fs(5.8), borderaxespad=0.2)

panel(axb, "b", "post hoc, OP1")
o = pd.read_csv(POSTHOC / "posthoc_oscillation_vs_Nd.csv")
exps = [("IF_plateau_0.5to1.5E1_on", "plateau onset", MODEL["full"], "o"), ("IF_step_P1dB_1.00to1.05E1", "+5% step", INK, "s"),
        ("IF_step_P1dB_1.05to1.00E1", "−5% step", MUTED, "s"), ("IF_plateau_0.5to1.5E1_off_recovery", "plateau recovery", INK2, "^"),
        ("LO_step_+1pct", "+1% LO step", MUTED, "D")]
for key, lab, col, mk in exps:
    g = o[o.experiment == key].sort_values("Nd")
    axb.plot(g.Nd, g.osc_period_us, color=col, marker=mk, ms=3.4, mec="white", mew=0.4, lw=1.1 if key.endswith("_on") else 0.8, label=lab)
    rec("b", f"osc period {lab}", g.Nd.values, g.osc_period_us.values, source=relpath(POSTHOC / "posthoc_oscillation_vs_Nd.csv"))
axb.set_xscale("log")
axb.set_xticks([1501, 4001, 8001])
axb.set_xticklabels(["1501", "4001", "8001"])
axb.xaxis.set_minor_formatter(plt.NullFormatter())
axb.set_xlim(1200, 10000)
axb.set_ylim(0, 4.6)
axb.set_xlabel(ND)
axb.set_ylabel("oscillation period (µs)")
axb.legend(loc="center right", fontsize=fs(5.6), borderaxespad=0.2, bbox_to_anchor=(1.0, 0.62))

save(fig, "figS1", supplement=True)
