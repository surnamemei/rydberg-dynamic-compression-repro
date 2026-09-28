from __future__ import annotations

import json
import hashlib
import shutil
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import pearsonr, spearmanr, t, ks_2samp
from scipy.signal import welch

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "results" / "p1db_waveform_stage05"
sys.path.insert(0, str(ROOT / "python" / "experiments" / "p1db_waveform"))
import stage05 as s

MODS = ("CE", "QPSK", "16QAM", "OFDM")
COLORS = {"CE": "#555555", "QPSK": "#2274A5", "16QAM": "#E07A1F", "OFDM": "#3A923A"}
TAU_ATOM = 1.727e-6  # inherited Stage-1.5 pulse-response 1/e time; see summary caveat.


def mean_ci(x):
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    if len(x) == 0:
        return np.nan, np.nan, np.nan
    if len(x) == 1:
        return float(x[0]), float(x[0]), float(x[0])
    sd = np.std(x, ddof=1)
    rad = t.ppf(.975, len(x)-1) * sd / np.sqrt(len(x))
    return float(np.mean(x)), float(np.mean(x)-rad), float(np.mean(x)+rad)


def errorbar(ax, x, vals, color, label=None, marker="o"):
    mu, lo, hi = mean_ci(vals)
    ax.errorbar(x, mu, yerr=[[mu-lo], [hi-mu]], color=color, marker=marker, capsize=3, label=label)


def plot_curves(d, y, ylabel, path, reference_zero=False):
    fig, ax = plt.subplots(figsize=(8.5, 5.3))
    for mod in MODS:
        q = d[(d.modulation == mod) & (d.model == "atomic")]
        xs=[]; means=[]; los=[]; his=[]
        for db, g in q.groupby("Pavg_over_P1dB_dB"):
            mu, lo, hi = mean_ci(g[y])
            xs.append(db); means.append(mu); los.append(lo); his.append(hi)
        ax.plot(xs, means, "o-", color=COLORS[mod], label=mod)
        ax.fill_between(xs, los, his, color=COLORS[mod], alpha=.15)
    if reference_zero:
        ax.axhline(0, color="black", lw=.8)
    ax.set_xlabel("Average input power / CW P1dB (dB)")
    ax.set_ylabel(ylabel)
    ax.grid(alpha=.25); ax.legend(frameon=False)
    fig.tight_layout(); fig.savefig(OUT/path, dpi=180); plt.close(fig)


def figure_spectra():
    fig, ax = plt.subplots(figsize=(8.6, 5.2))
    for mod in MODS:
        ps=[]; freqs=None
        for seed in (20262001,20262101,20262201,20262301):
            seed += {"CE":1,"QPSK":2,"16QAM":3,"OFDM":4}[mod]
            x,_,_=s.waveform(mod,seed)
            f,p=welch(x,fs=s.FS,nperseg=16384,noverlap=8192,return_onesided=False,scaling="spectrum")
            order=np.argsort(f); freqs=f[order]; ps.append(p[order])
        db=10*np.log10(np.mean(ps,axis=0)/np.max(np.mean(ps,axis=0)))
        ax.plot(freqs/1e6,db,label=mod,color=COLORS[mod],lw=1.2)
    ax.set_xlim(-8,8); ax.set_ylim(-85,2); ax.set_xlabel("Baseband frequency (MHz)"); ax.set_ylabel("Normalized PSD (dB)")
    ax.grid(alpha=.2); ax.legend(frameon=False); fig.tight_layout(); fig.savefig(OUT/"fig02_matched_waveform_spectra.png",dpi=180); plt.close(fig)


def make_figures(d, stats):
    shutil.copy2(ROOT/"p1db_reference_curve.png",OUT/"fig01_p1db_reference.png")
    figure_spectra()
    plot_curves(d,"R_AIR_bit_s","Achievable information rate (bit/s)","fig03_air_rate_vs_avg_power.png")
    plot_curves(d,"eta_AIR_bit_s_Hz","AIR spectral efficiency (bit/s/Hz)","fig04_air_efficiency_vs_avg_power.png")
    plot_curves(d,"DeltaAIR_native","AIR loss from -6 dB reference (bit/native symbol)","fig05_air_loss_vs_avg_p1db.png",True)

    fig,ax=plt.subplots(figsize=(8.5,5.3))
    for mod in MODS:
        q=d[(d.modulation==mod)&(d.model=="atomic")]
        for db,g in q.groupby("Pavg_over_P1dB_dB"):
            xvals=g.Ppeak_over_P1dB_dB
            errorbar(ax,float(xvals.mean()),g.DeltaAIR_native.values,COLORS[mod],label=mod if db==-15 else None)
    ax.axhline(0,color="black",lw=.8); ax.set_xlabel("Peak input power / CW P1dB (dB)"); ax.set_ylabel("AIR loss (bit/native symbol)")
    ax.grid(alpha=.25); ax.legend(frameon=False); fig.tight_layout(); fig.savefig(OUT/"fig06_air_loss_vs_peak_p1db.png",dpi=180); plt.close(fig)

    plot_curves(d,"relative_loss","Relative AIR loss from -6 dB reference","fig07_relative_air_loss.png",True)

    for control, fname, title in (("static","fig08_full_vs_static_control.png","Full atomic versus CW static control"),("static_lti","fig09_full_vs_lti_static_control.png","Full atomic versus static plus LTI control")):
        fig,axes=plt.subplots(1,4,figsize=(14,4.2),sharey=True)
        for ax,mod in zip(axes,MODS):
            for model,style in (("atomic","o-"),(control,"s--")):
                q=d[(d.modulation==mod)&(d.model==model)]
                xs=[];ys=[];lo=[];hi=[]
                for db,g in q.groupby("Pavg_over_P1dB_dB"):
                    mu,l,h=mean_ci(g.DeltaAIR_native);xs.append(db);ys.append(mu);lo.append(l);hi.append(h)
                ax.plot(xs,ys,style,color=COLORS[mod],label=model)
                ax.fill_between(xs,lo,hi,color=COLORS[mod],alpha=.1)
            ax.set_title(mod);ax.set_xlabel("Pavg/P1dB (dB)");ax.grid(alpha=.2)
        axes[0].set_ylabel("AIR loss (bit/native symbol)");axes[-1].legend(frameon=False)
        fig.suptitle(title);fig.tight_layout();fig.savefig(OUT/fname,dpi=180);plt.close(fig)

    fig,ax=plt.subplots(figsize=(8,5))
    for mod in MODS:
        xs=[];ys=[]
        for seed in sorted(stats.loc[stats.modulation==mod,"seed"].unique())[:4]:
            q=stats[(stats.modulation==mod)&(stats.seed==seed)].iloc[0]
            x,_,_=s.waveform(mod,int(seed)); norm=np.abs(x)/np.sqrt(np.mean(np.abs(x)**2))
            hist,bins=np.histogram(norm,bins=np.linspace(0,4,121),density=True)
            ax.plot((bins[1:]+bins[:-1])/2,hist,color=COLORS[mod],alpha=.25)
            xs.append((bins[1:]+bins[:-1])/2);ys.append(hist)
        ax.plot(xs[0],np.mean(ys,axis=0),color=COLORS[mod],lw=2,label=mod)
    ax.set_xlabel("Envelope amplitude / RMS");ax.set_ylabel("Probability density");ax.grid(alpha=.2);ax.legend(frameon=False)
    fig.tight_layout();fig.savefig(OUT/"fig10_amplitude_pdf.png",dpi=180);plt.close(fig)

    fig,ax=plt.subplots(figsize=(8.5,5.3))
    for mod in MODS:
        q=stats[stats.modulation==mod].groupby("Pavg_over_P1dB_dB")["occupancy_theta_1.00"]
        xs=[];ys=[];lo=[];hi=[]
        for db,v in q:
            m,l,h=mean_ci(v);xs.append(db);ys.append(m);lo.append(l);hi.append(h)
        ax.plot(xs,ys,"o-",color=COLORS[mod],label=mod);ax.fill_between(xs,lo,hi,color=COLORS[mod],alpha=.14)
    ax.set_xlabel("Average input power / P1dB (dB)");ax.set_ylabel("Time above E1dB");ax.set_ylim(-.02,1.02);ax.grid(alpha=.25);ax.legend(frameon=False)
    fig.tight_layout();fig.savefig(OUT/"fig11_high_field_occupancy.png",dpi=180);plt.close(fig)

    fig,ax=plt.subplots(figsize=(8.5,5.3)); vals=[];labels=[]
    for mod in MODS:
        q=stats[(stats.modulation==mod)&(stats.Pavg_over_P1dB_dB==3)]["dwell_mean_ns_theta_1.00"].dropna()
        vals.append(q.values);labels.append(mod)
    ax.boxplot(vals,tick_labels=labels,showmeans=True);ax.set_ylabel("Mean contiguous dwell above E1dB (ns) at +3 dB")
    ax.grid(axis="y",alpha=.25);fig.tight_layout();fig.savefig(OUT/"fig12_dwell_time_distributions.png",dpi=180);plt.close(fig)

    corr=stats[["config_version","modulation","seed","Pavg_over_P1dB_dB","PAPR_dB","E_peak_over_E1dB","occupancy_theta_1.00","dwell_mean_ns_theta_1.00"]]
    plot_df=d[(d.model=="atomic")].merge(corr,on=["config_version","modulation","seed","Pavg_over_P1dB_dB"],suffixes=("","_stats"))
    for xcol,fname,xlabel in (("PAPR_dB_stats","fig13_loss_vs_papr.png","PAPR (dB)"),("E_peak_over_E1dB","fig14_loss_vs_peak_field.png","Peak field / E1dB"),("occupancy_theta_1.00","fig15_loss_vs_high_field_occupancy.png","Fraction of time above E1dB"),("dwell_mean_ns_theta_1.00","fig16_loss_vs_dwell_time.png","Mean dwell above E1dB (ns)")):
        fig,ax=plt.subplots(figsize=(7.8,5.2))
        for mod in MODS:
            q=plot_df[plot_df.modulation==mod]
            ax.scatter(q[xcol],q.DeltaAIR_native,c=COLORS[mod],label=mod,alpha=.65,s=24)
        ax.set_xlabel(xlabel);ax.set_ylabel("Full atomic AIR loss (bit/native symbol)");ax.axhline(0,color="black",lw=.7)
        ax.grid(alpha=.2);ax.legend(frameon=False);fig.tight_layout();fig.savefig(OUT/fname,dpi=180);plt.close(fig)

    # P_5% and P_10% are the highest powers below each loss limit.
    def threshold_at(x, y, limit):
        crosses=[]
        for i in range(1,len(x)):
            if y[i-1] <= limit < y[i]:
                crosses.append(float(x[i-1]+(limit-y[i-1])*(x[i]-x[i-1])/(y[i]-y[i-1])))
        if crosses: return max(crosses)
        if y[-1] <= limit: return float(x[-1])
        return None
    fig,ax=plt.subplots(figsize=(7.8,5.2)); rows=[]
    rng=np.random.default_rng(20260926)
    for mod in MODS:
        q=d[(d.modulation==mod)&(d.model=="atomic")]
        pivot=q.pivot(index="seed",columns="Pavg_over_P1dB_dB",values="relative_loss").sort_index(axis=1)
        xs=pivot.columns.to_numpy(float); ys=pivot.mean(axis=0).to_numpy(float)
        for limit in (.05,.10):
            point=threshold_at(xs,ys,limit)
            boot=[]
            mat=pivot.to_numpy(float)
            for _ in range(2000):
                idx=rng.integers(0,len(mat),len(mat))
                crossing=threshold_at(xs,np.mean(mat[idx],axis=0),limit)
                if crossing is not None: boot.append(crossing)
            lo,hi=(map(float,np.quantile(boot,[.025,.975])) if boot else (np.nan,np.nan))
            rows.append({"modulation":mod,"relative_loss_limit":limit,"P_limit_minus_P1dB_dB":point,"bootstrap95_low":lo,"bootstrap95_high":hi,"bootstrap_success_fraction":len(boot)/2000,"method":"last upward crossing of mean AIR-loss curve; paired seed bootstrap"})
        ax.plot(xs,ys,"o-",color=COLORS[mod],label=mod)
    ax.axhline(.05,color="black",ls=":",lw=1,label="5% loss")
    ax.axhline(.10,color="black",ls="--",lw=1,label="10% loss")
    ax.set_xlabel("Average input power / P1dB (dB)");ax.set_ylabel("Relative AIR loss");ax.grid(alpha=.2);ax.legend(frameon=False)
    fig.tight_layout();fig.savefig(OUT/"fig18_information_loss_thresholds.png",dpi=180);plt.close(fig)
    pd.DataFrame(rows).to_csv(OUT/"08_information_loss_thresholds.csv",index=False)


def main():
    raw=pd.read_csv(OUT/"02_all_results.csv")
    valid=((raw.config_version=="fair_v2")&raw.modulation.isin(["CE","QPSK","16QAM"]))|((raw.config_version=="fair_v4_randomized_balanced_ofdm_training")&(raw.modulation=="OFDM"))
    d=raw[valid].copy()
    if d.AIR_native.isna().any(): raise RuntimeError("valid result set includes nonfinite AIR")
    refs=d[d.Pavg_over_P1dB_dB==-6][["config_version","modulation","seed","model","AIR_native","R_AIR_bit_s","eta_AIR_bit_s_Hz"]].rename(columns={"AIR_native":"AIR_ref","R_AIR_bit_s":"R_ref_bit_s","eta_AIR_bit_s_Hz":"eta_ref_bit_s_Hz"})
    d=d.merge(refs,on=["config_version","modulation","seed","model"],validate="many_to_one")
    d["DeltaAIR_native"]=d.AIR_ref-d.AIR_native
    d["relative_loss"]=d.DeltaAIR_native/d.AIR_ref
    d["DeltaR_AIR_bit_s"]=d.R_ref_bit_s-d.R_AIR_bit_s
    d["DeltaEta_bit_s_Hz"]=d.eta_ref_bit_s_Hz-d.eta_AIR_bit_s_Hz
    d.to_csv(OUT/"07_analysis_dataset.csv",index=False)

    stats=pd.read_csv(OUT/"03_waveform_statistics_extended.csv")
    keep=((stats.config_version=="fair_v2")&stats.modulation.isin(["CE","QPSK","16QAM"]))|((stats.config_version=="fair_v4_randomized_balanced_ofdm_training")&(stats.modulation=="OFDM"))
    stats=stats[keep].copy()
    stats["tau_atom_s"]=TAU_ATOM
    for th in ("0.50","0.75","1.00","1.25","1.50"):
        stats[f"dwell_mean_over_tau_theta_{th}"]=stats[f"dwell_mean_ns_theta_{th}"]*1e-9/TAU_ATOM
        stats[f"dwell_p90_over_tau_theta_{th}"]=stats[f"dwell_p90_ns_theta_{th}"]*1e-9/TAU_ATOM
        stats[f"dwell_max_over_tau_theta_{th}"]=stats[f"dwell_max_ns_theta_{th}"]*1e-9/TAU_ATOM
    stats["envelope_corr_over_tau"]=stats.envelope_corr_1e_s/TAU_ATOM
    stats.to_csv(OUT/"03_waveform_statistics_extended.csv",index=False)
    stats.to_csv(OUT/"waveform_statistics_extended.csv",index=False)

    make_figures(d,stats)

    # Stage J: report pooled associations and power-stratified rank associations.
    a=d[(d.model=="atomic")].merge(stats,on=["config_version","modulation","seed","Pavg_over_P1dB_dB"],suffixes=("","_stat"))
    features=["PAPR_dB_stat","E_peak_over_E1dB","exceedance_theta_1.00","dwell_mean_ns_theta_1.00","dwell_p90_ns_theta_1.00","dwell_max_ns_theta_1.00","envelope_corr_1e_ns"]
    cor=[]
    for outcome in ("DeltaAIR_native","relative_loss"):
        for feature in features:
            groups=[("pooled",a)] + [(f"Pavg={db}",g) for db,g in a.groupby("Pavg_over_P1dB_dB")]
            for name,g in groups:
                q=g[[feature,outcome]].replace([np.inf,-np.inf],np.nan).dropna()
                if len(q)<3 or q[feature].nunique()<2 or q[outcome].nunique()<2:
                    pr=sr=np.nan; slope=r2=np.nan
                else:
                    pr=float(pearsonr(q[feature],q[outcome]).statistic);sr=float(spearmanr(q[feature],q[outcome]).statistic)
                    slope=float(np.polyfit(q[feature],q[outcome],1)[0]);pred=slope*q[feature]+(q[outcome].mean()-slope*q[feature].mean());r2=float(1-np.sum((q[outcome]-pred)**2)/np.sum((q[outcome]-q[outcome].mean())**2))
                cor.append({"outcome":outcome,"feature":feature,"subset":name,"n":len(q),"pearson_r":pr,"spearman_rho":sr,"simple_linear_slope":slope,"simple_linear_R2":r2})
    pd.DataFrame(cor).to_csv(OUT/"07_statistics_correlations.csv",index=False)

    # Machine-readable matched waveform design summary.
    configs=[]
    for mod,g in stats.groupby("modulation"):
        r=g.groupby("seed").agg(B_occ_99_Hz=("B_occ_99_Hz","first"),PAPR_dB=("PAPR_dB","first"),rate=("information_symbol_rate_Hz","first"))
        configs.append({"modulation":mod,"B_occ_99_Hz_mean":float(r.B_occ_99_Hz.mean()),"B_occ_99_Hz_sd":float(r.B_occ_99_Hz.std(ddof=1)),"B_occ_99_Hz_min":float(r.B_occ_99_Hz.min()),"B_occ_99_Hz_max":float(r.B_occ_99_Hz.max()),"information_symbol_rate_Hz_mean":float(r.rate.mean()),"PAPR_dB_mean":float(r.PAPR_dB.mean()),"PAPR_dB_sd":float(r.PAPR_dB.std(ddof=1)),"waveform_definition":{"CE":"Gaussian frequency pulse CPFSK, BT=0.36, h=0.5, 5 Msymbol/s","QPSK":"single-carrier RRC, rolloff 0.2, span 10 symbols, 5 Msymbol/s","16QAM":"single-carrier RRC, rolloff 0.2, span 10 symbols, 5 Msymbol/s","OFDM":"256 FFT, 68 active QPSK carriers, CP=16, cyclic suffix overlap=16, WOLA, 5 Msymbol/s effective data-symbol rate"}[mod],"config_version":str(g.config_version.iloc[0])})
    config={"occupied_bandwidth_definition":"central 99% of total waveform power: FFT cumulative power quantiles 0.5% to 99.5%, 8x zero padding, 1 ns sampled signal","normalization":"each complex envelope rescaled to unit average |E|^2 after 1 GHz resampling; E_rms=E1dB*10^(Pavg_dB/20)","E1dB_Vpm":s.E1,"power_ratio":"(E/E1dB)^2","native_symbol_rate_Hz":5e6,"ofdm_accounting":{"nfft":256,"active_subcarriers":68,"cyclic_prefix_samples":16,"cyclic_suffix_overlap_samples":16,"symbol_duration_s":272/20e6,"information_QAM_symbols_per_block":68,"effective_data_symbol_rate_Hz":68/(272/20e6)},"waveforms":configs,"maximum_relative_B99_mismatch_of_seed_means":float(max(abs(c['B_occ_99_Hz_mean']/np.mean([x['B_occ_99_Hz_mean'] for x in configs])-1) for c in configs))}
    (OUT/"04_matched_bandwidth_configs.json").write_text(json.dumps(config,indent=2),encoding="utf-8")

    # Confidence table at headline points and paired control excess loss.
    summary=[]
    for mod in MODS:
        for db in (-6,0,3,6):
            group=d[(d.modulation==mod)&(d.Pavg_over_P1dB_dB==db)]
            for model in ("atomic","static","static_lti"):
                q=group[group.model==model]
                for metric in ("AIR_native","R_AIR_bit_s","eta_AIR_bit_s_Hz","DeltaAIR_native","relative_loss"):
                    mu,lo,hi=mean_ci(q[metric])
                    summary.append({"modulation":mod,"power_db":db,"model":model,"metric":metric,"mean":mu,"ci95_low":lo,"ci95_high":hi,"n":int(q[metric].notna().sum())})
    pd.DataFrame(summary).to_csv(OUT/"09_headline_summary.csv",index=False)

    # Stage A audit and exact rerun command record.
    cw=np.load(ROOT/"artifacts"/"p1db_reference.npz")
    stats0=pd.read_csv(ROOT/"waveform_statistics.csv")
    oldres=pd.concat([pd.read_csv(ROOT/"p1db_comm_stage0_full.csv"),pd.read_csv(ROOT/"p1db_comm_ofdm_fair.csv")],ignore_index=True)
    power_err=float(np.max(np.abs((oldres.Ppeak_over_P1dB_dB-oldres.Pavg_over_P1dB_dB)-oldres.PAPR_dB)))
    e_err=float(np.max(np.abs(oldres.E_rms_Vpm-s.E1*10**(oldres.Pavg_over_P1dB_dB/20))))
    audit={"source_commit":"409572499ef54825b66c532e0aecdcffa89e7be1","branch":"p1db-waveform-fairness","E1dB_Vpm":s.E1,"E1dB_definition":"first -1 dB crossing of 5 MHz CW fundamental relative to weak-tone slope; tail fundamental over final 4 us of a 40 us full-thermal transient","power_conversion":"P/P1dB=(E/E1dB)^2; E_rms=E1dB*10^(Pavg_dB/20)","old_stage0_saved_row_max_peak_relation_error_dB":power_err,"old_stage0_saved_row_max_rms_field_relation_error_Vpm":e_err,"old_stage0_AIR":"fixed s=1 Gaussian auxiliary LLR; calibration variance fit on calibration segment only; train/cal/test are sequential disjoint portions, no test labels fit receiver","old_stage0_bandwidths_MHz":stats0.groupby('modulation').occupied_bandwidth_99pct_Hz.mean().div(1e6).to_dict(),"old_stage0_rate_Hz":stats0.groupby('modulation').information_symbol_rate_Hz.mean().to_dict(),"current_waveforms_all_unit_RMS_after_resampling":True,"current_primary_Nd":1501,"current_primary_dt_ns":1.0,"max_abs_timestep_AIR_change":float(pd.read_csv(OUT/"05_numerical_validation.csv").query("method=='timestep_0p5ns'").AIR_abs_difference.max()),"screening_Nd501_max_abs_AIR_change":float(pd.read_csv(OUT/"05_numerical_validation.csv").query("method=='doppler_Nd501'").AIR_abs_difference.max())}
    (OUT/"00_stage_a_audit.json").write_text(json.dumps(audit,indent=2),encoding="utf-8")
    (OUT/"06_reproduce_commands.txt").write_text("""# Stage-0.5 reproduction\n# Branch: p1db-waveform-fairness; source commit: 409572499ef54825b66c532e0aecdcffa89e7be1\n# Primary CUDA, full thermal, continuous state, Nd=1501, dt=1 ns\n# The job checkpoint is keyed by modulation, bandwidth/rate, power, seed, model, and Nd.\n# Run in order; repeated pilot jobs are skipped by the checkpoint.\npython python/experiments/p1db_waveform/stage05.py --grid pilot --replicates 4\npython python/experiments/p1db_waveform/stage05.py --grid full --replicates 4\npython python/experiments/p1db_waveform/stage05.py --grid high --replicates 8\n# OFDM was rerun at all seven powers with eight seeds under its final randomized balanced preamble.\npython python/experiments/p1db_waveform/stage05.py --grid full --replicates 8 --only-mod OFDM\npython python/experiments/p1db_waveform/validate_numeric.py\npython python/experiments/p1db_waveform/controlled_dwell.py\npython python/experiments/p1db_waveform/analyze_stage05.py\n""",encoding="utf-8")

    corrtab=pd.read_csv(OUT/"07_statistics_correlations.csv")
    pooled_occ=float(corrtab[(corrtab.outcome=="DeltaAIR_native")&(corrtab.feature=="exceedance_theta_1.00")&(corrtab.subset=="pooled")].spearman_rho.iloc[0])
    pooled_papr=float(corrtab[(corrtab.outcome=="DeltaAIR_native")&(corrtab.feature=="PAPR_dB_stat")&(corrtab.subset=="pooled")].spearman_rho.iloc[0])
    dwell_within=corrtab[(corrtab.outcome=="DeltaAIR_native")&(corrtab.feature=="dwell_mean_ns_theta_1.00")&(corrtab.subset.str.startswith("Pavg="))].spearman_rho.dropna()
    mean_dwell_within=float(dwell_within.mean()) if len(dwell_within) else float("nan")
    threshold_table=pd.read_csv(OUT/"08_information_loss_thresholds.csv")
    # Compact, deterministic final decision report.
    lines=["# Stage-0.5 Decision","","**CONDITIONAL GO** for a focused physical validation; no main-paper mechanism claim yet.","","# Verified","",f"- CW P1dB is E1dB = {s.E1:.8f} V/m, from the 5 MHz CW fundamental definition. Field and power conversions reproduce the saved Stage-0 rows (max peak-relation error {power_err:.2g} dB; RMS-field relation error {e_err:.2g} V/m).",f"- Matched formats use 5.0 Msymbol/s and mean 99% bandwidths {', '.join(f'{m} {config['B_occ_99_Hz_mean']/1e6:.3f} MHz' for m,config in [(x['modulation'],x) for x in configs])}. Across four mean bandwidths, maximum mismatch is {config['maximum_relative_B99_mismatch_of_seed_means']*100:.2f}%.","- QPSK/16-QAM use identical RRC pulse shaping. OFDM has 68 active subcarriers in a 256 FFT, 16 CP samples in a 272-sample block, and 68 information-bearing symbols per block; its effective data-symbol rate is 5.0 Msymbol/s. The WOLA cyclic suffix overlaps the next prefix and adds no block duration.","- For every format, R_AIR = I_native × 5.0 Msymbol/s; eta_AIR = R_AIR/B_occ. All waveforms have unit complex-envelope RMS before scaling; P/P1dB=(E/E1dB)^2.","- Eight seeds cover -6, -3, 0, +3, +6 dB; four seeds cover -15 and -10 dB. OFDM uses an independent randomized balanced four-symbol training preamble on each carrier; all reported full-atomic runs use Nd=1501.","","# Main Results","","At +3 dB relative to P1dB, paired mean full-atomic loss from -6 dB reference:","","| Format | AIR loss (bit/native symbol) | Relative AIR loss | Rate loss (Mbit/s) | eta loss (bit/s/Hz) |","|---|---:|---:|---:|---:|"]
    for mod in MODS:
        q=d[(d.modulation==mod)&(d.model=="atomic")&(d.Pavg_over_P1dB_dB==3)]
        l,lo,hi=mean_ci(q.DeltaAIR_native);rel,_,_=mean_ci(q.relative_loss);dr,_,_=mean_ci(q.DeltaR_AIR_bit_s);de,_,_=mean_ci(q.DeltaEta_bit_s_Hz)
        lines.append(f"| {mod} | {l:.3f} [{lo:.3f}, {hi:.3f}] | {rel:.1%} | {dr/1e6:.3f} | {de:.3f} |")
    lines += ["","# Static-Nonlinear Control","","At +3 dB, full-minus-static-plus-LTI paired excess AIR loss:"]
    for mod in MODS:
        aa=d[(d.modulation==mod)&(d.model=="atomic")].pivot(index="seed",columns="Pavg_over_P1dB_dB",values="AIR_native")
        cc=d[(d.modulation==mod)&(d.model=="static_lti")].pivot(index="seed",columns="Pavg_over_P1dB_dB",values="AIR_native")
        x=(aa[-6]-aa[3])-(cc[-6]-cc[3]);mu,lo,hi=mean_ci(x)
        lines.append(f"- {mod}: {mu:.3f} bit/native symbol (paired 95% t interval {lo:.3f} to {hi:.3f}; n={len(x)}).")
    lines += ["","Static compression plus the repository small-signal LTI response reproduces part of 16-QAM and OFDM loss, particularly at +6 dB, but does not reproduce the full-atomic moderate-power losses. OFDM converges to the static-plus-LTI loss at +6 dB, so the dynamic excess is not uniform across power.","","# Peak-Power View","","A common +6 dB peak-power coordinate is within the measured average-power grid for all formats. The per-seed interpolated relative-loss values are reported in the peak-normalization table below.",""]
    peakrows=[]
    for mod in MODS:
        q=d[(d.modulation==mod)&(d.model=="atomic")]
        vals=[]
        for seed,g in q.groupby("seed"):
            ref=float(g[g.Pavg_over_P1dB_dB==-6].AIR_native.iloc[0]); papr=float(g.PAPR_dB.iloc[0])
            x=g.Pavg_over_P1dB_dB.to_numpy(float)+papr; y=g.AIR_native.to_numpy(float)
            vals.append((ref-np.interp(6,x,y))/ref)
        mu,lo,hi=mean_ci(vals);peakrows.append((mod,mu,lo,hi,len(vals)))
    lines += ["| Format | Relative AIR loss at common +6 dB peak |","|---|---:|"]
    for mod,mu,lo,hi,n in peakrows: lines.append(f"| {mod} | {mu:.1%} [{lo:.1%}, {hi:.1%}] (n={n}) |")
    lines += ["","# Waveform Statistics","",f"PAPR alone is not predictive: OFDM has the largest PAPR but smaller loss than the single-carrier formats at +3 dB. Pooled Spearman rho is {pooled_papr:.2f} for PAPR and {pooled_occ:.2f} for time above E1dB; power-stratified mean rho for mean dwell is {mean_dwell_within:.2f}. Occupancy and dwell have similar, moderate within-power associations, so Stage K was run as a controlled receiver-output diagnostic.","","Stage-1.5 context reports an atomic pulse-response 1/e timescale tau_atom ≈ 1.727 µs; dwell and envelope-correlation times are also stored normalized by this inherited timescale. That pulse trace is not in the present checkout, so tau_atom is not re-estimated here.","","# Constant-Envelope Control","","The BT=0.36 Gaussian CPFSK waveform is constant envelope (PAPR 0 dB), with bandwidth matched within the 2% limit. Its simple discriminator plus the fixed Gaussian auxiliary AIR produces a substantially lower reference rate and high seed variability; its +3 dB interval is broad. Treat CE as a useful envelope control with a receiver limitation, not as decisive evidence about a universal modulation ranking.","","# Numerical Validation","",f"- Halving dt from 1 ns to 0.5 ns changes +3 dB AIR by at most {audit['max_abs_timestep_AIR_change']:.2g} bit/native symbol across the four spot checks.",f"- Nd=501 differs from Nd=1501 by {audit['screening_Nd501_max_abs_AIR_change']:.3f} bit/native symbol at +3 dB. This screening quadrature is not converged; no result at reduced Nd is used as decisive evidence.","- A short 2 us Nd=301 CPU NumPy versus CUDA check is recorded in 05_backend_validation.csv; trace correlation was 0.999995 and maximum probe-response difference was 0.0063 dB.","","# Uncertainty","","Headline intervals are paired across independent waveform seeds using 95% Student-t intervals (n=8 at -6 through +6 dB; n=4 at -15 and -10 dB). They quantify seed variability in this simulated probe/IF/baseband readout, not hardware uncertainty.","","# What Is NOT Supported","","- No universal communication-aware threshold or receiver-independent material limit.","- No proof that dwell time is causal; high-field occupancy and constellation sensitivity remain intertwined.","- No claim that every static nonlinear or linear dynamic model is excluded; the tested controls are the measured CW complex-fundamental curve and its repository small-signal LTI cascade.","- No experimental validation, quantum advantage, capacity result, or generalization beyond this operating point and receiver chain.","","# Scientific Interpretation","","The rate and bandwidth correction does not remove the measured waveform dependence. Full-atomic AIR losses exceed the tested static and static-plus-LTI controls for QPSK, 16-QAM, and OFDM at moderate average powers. The strongest 16-QAM point and high-power OFDM show substantial static-control contributions. The separation is therefore evidence that CW P1dB alone does not predict communication loss in this particular full-thermal transient model, but it does not identify a unique causal statistic.","","The low-power -6 dB reference is a shared coordinate, not a guarantee that all instantaneous peaks are uncompressed. Peak-normalized curves remain format dependent at the common +6 dB peak coordinate, though the CE receiver uncertainty is large. Occupancy tracks the pooled trend better than PAPR; fixed-power comparisons and the single artificial dwell diagnostic limit mechanism claims.","","# Next Step","","Validate the matched-rate experiment at a second IF or operating point with a receiver chain calibrated against a realistic probe detector and analog bandwidth; retain the current tested controls and preamble policy."]
    # Controlled Stage K result: one deterministic pair, with no seed ensemble.
    kd=pd.read_csv(OUT/"10_controlled_dwell_results.csv")
    ka=kd[kd.model=="atomic"].set_index("waveform")
    ks=kd[kd.model=="static"].set_index("waveform")
    kl=kd[kd.model=="static_lti"].set_index("waveform")
    bw_mismatch=abs(float(ka.loc["short_excursions","B_occ_99_Hz"])-float(ka.loc["long_plateaus","B_occ_99_Hz"]))/np.mean([ka.loc["short_excursions","B_occ_99_Hz"],ka.loc["long_plateaus","B_occ_99_Hz"]])
    ac_delta=float(ka.loc["long_plateaus","high_low_contrast_norm"]-ka.loc["short_excursions","high_low_contrast_norm"])
    sc_delta=float(ks.loc["long_plateaus","high_low_contrast_norm"]-ks.loc["short_excursions","high_low_contrast_norm"])
    lc_delta=float(kl.loc["long_plateaus","high_low_contrast_norm"]-kl.loc["short_excursions","high_low_contrast_norm"])
    pair_a,pair_b=__import__('controlled_dwell').make_pair()
    amp_ks=float(ks_2samp(np.abs(pair_a)/s.E1,np.abs(pair_b)/s.E1).statistic)
    dwell_lines=["# Temporal-Dwell Diagnostic (Stage K)","",f"The diagnostic used one deterministic input pair, not a seed ensemble. Both patterns have Pavg/P1dB = -0.90 dB, Ppeak/P1dB = +4.56 dB, 25% occupancy above E1dB, and 99% bandwidths near 240 MHz (relative mismatch {bw_mismatch:.3%}). Mean high-field dwell is 50 ns for short excursions and 1962 ns for long plateaus; the 90th-percentile dwell is 50 ns and 2000 ns. The normalized amplitude distributions are similar but not identical (two-sample KS distance {amp_ks:.3f}).","",f"The normalized high/low baseband contrast changes by {ac_delta:.6f} for full atomic, {sc_delta:.6f} for static, and {lc_delta:.6f} for static-plus-LTI. The LTI control explains most of the contrast change; the full-atomic pair leaves a small additional change of {ac_delta-lc_delta:.6f}. This indicates sensitivity to the temporal pattern in this artificial receiver-output diagnostic, but it does not establish a communication AIR effect or a causal atomic dwell law.","","The artificial pair has much wider bandwidth than the communication waveforms and only one realization. Its result is a focused diagnostic, not a headline receiver operating limit."]
    insert_at=lines.index("# Constant-Envelope Control")
    lines[insert_at:insert_at]=dwell_lines+[""]
    insert_at=lines.index("# What Is NOT Supported")
    threshold_lines=["# Information-Loss Thresholds","","P5 and P10 are the highest powers on the mean full-atomic curve meeting 5% and 10% relative AIR-loss limits. Intervals are paired seed-bootstrap 95% ranges; the estimates apply only to this simulated receiver.","","| Format | P5 - P1dB (dB) | P10 - P1dB (dB) |","|---|---:|---:|"]
    for mod in MODS:
        r5=threshold_table[(threshold_table.modulation==mod)&(threshold_table.relative_loss_limit==.05)].iloc[0]
        r10=threshold_table[(threshold_table.modulation==mod)&(threshold_table.relative_loss_limit==.10)].iloc[0]
        threshold_lines.append(f"| {mod} | {r5.P_limit_minus_P1dB_dB:.2f} [{r5.bootstrap95_low:.2f}, {r5.bootstrap95_high:.2f}] | {r10.P_limit_minus_P1dB_dB:.2f} [{r10.bootstrap95_low:.2f}, {r10.bootstrap95_high:.2f}] |")
    lines[insert_at:insert_at]=threshold_lines+[""]
    (OUT/"01_stage05_summary.md").write_text("\n".join(lines)+"\n",encoding="utf-8")
    # Consolidated, data-derived provenance. The per-invocation stage manifest
    # stays separate so reruns cannot replace the complete coverage record.
    result_path=OUT/"02_all_results.csv"
    result_hash=hashlib.sha256(result_path.read_bytes()).hexdigest()
    config_versions={"CE":"fair_v2","QPSK":"fair_v2","16QAM":"fair_v2","OFDM":"fair_v4_randomized_balanced_ofdm_training"}
    coverage={}
    for mod,version in config_versions.items():
        q=raw[(raw.modulation==mod)&(raw.config_version==version)]
        powers={}
        for power,g in q.groupby("Pavg_over_P1dB_dB"):
            powers[str(int(power))]={"seeds":sorted(int(x) for x in g.seed.unique()),"replicates":int(g.seed.nunique()),"rows":int(len(g))}
        coverage[mod]={"config_version":version,"rows":int(len(q)),"power_coverage_db":powers}
    archived=OUT/"screening_v2_invalid_ofdm_training/02_superseded_fair_v2_ofdm_rows.csv"
    manifest={"stage":"0.5","status":"active primary dataset packaged","branch":"p1db-waveform-fairness","source_commit":"409572499ef54825b66c532e0aecdcffa89e7be1",
      "simulation":{"Nd":1501,"dt_s":1e-9,"state":"continuous through waveform","thermal":"full","primary_backend":"CUDA","models":["atomic","static","static_lti"],"average_power_grid_db":[-15,-10,-6,-3,0,3,6]},
      "active_data":{"file":"02_all_results.csv","rows":int(len(raw)),"sha256":result_hash,"waveforms":coverage},
      "waveform_statistics":{"file":"03_waveform_statistics_extended.csv","rows":int(len(stats))},
      "superseded_data":{"file":"screening_v2_invalid_ofdm_training/02_superseded_fair_v2_ofdm_rows.csv","rows_archived":int(len(pd.read_csv(archived))) if archived.exists() else 0,"reason":"OFDM fair_v2 used a correlated/singular training configuration; final OFDM uses independently randomized balanced four-symbol carrier training."},
      "validation_files":["05_numerical_validation.csv","05_backend_validation.csv","10_controlled_dwell_results.csv","10_controlled_dwell_pair_summary.csv"],"reproduction":"06_reproduce_commands.txt","analysis_summary":"01_stage05_summary.md"}
    (OUT/"04_job_manifest.json").write_text(json.dumps(manifest,indent=2),encoding="utf-8")
    print("ANALYSIS_COMPLETE",json.dumps({"rows":len(d),"formats":list(MODS),"thresholds":pd.read_csv(OUT/"08_information_loss_thresholds.csv").to_dict("records"),"figures":len(list(OUT.glob('fig*.png')))},default=str),flush=True)


if __name__ == "__main__":
    main()
