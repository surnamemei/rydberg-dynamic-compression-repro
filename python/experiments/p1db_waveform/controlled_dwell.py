from __future__ import annotations

import csv
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT/"python"/"experiments"/"p1db_waveform"))
import stage05 as s

OUT=ROOT/"results"/"p1db_waveform_stage05"
N=204000
HIGH=1.5*s.E1
LOW=.5*s.E1


def make_pair():
    # Shared two-level base pattern: 25% high, 75% low before the common quadrature dither.
    # Short pattern: 50 ns high / 150 ns low. Long pattern: 2 us high / 6 us low.
    a=np.tile(np.r_[np.full(50,HIGH),np.full(150,LOW)],N//200)
    b=np.tile(np.r_[np.full(2000,HIGH),np.full(6000,LOW)],25)
    b=np.r_[b,np.full(1000,HIGH),np.full(3000,LOW)]
    assert len(a)==N and len(b)==N
    # A shared period-200 ns quadrature multisine supplies a common occupied band
    # while remaining below the high/low threshold boundaries.
    rng=np.random.default_rng(20260926)
    X=np.zeros(200//2+1,complex)
    amps=np.ones(24); amps[-2:]=2.5
    X[1:25]=amps*np.exp(1j*rng.uniform(0,2*np.pi,24))
    d=np.fft.irfft(X,n=200)
    d=d/np.std(d)*(.25*s.E1)
    shift=25-int(np.argmax(np.abs(d)))
    d=np.roll(d,shift)
    dq=np.tile(d,N//200)
    return a.astype(complex)+1j*dq,b.astype(complex)+1j*dq


def occupied99(x):
    return s.occupied99(x)


def dwell(values):
    mask=np.abs(values)>s.E1*(1+1e-12)
    edge=np.diff(np.r_[False,mask,False].astype(np.int8))
    starts,ends=np.flatnonzero(edge==1),np.flatnonzero(edge==-1)
    return (ends-starts)*1e-9


def output_bb(intensity):
    t=np.arange(len(intensity))/s.FS
    return s.old.sosfilt(s.old.sos,(intensity-s.old.dc)*np.exp(-2j*np.pi*s.old.IF*t))[::s.UPSAMPLE]


def main():
    previous=OUT/"10_controlled_dwell_results.csv"
    if previous.exists():
        old=__import__('pandas').read_csv(previous)
        bw=old.groupby("waveform").B_occ_99_Hz.first()
        mismatch=abs(float(bw.iloc[0])-float(bw.iloc[1]))/np.mean(bw.to_numpy())
        if mismatch>.02:
            archive=OUT/"10_controlled_dwell_unmatched_spectrum.csv"
            serial=2
            while archive.exists():
                archive=OUT/f"10_controlled_dwell_unmatched_spectrum_{serial}.csv";serial+=1
            previous.replace(archive)
    pair_a,pair_b=make_pair();waves={"short_excursions":pair_a,"long_plateaus":pair_b}
    dt=1e-9;t=np.arange(N)*dt
    details=[];traces={}
    for name,unit in waves.items():
        inenv=unit
        intensity={
            "atomic":s.atomic_intensity(inenv,nd=1501),
            "static":s.static_intensity(inenv),
            "static_lti":s.static_lti_intensity(inenv),
        }
        traces[name]={"input":np.abs(inenv)[::s.UPSAMPLE]/s.E1}
        for model,y in intensity.items():
            bb=output_bb(y);traces[name][model]=np.abs(bb)/s.E1
            high=(np.abs(inenv)[::s.UPSAMPLE]>s.E1)
            low=~high
            # Flat-top summaries omit 100 ns after each transition for the long case;
            # short pulses have no flat interior, so all samples are summarized.
            vals={"bb_rms_norm":float(np.sqrt(np.mean(np.abs(bb)**2))/s.E1),"bb_mean_abs_norm":float(np.mean(np.abs(bb))/s.E1),"high_mean_abs_norm":float(np.mean(np.abs(bb[high]))/s.E1),"low_mean_abs_norm":float(np.mean(np.abs(bb[low]))/s.E1),"high_low_contrast_norm":float((np.mean(np.abs(bb[high]))-np.mean(np.abs(bb[low])))/s.E1)}
            details.append({"waveform":name,"model":model,"Nd":1501 if model=="atomic" else "control","dt_ns":1.0,"duration_us":N/1000,"E_rms_Vpm":float(np.sqrt(np.mean(np.abs(inenv)**2))),"E_peak_Vpm":float(np.max(np.abs(inenv))),"Pavg_over_P1dB_dB":float(10*np.log10(np.mean(np.abs(inenv)**2)/s.E1**2)),"Ppeak_over_P1dB_dB":float(20*np.log10(np.max(np.abs(inenv))/s.E1)),"PAPR_dB":float(10*np.log10(np.max(np.abs(inenv)**2)/np.mean(np.abs(inenv)**2))),"B_occ_99_Hz":occupied99(inenv),"high_occupancy":float(np.mean(np.abs(inenv)>s.E1*(1+1e-12))),"high_dwell_mean_ns":float(np.mean(dwell(inenv))*1e9),"high_dwell_p90_ns":float(np.quantile(dwell(inenv),.9)*1e9),"high_dwell_max_ns":float(np.max(dwell(inenv))*1e9),**vals})
    df=__import__('pandas').DataFrame(details)
    config={"duration_s":N/1e9,"high_level_before_dither_Vpm":HIGH,"low_level_before_dither_Vpm":LOW,"short_pattern":{"high_ns":50,"low_ns":150,"period_ns":200},"long_pattern":{"high_us":2,"low_us":6,"period_us":8,"complete_periods":25,"remainder_high_us":1,"remainder_low_us":3},"shared_quadrature_dither":{"period_ns":200,"fourier_bins":list(range(1,25)),"outer_bin_amplitude_scale":2.5,"rms_fraction_of_E1dB":0.25,"random_seed":20260926,"roll_max_to_sample":25},"threshold_E1dB_Vpm":s.E1,"same_average_power_peak_and_B99_target":True,"full_atomic_Nd":1501,"dt_s":1e-9}
    (OUT/"10_controlled_dwell_config.json").write_text(__import__('json').dumps(config,indent=2),encoding="utf-8")
    df.to_csv(OUT/"10_controlled_dwell_results.csv",index=False)
    # Same-input controls make a waveform-pair delta directly interpretable.
    wide=df.pivot(index="model",columns="waveform",values=["bb_rms_norm","high_low_contrast_norm"])
    wide.columns=[f"{a}_{b}" for a,b in wide.columns]
    wide.reset_index().to_csv(OUT/"10_controlled_dwell_pair_summary.csv",index=False)
    fig,axes=plt.subplots(2,1,figsize=(10,7),sharex=True)
    time_us=np.arange(N//s.UPSAMPLE)[::2]*0.1
    for ax,name in zip(axes,("short_excursions","long_plateaus")):
        ax.plot(time_us,traces[name]["input"][::2],color="black",lw=1.5,label="input |E|/E1dB")
        for model,style in (("atomic","-"),("static","--"),("static_lti",":")):
            ax.plot(time_us,traces[name][model][::2],style,lw=1.0,label=f"{model} |baseband output|/E1dB")
        ax.set_xlim(0,20);ax.set_ylabel(name.replace("_"," "));ax.grid(alpha=.2)
    axes[-1].set_xlabel("Time (µs)");axes[0].legend(ncol=2,fontsize=8,frameon=False)
    fig.suptitle("Controlled dwell diagnostic: same RMS, peak, and 99% bandwidth")
    fig.tight_layout();fig.savefig(OUT/"fig17_controlled_dwell_experiment.png",dpi=180);plt.close(fig)
    print(df.to_string(index=False),flush=True)


if __name__=="__main__":main()





