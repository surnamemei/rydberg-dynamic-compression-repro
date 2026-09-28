from __future__ import annotations
from pathlib import Path
import os,sys,csv,json,argparse,time
import numpy as np
from scipy.signal import resample_poly,fftconvolve,butter,sosfilt
from scipy.special import logsumexp
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/"python"))
_CUDA_DLL_DIR = os.environ.get("RYDBERG_CUDA_DLL_DIR", str(ROOT / ".cuda-env" / "Library" / "bin"))  # Windows CUDA runtime DLLs
if hasattr(os, "add_dll_directory") and os.path.isdir(_CUDA_DLL_DIR):  # Linux finds the runtime via the rpath of cu_ryd.so
    os.add_dll_directory(_CUDA_DLL_DIR)
from utils.transient_quantum import SimConfig,TransientQuantumSimulator,configure_raqr
from utils.sim_SingleCarrier import SingleCarrier
CW=np.load(ROOT/"artifacts"/"p1db_reference.npz")
E1=float(CW["e1db"]);AMPS=CW["amps"];HARM=CW["harm"]
FS=1e9;FS20=20e6;IF=5e6;ND=1501;NS20=2880;UPSAMPLE=50
SIGMA=float(np.load(ROOT/"artifacts"/"stage1_5_frozen.npz")["noise_sigma"])
raqr=configure_raqr("Transit");sim=TransientQuantumSimulator()
sos=butter(8,3e6,fs=FS,output="sos")
dc=float(CW["means"][0])
LTI=np.load(ROOT/"artifacts"/"p1db_lti_kernel.npz")
HR,HI=LTI["hR"],LTI["hI"]
QPSK=np.exp(1j*(np.pi/4+np.pi/2*np.arange(4)))
grid=np.array([-3,-1,1,3])
QAM=(grid[:,None]+1j*grid[None,:]).ravel()/np.sqrt(10)
BITS4=np.array([[0,0],[0,1],[1,1],[1,0]],dtype=int)
BITS_QPSK=BITS4
BITS_QAM=np.array([np.r_[BITS4[i],BITS4[j]] for i in range(4) for j in range(4)])
RRC=SingleCarrier._rrcosdesign(None,beta=.2,span=10,sps=4)
ACTIVE=np.r_[np.arange(1,10),np.arange(64-9,64)]
def waveform(mod,seed):
    rng=np.random.default_rng(seed)
    if mod in ("QPSK","16QAM"):
        alphabet=QPSK if mod=="QPSK" else QAM
        labels=rng.integers(0,len(alphabet),720)
        up=np.zeros(NS20,dtype=complex);up[::4]=alphabet[labels]
        tx20=fftconvolve(up,RRC,mode="same")
        shape=(720,)
    else:
        alphabet=QPSK
        labels=rng.integers(0,4,(36,18))
        freq=np.zeros((36,64),dtype=complex)
        freq[:,ACTIVE]=alphabet[labels]
        timed=np.fft.ifft(freq,axis=1)*np.sqrt(64)
        # Weighted overlap-add OFDM: the taper is confined to CP and cyclic suffix.
        window=np.ones(88)
        ramp=np.sin(np.pi/2*(np.arange(8)+.5)/8)**2
        window[:8]=ramp
        window[-8:]=ramp[::-1]
        tx20=np.zeros(36*80+8,dtype=complex)
        for j in range(36):
            block=np.r_[timed[j,-16:],timed[j],timed[j,:8]]
            tx20[j*80:j*80+88]+=block*window
        tx20=tx20[:NS20]
        shape=(36,18)
    tx=resample_poly(tx20,UPSAMPLE,1)
    tx=tx/np.sqrt(np.mean(np.abs(tx)**2))
    assert len(tx)==NS20*UPSAMPLE
    return tx,labels,alphabet,shape
def extract(intensity,mod):
    n=len(intensity);t=np.arange(n)/FS
    bb=sosfilt(sos,(intensity-dc)*np.exp(-2j*np.pi*IF*t))[::UPSAMPLE]
    if mod=="OFDM":
        frame=bb.reshape(36,80)[:,16:]
        return (np.fft.fft(frame,axis=1)/np.sqrt(64))[:,ACTIVE]
    matched=fftconvolve(bb,RRC,mode="same")
    return matched[6::4]
def static_intensity(env):
    a=np.abs(env)
    c=np.interp(a,AMPS,HARM.real)+1j*np.interp(a,AMPS,HARM.imag)
    out=np.zeros(len(env),dtype=complex)
    np.divide(c*env,a,out=out,where=a>0)
    t=np.arange(len(env))/FS
    return dc+np.real(out*np.exp(2j*np.pi*IF*t))
def linear_intensity(env):
    t=np.arange(len(env))/FS
    field=env*np.exp(2j*np.pi*IF*t)
    return dc+fftconvolve(field.real,HR)[:len(env)]+fftconvolve(field.imag,HI)[:len(env)]
def atomic_intensity(env):
    n=len(env);t=np.arange(n)/FS
    cfg=SimConfig(300,n,1/FS,n/FS,t)
    res=sim.run(raqr,cfg,raqr.A_LO+env*np.exp(2j*np.pi*IF*t),Nd=ND,device="cuda")
    return 10**(res.probeResponse/10)
def fit_gain(x,y):
    xm=np.mean(x);ym=np.mean(y)
    g=np.vdot(x-xm,y-ym)/np.vdot(x-xm,x-xm)
    return g,ym-g*xm
def ll(z,alphabet,var):
    return -np.abs(z[...,None]-alphabet[None,...])**2/var
def info(z,labels,alphabet,var):
    flat=z.ravel();truth=labels.ravel();logq=ll(flat,alphabet,var)
    density=np.log2(len(alphabet))+(logq[np.arange(len(flat)),truth]-logsumexp(logq,axis=1))/np.log(2)
    return float(np.mean(density))
def evaluate(raw,labels,alphabet,mod):
    if mod!="OFDM": labels=labels[:len(raw)]
    x=alphabet[labels]
    if mod=="OFDM":
        g=np.empty(18,dtype=complex);b=np.empty(18,dtype=complex)
        for j in range(18):g[j],b[j]=fit_gain(x[:12,j],raw[:12,j])
        z=(raw-b)/g
        cal=z[12:24];labcal=labels[12:24]
        test=z[24:];labtest=labels[24:]
    else:
        g,b=fit_gain(x[:240],raw[:240])
        z=(raw-b)/g
        cal=z[240:480];labcal=labels[240:480]
        test=z[480:];labtest=labels[480:]
    base=max(float(np.mean(np.abs(cal-alphabet[labcal])**2)),1e-12)
    var,s=base,1.0
    actual=info(test,labtest,alphabet,var)
    pred=np.argmin(np.abs(test.ravel()[:,None]-alphabet[None,:])**2,axis=1)
    truth=labtest.ravel()
    bits=BITS_QAM if mod=="16QAM" else BITS_QPSK
    return {"AIR":actual,"SER":float(np.mean(pred!=truth)),"BER":float(np.mean(bits[pred]!=bits[truth])),"EVM":float(np.sqrt(np.mean(np.abs(test-alphabet[labtest])**2))),"gmi_s":float(s),"cal_var":float(var),"gain_mag":float(np.mean(np.abs(g)))}
def run(mod,seed,points):
    tx,labels,alphabet,shape=waveform(mod,seed)
    papr=float(np.max(np.abs(tx)**2)/np.mean(np.abs(tx)**2))
    stats={"modulation":mod,"seed":seed,"rms_normalized":1.0,"peak_normalized":float(np.max(np.abs(tx))),"PAPR_dB":10*np.log10(papr),"PAPR_linear":papr,"sample_rate_Hz":FS20,"information_symbol_rate_Hz":5e6 if mod!="OFDM" else 18/(80/FS20),"occupied_bandwidth_nominal_Hz":6e6 if mod!="OFDM" else 18*(FS20/64),"n_information_symbols":int(np.size(labels)),"constellation_size":len(alphabet)}
    rows=[]
    rng=np.random.default_rng(seed+9991)
    noise=rng.normal(0,SIGMA,len(tx))
    for db in points:
        E=E1*10**(db/20)
        env=E*tx
        start=time.perf_counter()
        atom=atomic_intensity(env)
        stat=static_intensity(env)
        linear=linear_intensity(env)
        if not np.isfinite(atom).all():raise RuntimeError(f"Nonfinite atomic output: {mod} {db}")
        for model,intensity in (("atomic",atom),("static",stat),("linear",linear)):
            rx=extract(intensity+noise,mod)
            result=evaluate(rx,labels,alphabet,mod)
            row={"modulation":mod,"seed":seed,"model":model,"Pavg_over_P1dB_dB":db,"Ppeak_over_P1dB_dB":db+stats["PAPR_dB"],"E_rms_Vpm":E,"E_peak_Vpm":E*stats["peak_normalized"],"PAPR_dB":stats["PAPR_dB"],"n_test":int(len(labels)-480 if mod!="OFDM" else 12*18),**result}
            rows.append(row)
        print("CASE",mod,seed,db,"atomic",rows[-3]["AIR"],"static",rows[-2]["AIR"],"linear",rows[-1]["AIR"],"seconds",time.perf_counter()-start,flush=True)
    return stats,rows
if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("--full",action="store_true");p.add_argument("--replicates",type=int,default=1)
    args=p.parse_args();points=[-15,-6,0,3] if not args.full else [-15,-10,-6,-3,0,3,6]
    rows=[];stats=[]
    for mod in ("QPSK","16QAM","OFDM"):
        for rep in range(args.replicates):
            st,rr=run(mod,20261001+100*rep+{"QPSK":1,"16QAM":2,"OFDM":3}[mod],points)
            stats.append(st);rows.extend(rr)
    with open(ROOT/("p1db_comm_stage0_pilot.csv" if not args.full else "p1db_comm_stage0_full.csv"),"w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    with open(ROOT/("waveform_statistics_pilot.csv" if not args.full else "waveform_statistics_full.csv"),"w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(stats[0]));w.writeheader();w.writerows(stats)









