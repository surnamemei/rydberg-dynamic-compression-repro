from pathlib import Path
import os,sys,json,time
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/"python"))
_CUDA_DLL_DIR = os.environ.get("RYDBERG_CUDA_DLL_DIR", str(ROOT / ".cuda-env" / "Library" / "bin"))  # Windows CUDA runtime DLLs
if hasattr(os, "add_dll_directory") and os.path.isdir(_CUDA_DLL_DIR):  # Linux finds the runtime via the rpath of cu_ryd.so
    os.add_dll_directory(_CUDA_DLL_DIR)
from utils.transient_quantum import SimConfig,TransientQuantumSimulator,configure_raqr
amps=np.array([0,.005,.01,.02,.03,.04,.05,.06,.07,.08,.09,.10,.12,.15,.20,.25,.30,.40,.50,.60,.70])
dt=1e-9;n=40000;t=np.arange(n)*dt;f=5e6
raqr=configure_raqr("Transit");cfg=SimConfig(300,n,dt,n*dt,t);sim=TransientQuantumSimulator()
harm=[];means=[]
for a in amps:
    field=raqr.A_LO+a*np.exp(2j*np.pi*f*t)
    start=time.perf_counter()
    res=sim.run(raqr,cfg,field,Nd=1501,device="cuda")
    intensity=10**(res.probeResponse/10)
    tail=intensity[-4000:]
    c=2*np.mean(tail*np.exp(-2j*np.pi*f*t[-4000:]))
    harm.append(c);means.append(tail.mean())
    print("CW",a,abs(c),np.angle(c),time.perf_counter()-start,flush=True)
harm=np.array(harm);means=np.array(means)
slope=np.mean(np.abs(harm[1:3])/amps[1:3])
gain=np.full(len(amps),np.nan);gain[1:]=20*np.log10(np.abs(harm[1:])/(slope*amps[1:]))
j=np.flatnonzero(gain[1:]<=-1)[0]+1
x=np.interp(-1,gain[j-1:j+1][::-1],amps[j-1:j+1][::-1])
print("P1DB",json.dumps({"E_1dB_Vpm":float(x),"P1dB_relative_field_squared":float(x*x),"IF_Hz":f,"LO_Vpm":raqr.A_LO,"slope_intensity_per_Vpm":float(slope),"gain_cross_pair":gain[j-1:j+1].tolist()}),flush=True)
np.savez_compressed(ROOT/"artifacts"/"p1db_reference.npz",amps=amps,harm=harm,means=means,gain_db=gain,e1db=x,slope=slope)
fig,ax=plt.subplots(1,2,figsize=(9,4))
ax[0].plot(amps[1:],np.abs(harm[1:]),"o-",label="CW fundamental")
ax[0].plot(amps[1:],slope*amps[1:],"--",label="small-signal extrapolation")
ax[0].axvline(x,color="k",ls=":",label=f"E1dB={x:.4f} V/m")
ax[0].set_xlabel("CW field amplitude (V/m)");ax[0].set_ylabel("probe intensity fundamental amplitude");ax[0].legend()
ax[1].plot(amps[1:],gain[1:],"o-");ax[1].axhline(-1,color="k",ls="--");ax[1].axvline(x,color="k",ls=":")
ax[1].set_xlabel("CW field amplitude (V/m)");ax[1].set_ylabel("fundamental gain relative to weak tone (dB)")
for a in ax:a.grid(alpha=.2)
fig.suptitle("Transit RAQR, LO 0.5 V/m, IF 5 MHz, Nd=1501, 40 us CW")
fig.tight_layout();fig.savefig(ROOT/"p1db_reference_curve.png",dpi=180)


