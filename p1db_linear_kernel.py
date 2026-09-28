from pathlib import Path
import os,sys
import numpy as np
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/"python"))
_CUDA_DLL_DIR = os.environ.get("RYDBERG_CUDA_DLL_DIR", str(ROOT / ".cuda-env" / "Library" / "bin"))  # Windows CUDA runtime DLLs
if hasattr(os, "add_dll_directory") and os.path.isdir(_CUDA_DLL_DIR):  # Linux finds the runtime via the rpath of cu_ryd.so
    os.add_dll_directory(_CUDA_DLL_DIR)
from utils.transient_quantum import SimConfig,TransientQuantumSimulator,configure_raqr
n=15000;pulse=1000;amp=2.0;dt=1e-9;t=np.arange(n)*dt
raqr=configure_raqr("Transit");cfg=SimConfig(300,n,dt,n*dt,t);sim=TransientQuantumSimulator()
def response(spike):
 field=np.full(n,raqr.A_LO,dtype=complex);field[pulse]+=spike
 out=sim.run(raqr,cfg,field,Nd=1501,device="cuda")
 return 10**(out.probeResponse/10)
hp=response(amp);hm=response(-amp);ip=response(1j*amp);im=response(-1j*amp)
hR=(hp[pulse:]-hm[pulse:])/(2*amp)
hI=(ip[pulse:]-im[pulse:])/(2*amp)
np.savez_compressed(ROOT/"artifacts"/"p1db_lti_kernel.npz",hR=hR,hI=hI,pulse_amplitude_Vpm=amp,dt_s=dt,LO_Vpm=raqr.A_LO,Nd=1501)
print("LTI_KERNEL",len(hR),float(np.max(np.abs(hR))),float(np.max(np.abs(hI))),flush=True)

