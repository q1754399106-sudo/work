import numpy as np
from coherent_filter import suppress_coherent,PARAMETERS
MULTIBAND_PARAMETERS=dict(band_edges_relative_to_observed_f90=[0,.4,.8,1.3,2.0],trace_rms_preconditioning='For prediction only; restore physical amplitudes before subtraction',base_parameters=PARAMETERS)
def suppress_multiband(D,dt,x):
 D=np.asarray(D,dtype=float);n,nx=D.shape;N=2**int(np.ceil(np.log2(2*n)));f=np.fft.rfftfreq(N,dt);F=np.fft.rfft(D,n=N,axis=0);live=np.sum(D*D,axis=0)>0
 if live.sum()<8:return D.copy(),np.zeros_like(D),dict(bands=[],reason='Insufficient live traces')
 power=np.mean(np.abs(F[:,live])**2,axis=1);sel=(f>=2)&(f<=min(.45/dt,120));cum=np.cumsum(np.where(sel,power,0));f90=float(f[np.searchsorted(cum,.9*cum[-1])]) if cum[-1]>0 else 30
 edges=f90*np.array(MULTIBAND_PARAMETERS['band_edges_relative_to_observed_f90']);edges[-1]=min(edges[-1],.48/dt)
 # Smooth nested low-passes; differences form an exact partition of unity.
 low=[]
 for cutoff in edges[1:]:
  w=.5*(1+np.cos(np.pi*np.clip((f-.85*cutoff)/(.30*cutoff),0,1)));w[f>=1.15*cutoff]=0;low.append(w)
 masks=[low[0]]+[low[j]-low[j-1] for j in range(1,len(low))]+[1-low[-1]]
 assert np.max(np.abs(np.sum(masks,axis=0)-1))<1e-12
 removed=np.zeros_like(D);diags=[]
 scales=np.sqrt(np.mean(D*D,axis=0));scales=np.where(live,np.maximum(scales,np.median(scales[live])*.05),1)
 for j,mask in enumerate(masks):
  if j==len(masks)-1:continue # Preserve the highest-frequency remainder.
  band=np.fft.irfft(F*mask[:,None],n=N,axis=0)[:n];normalized=band/scales[None,:]
  _,model,diag=suppress_coherent(normalized,dt,x);model*=scales[None,:];removed+=model;diags.append(dict(band=j+1,upper_edge_hz=float(edges[j+1]),diagnostics=diag))
 out=D-removed;assert np.max(np.abs(out+removed-D))<1e-12*max(1,np.max(np.abs(D)));assert np.all(out[:,~live]==D[:,~live])
 return out,removed,dict(parameters=MULTIBAND_PARAMETERS,f90_hz=f90,bands=diags,removed_norm_fraction=float(np.linalg.norm(removed)/np.linalg.norm(D)))
