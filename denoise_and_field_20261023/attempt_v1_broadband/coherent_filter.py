"""Conservative long-aperture linear-coherence subtraction, data only.

Research prototype, not an implementation of a published named algorithm.
No clean data, true Q, event times or synthetic interference trajectories.
"""
import numpy as np
PARAMETERS=dict(minimum_apparent_velocity_m_s=400.,slope_count=1201,minimum_peak_semblance=.65,gate_onset=.45,maximum_components=3,subtraction_fraction=.60,scan_dt_min_s=.004,minimum_aperture_coverage=.70)

def _smooth(a,n):return np.convolve(a,np.ones(n)/n,mode='same')

def suppress_coherent(D,dt,x):
 D=np.asarray(D,dtype=float);x=np.asarray(x,dtype=float).reshape(-1)
 assert D.ndim==2 and D.shape[1]==len(x) and np.isfinite(D).all() and np.all(np.diff(x)>0)
 nt,nx=D.shape;t=np.arange(nt)*dt;live=np.sqrt(np.mean(D*D,axis=0))>max(np.max(np.abs(D))*1e-12,np.finfo(float).tiny)
 if live.sum()<8 or not np.any(D):return D.copy(),np.zeros_like(D),dict(reason='Too few live traces or zero data',components=[])
 stride=max(1,int(round(PARAMETERS['scan_dt_min_s']/dt)));ts=t[::stride];dts=dt*stride
 # Anti-alias the scan only; output is the unchanged input minus its model.
 nfft=2**int(np.ceil(np.log2(nt*2)));freq=np.fft.rfftfreq(nfft,dt);F=np.fft.rfft(D,n=nfft,axis=0)
 ps=np.mean(np.abs(F[:,live])**2,axis=1);sel=(freq>=2)&(freq<=min(.45/dt,120));cum=np.cumsum(np.where(sel,ps,0));f90=float(freq[np.searchsorted(cum,.9*cum[-1])]) if cum[-1]>0 else 30.
 cutoff=min(.40/dts,max(15.,2*f90));taper=np.ones_like(freq);a=(freq-cutoff*.8)/max(.2*cutoff,1e-9);taper[freq>cutoff*.8]=.5*(1+np.cos(np.pi*np.clip(a[freq>cutoff*.8],0,1)));taper[freq>=cutoff]=0
 scan=np.fft.irfft(F*taper[:,None],n=nfft,axis=0)[:nt:stride]
 xc=float(np.mean(x[live]));xx=x-xc;aperture=float(np.ptp(x[live]));pmin=2/max(f90*aperture,1)
 slopes=np.linspace(-1/PARAMETERS['minimum_apparent_velocity_m_s'],1/PARAMETERS['minimum_apparent_velocity_m_s'],PARAMETERS['slope_count']);eligible=np.abs(slopes)>=pmin
 width=max(5,int(round(max(.024,1.5/max(f90,5))/dts)));width+=1-width%2
 removed=np.zeros_like(D);work=scan.copy();excluded=np.zeros((len(slopes),len(ts)),bool);components=[]
 def aligned(y,p):
  a=np.zeros_like(y);valid=np.zeros_like(y,dtype=bool)
  for j in np.where(live)[0]:
   at=ts+p*xx[j];valid[:,j]=(at>=ts[0])&(at<=ts[-1]);a[:,j]=np.interp(at,ts,y[:,j],left=0,right=0)
  count=valid.sum(axis=1);total=a.sum(axis=1);energy=np.sum(a*a,axis=1)
  num=_smooth(total*total,width);den=_smooth(count*energy,width)
  sem=np.divide(num,den,out=np.zeros_like(num),where=den>0);coverage=count/live.sum()
  return a,count,sem,coverage,_smooth(energy/np.maximum(count,1),width)
 for iteration in range(PARAMETERS['maximum_components']):
  best=None
  for ip,p in enumerate(slopes):
   if not eligible[ip]:continue
   aa,count,sem,coverage,power=aligned(work,p)
   ok=(coverage>=PARAMETERS['minimum_aperture_coverage'])&(sem>=PARAMETERS['minimum_peak_semblance'])&~excluded[ip]
   # Coherence plus actual energy, not coherence of a numerically quiet tail.
   score=np.where(ok,sem**2*power,0);k=int(np.argmax(score));value=float(score[k])
   if value>0 and (best is None or value>best[0]):best=(value,ip,k,aa,count,sem,coverage,power)
  if best is None:break
  _,ip,k,aa,count,sem,coverage,power=best;p=float(slopes[ip]);tau=float(ts[k])
  gate=np.clip((sem-PARAMETERS['gate_onset'])/(PARAMETERS['minimum_peak_semblance']-PARAMETERS['gate_onset']),0,1)
  # Only the detected packet, not every event sharing its slope.
  half=max(.06,3/max(f90,10));window=np.exp(-.5*((ts-tau)/(half/2))**2);gate*=window
  stacked=np.sum(aa,axis=1)/np.maximum(count,1);wave=stacked*gate
  model=np.zeros_like(D)
  for j in np.where(live)[0]:model[:,j]=np.interp(t-p*xx[j],ts,wave,left=0,right=0)
  model*=PARAMETERS['subtraction_fraction'];removed+=model
  # The same packet cannot be repeatedly removed up to 94% by 3 iterations.
  radius=max(3*(slopes[1]-slopes[0]),2/max(f90*aperture,1));excluded|=(np.abs(slopes[:,None]-p)<=radius)&(np.abs(ts[None,:]-tau)<=2*half)
  work=scan-removed[::stride]
  components.append(dict(slope_s_per_m=p,apparent_velocity_m_s=1/abs(p),time_at_aperture_center_s=tau,semblance=float(sem[k]),aperture_center_m=xc,packet_half_width_s=half))
 out=D-removed;assert np.max(np.abs(out+removed-D))<=1e-12*max(1,np.max(np.abs(D)))
 assert np.all(out[:,~live]==D[:,~live])
 return out,removed,dict(parameters=PARAMETERS,components=components,live_traces=int(live.sum()),f90_hz=f90,scan_dt_s=dts,removed_norm_fraction=float(np.linalg.norm(removed)/np.linalg.norm(D)),purpose='Estimate long straight coherent packets. May remove true linear events; inspection required, no default promotion without clean-signal gate.')
