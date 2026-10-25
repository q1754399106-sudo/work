from pathlib import Path
import numpy as np,json,hashlib,shutil,datetime
ROOT=Path('C:\\Users\\qq175\\Documents\\ChatGPT\\研究\\q_identifiable_cluster_20261023\\generated\\FreshThinLowFrequency');BASE=ROOT.parent
OLD=BASE/'unavailable_original_source'
for folder in [ROOT/'inputs',ROOT/'truth',ROOT/'analysis',ROOT/'figures']:folder.mkdir(exist_ok=True)
if (OLD/'core_p1').exists():
 shutil.copytree(OLD/'core_p1',ROOT/'core_p1',dirs_exist_ok=True)
 shutil.copy2(OLD/'scripts/run_local.py',ROOT/'scripts/run_local.py')
 hashes={p.relative_to(OLD/'core_p1').as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in (OLD/'core_p1').rglob('*.m')}
else:
 # The delivered experiment already includes the complete frozen solver.
 hashes=json.loads((ROOT/'FROZEN_SOLVER_HASHES.json').read_text(encoding='utf-8'))
assert all(hashlib.sha256((ROOT/'core_p1'/p).read_bytes()).hexdigest()==h for p,h in hashes.items())
(ROOT/'FROZEN_SOLVER_HASHES.json').write_text(json.dumps(hashes,indent=2),encoding='utf-8')
plan=dict(name='ComplexFieldLike_OneGeology',registered_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),seed=20262901,nx=99,dx_m=30,nz=441,dz_m=5,nt=851,dt_s=.002,offset_m=300,
 geology='60 irregular beds, folding, two faults, channel erosion, pinchout wedge, gas/low-Q lens, smooth heterogeneity',
 wavelet='Spatially varying 36 Hz +/-7% Ricker, phase 24 +/-7 degrees, delayed source ghost -0.16 at 12 ms and +0.04 at 25 ms',
 forward='Column-wise 2D-property quasi-1D convolution; cumulative vertical dt/Q multiplied by straight/RMS finite-offset path factor; constant-Q amplitude/dispersion. Kinematic diffraction/multiples. Not 2D elastic/viscoelastic full-wave simulation.',
 calibration='Independent 18 dB known-four-layer calibration, using spatial-average interval properties and the same mixed-phase source at x=1470 m. True Q not passed to solvers. One shared estimated calibration transferred to all lateral positions.',
 observations=[dict(name='Moderate',snr_db=14,seed=20262911,coherent_energy_fraction=.18,multiple_scale=1.,static_std_ms=1.5),dict(name='Strong',snr_db=8,seed=20262912,coherent_energy_fraction=.30,multiple_scale=1.5,static_std_ms=3.)],
 noise='AR(1)=0.55 temporally colored/spatially correlated Gaussian noise, ground roll, linear interference, impulsive bursts, two dead traces (1-based 23 and 78)',
 scoring='Full and live-trace relative waveform error; waveform correlation; deep window t>=1.05s; low-energy truth-local-RMS 20th-45th percentile region (evaluation only). No individual amplitude claim for overlapping dense reflectors.',
 methods=['FixedSplit','FixedSplitSharedCal','FullCalRefit','PrunedFullCalRefit'],solver_policy='All original production MATLAB files frozen; process independent adjacent blocks of 11 Section traces, same 450 m offset; reuse prepared calibration/dictionary. No new thresholds selected on this experiment.')
(ROOT/'PLAN.json').write_text(json.dumps(plan,ensure_ascii=False,indent=2),encoding='utf-8')
rng=np.random.default_rng(plan['seed']);x=np.arange(plan['nx'])*plan['dx_m'];z=np.arange(plan['nz'])*plan['dz_m'];t=np.arange(plan['nt'])*plan['dt_s'];X,Z=np.meshgrid(x,z)
fold=45*np.sin(2*np.pi*x/2940)+18*np.sin(4*np.pi*x/2940+.4)
throw1=155*(x>=1260);throw2=-85*(x>=2280)
strata=Z-fold[None,:]-(throw1[None,:]+throw2[None,:])*np.clip((Z-300)/900,0,1)
horizons=330+np.r_[0,np.cumsum(rng.uniform(14,29,49))]
layer=np.searchsorted(horizons,strata);nlay=len(horizons)+1
av=rng.normal(0,48,nlay);ar=rng.normal(0,.025,nlay);aq=rng.normal(0,.12,nlay)
av[0]=ar[0]=aq[0]=0
base_depth=np.maximum(strata-330,0)
vp=1900+.60*base_depth+av[layer]*(1+.15*np.sin(X/450+layer*.3))
rho=1.95+.00024*base_depth+ar[layer]
Q=(45+(0.025)*base_depth)*np.exp(aq[layer]+.10*np.sin(X/700))
for depth,contrast in [(650,100),(1060,140),(1470,130)]:
 vp+=contrast*(strata>depth);rho+=.035*(strata>depth);Q+=6*(strata>depth)
channel=((X-1710)/360)**2+((strata-815)/110)**2<1
vp=np.where(channel,2220+.10*(strata-815),vp);rho=np.where(channel,2.05,rho);Q=np.where(channel,40,Q)
wedge_h=42*np.clip(1-(x-180)/2190,0,1);wedge=(strata>=1180)&(strata<1180+wedge_h[None,:])
vp+=180*wedge;rho+=.075*wedge;Q*=np.where(wedge,.82,1)
gas=((X-2080)/330)**2+((strata-1215)/85)**2<1
vp-=220*gas;rho-=.10*gas;Q*=np.where(gas,.52,1)
hetero=rng.normal(size=vp.shape)
for _ in range(3):hetero=(np.roll(hetero,1,0)+2*hetero+np.roll(hetero,-1,0))/4;hetero=(np.roll(hetero,1,1)+2*hetero+np.roll(hetero,-1,1))/4
mask=strata>330;vp+=12*hetero*mask;rho+=.006*hetero*mask
vp=np.clip(vp,1750,3650);rho=np.clip(rho,1.85,2.65);Q=np.clip(Q,18,160)
vs=vp*(.46+.015*np.sin(strata/370))*np.where(gas,1.04,1)
vs=np.where(strata<=330,874.,vs)
dz=plan['dz_m'];nt=plan['nt'];N=4096;freq=np.arange(N//2+1)/(N*plan['dt_s']);b=-np.pi*freq+2j*freq*np.log(np.maximum(freq,np.finfo(float).tiny)/100);b[0]=0
E=np.exp(-2j*np.pi*freq[:,None]*t[None,:])
def attenuate(s,a):
 Y=np.sum(np.exp(b[:,None]*a[None,:])*E*s[None,:],axis=1);Y[-1]=Y[-1].real
 return np.fft.irfft(Y,n=N)[:nt]
def phase_source(s,phase):
 f=np.fft.rfft(s,n=N);rot=np.exp(-1j*np.deg2rad(phase))*np.ones_like(f);rot[0]=1;rot[-1]=np.cos(np.deg2rad(phase));w=np.fft.irfft(f*rot,n=N)[:nt]
 return w-.16*np.interp(t-.012,t,w,left=0,right=0)+.04*np.interp(t-.025,t,w,left=0,right=0)
def acquisition(s,gain,static_ms):return gain*np.interp(t-static_ms*.001,t,s,left=0,right=0)
primary=np.zeros((nt,len(x)));atten=np.zeros_like(primary);multiple=np.zeros_like(primary);zero=np.zeros_like(primary);cumtime=np.zeros_like(vp);timeQ=np.zeros_like(primary)
for j,xx in enumerate(x):
 travel=np.r_[0,np.cumsum(2*dz/((vp[:-1,j]+vp[1:,j])/2))];cumtime[:,j]=travel
 tau=np.diff(travel);vr=.5*(vp[:-1,j]+vp[1:,j]);rr=.5*(rho[:-1,j]+rho[1:,j]);vsr=.5*(vs[:-1,j]+vs[1:,j])
 A=.5*(np.diff(vp[:,j])/vr+np.diff(rho[:,j])/rr)
 B=.5*np.diff(vp[:,j])/vr-2*(vsr/vr)**2*(np.diff(rho[:,j])/rr+2*np.diff(vs[:,j])/vsr)
 events=travel[1:];vrms=np.sqrt(np.cumsum(vr**2*tau)/events)
 sin2=plan['offset_m']**2/(plan['offset_m']**2+(vrms*events)**2);coeff=A+B*sin2
 fpeak=30*(1+.07*np.sin(2*np.pi*xx/2940+.3));phase=35+7*np.sin(xx/560)
 zz=np.pi*fpeak*(t[:,None]-events[None,:]);W=(1-2*zz**2)*np.exp(-zz**2);clean=W@coeff;clean0=W@A
 for center,time0,amp in [(1260,.79,.018),(2280,1.02,.014),(1710,.72,.010)]:
  tt=np.sqrt(time0*time0+4*(xx-center)**2/2400**2);u=np.pi*fpeak*(t-tt);clean+=amp/(1+abs(xx-center)/300)*(1-2*u*u)*np.exp(-u*u)
 clean=phase_source(clean,phase);clean0=phase_source(clean0,phase)
 aqcum=np.r_[0,np.cumsum(tau/np.sqrt(Q[:-1,j]*Q[1:,j]))];a0=np.interp(t,travel,aqcum)
 vtime=np.interp(t,events,vrms,left=vp[0,j],right=vrms[-1]);factor=np.sqrt(1+(plan['offset_m']/np.maximum(vtime*t,vp[0,j]*.02))**2);a=a0*factor
 primary[:,j]=clean;zero[:,j]=clean0;atten[:,j]=attenuate(clean,a);timeQ[:,j]=np.divide(t,a0,out=np.repeat(Q[0,j],nt),where=a0>0)
 dirty=.14*np.interp(t-(.18+.025*np.sin(xx/600)),t,clean,left=0,right=0)+.08*np.interp(t-.36,t,clean,left=0,right=0)
 multiple[:,j]=attenuate(dirty,a)
 if (j+1)%11==0:print('Generated forward traces',j+1,'/',len(x),flush=True)
# Independent calibration: a well-based coarse 4-layer geometry, same source.
edges=np.array([0,330,780,1260,2205]);cvp=[];cvs=[];crho=[];cq=[]
for lo,hi in zip(edges[:-1],edges[1:]):
 k=(z>=lo)&(z<hi);cvp.append(float(np.mean(vp[k,:])));cvs.append(float(np.mean(vs[k,:])));crho.append(float(np.mean(rho[k,:])));cq.append(float(1/np.mean(1/Q[k,:])))
cvp=np.array(cvp);ch=np.diff(edges)[:3];ctime=np.cumsum(2*ch/cvp[:3]);coff=np.arange(11)*150
public=dict(dt=plan['dt_s'],nt=nt,t=t.tolist(),offset=coff.tolist(),vp=cvp.tolist(),h=ch.tolist(),fref=100,nfft=N,lambda_x=.001,lambda_fixed=.001,lambda_grid=np.logspace(-7,0,41).tolist(),cov_floor=1e-6,nominal_frequency=30)
def exact_attenuation(offset):
 bounds=np.r_[0,np.cumsum(2*ch/cvp[:3]),np.inf];tau=np.maximum(0,np.minimum(t[:,None],bounds[1:])-bounds[:-1]);active=tau>0
 vmax=np.max(np.where(active,cvp,0),axis=1);vmax=np.where(vmax>0,vmax,cvp[0]);lo=np.zeros(nt);hi=(1-1e-12)/vmax
 for _ in range(70):
  pm=(lo+hi)/2;xx=np.sum(tau*cvp**2*pm[:,None]/np.sqrt(np.maximum(1-(pm[:,None]*cvp)**2,np.finfo(float).tiny)),axis=1);take=xx<offset;lo=np.where(take,pm,lo);hi=np.where(take,hi,pm)
 p=(lo+hi)/2
 if offset==0:p[:]=0
 a=np.sum(tau/np.array(cq)/np.sqrt(np.maximum(1-(p[:,None]*cvp)**2,np.finfo(float).tiny)),axis=1);a[t<=0]=offset/(cvp[0]*cq[0]);return a,p
cal=np.zeros((nt,11));rngcal=np.random.default_rng(20262902)
for j,offset in enumerate(coff):
 a,p=exact_attenuation(offset);pp=np.interp(ctime,t,p);am=np.array([.10,-.045,.075])+np.array([-.045,.025,-.030])*(pp*cvp[:3])**2
 u=np.pi*30*(t[:,None]-ctime);s=((1-2*u*u)*np.exp(-u*u))@am;s=phase_source(s,35+7*np.sin(1470/560));s=attenuate(s,a)
 noise=rngcal.normal(size=nt);cal[:,j]=s+noise/np.sqrt(np.mean(noise**2))*np.sqrt(np.mean(s*s))*10**(-18/20)
live=np.ones(len(x),bool);live[[22,77]]=False
datasets=[];truth_data=dict(x=x,z=z,t=t,vp=vp,vs=vs,rho=rho,Q=Q,stratigraphic_depth=strata,channel_mask=channel,wedge_mask=wedge,gas_mask=gas,primary_unattenuated=primary,primary_attenuated=atten,multiples_attenuated=multiple,zero_angle_primary=zero,vertical_time=cumtime,effective_Q=timeQ,live=live,calibration_true_Q=np.array(cq))
for cfg in plan['observations']:
 r=np.random.default_rng(cfg['seed']);gain=np.exp(.09*r.normal(size=len(x)));stat=cfg['static_std_ms']*r.normal(size=len(x));target=np.column_stack([acquisition(primary[:,j],gain[j],stat[j]) for j in range(len(x))]);base=np.column_stack([acquisition(atten[:,j]+cfg['multiple_scale']*multiple[:,j],gain[j],stat[j]) for j in range(len(x))])
 noise=r.normal(size=base.shape)
 for k in range(1,nt):noise[k,:]=.55*noise[k-1,:]+np.sqrt(1-.55**2)*noise[k,:]
 noise=(np.roll(noise,1,axis=1)+2*noise+np.roll(noise,-1,axis=1))/4
 rms=np.sqrt(np.mean(np.column_stack([acquisition(atten[:,j],gain[j],stat[j]) for j in range(len(x))])**2,axis=0));noise=noise/np.sqrt(np.mean(noise*noise,axis=0))*rms[None,:]*10**(-cfg['snr_db']/20)
 coherent=np.zeros_like(base)
 for j,xx in enumerate(x):
  u=t-(.18+abs(xx-420)/1850);ground=np.sin(2*np.pi*(10*u+8*u*u))*np.exp(-(u/.095)**2)
  v=t-(.49+.00010*xx);linear=np.cos(2*np.pi*29*v)*np.exp(-(v/.045)**2);coherent[:,j]=ground+.35*linear
 coherent*=cfg['coherent_energy_fraction']*np.linalg.norm(base)/np.linalg.norm(coherent)
 bursts=np.zeros_like(base)
 for j in [16,51,82]:
  idx=r.integers(170,nt-50,size=4);bursts[idx,j]=r.choice([-1.,1.],4)*6*rms[j]
 observed=base+noise+coherent+bursts;observed[:,~live]=0
 truth_data[cfg['name']+'_target']=target;truth_data[cfg['name']+'_noise']=noise;truth_data[cfg['name']+'_coherent']=coherent;truth_data[cfg['name']+'_gain']=gain;truth_data[cfg['name']+'_statics_ms']=stat
 datasets.append(dict(name=cfg['name'],snr_db=cfg['snr_db'],D=observed.tolist()))
np.savez_compressed(ROOT/'truth/model_truth.npz',**truth_data)
payload=dict(public=public,calibrationD=cal.tolist(),section_offset_m=plan['offset_m'],datasets=datasets)
(ROOT/'core_p1/inputs/complex_input.json').write_text(json.dumps(payload,separators=(',',':')),encoding='utf-8')
(ROOT/'inputs/calibration_public.json').write_text(json.dumps(public,indent=2),encoding='utf-8')
(ROOT/'truth/calibration_truth.json').write_text(json.dumps(dict(vp=cvp.tolist(),vs=cvs,rho=crho,Q=cq,event_t=ctime.tolist(),source_phase_deg=24+7*np.sin(1470/560)),indent=2),encoding='utf-8')
print('Model and two observations ready. No truth array supplied in solver input.',flush=True)
