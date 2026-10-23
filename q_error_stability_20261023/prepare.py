from pathlib import Path
from datetime import datetime,timezone
import json,hashlib,shutil
import numpy as np
R=Path(__file__).resolve().parent;BASE=R.parent
for name in ['inputs','source','analysis','figures']: (R/name).mkdir(exist_ok=True)
src=BASE/'denoise_and_field_20261023/core_p1/matlab_src'
hashes={}
for p in src.glob('*.m'):shutil.copy2(p,R/'source'/p.name);hashes[p.name]=hashlib.sha256(p.read_bytes()).hexdigest()
(R/'FROZEN_CORE_HASHES.json').write_text(json.dumps(hashes,indent=2),encoding='utf-8')
errors=[-.5,-.4,-.3,-.2,-.1,0,.1,.2,.3,.4,.5]
scenarios=[dict(name=f'bias_{int(round(100*e)):+d}',error_percent=100*e,factors=[1+e]*4) for e in errors]
scenarios += [dict(name='shallow_low_deep_high',error_percent=None,factors=[.7,.9,1.1,1.3]),dict(name='shallow_high_deep_low',error_percent=None,factors=[1.3,1.1,.9,.7])]
plan=dict(registered_utc=datetime.now(timezone.utc).isoformat(),task='Synthetic Q-error stability; stop all real-data experiments',
    primary='Known four-layer Q and source; native reflection dictionaries, no field waveform adapter',
    definition='Uniform Q error = 100*(Q_used/Q_true-1) on all four intervals; additional two depth-dependent error patterns',
    primary_true_Q=[40.,60.,85.,110.],primary_frequency_hz=32.,primary_dt_s=.004,primary_nt=451,primary_traces=11,offset_m=450.,
    primary_scenarios=scenarios,noise_conditions=[dict(name='Moderate',snr_db=14,seed=2026102301),dict(name='Strong',snr_db=8,seed=2026102302)],
    complex='Frozen old complex-model observations, 11 predetermined lateral traces. Seven uniform perturbations around the estimated layered-Q profile, not true spatial Q.',
    complex_errors_percent=[-50,-30,-10,0,10,30,50],
    methods=['FixedSplitSharedQ','FullCalRefit','PrunedFullCalRefit'],
    policy='Same imposed Q and source for all methods, original inverse cores unchanged. Freeze geometry, wavelet, regularization and noise estimates across Q. No recalibration or retuning.',
    noise_gain='Primary only: norm(F_Q(D_noisy)-F_Q(D_clean))/norm(known added noise), using identical noisy-data sigma for both solves. Nonlinear noise response, not simple linear gain.',
    weak_scoring='Known weak-event template projection for primary; old fixed low-energy mask for complex. Truth is evaluation/design only.',
    practical_gate='Uniform cases only: total error <=1.10*zero-bias error AND weak error <=1.15*zero-bias weak error AND noise response gain <=2*zero-bias gain. Descriptive chosen tolerance, not a confidence limit.',
    scope='Two synthetic conditions on one controlled geometry plus a fixed complex-model subset. Does not prove universal robustness or independent geology generalization.')
(R/'PLAN.json').write_text(json.dumps(plan,ensure_ascii=False,indent=2),encoding='utf-8')
dt=.004;nt=451;t=np.arange(nt)*dt;nx=11;xp=np.arange(nx)*30.
c=dict(dt=dt,nt=nt,t=t.tolist(),nfft=2048,fref=100.,lambda_x=.001,lambda_fixed=.001,vp=[1800.,2300.,2800.,3300.],h=[230.,300.,350.],Q=plan['primary_true_Q'],fpeak=32.)
event_times=np.array([.24,.36,.49,.65,.84,1.04,1.23,1.42,1.60]);amps=np.array([.10,.015,-.05,.008,.08,-.009,.04,.006,-.025]);weak=[1,3,5,7]
ref=np.zeros((nt,nx));event_idx=[]
for j in range(nx):
    times=event_times+.008*np.sin(j*.45+np.arange(len(amps))*.9)
    ids=np.round(times/dt).astype(int);av=amps*(1+.10*np.sin(j*.6+np.arange(len(amps))))
    ref[ids,j]=av;event_idx.append(ids.tolist())
rng=np.random.default_rng(20261023);noise=[]
for ns in plan['noise_conditions']:
    rng=np.random.default_rng(ns['seed']);a=rng.standard_normal((nt,nx))
    for i in range(1,nt):a[i]+=.45*a[i-1]
    a=.9*a+.1*a.mean(axis=1,keepdims=True);noise.append(dict(**ns,unit_noise=a.tolist()))
primary=dict(config=c,reflectivity=ref.tolist(),event_indices_zero_based=event_idx,event_amplitudes=amps.tolist(),weak_event_indices_zero_based=weak,x=xp.tolist(),noise=noise,scenarios=scenarios)
(R/'inputs/primary.json').write_text(json.dumps(primary),encoding='utf-8')

old=BASE/'complex_fieldlike_20261023';run=Path(json.loads((old/'run.json').read_text(encoding='utf-8'))['run'])
cal=json.loads((BASE/'denoise_and_field_20261023/core_p1/inputs/shared_calibration.json').read_text(encoding='utf-8'))
shutil.copy2(BASE/'denoise_and_field_20261023/core_p1/inputs/shared_calibration.json',R/'inputs/shared_calibration.json')
truth=np.load(old/'truth/model_truth.npz');idx=np.round(np.linspace(0,98,11)).astype(int)
datasets=[]
for name in ['Moderate','Strong']:
    d=json.loads((run/'matlab'/f'{name}.json').read_text(encoding='utf-8'));mask=np.load(old/'analysis'/f'{name}_evaluation_masks.npz')['low_energy']
    datasets.append(dict(name=name,D=np.array(d['observed'])[:,idx].tolist(),target=truth[name+'_target'][:,idx].tolist(),low_energy=mask[:,idx].astype(int).tolist()))
complex_input=dict(public=cal['corrected'],estimate=cal['estimate'],offset_m=450.,x=truth['x'][idx].tolist(),trace_indices_one_based=(idx+1).tolist(),datasets=datasets,
    scenarios=[dict(name=f'calibration_bias_{e:+d}',error_percent=e,factors=[1+e/100]*4) for e in plan['complex_errors_percent']])
(R/'inputs/complex.json').write_text(json.dumps(complex_input),encoding='utf-8')
print(json.dumps(dict(core_files=len(hashes),primary_scenarios=len(scenarios),complex_scenarios=7,complex_traces_one_based=(idx+1).tolist())))
