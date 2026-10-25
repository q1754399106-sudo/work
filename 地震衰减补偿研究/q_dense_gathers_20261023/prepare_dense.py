from pathlib import Path
from datetime import datetime,timezone
import json,hashlib,shutil
import numpy as np
R=Path(__file__).resolve().parent;B=R.parent;OLD=B/'q_error_stability_20261023'
for d in ['inputs','source','results','figures','analysis']:(R/d).mkdir(exist_ok=True)
hashes=json.loads((OLD/'FROZEN_CORE_HASHES.json').read_text(encoding='utf-8'))
for p,h in hashes.items():
    assert hashlib.sha256((OLD/'source'/p).read_bytes()).hexdigest()==h
    shutil.copy2(OLD/'source'/p,R/'source'/p)
(R/'FROZEN_CORE_HASHES.json').write_text(json.dumps(hashes,indent=2),encoding='utf-8')
errors=[-30,-10,0,10,30]
plan=dict(created_utc=datetime.now(timezone.utc).isoformat(),purpose='Replace 11-trace figures with computed 99-trace gathers and legible exports',
    errors_percent=errors,methods=['FixedSplitSharedQ','FullCalRefit','PrunedFullCalRefit'],
    primary='New 99-trace synthesis on the same 0-300m aperture and continuous event-time/amplitude laws. New fixed noise draws; no interpolation of processed traces.',
    complex='All original 99 traces, original observations and targets unchanged, including dead traces 23 and 78. No spatial interpolation.',
    solver='Original 65 core files, contiguous 11-trace blocks. Same Q/wavelet/noise estimate across methods and frozen noise estimates across Q.',
    execution='MATLAB 4 computation threads; all five Q cases recomputed, including zero. No solver parameters tuned for denser figures.',
    display='Raster PNG >=400 dpi plus vector PDF; large comparisons split into separate bias pages. Shared color limits per comparison, with declared common clipping for detail views; no tracewise gains.',
    comparability='Primary is a newly sampled synthetic observation, so old 11-trace numbers are not directly interchangeable. Complex full-gather statistics differ from old subsample statistics.')
(R/'PLAN.json').write_text(json.dumps(plan,ensure_ascii=False,indent=2),encoding='utf-8')
old=json.loads((OLD/'inputs/primary.json').read_text(encoding='utf-8'));nx=99;nt=old['config']['nt'];dt=old['config']['dt'];u=np.linspace(0,10,nx)
times=np.array([.24,.36,.49,.65,.84,1.04,1.23,1.42,1.60]);amps=np.array(old['event_amplitudes']);reflectivity=np.zeros((nt,nx));ids=[]
for j,s in enumerate(u):
    ti=times+.008*np.sin(s*.45+np.arange(len(amps))*.9);ii=np.round(ti/dt).astype(int);reflectivity[ii,j]=amps*(1+.10*np.sin(s*.6+np.arange(len(amps))));ids.append(ii.tolist())
noise=[]
for k,cond in enumerate(['Moderate','Strong']):
    seed=2026102311+k;rng=np.random.default_rng(seed);a=rng.standard_normal((nt,nx))
    for i in range(1,nt):a[i]+=.45*a[i-1]
    a=.9*a+.1*a.mean(axis=1,keepdims=True)
    noise.append(dict(name=cond,snr_db=14 if k==0 else 8,seed=seed,unit_noise=a.tolist()))
primary=dict(config=old['config'],reflectivity=reflectivity.tolist(),event_indices_zero_based=ids,x=np.linspace(0,300,nx).tolist(),noise=noise,errors_percent=errors)
(R/'inputs/primary.json').write_text(json.dumps(primary),encoding='utf-8')
oldci=json.loads((OLD/'inputs/complex.json').read_text(encoding='utf-8'));src=B/'complex_fieldlike_20261023';run=Path(json.loads((src/'run.json').read_text(encoding='utf-8'))['run']);truth=np.load(src/'truth/model_truth.npz')
ds=[];origin=[]
for cond in ['Moderate','Strong']:
    p=run/'matlab'/f'{cond}.json';d=json.loads(p.read_text(encoding='utf-8'));D=np.array(d['observed']);assert D.shape==(851,99)
    # Verify the earlier 11-trace subset really came from these exact observations.
    idx=np.array(oldci['trace_indices_one_based'])-1;subset=np.array(next(a for a in oldci['datasets'] if a['name']==cond)['D']);assert np.array_equal(D[:,idx],subset)
    ds.append(dict(name=cond,D=D.tolist(),target=truth[cond+'_target'].tolist()))
    origin.append(dict(file=str(p),sha256=hashlib.sha256(p.read_bytes()).hexdigest(),observation_sha256=hashlib.sha256(D.tobytes()).hexdigest()))
complex_input=dict(public=oldci['public'],estimate=oldci['estimate'],datasets=ds,x=truth['x'].tolist(),errors_percent=errors,live=truth['live'].astype(int).tolist())
(R/'inputs/complex.json').write_text(json.dumps(complex_input),encoding='utf-8');(R/'INPUT_PROVENANCE.json').write_text(json.dumps(origin,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps(dict(primary_shape=[nt,nx],complex_shape=[851,99],Q_errors=errors,cores=len(hashes))))
