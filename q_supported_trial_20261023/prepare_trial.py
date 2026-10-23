from pathlib import Path
from datetime import datetime,timezone
import csv,json,hashlib,shutil,struct
import numpy as np
R=Path(__file__).resolve().parent;SRC=R.parent/'q_reference_data_20261023';CORE=R.parent/'denoise_and_field_20261023/core_p1/matlab_src'
for d in ['inputs','source','analysis','figures','raw']: (R/d).mkdir(exist_ok=True)
P=dict(registered_utc=datetime.now(timezone.utc).isoformat(),experiment='External-reference-Q real-waveform adaptation and controlled added attenuation',
    source='https://zenodo.org/records/6932666',Q_selection='Published alpha=0 Z-component results, not re-estimated from trial windows',
    waveform_window_relative_P_s=[-4,12],background_relative_P_s=[-24,-8],dt_s=.05,
    shared_bandpass_hz=[.2,.4,6.5,8.0],evaluated_relative_P_s=[-1,7],fref_hz=1.,
    adaptation='Stationary single-path tstar operator; offset stage is identity; optimized dictionary W=I estimates an unknown band-limited source pulse, not layered reflectivity or Ricker atoms.',
    fixed_lambda=1e-3,optimized_penalty='Unchanged core: .5*sqrt(2*log(nt))*sigma, support refit and group pruning',
    full_path_tstar_factors=[.75,1.,1.25],controlled_extra_tstar_s=[.03,.06,.12],controlled_added_noise_sigma_factors=[0.,1.],
    seed=20261023,development_record='2016-02-10T210245/LWLI',evaluation_records='Remaining 3 records; two earthquakes total, not four independent events',
    rule='No parameter tuning based on controlled-test errors. Published Q is an estimated reference, not physical ground truth.',
    controlled_reference='Band-limited real observed P-wave window, already attenuated by nature; recovery only reverses newly added attenuation.',
    full_path_quality='Report replay residual, amplitude gain, background response and tstar sensitivity, not truth error or measured SNR improvement.')
(R/'PLAN.json').write_text(json.dumps(P,ensure_ascii=False,indent=2),encoding='utf-8')
files=['attenuation_matrix.m','ridge_operator.m','fixed_v1_apply.m','reflectivity_refit_prepared.m','reflectivity_prune_prepared.m','refit_lasso.m','uq_writejson.m']
hashes={}
for f in files:shutil.copy2(CORE/f,R/'source'/f);hashes[f]=hashlib.sha256((CORE/f).read_bytes()).hexdigest()
(R/'FROZEN_CORE_HASHES.json').write_text(json.dumps(hashes,indent=2),encoding='utf-8')
meta=json.loads((SRC/'analysis/VERIFIED_SUMMARY.json').read_text(encoding='utf-8'));items=[];dt=.05;tt=np.arange(321)*dt
def band(a,dt):
    n=2**int(np.ceil(np.log2(2*len(a))));f=np.fft.rfftfreq(n,dt)
    rise=lambda z:.5-.5*np.cos(np.pi*np.clip(z,0,1))
    m=rise((f-.2)/.2)*(1-rise((f-6.5)/1.5))
    return np.fft.irfft(np.fft.rfft(a,n=n)*m,n=n)[:len(a)]
for idx,row in enumerate(meta['sumatra']['sample_waveforms']):
    p=SRC/row['file'];b=p.read_bytes();a=np.frombuffer(b,dtype='<f4',count=struct.unpack_from('<i',b,316)[0],offset=632).astype(float)
    dt0=struct.unpack_from('<f',b,0)[0];begin=struct.unpack_from('<f',b,20)[0];pick=struct.unpack_from('<f',b,32)[0];tax=begin+np.arange(len(a))*dt0
    assert abs(pick/row['author_constant_Q']-row['author_tstar_s'])<1e-5
    filtered=band(a,dt0);signal=np.interp(pick-4+tt,tax,filtered);noise=np.interp(pick-24+tt,tax,filtered)
    sigma=np.median(abs(noise-np.median(noise)))/.6744897501960817
    scale=max(abs(signal).max(),1e-30);signal/=scale;noise/=scale;sigma=max(sigma/scale,1e-6)
    items.append(dict(name=row['event']+'_'+row['station'],station=row['station'],event=row['event'],Q=row['author_constant_Q'],tstar=row['author_tstar_s'],P_pick_s=pick,dt_original_s=dt0,
        D=signal.tolist(),noise=noise.tolist(),sigma=float(sigma),scale=scale,source_filename=p.name,source_sha256=hashlib.sha256(b).hexdigest(),split='development' if idx==0 else 'evaluation'))
    shutil.copy2(p,R/'raw'/p.name)
dev=np.array(items[0]['D']);freq=np.fft.rfftfreq(512,dt);pow=abs(np.fft.rfft(dev,n=512))**2;pow[(freq<.4)|(freq>6.5)]=0
f50=float(freq[np.searchsorted(np.cumsum(pow),.5*sum(pow))]);f50=max(.5,f50)
data=dict(dt=dt,nt=321,nfft=1024,fref=1.,cluster_frequency_hz=f50,items=items,plan=P)
(R/'inputs/trial.json').write_text(json.dumps(data,ensure_ascii=False),encoding='utf-8')
shutil.copy2(SRC/'sumatra/results.csv',R/'raw/sumatra_author_results.csv')
(R/'inputs/geometry_and_Q.json').write_text(json.dumps([{k:v for k,v in a.items() if k not in ['D','noise']} for a in items],indent=2),encoding='utf-8')
print(json.dumps(dict(records=len(items),cluster_frequency=f50,noise_sigma=[a['sigma'] for a in items],source_files_unchanged=True)))
