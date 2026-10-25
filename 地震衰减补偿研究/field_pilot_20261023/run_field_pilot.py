"""Frozen adjacent-shot field pilot. No raw file edits; no assumed physical Q."""
from pathlib import Path
import csv, hashlib, json, struct, sys, time
from datetime import datetime, timezone
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parent
OLD=ROOT.parent/'denoise_and_field_20261023'
for d in ['data','analysis','figures']: (ROOT/d).mkdir(exist_ok=True)
sys.path.insert(0,str(OLD/'scripts' if (OLD/'scripts').exists() else ROOT/'scripts'))
from fk_baseline import soft_fk
plt.rcParams.update({'font.sans-serif':['Microsoft YaHei','DejaVu Sans'],'axes.unicode_minus':False,'font.size':10})
PLAN=dict(created_utc=datetime.now(timezone.utc).isoformat(),development_records=[111,112,113,114],held_out_records=[115,116,117,118],
    split_note='Adjacent records within one survey, not independent surveys. 111 previously seen; 112-118 new this run. Held-out waveforms are not parsed until parameters and event windows frozen.',
    method='Local overlap-add directional f-k, split left/right spread, soft slow-wave time gate and central-offset protection',
    local_window_samples=512,local_window_traces=24,hop_samples=128,hop_traces=6,
    noise_velocity_full_m_s=1000.,noise_velocity_transition_m_s=2200.,
    noise_frequency_full_relative_f90=.30,noise_frequency_zero_relative_f90=.70,
    protected_central_offset_m=250.,maximum_time_s=3.,
    candidate_strengths=[.20,.35,.50],control_maximum_relative_change=.02,
    spectral_ratio_band_hz=[18,55],event_search_s=[.65,3.6],event_separation_s=.35,
    q_window_s=.192,q_minimum_r2=.70,q_negative_slope_fraction_minimum=.80,q_time_shift_relative_spread_maximum=.30,
    q_cross_shot_relative_mad_maximum=.35,
    no_physical_Q_claim='Event matching, source/receiver effects and propagation paths remain unvalidated; apparent Q is diagnostic only.',
    promotion='Do not replace production entry. Do not run field Q compensation unless event, spectral stability and propagation-calibration conditions pass.')
if (ROOT/'PLAN.json').exists():
    previous=json.loads((ROOT/'PLAN.json').read_text(encoding='utf-8'))
    assert all(previous[k]==v for k,v in PLAN.items() if k!='created_utc')
    PLAN=previous
else:
    (ROOT/'PLAN.json').write_text(json.dumps(PLAN,ensure_ascii=False,indent=2),encoding='utf-8')
print('Plan written before development processing',flush=True)
source=OLD/'field_data/L23535.SGY' if (OLD/'field_data/L23535.SGY').exists() else ROOT/'raw/L23535.SGY'
raw=source.read_bytes();sourcehash=hashlib.sha256(raw).hexdigest()
assert sourcehash=='a07d59ee732100b36a985f83ed3222572d5b2e45aee966ee639ac0373cc08487'
dt=.002;nt=3000;step=12240;t=np.arange(nt)*dt;x=(-5225+110*np.arange(96))*.3048;dx=float(x[1]-x[0])
headerrecords={}
for j in range((len(raw)-3600)//step):
    pos=3600+j*step;rec,chan=struct.unpack_from('>ii',raw,pos+8)
    if chan<=96:headerrecords.setdefault(rec,[]).append((chan,pos))
def read_record(rec):
    entries=sorted(headerrecords[rec]);assert [e[0] for e in entries]==list(range(1,97))
    a=[]
    for chan,pos in entries:
        w=np.frombuffer(raw,dtype='>u4',count=nt,offset=pos+240).astype(np.uint32)
        v=np.where(w>>31,-1.,1.)*((w&0xffffff)/16777216.)*np.power(16.,((w>>24)&127).astype(int)-64)
        v[w==0]=0;a.append(v)
    D=np.column_stack(a);assert np.isfinite(D).all()
    return D

def observed_f90(D):
    f=np.fft.rfftfreq(4096,dt);p=np.mean(abs(np.fft.rfft(D,n=4096,axis=0))**2,axis=1)
    p[(f<8)|(f>100)]=0;return float(f[np.searchsorted(np.cumsum(p),.9*sum(p))])

def cos_rise(q):return .5-.5*np.cos(np.pi*np.clip(q,0,1))

def local_model(D,f90):
    """Linear candidate with development-frozen masks; preserve zero traces."""
    tw,xw=PLAN['local_window_samples'],PLAN['local_window_traces'];ht,hx=128,6
    padt=tw;padx=xw
    result=np.zeros_like(D)
    f=np.fft.rfftfreq(2*tw,dt);k=np.fft.fftfreq(2*xw,dx)
    p=np.divide(abs(k)[None,:],f[:,None],out=np.zeros((len(f),len(k))),where=f[:,None]>0)
    cone=cos_rise((p-1/2200)/(1/1000-1/2200))
    low=1-cos_rise((f-.30*f90)/(.40*f90));mask=cone*low[:,None];mask[0]=0
    # Spatial alias region is excluded, rather than folded into a false velocity.
    mask[f>1000/(2*dx),:]=0
    win=np.sqrt(np.hanning(tw)[:,None]*np.hanning(xw)[None,:])
    for inds in [np.where(x<0)[0],np.where(x>0)[0]]:
        a=np.pad(D[:,inds],((padt,padt),(padx,padx)))
        estimate=np.zeros_like(a);weights=np.zeros_like(a)
        for ti in range(0,a.shape[0]-tw+1,ht):
            for xi in range(0,a.shape[1]-xw+1,hx):
                z=a[ti:ti+tw,xi:xi+xw]*win
                S=np.fft.fft(np.fft.rfft(z,n=2*tw,axis=0),n=2*xw,axis=1)
                m=np.fft.irfft(np.fft.ifft(S*mask,axis=1)[:,:xw],n=2*tw,axis=0)[:tw]
                estimate[ti:ti+tw,xi:xi+xw]+=m*win;weights[ti:ti+tw,xi:xi+xw]+=win**2
        estimate=np.divide(estimate,weights,out=np.zeros_like(estimate),where=weights>1e-12)
        result[:,inds]=estimate[padt:padt+nt,padx:padx+len(inds)]
    ab=abs(x)[None,:];tt=t[:,None]
    gate=cos_rise((ab-250)/200)*cos_rise((tt-(ab/2200-.06))/.10)
    gate*=1-cos_rise((tt-(ab/400+.15))/.20)
    gate*=1-cos_rise((tt-2.7)/.3)
    result*=gate;result[:,np.sum(D**2,axis=0)==0]=0
    return result

def probe(v,t0,f0):
    tau=np.sqrt(t0*t0+(x/v)**2);u=np.pi*f0*(t[:,None]-tau[None,:])
    return (1-2*u*u)*np.exp(-u*u)

def slow_power(D,f90):
    S=np.fft.fft(np.fft.rfft(D,n=4096,axis=0),n=256,axis=1)
    f=np.fft.rfftfreq(4096,dt);k=np.fft.fftfreq(256,dx)
    mask=(f[:,None]>=4)&(f[:,None]<=min(.5*f90,1000/(2*dx)))&(abs(k)[None,:]>f[:,None]/1500)
    return float(np.sum(abs(S[mask])**2))

dev={r:read_record(r) for r in PLAN['development_records']}
f90=float(np.median([observed_f90(a) for a in dev.values()]))
models={r:local_model(a,f90) for r,a in dev.items()}
probes=[(v,t0,f,probe(v,t0,f)) for v in [2000.,3000.,4500.] for t0 in [.5,1.,1.8] for f in [20.,35.,55.]]
probe_models=[local_model(a,f90) for v,t0,f,a in probes]
selection=[]
for strength in PLAN['candidate_strengths']:
    changes=[float(strength*np.linalg.norm(m)/np.linalg.norm(q[-1])) for q,m in zip(probes,probe_models)]
    reductions=[1-slow_power(a-strength*models[r],f90)/slow_power(a,f90) for r,a in dev.items()]
    selection.append(dict(strength=strength,probe_maximum_relative_change=max(changes),development_slow_cone_power_reduction_mean=float(np.mean(reductions)),passed=max(changes)<=.02))
passing=[r for r in selection if r['passed']]
strength=max(passing,key=lambda z:z['development_slow_cone_power_reduction_mean'])['strength'] if passing else 0.
print('Development selection',selection,'chosen',strength,flush=True)

def bandpass(D,lo=15,hi=65):
    f=np.fft.rfftfreq(4096,dt);m=cos_rise((f-(lo-5))/5)*(1-cos_rise((f-hi)/8))
    return np.fft.irfft(np.fft.rfft(D,n=4096,axis=0)*m[:,None],n=4096,axis=0)[:nt]

vels=np.arange(1800.,5001.,100.)
t0s=np.arange(.65,3.601,.008)
def semblance(D):
    A=bandpass(D);inds=np.where((abs(x)>=280)&(abs(x)<=1100))[0]
    # Trace scaling for picking only, never used for amplitude spectra or output.
    A=A/np.maximum(np.sqrt(np.mean(A*A,axis=0)),np.finfo(float).tiny)[None,:]
    scores=np.zeros((len(vels),len(t0s)))
    for iv,v in enumerate(vels):
        vals=[]
        for lag in np.arange(-.016,.0161,.004):
            at=np.sqrt(t0s[:,None]**2+(x[inds][None,:]/v)**2)+lag
            z=np.column_stack([np.interp(at[:,j],t,A[:,i]) for j,i in enumerate(inds)])
            vals.append((np.sum(z,axis=1)**2,len(inds)*np.sum(z*z,axis=1)))
        num=sum(a for a,b in vals);den=sum(b for a,b in vals)
        scores[iv]=np.divide(num,den,out=np.zeros_like(num),where=den>0)
    return scores

devsem=np.mean([semblance(a-strength*models[r]) for r,a in dev.items()],axis=0)
best=devsem.max(axis=0);candidates=[]
for i in np.argsort(best)[::-1]:
    ti=float(t0s[i]);vi=int(np.argmax(devsem[:,i]))
    if all(abs(ti-e['t0_s'])>=.35 for e in candidates):candidates.append(dict(t0_s=ti,rms_velocity_m_s=float(vels[vi]),development_mean_semblance=float(best[i])))
    if len(candidates)==5:break
candidates.sort(key=lambda q:q['t0_s'])
# Freeze two strongest separated candidates; automated hypotheses, not interpreted horizons.
events=sorted(sorted(candidates,key=lambda q:q['development_mean_semblance'],reverse=True)[:2],key=lambda q:q['t0_s'])
frozen=dict(observed_f90_hz=f90,subtraction_strength=strength,selection=selection,event_hypotheses=events,
            qualification='Artificial reflection probes test selected velocity/frequency families only. They do not certify true-wave preservation.')
frozen_text=json.dumps(frozen,indent=2)
if (ROOT/'FROZEN_PARAMETERS.json').exists():
    assert (ROOT/'FROZEN_PARAMETERS.json').read_text(encoding='utf-8')==frozen_text,'Frozen parameters may not change after held-out inspection'
else:
    (ROOT/'FROZEN_PARAMETERS.json').write_text(frozen_text,encoding='utf-8')
frozen_hash=hashlib.sha256((ROOT/'FROZEN_PARAMETERS.json').read_bytes()).hexdigest()
print('Frozen before held-out parsing',frozen_hash,'events',events,flush=True)

all_data=dict(dev)
all_data.update({r:read_record(r) for r in PLAN['held_out_records']})
metrics=[];qrows=[];qfits=[];controls=[];nmos={};spectra=[]

def event_spectrum(D,j,event,shift=0.,win_s=.192):
    center=np.sqrt(event['t0_s']**2+(x[j]/event['rms_velocity_m_s'])**2)+shift
    n=int(round(win_s/dt));n+=n%2;at=center+(np.arange(n)-n//2)*dt
    a=np.interp(at,t,D[:,j],left=0,right=0);a-=a.mean()
    # Three deterministic sine tapers, power average; no zero padding treated as independent bins.
    tap=np.sqrt(2/(n+1))*np.sin(np.pi*np.arange(1,n+1)[:,None]*np.arange(1,4)[None,:]/(n+1))
    S=np.fft.rfft(a[:,None]*tap,axis=0)
    return np.fft.rfftfreq(n,dt),np.sqrt(np.mean(abs(S)**2,axis=1)),center

def q_measure(D,rec,stage):
    for shift in [-.016,0.,.016]:
        for j in np.where((abs(x)>=280)&(abs(x)<=1000))[0]:
            f,a,ta=event_spectrum(D,j,events[0],shift)
            _,b,tb=event_spectrum(D,j,events[1],shift)
            sel=(f>=18)&(f<=55);ff=f[sel];yy=np.log(np.maximum(b[sel],1e-30)/np.maximum(a[sel],1e-30))
            slope,inter=np.polyfit(ff,yy,1);pred=inter+slope*ff
            r2=1-float(np.sum((yy-pred)**2)/max(np.sum((yy-yy.mean())**2),1e-30))
            q=-np.pi*(tb-ta)/slope if slope<0 else None
            row=dict(record=rec,stage=stage,shift_s=shift,trace=j+1,offset_m=float(x[j]),slope_per_hz=float(slope),r2=r2,apparent_q=q,frequency_bins=int(sel.sum()))
            qrows.append(row)
            if shift==0 and j==34:spectra.append(dict(record=rec,stage=stage,f=ff.tolist(),log_ratio=yy.tolist(),fit=pred.tolist(),r2=r2))

def diagnostic(D,Y,M,rec):
    late=(t[:,None]>=1.5)&(t[:,None]<=3.5)&(abs(x)[None,:]<800)
    center=abs(x)<250
    return dict(record=rec,split='held_out' if rec in PLAN['held_out_records'] else 'development',
        removed_norm_fraction=float(np.linalg.norm(M)/np.linalg.norm(D)),
        slow_cone_power_reduction=1-slow_power(Y,f90)/slow_power(D,f90),
        central_offset_relative_change=float(np.linalg.norm(M[:,center])/max(np.linalg.norm(D[:,center]),1e-30)),
        late_near_offset_relative_change=float(np.linalg.norm(M[late])/max(np.linalg.norm(D[late]),1e-30)),
        conservation_max_abs=float(np.max(abs(D-Y-M))),
        conservation_relative_max=float(np.max(abs(D-Y-M))/max(np.max(abs(D)),np.finfo(float).tiny)))

for rec,D in all_data.items():
    M=strength*(models[rec] if rec in models else local_model(D,f90));Y=D-M
    B,BM,bd=soft_fk(D,dt,dx)
    met=diagnostic(D,Y,M,rec);met['baseline_removed_norm_fraction']=bd['removed_norm_fraction'];metrics.append(met)
    np.savez_compressed(ROOT/'data'/f'record_{rec}.npz',original=D,local=Y,removed=M,baseline=B,baseline_removed=BM,t=t,x=x,dt=dt)
    q_measure(D,rec,'raw');q_measure(Y,rec,'local')
    nmos[rec]=dict(raw=semblance(D),local=semblance(Y))
    print('Processed',rec,met,flush=True)
    if rec in [111,114,115,118]:
        arrays=[D,B,Y,M];gain=np.minimum(t,3.)**1.5;display=[a*gain[:,None] for a in arrays]
        lim=float(np.quantile(abs(display[0][t<=3.5]),.995));fig,axs=plt.subplots(1,4,figsize=(16,7),sharex=True,sharey=True,layout='constrained')
        for ax,a,title in zip(axs,display,['原始炮集','原全局软 f-k 基线','局部去干扰候选','局部候选减去分量']):
            im=ax.imshow(a,extent=[x[0]/1000,x[-1]/1000,t[-1],0],aspect='auto',cmap='RdBu_r',vmin=-lim,vmax=lim,interpolation='nearest');ax.set_ylim(3.5,0);ax.set_title(title);ax.set_xlabel('偏移距 / km')
        for ev in events:
            tau=np.sqrt(ev['t0_s']**2+(x/ev['rms_velocity_m_s'])**2)
            axs[0].plot(x/1000,tau,color='#222222',lw=.65,ls='--',alpha=.65)
        axs[0].set_ylabel('时间 / s');fig.colorbar(im,ax=axs,label='共用幅度色标（统一显示增益）',shrink=.85)
        fig.suptitle(f'Alaska 记录 {rec} / '+met['split']+' / 局部候选减去范数 '+f"{met['removed_norm_fraction']:.2%}")
        fig.supxlabel('统一 t^1.5 显示增益与原始面板 99.5% 分位色标；不参与处理。虚线为自动事件假设，未确认是同一地质反射。',fontsize=9)
        fig.savefig(ROOT/'figures'/f'record_{rec}.png',dpi=150);fig.savefig(ROOT/'figures'/f'record_{rec}.pdf');plt.close(fig)

# Frozen linear-filter injection check: exact difference between F(D+S) and F(D).
for (v,t0,f,P),PM in zip(probes,probe_models):
    controls.append(dict(velocity_m_s=v,t0_s=t0,frequency_hz=f,relative_change=float(strength*np.linalg.norm(PM)/np.linalg.norm(P))))
assert (ROOT/'FROZEN_PARAMETERS.json').read_bytes() and frozen_hash==hashlib.sha256((ROOT/'FROZEN_PARAMETERS.json').read_bytes()).hexdigest()
assert max(m['conservation_relative_max'] for m in metrics)<1e-12
assert sourcehash==hashlib.sha256(source.read_bytes()).hexdigest()

summaryQ=[]
for rec in all_data:
    for stage in ['raw','local']:
        z=[a for a in qrows if a['record']==rec and a['stage']==stage and a['shift_s']==0]
        neg=[a for a in z if a['slope_per_hz']<0]
        good=[a for a in neg if a['r2']>=.7]
        qs=[a['apparent_q'] for a in good]
        summaryQ.append(dict(record=rec,stage=stage,negative_slope_fraction=len(neg)/len(z),high_r2_negative_fraction=len(good)/len(z),median_r2=float(np.median([a['r2'] for a in z])),median_apparent_q=float(np.median(qs)) if qs else None))
heldQ=[a for a in summaryQ if a['stage']=='local' and a['record'] in PLAN['held_out_records']]
failure=[]
if min(e['development_mean_semblance'] for e in events)<.18:failure.append('At least one automatically picked event has low development semblance (<0.18).')
if any(a['negative_slope_fraction']<.8 for a in heldQ):failure.append('Negative spectral slope fraction below 80% on at least one held-out record.')
if any(a['high_r2_negative_fraction']<.5 for a in heldQ):failure.append('Fewer than half the trace measurements have negative slope and R2 >= 0.70 on at least one held-out record.')
qshifts=[]
for rec in PLAN['held_out_records']:
    med=[]
    for shift in [-.016,0.,.016]:
        vals=[a['apparent_q'] for a in qrows if a['record']==rec and a['stage']=='local' and a['shift_s']==shift and a['apparent_q'] is not None and a['r2']>=.7]
        med.append(float(np.median(vals)) if vals else None)
    spread=(max(med)-min(med))/np.median(med) if all(v is not None for v in med) else None
    qshifts.append(dict(record=rec,shift_medians=med,relative_spread=spread))
if any(a['relative_spread'] is None or a['relative_spread']>.3 for a in qshifts):failure.append('Spectral-window shift sensitivity fails the 30% spread threshold or lacks valid measurements.')
valid_shot_q=[a['median_apparent_q'] for a in heldQ if a['median_apparent_q'] is not None]
cross_shot_relative_mad=float(np.median(abs(np.array(valid_shot_q)-np.median(valid_shot_q)))/np.median(valid_shot_q)) if len(valid_shot_q)==len(heldQ) else None
if cross_shot_relative_mad is None or cross_shot_relative_mad>.35:failure.append('Cross-shot apparent Q relative MAD exceeds 35% or lacks valid measurements.')
failure.append('No independently interpreted reflection pair, source/receiver calibration or confirmed propagation model. Physical Q remains non-identifiable in this pilot.')

for file,rows in [('field_metrics.csv',metrics),('reflection_probe_controls.csv',controls),('spectral_ratio_measurements.csv',qrows),('q_summary.csv',summaryQ)]:
    with (ROOT/'analysis'/file).open('w',encoding='utf-8-sig',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

result=dict(source_sha256=sourcehash,frozen_parameters_sha256=frozen_hash,metrics=metrics,probe_controls=controls,event_hypotheses=events,
    apparent_Q_summary=summaryQ,window_sensitivity=qshifts,cross_shot_relative_mad=cross_shot_relative_mad,q_gate=dict(passed=False,reasons=failure),
    field_Q_compensation_run=False,reason='Q feasibility gate failed; no default or arbitrary Q inserted.',
    note='Slow-cone reduction is a processing metric, not measured SNR improvement; probes are limited synthetic controls; holdouts belong to the same line.')
(ROOT/'analysis/results.json').write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf-8')
(ROOT/'analysis/spectral_examples.json').write_text(json.dumps(spectra,indent=2),encoding='utf-8')
np.savez_compressed(ROOT/'analysis/semblance_scans.npz',velocities=vels,t0s=t0s,development_mean=devsem,**{f'{r}_{k}':a for r,z in nmos.items() for k,a in z.items()})

fig,axs=plt.subplots(1,3,figsize=(15,4.8),layout='constrained')
recs=[a['record'] for a in metrics];axs[0].bar(recs,[100*a['removed_norm_fraction'] for a in metrics],color=['#597b8c' if a['split']=='development' else '#d57732' for a in metrics]);axs[0].set(xlabel='记录号',ylabel='减去输入范数 / %',title='开发炮：蓝；留出炮：橙')
axs[1].bar(recs,[100*a['slow_cone_power_reduction'] for a in metrics],color='#087f8c');axs[1].set(xlabel='记录号',ylabel='慢速低频谱锥功率变化 / %',title='正值表示该谱锥功率降低，并非 SNR')
axs[2].bar(range(len(controls)),[100*a['relative_change'] for a in controls],color='#597b8c');axs[2].axhline(2,color='#b33e35',ls='--',label='2% 控制门槛');axs[2].legend();axs[2].set(xlabel='人工反射探针序号（27 个）',ylabel='相对波形变化 / %',title='仅覆盖指定速度 / 频率族')
fig.savefig(ROOT/'figures/metrics.png',dpi=160);fig.savefig(ROOT/'figures/metrics.pdf');plt.close(fig)

fig,axs=plt.subplots(1,3,figsize=(15,5),layout='constrained')
im=axs[0].imshow(devsem,extent=[t0s[0],t0s[-1],vels[-1],vels[0]],aspect='auto',vmin=0,vmax=max(.3,devsem.max()),cmap='viridis');axs[0].set(xlabel='零偏移时间假设 / s',ylabel='均方根速度假设 / (m/s)',title='开发炮平均速度相似性');fig.colorbar(im,ax=axs[0],label='semblance')
for e in events:axs[0].plot(e['t0_s'],e['rms_velocity_m_s'],'rx')
for z in summaryQ:
    if z['stage']=='local':axs[1].plot(z['record'],z['negative_slope_fraction'],'o',color='#087f8c');axs[1].plot(z['record'],z['high_r2_negative_fraction'],'s',color='#d57732')
axs[1].plot([],[],'o',color='#087f8c',label='负斜率比例');axs[1].plot([],[],'s',color='#d57732',label='负斜率且 R²≥0.70');axs[1].axhline(.8,color='gray',ls='--');axs[1].set(xlabel='记录号',ylabel='测量比例',ylim=(0,1),title='光谱比一致性');axs[1].legend(fontsize=9)
example=next(z for z in spectra if z['record']==115 and z['stage']=='local');axs[2].plot(example['f'],example['log_ratio'],'o-',label='测量');axs[2].plot(example['f'],example['fit'],'--',label='直线拟合');axs[2].set(xlabel='频率 / Hz',ylabel='ln(晚事件 / 早事件振幅谱)',title=f"留出 115 炮 / 35 道 / R²={example['r2']:.2f}");axs[2].legend()
fig.suptitle('表观 Q 可估计性检查：自动事件假设不等于确认反射层位')
fig.savefig(ROOT/'figures/q_feasibility.png',dpi=160);fig.savefig(ROOT/'figures/q_feasibility.pdf');plt.close(fig)
print(json.dumps(dict(output=str(ROOT),chosen_strength=strength,f90=f90,q_gate=result['q_gate']),ensure_ascii=True),flush=True)
