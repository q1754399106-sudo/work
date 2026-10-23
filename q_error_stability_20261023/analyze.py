from pathlib import Path
import csv,hashlib,json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
R=Path(__file__).resolve().parent
plt.rcParams.update({'font.sans-serif':['Microsoft YaHei','DejaVu Sans'],'axes.unicode_minus':False,'font.size':10})
methods=['fixed','unpruned','optimized'];labels=['最初分步法（共享 Q / 子波）','未剪枝重拟合','优化剪枝法'];colors=['#d47832','#547daf','#087f8c']
metrics=[];runs={};references={};checks={};gates=[]
pi=json.loads((R/'inputs/primary.json').read_text(encoding='utf-8'));ci=json.loads((R/'inputs/complex.json').read_text(encoding='utf-8'))
def arr(v):return np.asarray(v,dtype=float)
def weak_coeff(y):
    out=[];true=[];ref=arr(pi['reflectivity']);dt=pi['config']['dt'];f=pi['config']['fpeak']
    for j,events in enumerate(pi['event_indices_zero_based']):
        for e in pi['weak_event_indices_zero_based']:
            k=events[e];idx=np.arange(max(0,k-10),min(len(y),k+11));z=np.pi*f*((idx-k)*dt);w=(1-2*z*z)*np.exp(-z*z)
            out.append(np.dot(y[idx,j],w)/np.dot(w,w));true.append(ref[k,j])
    return np.array(out),np.array(true)
def base_metric(out,target,t,low=None):
    err=out-target;deep=t>=1.05
    m=dict(relative_error=float(np.linalg.norm(err)/np.linalg.norm(target)),deep_error=float(np.linalg.norm(err[deep])/np.linalg.norm(target[deep])),correlation=float(np.corrcoef(out.ravel(),target.ravel())[0,1]),output_amplitude_ratio=float(np.linalg.norm(out)/np.linalg.norm(target)))
    if low is not None:m['weak_error']=float(np.linalg.norm(err[low])/np.linalg.norm(target[low]))
    return m

for protocol in ['primary','complex']:
    folder=R/('matlab_'+protocol);check=json.loads((folder/'checks.json').read_text(encoding='utf-8'));checks[protocol]=check
    assert check['complete'] and check['optimization_flags']==0
    assert int((R/f'{protocol}_exit_code.txt').read_text(encoding='utf-8-sig'))==0
    ref=json.loads((folder/'reference.json').read_text(encoding='utf-8'));references[protocol]=ref;t=arr(ref['t']).ravel();runs[protocol]={}
    for p in folder.glob('*.json'):
        if p.name in ['checks.json','progress.json','reference.json']:continue
        d=json.loads(p.read_text(encoding='utf-8'));scenario=d['scenario'];error=scenario['error_percent'];uniform=isinstance(error,(int,float))
        runs[protocol][scenario['name']]=d
        for record,reference in zip(d['records'],ref['datasets']):
            assert record['name']==reference['name'];target=arr(reference['target']);low=arr(reference['low_energy']).astype(bool) if protocol=='complex' else None
            for method in methods:
                out=arr(record[method]);assert np.isfinite(out).all()
                row=dict(protocol=protocol,condition=record['name'],scenario=scenario['name'],uniform=uniform,Q_error_percent=float(error) if uniform else None,method=method,**base_metric(out,target,t,low))
                if protocol=='primary':
                    est,true=weak_coeff(out);cleanout=arr(record[method+'_clean']);n=arr(reference['noise'])
                    row.update(weak_error=float(np.linalg.norm(est-true)/np.linalg.norm(true)),weak_amplitude_ratio_median=float(np.median(est/true)),noise_response_gain=float(np.linalg.norm(out-cleanout)/np.linalg.norm(n)),clean_distortion=float(np.linalg.norm(cleanout-target)/np.linalg.norm(target)))
                else:
                    early=(t>=.02)&(t<=.18);row.update(weak_amplitude_ratio_median=None,noise_response_gain=None,clean_distortion=None,early_false_signal_ratio=float(np.linalg.norm(out[early])/np.linalg.norm(target)))
                metrics.append(row)

# Derived changes always refer to each method's own zero-perturbation output.
for row in metrics:
    baseline=next(a for a in metrics if a['protocol']==row['protocol'] and a['condition']==row['condition'] and a['method']==row['method'] and a['uniform'] and a['Q_error_percent']==0)
    row['error_ratio_to_zero']=row['relative_error']/baseline['relative_error'];row['weak_error_ratio_to_zero']=row['weak_error']/baseline['weak_error']
    d=runs[row['protocol']][row['scenario']];b=next(z for z in runs[row['protocol']].values() if isinstance(z['scenario']['error_percent'],(int,float)) and z['scenario']['error_percent']==0)
    y=arr(next(z for z in d['records'] if z['name']==row['condition'])[row['method']]);yb=arr(next(z for z in b['records'] if z['name']==row['condition'])[row['method']])
    row['output_change_from_zero']=float(np.linalg.norm(y-yb)/max(np.linalg.norm(yb),1e-30))
    if row['protocol']=='primary':
        row['noise_gain_ratio_to_zero']=row['noise_response_gain']/max(baseline['noise_response_gain'],1e-30)
        row['practical_gate_passed']=bool(row['error_ratio_to_zero']<=1.10 and row['weak_error_ratio_to_zero']<=1.15 and row['noise_gain_ratio_to_zero']<=2)
    else:row['noise_gain_ratio_to_zero']=None;row['practical_gate_passed']=None

for cond in ['Moderate','Strong']:
    for method in methods:
        z=sorted([a for a in metrics if a['protocol']=='primary' and a['condition']==cond and a['method']==method and a['uniform']],key=lambda a:a['Q_error_percent'])
        mid=next(i for i,a in enumerate(z) if a['Q_error_percent']==0);lo=hi=mid
        while lo>0 and z[lo-1]['practical_gate_passed']:lo-=1
        while hi+1<len(z) and z[hi+1]['practical_gate_passed']:hi+=1
        gates.append(dict(condition=cond,method=method,contiguous_tested_interval_percent=[z[lo]['Q_error_percent'],z[hi]['Q_error_percent']],left_reaches_test_limit=lo==0,right_reaches_test_limit=hi==len(z)-1,definition='Overall error <=110%, weak amplitude error <=115%, nonlinear noise gain <=200% of own zero-bias baseline. Tested grid only, not a confidence interval.'))
keys=[]
for row in metrics:
    for k in row:
        if k not in keys:keys.append(k)
with (R/'analysis/metrics.csv').open('w',encoding='utf-8-sig',newline='') as f:w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows(metrics)

for protocol in ['primary','complex']:
    fig,axs=plt.subplots(2,3,figsize=(16,8),layout='constrained')
    metrics_used=[('relative_error','整体相对波形误差'),('weak_error','弱事件振幅误差' if protocol=='primary' else '固定低能量区波形误差'),('noise_response_gain','非线性噪声响应增益') if protocol=='primary' else ('early_false_signal_ratio','早时伪信号 / 参考范数')]
    for i,condition in enumerate(['Moderate','Strong']):
        for j,(key,ylabel) in enumerate(metrics_used):
            ax=axs[i,j]
            for method,label,color in zip(methods,labels,colors):
                z=sorted([a for a in metrics if a['protocol']==protocol and a['condition']==condition and a['method']==method and a['uniform']],key=lambda a:a['Q_error_percent'])
                ax.plot([a['Q_error_percent'] for a in z],[a[key] for a in z],'o-',label=label,color=color,ms=4)
            ax.axvline(0,color='gray',ls=':',lw=.8);ax.set(xlabel='Q 相对真值误差 / %' if protocol=='primary' else '相对已标定 Q 的扰动 / %',ylabel=ylabel,title=condition,ylim=(0,None));ax.grid(alpha=.2)
            if i==0 and j==0:ax.legend(fontsize=8)
    fig.suptitle('已知真 Q 的四层反射模型：只改变补偿 Q' if protocol=='primary' else '复杂模型固定 11 道：扰动既有标定 Q，不代表空间真 Q 误差')
    fig.supxlabel('同一观测、子波、几何、正则化和噪声估计；Q 低估意味着更强的吸收补偿。各面板纵轴独立，均从 0 开始。',fontsize=9)
    fig.savefig(R/'figures'/f'{protocol}_stability.png',dpi=155);fig.savefig(R/'figures'/f'{protocol}_stability.pdf');plt.close(fig)

fig,axs=plt.subplots(1,2,figsize=(13,5),layout='constrained')
for ax,condition in zip(axs,['Moderate','Strong']):
    for method,label,color in zip(methods,labels,colors):
        z=sorted([a for a in metrics if a['protocol']=='primary' and a['condition']==condition and a['method']==method and a['uniform']],key=lambda a:a['Q_error_percent'])
        ax.plot([a['Q_error_percent'] for a in z],[a['weak_amplitude_ratio_median'] for a in z],'o-',label=label,color=color)
    ax.axhline(1,color='gray',ls='--',label='参考振幅比 1');ax.set(xlabel='Q 相对真值误差 / %',ylabel='弱反射投影振幅比中位数',title=condition);ax.legend(fontsize=8);ax.grid(alpha=.2)
fig.suptitle('弱反射恢复：保留符号的事件模板投影，不逐道归一化')
fig.savefig(R/'figures/weak_amplitudes.png',dpi=160);fig.savefig(R/'figures/weak_amplitudes.pdf');plt.close(fig)

fig,axs=plt.subplots(1,2,figsize=(13,5),layout='constrained')
for ax,condition in zip(axs,['Moderate','Strong']):
    for method,label,color in zip(methods,labels,colors):
        z=sorted([a for a in metrics if a['protocol']=='primary' and a['condition']==condition and a['method']==method and a['uniform']],key=lambda a:a['Q_error_percent'])
        ax.plot([a['Q_error_percent'] for a in z],[a['error_ratio_to_zero'] for a in z],'o-',label=label,color=color)
    ax.axhline(1.1,color='gray',ls='--',label='整体误差 +10%');ax.set(xlabel='Q 相对真值误差 / %',ylabel='整体误差 / 本方法零误差基线',title=condition,ylim=(0,None));ax.legend(fontsize=8);ax.grid(alpha=.2)
fig.suptitle('相对于各自基线的敏感性：曲线平不等于绝对误差低')
fig.savefig(R/'figures/relative_sensitivity.png',dpi=160);fig.savefig(R/'figures/relative_sensitivity.pdf');plt.close(fig)

for protocol in ['primary','complex']:
    ref=references[protocol];t=arr(ref['t']).ravel();rr=next(a for a in ref['datasets'] if a['name']=='Strong');x=arr(pi['x'] if protocol=='primary' else ci['x'])/1000
    views=[]
    for error in [-50,0,50]:
        d=next(a for a in runs[protocol].values() if isinstance(a['scenario']['error_percent'],(int,float)) and a['scenario']['error_percent']==error)
        r=next(a for a in d['records'] if a['name']=='Strong');views.append((error,r))
    arrays=[arr(rr['target']),arr(rr['D'])]+[arr(a[m]) for e,a in views for m in ['fixed','optimized']];limit=max(float(abs(a).max()) for a in arrays)
    fig,axs=plt.subplots(2,4,figsize=(17,9),layout='constrained',sharex=True,sharey=True)
    panels=[(axs[0,0],arr(rr['target']),'纯净参考'),(axs[1,0],arr(rr['D']),'同一观测')]
    for col,(error,r) in enumerate(views,1):panels.extend([(axs[0,col],arr(r['fixed']),f'分步法 / Q {error:+d}%'),(axs[1,col],arr(r['optimized']),f'优化剪枝 / Q {error:+d}%')])
    for ax,y,title in panels:
        im=ax.pcolormesh(x,t,y,shading='nearest',cmap='RdBu_r',vmin=-limit,vmax=limit,rasterized=True);ax.set_title(title);ax.set_xlabel('横向位置 / km');ax.invert_yaxis()
    # Shared axes already propagate inversion; set once to a deterministic order.
    axs[0,0].set_ylim(t[-1],t[0]);axs[0,0].set_ylabel('时间 / s');axs[1,0].set_ylabel('时间 / s');fig.colorbar(im,ax=axs,label='八面板共同原幅度色标',shrink=.85)
    fig.suptitle(('真 Q 误差' if protocol=='primary' else '相对标定 Q 扰动')+' / Strong：低估、基线、高估的波形对比')
    fig.supxlabel('共同绝对最大幅度色标，不截幅、不逐道归一化；复杂模型只展示预先固定的 11 道，未重新生成观测。',fontsize=9)
    fig.savefig(R/'figures'/f'{protocol}_gathers.png',dpi=150);fig.savefig(R/'figures'/f'{protocol}_gathers.pdf');plt.close(fig)

# Depth-dependent error patterns get a separate table, not a fake uniform-error x coordinate.
mixed=[a for a in metrics if a['protocol']=='primary' and not a['uniform']]
with (R/'analysis/depth_pattern_metrics.csv').open('w',encoding='utf-8-sig',newline='') as f:w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows(mixed)
summary=dict(checks=checks,metrics=metrics,practical_intervals=gates,source_core_files_unchanged=65,
    distinctions='Primary true-Q errors versus complex calibrated-Q perturbations are separate. Practical intervals use chosen thresholds on a finite grid, not confidence limits.',
    no_real_data=True)
assert all(hashlib.sha256((R/'source'/p).read_bytes()).hexdigest()==h for p,h in json.loads((R/'FROZEN_CORE_HASHES.json').read_text(encoding='utf-8')).items())
(R/'analysis/summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps(dict(checks=checks,practical_intervals=gates,rows=len(metrics)),ensure_ascii=True,indent=2))
