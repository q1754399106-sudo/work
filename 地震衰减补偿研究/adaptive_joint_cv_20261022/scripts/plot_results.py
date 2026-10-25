from pathlib import Path
import json,hashlib,html,sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1]
def run(pointer):return Path(json.loads((ROOT/pointer).read_text(encoding='utf-8'))['run'])
OUT=ROOT/'figures';OUT.mkdir(exist_ok=True);(OUT/'gathers').mkdir(exist_ok=True);(OUT/'wiggles').mkdir(exist_ok=True)
records=[]
manifest_path=OUT/'figure_manifest.json';cached={}
if manifest_path.exists():cached={r['case']:r for r in json.loads(manifest_path.read_text(encoding='utf-8'))}
sources=[('development','pruning_run.json'),('validation','pruned_validation_run.json')]
if '--development-only' in sys.argv:sources=sources[:1]
else:
    check=json.loads((run('pruned_validation_run.json')/'matlab/checks.json').read_text(encoding='utf-8'));assert check['complete'] and check['cases']==75
for label,pointer in sources:
    for path in sorted((run(pointer)/'matlab/cases').glob('*.json')):
        d=json.loads(path.read_text(encoding='utf-8'));meta=d['meta'];case=meta['id']
        prior=cached.get(case);source_hash=hashlib.sha256(path.read_bytes()).hexdigest()
        if prior and prior['source_sha256']==source_hash and (ROOT/prior['png']).exists() and hashlib.sha256((ROOT/prior['png']).read_bytes()).hexdigest()==prior['image_sha256']:
            records.append(prior)
            if len(records)%20==0:print(f'Reused {len(records)} verified gathers',flush=True)
            continue
        t=np.asarray(d['t']).reshape(-1);truth=np.asarray(d['truth']);observed=np.asarray(d['observed']);fixed=np.asarray(d['fixed']);base=np.asarray(d['baseline']);new=np.asarray(d['optimized'])
        arrays=[truth,observed,fixed,base,new,fixed-truth,base-truth,new-truth]
        amp=max(float(np.max(np.abs(a))) for a in arrays[:5]);err=max(float(np.max(np.abs(a))) for a in arrays[5:])
        titles=['Reference','Observed','Original FixedSplit','Previous method','Optimized','FixedSplit error','Previous error','Optimized error']
        fig=plt.figure(figsize=(15,8),layout='constrained');gs=fig.add_gridspec(2,5,width_ratios=[1,1,1,1,.045])
        axes=[];ampim=errim=None
        for k,(a,title) in enumerate(zip(arrays,titles)):
            ax=fig.add_subplot(gs[k//4,k%4]);limit=amp if k<5 else err
            im=ax.imshow(a,extent=[.5,a.shape[1]+.5,t[-1],t[0]],aspect='auto',interpolation='nearest',cmap='RdBu_r',vmin=-limit,vmax=limit)
            if k<5:ampim=im
            else:errim=im
            ax.set_title(title,fontsize=10);ax.set_xlabel('Trace');
            if k%4==0:ax.set_ylabel('Time (s)')
            axes.append(ax)
        fig.colorbar(ampim,cax=fig.add_subplot(gs[0,4]),label='Amplitude: first five panels')
        fig.colorbar(errim,cax=fig.add_subplot(gs[1,4]),label='Error: last three panels')
        fig.suptitle(case,fontsize=11)
        dest=OUT/'gathers'/(case+'.png');fig.savefig(dest,dpi=150);plt.close(fig)
        predetermined=(meta['scenario'] in ['WeakSection','SmoothAVO'] and meta['snr_db']==12 and meta['seed'] in [20261411,20262291,20262311] and meta.get('suite','normal')=='normal')
        if predetermined:
            fig,axs=plt.subplots(1,4,figsize=(12,6),sharex=True,sharey=True,layout='constrained');gain=.42/max(float(np.max(np.abs(a))) for a in [truth,fixed,base,new])
            for ax,a,title in zip(axs,[truth,fixed,base,new],['Reference','Original FixedSplit','Previous method','Optimized']):
                for j in range(a.shape[1]):
                    center=j+1;line=center+gain*a[:,j];ax.plot(line,t,color='black',lw=.55);ax.fill_betweenx(t,center,line,where=a[:,j]>0,color='black',alpha=.65)
                ax.set_title(title,fontsize=10);ax.set_xlabel('Trace');ax.set_xlim(0,a.shape[1]+1);ax.set_ylim(t[-1],t[0])
            axs[0].set_ylabel('Time (s)');fig.suptitle(case,fontsize=10)
            fig.supxlabel(f'Common gain for every trace and method: {gain:.5g}; no per-trace normalization',fontsize=8)
            fig.savefig(OUT/'wiggles'/(case+'.png'),dpi=170);fig.savefig(OUT/'wiggles'/(case+'.pdf'));plt.close(fig)
        records.append({'case':case,'split':label,'meta':meta,'source':str(path),'source_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'png':dest.relative_to(ROOT).as_posix(),'amplitude_limit':amp,'error_limit':err,'predetermined_wiggle':predetermined,'image_sha256':hashlib.sha256(dest.read_bytes()).hexdigest()})
        if len(records)%20==0:print(f'Plotted {len(records)} gathers',flush=True)
(OUT/'figure_manifest.json').write_text(json.dumps(records,ensure_ascii=False,indent=2),encoding='utf-8')
if '--development-only' in sys.argv:
    print('Development figures complete; full gallery pending fresh validation.',flush=True);raise SystemExit(0)
s=json.loads((ROOT/'analysis/final_summary.json').read_text(encoding='utf-8'))
fig,axs=plt.subplots(1,2,figsize=(11,4.5),layout='constrained');methods=['FixedSplit','DecoupledWeak','PrunedFullCalRefit'];colors=['#777777','#0072B2','#D55E00'];labels=['Original','Previous','Optimized'];groups=['development','validation','colored_stress'];x=np.arange(3)
for ax,key,title in zip(axs,['RelativeError','WeakEventError'],['Waveform error','Weak reflection amplitude error']):
    for j,(m,color,name) in enumerate(zip(methods,colors,labels)):
        values=[s[g]['means'][m][key] for g in groups];ax.bar(x+(j-1)*.24,values,width=.23,color=color,label=name)
    ax.set_xticks(x,['Development\n100 cases','Fresh geology\n60 cases','Fresh colored noise\n15 cases']);ax.set_title(title);ax.set_ylabel('Mean error (lower is better)');ax.set_ylim(bottom=0);ax.grid(axis='y',alpha=.2);ax.spines[['top','right']].set_visible(False)
axs[0].legend(fontsize=9);fig.supxlabel('Waveform uses all cases; weak-event metric uses 20 / 12 / 3 WeakSection cases. Means without uncertainty intervals.',fontsize=8)
fig.savefig(OUT/'overall_comparison.png',dpi=180);fig.savefig(OUT/'overall_comparison.pdf');plt.close(fig)
rows='\n'.join('<option value="'+html.escape(r['png'])+'">'+html.escape(r['case'])+'</option>' for r in records)
page='''<!doctype html><html lang="zh-CN"><meta charset="utf-8"><title>统一优化道集对比</title><style>body{font:16px system-ui;max-width:1500px;margin:24px auto;padding:0 20px;background:#f7f8fa}select{max-width:100%;padding:10px}img{width:100%;background:white;margin-top:20px}p{line-height:1.7}</style><h1>统一优化道集对比</h1><p>110组开发资料（含10组有色噪声）、60组全新模型主验证、15组全新有色噪声压力测试。所有案例均展示，未按结果挑选。前五幅共用振幅色标，后三幅共用误差色标；不做逐道归一化。</p><select id="cases">'''+rows+'''</select><img id="gather" alt="道集与误差对比图"><script>const s=document.getElementById('cases'),im=document.getElementById('gather');s.onchange=()=>im.src=s.value;s.onchange();</script></html>'''
(ROOT/'index.html').write_text(page,encoding='utf-8')
print(f'Completed {len(records)} gathers and {sum(r["predetermined_wiggle"] for r in records)} predetermined wiggle plots',flush=True)
