"""Reanalyse existing paired MATLAB results; never tune or modify source data."""
from pathlib import Path
import argparse, json, hashlib, shutil, sys, platform, itertools
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.lines import Line2D
from PIL import Image, ImageOps, ImageDraw

ROOT = Path(__file__).resolve().parent
MODELS = {'SoftLowQ':'强吸收软地层 / Soft low Q', 'VelocityReversal':'速度反转 / Velocity reversal', 'ThinHighFrequency':'薄层高频 / Thin high frequency'}
SCENES = {'SmoothAVO':'连续AVO / Smooth AVO','FaultSection':'断层 / Fault','PinchoutSection':'尖灭 / Pinchout','WeakSection':'弱反射 / Weak','CoherentNoiseAVO':'相干干扰 / Coherent noise'}
LABELS = {'Direct':'一步法 / Direct','Split':'常规分步 / Split','Frozen':'Frozen分步 / Frozen split','Target08':'Target08分步 / Target08 split'}
SHORT = ['一步\nDirect','常规分步\nSplit','Frozen','Target08']
MAP = {'Direct':('Baseline','Direct'), 'Split':('Baseline','Split'), 'Frozen':('Frozen','Split'), 'Target08':('Tuned','Split')}
COLORS = {'Direct':'#333333','Split':'#E69F00','Frozen':'#0072B2','Target08':'#009E73'}
MARKERS = {'Direct':'o','Split':'^','Frozen':'s','Target08':'D'}
STYLES = {'Direct':':','Split':'--','Frozen':'-','Target08':'-.'}
METRICS = ['RelativeError','EventAmplitudeError','AVOGradientError','WeakEventError','FaultNeighborhoodError','PinchoutTailError']
KEY = ['Model','Suite','Scenario','Seed']

def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p): return json.loads(p.read_text(encoding='utf-8-sig'))
def dump(p,x): p.write_text(json.dumps(x,ensure_ascii=False,indent=2,allow_nan=False),encoding='utf-8')
def table(frame):
    def fmt(x):
        if isinstance(x,(float,np.floating)): return '—' if not np.isfinite(x) else f'{x:.6f}'
        return str(x)
    return '\n'.join(['| '+' | '.join(map(str,frame.columns))+' |','| '+' | '.join(['---']*len(frame.columns))+' |']+['| '+' | '.join(fmt(v) for v in row)+' |' for row in frame.itertuples(index=False,name=None)])

def generate(source):
    assert source.resolve()!=ROOT.resolve()
    assert not (ROOT/'REPORT_ZH_EN.md').exists(), 'Use a fresh output directory; preserve delivered report.'
    inputs=ROOT/'inputs'; inputs.mkdir(exist_ok=True); data=ROOT/'data'; data.mkdir(exist_ok=True); figs=ROOT/'figures'; figs.mkdir(exist_ok=True)
    source_files=['results/matlab_test/metrics.csv','results/matlab_test/spec_used.json','results/matlab_test/connection.json','results/matlab_test/run_status.json','results/matlab_test/selected_before_test.json','results/matlab_test/selection_evidence.json','results/matlab_search2/selected_presets.json','inputs/stability_spec.json','matlab_src/run_unknown_q_test.m','matlab_src/uq_build_prior.m','matlab_src/uq_select.m','matlab_src/context_v2.m','matlab_src/build_structure_models.m','matlab_src/score_structure.m','matlab_src/score_stability.m','matlab_src/solve_v2.m']
    source_files += [f'results/matlab_calibrate/estimate_{m}.json' for m in MODELS]
    source_files += [f'results/matlab_test/example_{m}_{scene}.json' for m in MODELS for scene in ['SmoothAVO','CoherentNoiseAVO']]
    provenance=[]
    for rel in source_files:
        p=source/rel; dest=inputs/rel; dest.parent.mkdir(exist_ok=True,parents=True); shutil.copy2(p,dest)
        provenance.append({'source':str(p),'snapshot':str(dest.relative_to(ROOT)),'sha256':sha(p),'bytes':p.stat().st_size})
    def load(rel): return read(inputs/rel)
    raw=pd.read_csv(inputs/'results/matlab_test/metrics.csv'); spec=load('results/matlab_test/spec_used.json'); status=load('results/matlab_test/run_status.json'); connection=load('results/matlab_test/connection.json')
    assert len(raw)==756 and status['status']=='COMPLETED' and status['inputs']==126
    assert connection['execution']=='PASS' and all(connection['preserved'].values())
    assert not raw.duplicated(KEY+['Method','Route']).any()
    selected=load('results/matlab_search2/selected_presets.json')
    assert selected==load('results/matlab_test/selected_before_test.json')
    assert all(x['Preset']['name']=='Target08' for x in selected)
    pieces=[]
    for method,(old,route) in MAP.items():
        a=raw[(raw.Method==old)&(raw.Route==route)].copy(); assert len(a)==126
        a['ComparisonMethod']=method; pieces.append(a)
    d=pd.concat(pieces,ignore_index=True)
    assert len(d)==504 and not d.duplicated(KEY+['ComparisonMethod']).any()
    assert d.groupby(KEY).ComparisonMethod.nunique().eq(4).all()
    for field in ['EffectiveQ','EstimatedFrequency','EstimatedSigmaRatio']:
        assert d.groupby(KEY)[field].nunique().eq(1).all(),field
    assert (d.FinalFlag==0).all() and d.NormalResidual.max()<2e-7
    assert set(d.Seed)==set(spec['test_seeds'])
    assert not set(d.Seed)&(set(spec['tuning_seeds'])|{spec['calibration_seed']})
    assert set(d[d.ComparisonMethod=='Frozen'].TargetMultiplier)=={1.0}
    assert set(d[d.ComparisonMethod=='Target08'].TargetMultiplier)=={0.8}
    assert set(d[d.ComparisonMethod.isin(['Frozen','Target08'])].Gamma)=={10}
    assert set(d.PilotMultiplier)=={4}
    # Identical selected hyperparameters within each direct/split pair.
    for method,g in raw.groupby('Method'):
        for field in ['Scale','TargetMismatch','Gamma','PilotMultiplier','TargetMultiplier','EstimatedSigmaRatio']:
            assert g.groupby(KEY)[field].nunique().eq(1).all(),(method,field)
    d.to_csv(data/'four_methods_all_504_rows.csv',index=False)
    raw.to_csv(data/'six_routes_all_756_rows.csv',index=False)
    seedmeans=d.groupby(['Suite','Model','ComparisonMethod','Seed'],sort=False)[METRICS].mean().reset_index()
    seedmeans.to_csv(data/'seed_means.csv',index=False)
    summary=seedmeans.groupby(['Suite','Model','ComparisonMethod']).RelativeError.agg(['mean','std']).reset_index()
    tail=d.groupby(['Suite','Model','ComparisonMethod']).RelativeError.agg(worst='max',above_one=lambda a:int((a>1).sum()),above_two=lambda a:int((a>2).sum()),n='count').reset_index()
    summary=summary.merge(tail,on=['Suite','Model','ComparisonMethod']); summary.to_csv(data/'overall_summary.csv',index=False)
    cells=d.groupby(['Suite','Model','Scenario','ComparisonMethod'])[METRICS].agg(['mean','std','max']); cells.to_csv(data/'scene_summary.csv')
    paired=[]; pair_summary=[]
    pivot=d.pivot(index=KEY,columns='ComparisonMethod',values='RelativeError')
    for ref,method in itertools.combinations(LABELS,2):
        z=pivot[[ref,method]].copy(); z['Delta']=z[method]-z[ref]; z['Reference']=ref; z['Compared']=method
        paired.append(z.reset_index().rename(columns={ref:'ReferenceError',method:'ComparedError'}))
        for (model,suite),g in z.reset_index().groupby(['Model','Suite']):
            pair_summary.append(dict(Model=model,Suite=suite,Reference=ref,Compared=method,ReferenceMean=g[ref].mean(),ComparedMean=g[method].mean(),MeanDelta=g.Delta.mean(),Reduction=1-g[method].mean()/g[ref].mean(),LowerFraction=(g.Delta<0).mean(),Pairs=len(g)))
    pd.concat(paired,ignore_index=True).to_csv(data/'all_pairwise_deltas.csv',index=False)
    effects=pd.DataFrame(pair_summary); effects.to_csv(data/'pairwise_summary.csv',index=False)
    d.sort_values('RelativeError',ascending=False).to_csv(data/'all_rows_worst_first.csv',index=False)
    d[d.ScaleEndpoint.astype(bool)].to_csv(data/'parameter_endpoints.csv',index=False)
    # Same-H control: the original direct counterparts remain visible.
    controls=[]
    for (model,suite,method),g in raw.groupby(['Model','Suite','Method']):
        p=g.pivot(index=['Scenario','Seed'],columns='Route',values='RelativeError')
        controls.append(dict(Model=model,Suite=suite,Method=method,DirectMean=p.Direct.mean(),SplitMean=p.Split.mean(),SplitMinusDirect=(p.Split-p.Direct).mean(),MaxAbsPairedDifference=abs(p.Split-p.Direct).max()))
    controls=pd.DataFrame(controls); controls.to_csv(data/'same_prior_route_control.csv',index=False)
    font=Path('C:/Windows/Fonts/msyh.ttc'); assert font.exists(); font_manager.fontManager.addfont(str(font))
    plt.rcParams.update({'font.family':'Microsoft YaHei','font.size':10,'axes.unicode_minus':False,'pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False})
    manifest=[]
    def handles(): return [Line2D([0],[0],color=COLORS[m],marker=MARKERS[m],ls=STYLES[m],label=LABELS[m]) for m in LABELS]
    def save(fig,name,desc):
        legends=[t.get_text() for l in fig.legends for t in l.get_texts()]+[t.get_text() for ax in fig.axes if ax.get_legend() for t in ax.get_legend().get_texts()]
        assert legends,name
        fig.savefig(figs/(name+'.png'),dpi=170,facecolor='white'); fig.savefig(figs/(name+'.pdf'),facecolor='white')
        manifest.append(dict(name=name,description=desc,legends=legends,inches=fig.get_size_inches().tolist(),dpi=170))
        plt.close(fig)
    # All seed means are shown; SD uses six seed means, not scene rows.
    fig,axes=plt.subplots(2,3,figsize=(16,9),sharey=True)
    for ri,suite in enumerate(['Main','LowerSNR']):
        for col,model in enumerate(MODELS):
            ax=axes[ri,col]; g=summary[(summary.Suite==suite)&(summary.Model==model)].set_index('ComparisonMethod')
            for k,m in enumerate(LABELS):
                points=seedmeans[(seedmeans.Suite==suite)&(seedmeans.Model==model)&(seedmeans.ComparisonMethod==m)].RelativeError.to_numpy()
                ax.scatter(k+np.linspace(-.10,.10,len(points)),points,color=COLORS[m],s=16,alpha=.55)
                ax.errorbar(k,g.loc[m,'mean'],yerr=g.loc[m,'std'],color=COLORS[m],marker=MARKERS[m],capsize=4,ls='none',ms=7)
            ax.set(title=MODELS[model]+'\n'+('18dB；五情景 / Five scenes' if ri==0 else '12dB；两AVO情景 / Two AVO scenes'),xticks=range(4),xticklabels=SHORT,ylabel='相对波形误差 / Relative waveform error',ylim=(0,None)); ax.grid(axis='y',alpha=.2)
    fig.legend(handles=handles()+[Line2D([0],[0],marker='o',color='#999999',ls='none',label='小点：每种子情景均值 / Small dots: seed means')],loc='lower center',ncol=3,frameon=False,fontsize=9)
    fig.suptitle('四方法同输入、同估计Q对照 / Four methods, same observations and estimated Q\n大点±误差棒：六种子均值±样本SD / Large points ± bars: mean ± sample SD of six seed means',y=.985)
    fig.subplots_adjust(top=.82,bottom=.14,hspace=.65,wspace=.25); save(fig,'01_overall_four_methods','Both suites separately; exact same inputs; all six seed means and SD; not independent geological repetitions.')
    fig,axes=plt.subplots(1,3,figsize=(17,6),sharey=True)
    main=d[d.Suite=='Main']
    for ax,model in zip(axes,MODELS):
        for k,m in enumerate(LABELS):
            a=main[(main.Model==model)&(main.ComparisonMethod==m)].groupby('Scenario').RelativeError.agg(['mean','std']).reindex(SCENES)
            ax.errorbar(np.arange(5)+(k-1.5)*.045,a['mean'],yerr=a['std'],color=COLORS[m],marker=MARKERS[m],ls=STYLES[m],capsize=2,lw=1)
        ax.set(title=MODELS[model],xticks=range(5),xticklabels=[x.replace(' / ','\n') for x in SCENES.values()],ylabel='相对波形误差 / Relative error',ylim=(0,None)); ax.tick_params(axis='x',labelsize=9); ax.grid(alpha=.2)
    fig.legend(handles=handles(),loc='lower center',ncol=4,frameon=False,fontsize=9); fig.suptitle('18dB逐情景：六种子均值±SD / Scene results: mean ± SD of six seeds'); fig.subplots_adjust(top=.83,bottom=.23,wspace=.22); save(fig,'02_all_scenes','All five scenes, no selected favorable subset; SD across six noise seeds.')
    riskdefs=[('SmoothAVO','AVOGradientError','连续AVO梯度 / Smooth AVO gradient'),('CoherentNoiseAVO','AVOGradientError','干扰下AVO梯度 / Coherent AVO gradient'),('WeakSection','WeakEventError','弱反射振幅 / Weak amplitude'),('FaultSection','FaultNeighborhoodError','断层邻域波形 / Fault neighborhood'),('PinchoutSection','PinchoutTailError','尖灭尾部泄漏 / Pinchout leakage')]
    fig,axes=plt.subplots(2,3,figsize=(16,9))
    for ax,(scene,metric,title) in zip(axes.ravel(),riskdefs):
        for k,m in enumerate(LABELS):
            z=main[(main.Scenario==scene)&(main.ComparisonMethod==m)].groupby('Model')[metric].agg(['mean','std']).reindex(MODELS)
            ax.errorbar(np.arange(3)+(k-1.5)*.04,z['mean'],yerr=z['std'],color=COLORS[m],marker=MARKERS[m],ls=STYLES[m],capsize=2,lw=1)
        ax.set(title=title,xticks=range(3),xticklabels=['强吸收\nLow Q','速度反转\nReversal','薄层\nThin'],ylabel='误差指标 / Error metric'); ax.grid(alpha=.2)
    axes[1,2].axis('off'); axes[1,2].text(.02,.8,'18dB；六种子均值±样本SD\nMean ± sample SD of six seeds\n\n各指标定义不同，不跨面板比较高低\nMetrics differ; compare methods within a panel\n\n尖灭指标按整幅真值范数归一化\nLeakage normalized by full truth norm',va='top',fontsize=11)
    fig.legend(handles=handles(),loc='lower center',ncol=4,frameon=False,fontsize=9); fig.suptitle('保幅及局部代价 / Amplitude and local costs'); fig.subplots_adjust(top=.91,bottom=.13,hspace=.45,wspace=.3); save(fig,'03_amplitude_and_local_risks','Inapplicable metric values remain NaN; pinchout normalization is inherited from MATLAB and uses entire truth norm.')
    # Paired errors and tails: full samples, never truncated to convenient quantiles.
    fig,axes=plt.subplots(1,3,figsize=(16,6),sharey=True)
    p=main.pivot(index=['Model','Scenario','Seed'],columns='ComparisonMethod',values='RelativeError')
    for ax,model in zip(axes,MODELS):
        for k,m in enumerate(LABELS):
            points=p.loc[model,m].to_numpy(); ax.scatter(k+np.linspace(-.13,.13,len(points)),points,color=COLORS[m],marker=MARKERS[m],s=22,alpha=.7)
        ax.axhline(1,color='#999999',ls='--',lw=1); ax.axhline(2,color='#999999',ls=':',lw=1)
        ax.set(title=MODELS[model],xticks=range(4),xticklabels=SHORT,ylabel='道集相对波形误差 / Gather relative error',ylim=(0,None)); ax.grid(axis='y',alpha=.2)
    fig.legend(handles=handles()+[Line2D([0],[0],color='#999999',ls='--',label='误差=1 / Error=1'),Line2D([0],[0],color='#999999',ls=':',label='误差=2 / Error=2')],loc='lower center',ncol=3,frameon=False,fontsize=9); fig.suptitle('18dB全部90组观测：最差结果保留 / All 90 main observations; full tails retained\n每模型30点/方法；位置偏移只用于避免重叠 / 30 points per model/method; deterministic horizontal offsets'); fig.subplots_adjust(top=.81,bottom=.20,wspace=.22); save(fig,'04_all_errors_and_failures','All 30 scenario-seed outcomes per model; fixed deterministic x offsets; full y range includes worst output.')
    # Pre-existing example: fixed first seed and scenario, all methods from one JSON Z.
    waveform_audit=[]
    def getexample(model,scene):
        ex=load(f'results/matlab_test/example_{model}_{scene}.json'); z=np.array(ex['Z'],dtype=float); truth=np.array(ex['truth']); obs=np.array(ex['observed']); t=np.array(ex['t']).ravel()
        assert z.shape==(601,11,2,3) and truth.shape==obs.shape==(601,11)
        signals={'Direct':z[:,:,0,0],'Split':z[:,:,1,0],'Frozen':z[:,:,1,1],'Target08':z[:,:,1,2]}
        for m,s in signals.items():
            expected=float(d[(d.Model==model)&(d.Suite=='Main')&(d.Scenario==scene)&(d.Seed==spec['test_seeds'][0])&(d.ComparisonMethod==m)].RelativeError.iloc[0]); error=np.linalg.norm(s-truth)/np.linalg.norm(truth)
            assert abs(expected-error)<1e-12
            waveform_audit.append(dict(Model=model,Scenario=scene,Method=m,Recomputed=error,Recorded=expected,AbsoluteDifference=abs(error-expected)))
        return ex,t,truth,obs,signals
    fig,axes=plt.subplots(3,2,figsize=(16,11))
    for ri,model in enumerate(MODELS):
        ex,t,truth,obs,signals=getexample(model,'CoherentNoiseAVO'); freq=np.fft.rfftfreq(len(t),t[1]-t[0]); denom=np.max(abs(np.fft.rfft(truth[:,-1])))
        curves={'Truth':truth[:,-1],'Observed':obs[:,-1],**{m:s[:,-1] for m,s in signals.items()}}
        for m,s in curves.items():
            color=COLORS.get(m,{'Truth':'#D55E00','Observed':'#AAAAAA'}.get(m)); label=LABELS.get(m,{'Truth':'无吸收真值 / Truth','Observed':'观测 / Observed'}.get(m)); ls=STYLES.get(m,'-' if m=='Truth' else '--')
            axes[ri,0].plot(t,s,color=color,label=label,ls=ls,lw=1)
            axes[ri,1].plot(freq,abs(np.fft.rfft(s))/denom,color=color,label=label,ls=ls,lw=1)
        axes[ri,0].set(title=MODELS[model],xlabel='时间 / Time (s)',ylabel='绝对振幅 / Absolute amplitude')
        axes[ri,1].set(title='全部频段 / Full frequency band',xlabel='频率 / Frequency (Hz)',ylabel='同一真值谱峰归一化 / Same truth-peak normalization',xlim=(0,250))
        for ax in axes[ri]: ax.grid(alpha=.2)
    h,l=axes[0,0].get_legend_handles_labels(); fig.legend(h,l,loc='lower center',ncol=3,frameon=False,fontsize=9); fig.suptitle('固定首种子20261101、1500m相干干扰道 / First seed, 1500m coherent-noise trace\n四方法同一输入；不逐曲线归一化 / Same input; no per-curve gain'); fig.subplots_adjust(top=.87,bottom=.12,hspace=.7,wspace=.3); save(fig,'05_fixed_trace_and_spectrum','Fixed first held-out seed and far offset; absolute trace amplitudes; spectra divide by one shared truth peak per model; no frequency truncation.')
    for mi,model in enumerate(MODELS):
        ex,t,truth,obs,signals=getexample(model,'SmoothAVO'); matrices={'Truth':truth,'Observed':obs,**signals}; vmax=max(float(abs(s).max()) for s in matrices.values())
        fig,axes=plt.subplots(1,6,figsize=(18,7),sharey=True)
        for ax,(m,s) in zip(axes,matrices.items()):
            im=ax.imshow(s,extent=(-.5,10.5,t[-1],t[0]),aspect='auto',interpolation='nearest',cmap='RdBu_r',vmin=-vmax,vmax=vmax)
            ax.set(title=LABELS.get(m,{'Truth':'无吸收真值\nTruth','Observed':'含噪观测\nObserved'}.get(m)).replace(' / ','\n'),xlabel='炮检距 / Offset (m)',xticks=[0,5,10],xticklabels=['0','750','1500'])
        axes[0].set_ylabel('时间 / Time (s)'); cax=fig.add_axes([.925,.23,.012,.55]); fig.colorbar(im,cax=cax,label='统一绝对振幅 / Shared absolute amplitude')
        fig.legend(handles=[Line2D([0],[0],color='#2166ac',lw=4,label='负振幅 / Negative'),Line2D([0],[0],color='#b2182b',lw=4,label='正振幅 / Positive')],loc='lower center',ncol=2,frameon=False)
        fig.suptitle(MODELS[model]+'：固定首种子连续AVO道集 / First-seed smooth AVO gather\n同一色标，保留全部幅值；不逐道归一化 / Shared color scale, all values retained; no trace normalization'); fig.subplots_adjust(top=.80,bottom=.14,left=.055,right=.90,wspace=.25); save(fig,f'0{6+mi}_same_gather_{model}','All methods from same stored input; each model has a symmetric full-range color scale shared across truth, observation and all four results; no clipping.')
    pd.DataFrame(waveform_audit).to_csv(data/'stored_waveform_metric_recheck.csv',index=False)
    # The two extra direct controls isolate the role of the same prior from route effects.
    fig,axes=plt.subplots(1,3,figsize=(15,6),sharey=True)
    for ax,model in zip(axes,MODELS):
        z=controls[(controls.Model==model)&(controls.Suite=='Main')].set_index('Method').reindex(['Baseline','Frozen','Tuned'])
        ax.plot(range(3),z.DirectMean,'o-',color='#333333',label='相同约束一步 / Direct, matched prior')
        ax.plot(range(3),z.SplitMean,'s--',color='#0072B2',label='相同约束分步 / Split, matched prior')
        ax.set(title=MODELS[model],xticks=range(3),xticklabels=['无结构\nNo H','Frozen','Target08'],ylabel='平均相对波形误差 / Mean relative error',ylim=(0,None)); ax.grid(alpha=.2)
    h,l=axes[0].get_legend_handles_labels(); fig.legend(h,l,loc='lower center',ncol=2,frameon=False); fig.suptitle('相同约束下的一步/分步对照 / Direct versus split with matched priors\n18dB；五情景×六种子等权 / Five scenes × six seeds, equal weight'); fig.subplots_adjust(top=.80,bottom=.20,wspace=.22); save(fig,'09_same_prior_direct_split','Ancillary control uses all 756 original rows; same H, lambda and scale within each paired route; complete nonstationary G and split XA remain different.')
    dump(ROOT/'figure_manifest.json',{'figures':manifest,'software':{'python':platform.python_version(),'numpy':np.__version__,'pandas':pd.__version__,'matplotlib':matplotlib.__version__},'uncertainty':'sample SD of six noise seed means; same geology and fixed calibration Q, not geological replication','source_sha256':sha(source/'results/matlab_test/metrics.csv')})
    def wide(suite,metric='mean'):
        return summary[summary.Suite==suite].pivot(index='Model',columns='ComparisonMethod',values=metric).reindex(index=MODELS,columns=LABELS).rename(index=MODELS,columns=LABELS).reset_index().rename(columns={'Model':'模型 / Model'})
    scene_table=main.pivot_table(index=['Model','Scenario'],columns='ComparisonMethod',values='RelativeError',aggfunc='mean').reindex(columns=LABELS).reset_index(); scene_table.Model=scene_table.Model.map(MODELS); scene_table.Scenario=scene_table.Scenario.map(SCENES); scene_table=scene_table.rename(columns={'Model':'模型 / Model','Scenario':'情景 / Scene',**LABELS})
    localtables=[]
    for scene,metric,label in riskdefs:
        z=main[main.Scenario==scene].pivot_table(index='Model',columns='ComparisonMethod',values=metric,aggfunc='mean').reindex(index=MODELS,columns=LABELS).rename(index=MODELS,columns=LABELS).reset_index(); localtables+=['### '+label,'',table(z),'']
    tailreport=summary[summary.Suite=='Main'][['Model','ComparisonMethod','worst','above_one','above_two','n']].copy(); tailreport.Model=tailreport.Model.map(MODELS); tailreport.ComparisonMethod=tailreport.ComparisonMethod.map(LABELS); tailreport.columns=['模型 / Model','方法 / Method','最差误差 / Worst','误差>1 / Above 1','误差>2 / Above 2','观测数 / Inputs']
    route=controls[controls.Suite=='Main'].copy(); route.Model=route.Model.map(MODELS); route=route.drop(columns='Suite'); route=route.rename(columns={'Model':'模型 / Model','Method':'约束 / Prior','DirectMean':'一步均值 / Direct mean','SplitMean':'分步均值 / Split mean','SplitMinusDirect':'分步减一步 / Split−direct','MaxAbsPairedDifference':'最大配对绝对差 / Max paired gap'})
    report=['# 一步、常规分步、Frozen、Target08：同条件对比 / Four-method paired comparison','','## Material Passport','','- 状态：已完成已有MATLAB输出的配对重分析，未重新执行MATLAB。','- 原实验：ResearchMATLAB，R2024b Update 6；756条方法/路线记录、126组含噪观测，最终求解标志均为0。','- 本次四方法提取504行：主试验90组、12dB压力36组，每组四种方法均齐全；保留其余252行作相同先验一步法对照。','- 既有测试结果此前已经查看；本次没有新调参、没有删样本，不把这次整理称为新的独立验证。','','## 先看结果','','**没有一个版本在三套模型及所有局部指标上都最好。强吸收优先Frozen；速度反转和薄层模型的整体波形优先Target08。常规一步和分步结果非常接近。**','', '## 四种方法究竟比什么','','| 方法 / Method | 定义 / Definition | 结构γ / Gamma | 先导λ倍数 / Pilot multiplier | 噪声残差目标倍数 / Target multiplier |','| --- | --- | ---: | ---: | ---: |','| 一步 / Direct | 完整G一步反演，无多道结构约束 / Full G, no H | 0 | — | 1 |','| 常规分步 / Split | 炮检距后深度补偿，传播协方差白化，无H / Offset then depth, covariance whitening | 0 | — | 1 |','| Frozen分步 / Frozen split | 保守先导学H，炮检距域AVO保护 / Pilot H, AVO protection in offset domain | 10 | 4 | 1 |','| Target08分步 / Target08 split | 与Frozen同结构及先导，更积极残差目标 / Same H and pilot, lower residual target | 10 | 4 | 0.8 |','','常规分步是自动选第二步参数的无结构版本；第一步lambda_x仍固定。它不是最初两步均固定参数的旧实现，也不是论文完整叠前到成像流程。剖面域Frozen/Target08只用先导H，不施加AVO投影。','','## 公平性核对','','同一观测在run_unknown_q_test中生成一次，随后全部方法共用D。真实分层Q只用于合成与事后评分；四者共用独立校准所得有效Q、估计主频和逐道噪声估计。当前有效Q约27.75、67.74、49.34，估计主频39.00、76.52、65.64Hz。已知速度和层厚、受限子波族与独立校准道集仍是先验。','','逐道初始lambda由同一完整算子及噪声规则给出；每种预设再按其残差目标选择整体kappa。kappa依据一步路线残差选择，并给对应分步路线共用；不是两路线各自按真值寻优。结构惩罚与残差目标是有意保留的算法差异。Frozen/Target08的相同约束一步法附在后面，避免把结构收益全归因于分步。','','## 主试验：18dB / Main','','每模型5情景×6种子=30组观测。先对五情景等权，再对六种子平均；相对波形误差越小越好。定义为‖恢复−无吸收真值‖F/‖真值‖F。','',table(wide('Main')),'','![四方法整体比较](figures/01_overall_four_methods.png)','','## 每类情景都保留 / All scenes','',table(scene_table),'','![逐情景对比](figures/02_all_scenes.png)','','## 12dB压力 / Lower SNR','','只包含连续AVO、相干干扰AVO两情景；不与五情景主均值混比。Q点估计和已选预设固定。','',table(wide('LowerSNR')),'','## 最差结果与退化计数 / Tails','','误差>1表示误差平方和超过真值平方和，不等于所有反射都不可见。计数单位是情景/种子道集；不当作独立地质重复。','',table(tailreport),'','![所有主试验点](figures/04_all_errors_and_failures.png)','','## 保幅与局部结构 / Local metrics','','AVO梯度相对误差使用模板振幅对sin²θ回归；弱反射误差对应第二反射振幅；断层指标评价邻域第5–8道；尖灭评价第9–11道第二事件附近窗口，并按整幅真值范数归一化。数值尺度不同，不跨指标比大小。全部定义沿用MATLAB评分代码；评分真值不进入参数选择。','',*localtables,'![局部指标](figures/03_amplitude_and_local_risks.png)','','## 固定首种子波形与道集 / Fixed examples','','例子固定20261101；频谱使用1500m相干干扰道，0–250Hz全部保留，每模型各曲线除以同一个真值谱峰。不是按方法各自归一化。','', '![固定波形及频谱](figures/05_fixed_trace_and_spectrum.png)','']
    for i,m in enumerate(MODELS): report += [f'![{MODELS[m]}](figures/0{6+i}_same_gather_{m}.png)','']
    report+=['## 相同结构信息也给一步法 / Matched-prior control','','四组主比较的一步法无H，而Frozen/Target08带H，所以它们的差异包含结构先验与选参差异。下面让一步法使用对应相同H、lambda和kappa；仍保留G与XA的非稳态矩阵复合误差。','',table(route),'','![相同约束路线对比](figures/09_same_prior_direct_split.png)','','这些同先验路线差距很小。现有证据主要支持某些场景中结构约束和补偿强度的价值，不支持分步形式普遍优于完整算子一步法。','','## 解释及使用判断','','- 强吸收：Frozen主均值及最差结果更好；Target08降低残差目标后可能把模型误差与相干干扰当成可恢复成分，最差2.8558。','- 速度反转、薄层：Target08整体波形与AVO梯度更好；但相对Frozen，三模型弱反射振幅和尖灭泄漏均更差。','- 相对普通方法，多道约束的收益依模型和场景而变；不能把平均改善等同于每个事件改善。','- 不提供按真值选择每道最优方法的“混合结果”，那会泄露评分答案。当前程序也没有已验证的自动可靠回退机制。','','## 数值与数据审计','','- 四方法键一一齐全，126组观测每组均为四行；Q、主频、噪声估计逐键完全一致。','- 24个保存例子/方法波形误差用NumPy从601×11波形重新计算，与MATLAB记录差异<1e-12；这不是完整反演重新求解。','- 3次原始kappa选择触及下界；本次四组提取包含其中3行；最差一例目标失配4.414%，保留且未按0.5%容限排除。','- 原实验求解最大法方程相对残差约1e-7；求解收敛不等于恢复物理可信。','- 原始数据及来源代码快照、SHA256见inputs/与SOURCE_MANIFEST.json；重分析只写本独立目录。','','## 能支持到哪里','','这是三套有限合成模型、固定有效Q估计下的新噪声测试。没有重采样校准Q，没有任意未知子波、实际NMO/偏移或实测资料验证。仍采用Shuey两项、理想NMO和列时变吸收矩阵。图中SD表示六噪声种子条件变化，不是置信区间，也不含Q估计和地质不确定性。所有配对比例仅作描述统计，不作显著性推断。','','## 复核方式','','在具备numpy、pandas、matplotlib、Pillow的Python中执行本脚本，使用--source指定原未知Q工程、--output指定新目录；输出不能覆盖已交付报告。源数据只读。无需连接或启动MATLAB。所有PNG与PDF均有中英图例；不是期刊排版认证。']
    pooled=main.groupby('ComparisonMethod').RelativeError.mean().reindex(LABELS)
    pooled_table=pd.DataFrame({'方法 / Method':[LABELS[m] for m in LABELS],'三模型等权误差 / Equal-model mean':pooled.values})
    pooled_table.to_csv(data/'equal_model_main_summary.csv',index=False)
    at=report.index('## 每类情景都保留 / All scenes')
    report[at:at]=['## 三模型等权汇总 / Equal-model summary','',table(pooled_table),'','本次三模型等权时Frozen平均误差最低；Target08与普通分步的汇总均值几乎相同，但逐模型差异很大。这个汇总权重不代表实测地质的出现频率，不能替代分模型和局部指标。','']
    (ROOT/'REPORT_ZH_EN.md').write_text('\n'.join(report)+'\n',encoding='utf-8')
    # Honest image previews and metadata audit.
    thumbs=[]
    for entry in manifest:
        png=figs/(entry['name']+'.png'); pdf=figs/(entry['name']+'.pdf')
        with Image.open(png) as im:
            assert im.width>1000 and im.height>600; assert pdf.read_bytes()[:4]==b'%PDF'
            thumb=ImageOps.contain(im.convert('RGB'),(600,390)); tile=Image.new('RGB',(620,425),'white'); tile.paste(thumb,((620-thumb.width)//2,25)); ImageDraw.Draw(tile).text((10,5),entry['name'],fill='black'); thumbs.append(tile)
    contact=Image.new('RGB',(1240,425*((len(thumbs)+1)//2)),'#dddddd')
    for k,tile in enumerate(thumbs): contact.paste(tile,((k%2)*620,(k//2)*425))
    contact.save(ROOT/'CONTACT_SHEET.png')
    original_hashes_unchanged=all(sha(Path(p['source']))==p['sha256'] for p in provenance); assert original_hashes_unchanged
    dump(ROOT/'SOURCE_MANIFEST.json',provenance)
    audit=dict(status='PASS',raw_rows=756,observations=126,four_method_rows=504,main_observations=90,stress_observations=36,keys_complete=True,same_Q_frequency_noise_estimate=True,original_hashes_unchanged=original_hashes_unchanged,source_matlab_version=connection['version'],source_final_flags=int(raw.FinalFlag.sum()),max_normal_residual=float(raw.NormalResidual.max()),waveform_rechecks=len(waveform_audit),waveform_max_difference=max(x['AbsoluteDifference'] for x in waveform_audit),fourway_endpoint_rows=int(d.ScaleEndpoint.sum()),fourway_max_target_mismatch=float(d.TargetMismatch.max()),new_matlab_runs=0,figures=len(manifest))
    dump(ROOT/'ANALYSIS_CHECKS.json',audit)
    print(json.dumps(audit,ensure_ascii=False,indent=2)); print(wide('Main').to_string(index=False)); print(wide('LowerSNR').to_string(index=False))

if __name__=='__main__':
    parser=argparse.ArgumentParser(); parser.add_argument('--source',type=Path,default=ROOT.parent/'attenuation_unknown_q_tuning_20261018'); parser.add_argument('--output',type=Path); args=parser.parse_args()
    if args.output:
        args.output.mkdir(exist_ok=False,parents=True)
        shutil.copy2(Path(__file__),args.output/Path(__file__).name)
        if (ROOT/'PLAN_ZH.md').exists(): shutil.copy2(ROOT/'PLAN_ZH.md',args.output/'PLAN_ZH.md')
        ROOT=args.output.resolve()
    generate(args.source)
