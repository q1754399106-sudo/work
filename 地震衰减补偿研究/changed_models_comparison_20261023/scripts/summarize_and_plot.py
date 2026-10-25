from pathlib import Path
import json,csv,hashlib,html,math,zipfile
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]
RUN=Path(json.loads((ROOT/'run.json').read_text(encoding='utf-8'))['run'])
DATA=RUN/'matlab'
OUT=ROOT/'figures';OUT.mkdir(exist_ok=True)
(OUT/'gathers').mkdir(exist_ok=True);(OUT/'wiggles').mkdir(exist_ok=True)
checks=json.loads((DATA/'checks.json').read_text(encoding='utf-8'))
assert checks['complete'] and checks['cases']==100 and checks['rows']==500
assert checks['optimization_flags']==0 and checks['baseline_flags']==0
rows=list(csv.DictReader((DATA/'metrics.csv').open(encoding='utf-8-sig')))
methods=['FixedSplit','FixedSplitSharedCal','DecoupledWeak','FullCalRefit','PrunedFullCalRefit']
metrics=['RelativeError','EventAmplitudeError','WeakEventError','AVOGradientError','NoiseWindowRMS','FaultNeighborhoodError','PinchoutTailError']
plan=json.loads((ROOT/'PLAN.json').read_text(encoding='utf-8'))
models=[m['name'] for m in plan['models']]
legacy_rows=[dict(r) for r in rows]
weak_cache={};method_arrays=dict(FixedSplit='fixed',FixedSplitSharedCal='fixed_shared',DecoupledWeak='baseline',FullCalRefit='unpruned',PrunedFullCalRefit='optimized')
for row in rows:
 row['WeakEventErrorLegacy']=row['WeakEventError'];row['WeakIdentifiableTraces']=0
 if row['Scenario']!='WeakSection':continue
 case=row['CaseID']
 if case not in weak_cache:
  d=json.loads((DATA/'cases'/(case+'.json')).read_text(encoding='utf-8'));model=next(m for m in plan['models'] if m['name']==row['Model']);t=np.asarray(d['t']).reshape(-1)
  aa={k:[] for k in ['truth']+list(method_arrays.values())};valid=[];trueamp=None
  for j in range(11):
   times=np.asarray(model['event_t']).copy();times[1]+=.012*np.sin(j*np.pi/10);z=np.pi*model['fpeak']*(t[:,None]-times[None,:]);W=(1-2*z*z)*np.exp(-z*z);P=np.linalg.pinv(W,rcond=max(W.shape)*np.finfo(float).eps)
   valid.append(np.linalg.cond(W.T@W)<=1e6)
   # Section geometry has identical reflection coefficients on all traces.
   # Recover the true coefficient on the first full-rank trace, avoiding
   # arbitrary decomposition of coincident wavelets on the middle trace.
   if j==0:assert valid[-1];trueamp=float((P@np.asarray(d['truth'])[:,j])[1])
   aa['truth'].append(trueamp)
   for k in method_arrays.values():aa[k].append(float((P@np.asarray(d[k])[:,j])[1]))
  aa={k:np.asarray(v) for k,v in aa.items()};valid=np.asarray(valid);weak_cache[case]=(aa,valid)
 aa,valid=weak_cache[case];a=aa[method_arrays[row['Method']]];tr=aa['truth'];legacy=np.linalg.norm(a-tr)/np.linalg.norm(tr)
 assert abs(legacy-float(row['WeakEventErrorLegacy']))<1e-7*max(1,abs(legacy))
 row['WeakEventError']=float(np.linalg.norm((a-tr)[valid])/np.linalg.norm(tr[valid]));row['WeakIdentifiableTraces']=int(valid.sum())
with (ROOT/'metrics_identifiable.csv').open('w',encoding='utf-8-sig',newline='') as f:
 writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
def mean(sub,key):
 vals=[float(r[key]) for r in sub if math.isfinite(float(r[key]))]
 return sum(vals)/len(vals) if vals else None
def aggregate(sub):
 return {method:{key:mean([r for r in sub if r['Method']==method],key) for key in metrics} for method in methods}
def paired(sub,baseline,key='RelativeError'):
 a={r['CaseID']:float(r[key]) for r in sub if r['Method']==baseline and math.isfinite(float(r[key]))}
 b={r['CaseID']:float(r[key]) for r in sub if r['Method']=='PrunedFullCalRefit' and math.isfinite(float(r[key]))}
 cases=sorted(set(a)&set(b));bad=[dict(case=k,baseline=a[k],optimized=b[k]) for k in cases if b[k]>=a[k]]
 return dict(n=len(cases),wins=len(cases)-len(bad),non_wins=bad)
summary={'scope':'Four new geological configurations, 80 normal observations and 20 colored-noise observations; all frozen solver production .m files unchanged. Repeated scenes/seeds are not independent geological models.','checks':checks,'suites':{}}
summary['weak_metric_note']='Post hoc scoring correction: exclude weak amplitudes with cond(W.T W)>1e6, same identifiability criterion inherited by EventAmplitudeError. ThinPolarity WeakSection trace 6 has coincident events. Raw legacy scores remain untouched in metrics.csv; corrected scores in metrics_identifiable.csv. Waveform/AVO scores unchanged.'
summary['legacy_means']={s:aggregate([r for r in legacy_rows if r['Suite']==s]) for s in ['normal','colored']}
for suite in ['normal','colored']:
 sub=[r for r in rows if r['Suite']==suite]
 summary['suites'][suite]=dict(means=aggregate(sub),per_model={m:aggregate([r for r in sub if r['Model']==m]) for m in models},paired={b:{k:paired(sub,b,k) for k in ['RelativeError','WeakEventError','AVOGradientError']} for b in ['FixedSplit','FixedSplitSharedCal','DecoupledWeak']})
summary['per_scenario']={s:aggregate([r for r in rows if r['Scenario']==s and r['Suite']=='normal']) for s in plan['scenarios']}
summary['normal_weak_regression_models']=[m for m in models if summary['suites']['normal']['per_model'][m]['PrunedFullCalRefit']['WeakEventError']>summary['suites']['normal']['per_model'][m]['FixedSplit']['WeakEventError']]
pairrows={}
for row in rows:pairrows.setdefault(row['CaseID'],{})[row['Method']]=row
worst=max(pairrows,key=lambda k:float(pairrows[k]['PrunedFullCalRefit']['RelativeError'])/float(pairrows[k]['FixedSplit']['RelativeError']))
summary['worst_waveform_ratio_case']=worst
# Audit production source identity against the frozen original algorithm.
source=RUN/'source/core_p1'
expected=json.loads((ROOT/'FROZEN_SOLVER_HASHES.json').read_text(encoding='utf-8'))
assert all(hashlib.sha256((source/p).read_bytes()).hexdigest()==h for p,h in expected.items())
summary['unchanged_solver_files']=len(expected)
manifest=[];weak_amplitudes=[];weak_examples={}
paths=sorted((DATA/'cases').glob('*.json'));assert len(paths)==100
for path in paths:
 d=json.loads(path.read_text(encoding='utf-8'));meta=d['meta'];case=meta['id'];t=np.asarray(d['t']).reshape(-1)
 truth,observed,fixed,new=[np.asarray(d[k]) for k in ['truth','observed','fixed','optimized']]
 if meta['scenario']=='WeakSection':
  aa,valid=weak_cache[case];tr=aa['truth'];rec=dict(case=case,model=meta['model'],suite=meta['suite'],identifiable_traces=int(valid.sum()),true_amplitude_rms=float(np.sqrt(np.mean(tr[valid]**2))))
  for k in ['fixed','optimized']:
   rec[k+'_absolute_rmse']=float(np.sqrt(np.mean((aa[k]-tr)[valid]**2)))
  weak_amplitudes.append(rec)
  if meta['snr_db']==12 and meta['seed']==20262411 and meta['suite']=='normal':weak_examples[meta['model']]={k:np.where(valid,v,np.nan) for k,v in aa.items()}
 arrays=[truth,fixed,new,observed,fixed-truth,new-truth]
 amp=max(float(np.abs(a).max()) for a in arrays[:4]);err=max(float(np.abs(a).max()) for a in arrays[4:])
 selected=(meta['scenario'] in ['WeakSection','SmoothAVO'] and meta['snr_db']==12 and meta['seed']==20262411 and meta['suite']=='normal');dest=OUT/'gathers'/(case+'.png')
 if dest.exists() and (not selected or (OUT/'wiggles'/(case+'.png')).exists()):
  manifest.append(dict(case=case,meta=meta,source=str(path),source_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),png=dest.relative_to(ROOT).as_posix(),image_sha256=hashlib.sha256(dest.read_bytes()).hexdigest(),amplitude_limit=amp,error_limit=err,predetermined_wiggle=selected));continue
 fig=plt.figure(figsize=(12,8),layout='constrained');gs=fig.add_gridspec(2,4,width_ratios=[1,1,1,.045]);ampim=errim=None
 for k,(a,title) in enumerate(zip(arrays,['Reference','Original FixedSplit','Optimized','Observed','FixedSplit error','Optimized error'])):
  ax=fig.add_subplot(gs[k//3,k%3]);limit=amp if k<4 else err
  im=ax.imshow(a,extent=[.5,11.5,t[-1],t[0]],aspect='auto',cmap='RdBu_r',vmin=-limit,vmax=limit,interpolation='nearest')
  if k<4:ampim=im
  else:errim=im
  ax.set_title(title);ax.set_xlabel('Trace')
  if k%3==0:ax.set_ylabel('Time (s)')
 fig.colorbar(ampim,cax=fig.add_subplot(gs[0,3]),label='Amplitude (four gather panels)')
 fig.colorbar(errim,cax=fig.add_subplot(gs[1,3]),label='Error (two error panels)')
 fig.suptitle(case,fontsize=10)
 dest=OUT/'gathers'/(case+'.png');fig.savefig(dest,dpi=150);plt.close(fig)
 selected=(meta['scenario'] in ['WeakSection','SmoothAVO'] and meta['snr_db']==12 and meta['seed']==20262411 and meta['suite']=='normal')
 if selected:
  fig,axes=plt.subplots(1,3,figsize=(10,6),sharey=True,layout='constrained');gain=.42/max(float(np.abs(a).max()) for a in [truth,fixed,new])
  for ax,a,title in zip(axes,[truth,fixed,new],['Reference','Original FixedSplit','Optimized']):
   for j in range(11):
    line=j+1+gain*a[:,j];ax.plot(line,t,'k',lw=.55);ax.fill_betweenx(t,j+1,line,where=a[:,j]>0,color='k',alpha=.65)
   ax.set_title(title);ax.set_xlabel('Trace');ax.set_xlim(0,12);ax.set_ylim(t[-1],t[0])
  axes[0].set_ylabel('Time (s)');fig.suptitle(case,fontsize=10)
  fig.supxlabel(f'Common gain={gain:.5g}; no per-trace normalization',fontsize=9)
  fig.savefig(OUT/'wiggles'/(case+'.png'),dpi=160);fig.savefig(OUT/'wiggles'/(case+'.pdf'));plt.close(fig)
 manifest.append(dict(case=case,meta=meta,source=str(path),source_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),png=dest.relative_to(ROOT).as_posix(),image_sha256=hashlib.sha256(dest.read_bytes()).hexdigest(),amplitude_limit=amp,error_limit=err,predetermined_wiggle=selected))
 if len(manifest)%20==0:print('Plotted',len(manifest),flush=True)
(OUT/'figure_manifest.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2),encoding='utf-8')
summary['weak_absolute_amplitudes']=weak_amplitudes
fig,axes=plt.subplots(2,2,figsize=(11,7),layout='constrained')
for ax,model,label in zip(axes.flat,models,['Deep low Q','Thin polarity','Low contrast','High frequency']):
 aa=weak_examples[model]
 for k,color,marker,label2 in [('truth','black','o','Reference'),('fixed','#777777','s','Original'),('optimized','#D55E00','^','Optimized')]:ax.plot(np.arange(1,12),aa[k],color=color,marker=marker,label=label2,lw=1,ms=3)
 ax.axhline(0,color='k',lw=.6,alpha=.4);ax.set_title(label);ax.set_xlabel('Trace');ax.set_ylabel('Weak reflection coefficient');ax.ticklabel_format(axis='y',style='sci',scilimits=(-2,2));ax.grid(alpha=.15)
axes[0,0].legend(fontsize=8);fig.supxlabel('12 dB, first seed; evaluation only. Panel scales differ. Thin-polarity trace 6: coincident events, unidentifiable amplitudes shown as a gap.',fontsize=8)
fig.savefig(OUT/'weak_amplitudes.png',dpi=180);fig.savefig(OUT/'weak_amplitudes.pdf');plt.close(fig)
labels=['Deep low Q','Thin polarity','Low contrast','High frequency'];colors=['#777777','#0072B2','#D55E00'];names=['Original','Original, shared calibration','Optimized'];displaymethods=['FixedSplit','FixedSplitSharedCal','PrunedFullCalRefit']
fig,axes=plt.subplots(1,3,figsize=(15,4.8),layout='constrained');x=np.arange(4)
for ax,key,title in zip(axes,['RelativeError','WeakEventError','AVOGradientError'],['Waveform error','Identifiable weak reflection error','AVO gradient error']):
 for j,(method,color,label) in enumerate(zip(displaymethods,colors,names)):
  vals=[summary['suites']['normal']['per_model'][m][method][key] for m in models]
  ax.bar(x+(j-1)*.25,vals,width=.24,label=label,color=color)
 ax.set_xticks(x,labels,rotation=20,ha='right');ax.set_ylim(bottom=0);ax.set_title(title);ax.set_ylabel('Mean relative error (lower is better)');ax.grid(axis='y',alpha=.2)
 ax.spines[['right','top']].set_visible(False)
axes[0].legend(fontsize=8)
fig.supxlabel('Per model: waveform n=20, weak n=4, AVO n=8; weak scoring excludes unidentifiable ThinPolarity trace 6. No confidence intervals.',fontsize=9)
fig.savefig(OUT/'comparison.png',dpi=180);fig.savefig(OUT/'comparison.pdf');plt.close(fig)
fig,axes=plt.subplots(1,2,figsize=(10,5),layout='constrained')
modelcolors=['#0072B2','#D55E00','#009E73','#CC79A7']
for ax,suite,title in zip(axes,['normal','colored'],['Normal noise (80 observations)','Colored noise (20 observations)']):
 limit=0
 for model,label,color in zip(models,labels,modelcolors):
  cases=[p for p in pairrows.values() if p['FixedSplit']['Model']==model and p['FixedSplit']['Suite']==suite]
  xs=[float(p['FixedSplit']['RelativeError']) for p in cases];ys=[float(p['PrunedFullCalRefit']['RelativeError']) for p in cases]
  ax.scatter(xs,ys,label=label,color=color,alpha=.8,s=25,edgecolors='white',linewidths=.3);limit=max(limit,max(xs),max(ys))
 limit*=1.07;ax.plot([0,limit],[0,limit],'k--',lw=1);ax.set_xlim(0,limit);ax.set_ylim(0,limit);ax.set_aspect('equal');ax.set_title(title);ax.set_xlabel('Original waveform error');ax.set_ylabel('Optimized waveform error');ax.grid(alpha=.15)
axes[0].legend(fontsize=8);fig.supxlabel('Each point is one shared observation. Below diagonal: optimized has lower error. All cases shown.',fontsize=9)
fig.savefig(OUT/'paired_errors.png',dpi=180);fig.savefig(OUT/'paired_errors.pdf');plt.close(fig)
fig,axes=plt.subplots(1,2,figsize=(10,5),sharey=True,layout='constrained')
for model,label in zip(plan['models'],labels):
 depth=np.r_[0,np.cumsum(np.diff([0]+model['event_t'])*np.array(model['vp'][:3])/2)];bound=np.r_[depth,depth[-1]+200]
 for ax,key in zip(axes,['vp','Q']):
  ax.stairs(model[key],bound,orientation='horizontal',label=label,lw=1.6)
axes[0].invert_yaxis();axes[0].set_ylabel('Depth (m)');axes[0].set_xlabel('P velocity (m/s)');axes[1].set_xlabel('Q');axes[1].legend(fontsize=8);fig.suptitle('New geological profiles (bottom layer drawn for 200 m)')
fig.savefig(OUT/'new_models.png',dpi=180);fig.savefig(OUT/'new_models.pdf');plt.close(fig)
(ROOT/'summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2,allow_nan=False),encoding='utf-8')
options='\n'.join(f'<option value="{html.escape(r["png"])}">{html.escape(r["case"])}</option>' for r in manifest)
(ROOT/'index.html').write_text('''<!doctype html><html lang="zh-CN"><meta charset="utf-8"><title>更换模型后的分步法对比</title><style>body{font:16px system-ui;max-width:1400px;margin:24px auto;padding:0 20px;background:#f7f8fa}img{width:100%;background:white}select{max-width:100%;padding:10px}p{line-height:1.7}</style><h1>更换模型后的分步法对比</h1><p>4个新地质模型，80组普通噪声、20组有色噪声。所有结果完整展示；同一案例的4个道集共用振幅色标，两个误差面板共用误差色标，不做逐道归一化。新增同校准原分步法对照见汇总图和原始指标。</p><img src="figures/new_models.png" alt="新模型速度与Q剖面"><img src="figures/comparison.png" alt="逐模型误差比较"><p><a href="RESULTS_ZH.md">实验报告</a> · <a href="metrics.csv">完整指标</a></p><select id="cases">'''+options+'''</select><img id="gather" alt="参考道集、两种方法及误差对比"><script>const s=document.getElementById('cases'),im=document.getElementById('gather');s.onchange=()=>im.src=s.value;s.onchange();</script></html>''',encoding='utf-8')
import shutil
shutil.copy2(DATA/'metrics.csv',ROOT/'metrics.csv')
shutil.copy2(OUT/'gathers'/(worst+'.png'),OUT/'worst_relative_case.png')
lines=['# 更换模型后的对比结果','','## Material Passport','','- 类型：合成地震补偿实验；状态：实际执行完成。','- 问题：冻结优化方法在4个新地质配置中，是否仍优于最初分步法复现？',f'- MATLAB退出状态见 runs/{RUN.name}/execution.json；原始数据、源码快照、全部结果均保留。',f'- {len(expected)}个原有MATLAB源码文件哈希与上一轮完全一致；只增加实验入口和新模型计划。','','## 设计','','4个全新地质模型，每个覆盖5种场景、18/12 dB、2个新种子，合计80组普通噪声；另加20组12 dB、AR(1)=0.35有色噪声。独立校准18 dB。每个方法共用同一份观测和校准资料。','主对照FixedSplit保留原固定双岭参数1e-3；FixedSplitSharedCal额外使用与优化方法完全相同的完整校准结果。保留上一版和未剪枝结果供归因检查。','全部模型/种子/场景在运行前登记于PLAN.json，算法规则未根据本轮结果调整。','','| 模型 | Q | 反射时间/s | 主频/Hz |','|---|---|---|---|']
for m in plan['models']:lines.append(f'| {m["name"]} | {m["Q"]} | {m["event_t"]} | {m["fpeak"]} |')
lines+=['','## 弱反射评分修正说明（事后复核）','','薄层极性反转模型的WeakSection第6道，第二反射被场景起伏移至0.672秒，与第三反射完全重合。单道数据无法分别确定两者振幅；旧WeakEventError仍评分这一不可辨识分量，因此原始数值会误导比较。','以下表格与汇总图的弱反射误差，采用与原EventAmplitudeError相同的可辨识条件：cond(W转置W)<=1e6。对所有方法统一执行；薄层WeakSection的5组观测各排除第6道的单事件振幅评分，其余弱反射道保留。全部100组观测的完整11道波形、AVO及其他指标均未因此删除或修改。','这是运行后发现退化模板而进行的评分修正，并非运行前登记的新主指标。原评分原样保留在metrics.csv和summary.json的legacy_means；修正评分在metrics_identifiable.csv。修正未用于调参或重新筛选模型。','','| 原始未修正弱反射评分 | 原分步法 | 优化方法 |','|---|---:|---:|']
for suite,zh in [('normal','普通噪声'),('colored','有色噪声')]:
 a=summary['legacy_means'][suite];lines.append(f'| {zh} | {a["FixedSplit"]["WeakEventError"]:.6g} | {a["PrunedFullCalRefit"]["WeakEventError"]:.6g} |')
for suite,zh in [('normal','普通噪声80组'),('colored','有色噪声20组')]:
 ss=summary['suites'][suite];lines+=['',f'## {zh}','','| 方法 | 波形误差 | 弱反射误差 | AVO梯度误差 | 断层邻域误差 | 尖灭尾部误差 |','|---|---:|---:|---:|---:|---:|']
 for method in methods:
  a=ss['means'][method];lines.append('| '+method+' | '+' | '.join(f'{a[k]:.6g}' for k in ['RelativeError','WeakEventError','AVOGradientError','FaultNeighborhoodError','PinchoutTailError'])+' |')
 for baseline in ['FixedSplit','FixedSplitSharedCal']:
  pair=ss['paired'][baseline]['RelativeError'];lines+=['',f'相对{baseline}，优化方法波形误差降低：{pair["wins"]}/{pair["n"]}。']
  if pair['non_wins']:lines.append('未改善案例：'+json.dumps(pair['non_wins'],ensure_ascii=False))
lines+=['','## 逐模型普通噪声结果','','| 模型 | 原分步波形误差 | 优化波形误差 | 原分步弱反射误差 | 优化弱反射误差 |','|---|---:|---:|---:|---:|']
for m in models:
 a=summary['suites']['normal']['per_model'][m];lines.append(f'| {m} | {a["FixedSplit"]["RelativeError"]:.6g} | {a["PrunedFullCalRefit"]["RelativeError"]:.6g} | {a["FixedSplit"]["WeakEventError"]:.6g} | {a["PrunedFullCalRefit"]["WeakEventError"]:.6g} |')
lines+=['','## 分场景普通噪声结果','','| 场景 | 原分步波形误差 | 优化波形误差 |','|---|---:|---:|']
for s in plan['scenarios']:
 a=summary['per_scenario'][s];lines.append(f'| {s} | {a["FixedSplit"]["RelativeError"]:.6g} | {a["PrunedFullCalRefit"]["RelativeError"]:.6g} |')
lines+=['',f'按“优化波形误差/原分步波形误差”自动选择的最差案例：{worst}，图见figures/worst_relative_case.png；不是随机或代表性样本。']
lines+=['','## 弱反射绝对误差补充','','弱反射相对误差会受真值振幅大小影响，以下补充绝对系数RMSE；不替代主指标，也不据此排除不利结果。评分在真值子波模板上估计系数，模板只用于事后评价。','','| 模型 | 真弱反射系数RMS | 原分步绝对RMSE | 优化绝对RMSE |','|---|---:|---:|---:|']
for m in models:
 aa=[a for a in weak_amplitudes if a['model']==m and a['suite']=='normal'];lines.append('| '+m+' | '+' | '.join(f'{np.mean([a[k] for a in aa]):.6g}' for k in ['true_amplitude_rms','fixed_absolute_rmse','optimized_absolute_rmse'])+' |')
lines+=['','## 结论与退化记录','',f'普通噪声下，弱反射平均误差高于原分步法的模型数为{len(summary["normal_weak_regression_models"])}/4：'+', '.join(summary['normal_weak_regression_models'])+'。这些结果不支持“所有新模型的弱反射都改善”。','深层低Q叠加相干干扰的道集显示，约0.65秒附近非真实反射的干扰被明显放大。其真实反射位置的AVO指标仍可较小，因此需要结合完整波形和伪事件评估，不能仅以AVO指标判断结果。','本轮保持算法冻结，未根据退化案例调整阈值或重新运行选择有利结果。新增同校准分步法用于分辨校准与反演路线的作用。']
lines+=['','## 数值与适用范围','',f'优化阶段未收敛标志：{checks["optimization_flags"]}；最大振幅重估正规方程残差：{checks["max_normal"]:.6g}。','相对误差使用真值评分，真值和模型名称不传入反演算法。弱反射指标只在WeakSection计算，AVO指标只在两个AVO场景计算；CSV中的NaN表示不适用。尖灭振幅不可辨识处沿用原评分排除规则，整体波形仍完整评分。','本轮更换地层速度、密度、Q、反射时刻和主频，场景模板、零相位Ricker子波家族、四层/11道和正演物理保持不变。这里只是新地质配置验证，不能称为跨子波、跨正演体系或实测验证。80/20组观测来自4个地质配置，不能当作80/20个独立地质模型。','FixedSplit为按最初毕设思路实现的双岭分步法复现，不能等同于已逐行核对的毕设原始源码。','全部100幅道集共用每案例统一振幅/误差色标；预先固定展示每模型12dB第一个种子的SmoothAVO和WeakSection共8幅波形图，不按效果选图。','弱反射组筛选可能删除真实弱信号；改善百分比不等于振幅恢复比例。指标均值未附置信区间，不据此声称统计显著。','','复现命令：使用scripts/run_local.py --entry run_changed_models --label replay --timeout 5400 --run-pointer <新路径>；后处理使用绘图环境运行scripts/summarize_and_plot.py。']
(ROOT/'RESULTS_ZH.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
print(json.dumps({'checks':checks,'means':{s:summary['suites'][s]['means'] for s in ['normal','colored']}},ensure_ascii=False,indent=2),flush=True)
