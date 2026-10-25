from pathlib import Path
import json,csv,math,statistics,hashlib
ROOT=Path(__file__).resolve().parents[1]
def readrun(pointer):return Path(json.loads((ROOT/pointer).read_text(encoding='utf-8'))['run'])
def readrows(run):return list(csv.DictReader((run/'matlab/metrics.csv').open(encoding='utf-8-sig')))
def finite_mean(rows,key):
    values=[float(r[key]) for r in rows if math.isfinite(float(r[key]))]
    return statistics.mean(values) if values else None
dev=json.loads((ROOT/'analysis/refit_development_summary.json').read_text(encoding='utf-8'))
validation=readrun('validation_run.json');vcheck=json.loads((validation/'matlab/checks.json').read_text(encoding='utf-8'));assert vcheck['complete'] and vcheck['cases']==50
control=readrun('control_run.json');ccheck=json.loads((control/'matlab/checks.json').read_text(encoding='utf-8'));assert ccheck['complete'] and ccheck['cases']==60
api=readrun('api_run.json');api_check=json.loads((api/'matlab/checks.json').read_text(encoding='utf-8'));assert api_check['complete'] and api_check['api_replay_relative_error']<1e-8 and api_check['truth_input_rejected']
keys=list(dev['means']['FullCalRefit']);vrows=readrows(validation);cr=readrows(control);g=json.loads((ROOT/'core_p1/inputs/adaptive_cv_plan.json').read_text(encoding='utf-8'))['gates']
def summarize(rows):
    methods=sorted({r['Method'] for r in rows});models=sorted({r['Model'] for r in rows})
    means={m:{k:finite_mean([r for r in rows if r['Method']==m],k) for k in keys} for m in methods}
    permodel={model:{m:{k:finite_mean([r for r in rows if r['Method']==m and r['Model']==model],k) for k in keys} for m in methods} for model in models}
    return {'cases':len({r['CaseID'] for r in rows}),'means':means,'per_model':permodel}
main=summarize([r for r in vrows if r['Suite']=='normal']);stress=summarize([r for r in vrows if r['Suite']=='colored']);assert main['cases']==40 and stress['cases']==10
def gates(s):
    a=s['means']['FullCalRefit'];b=s['means']['DecoupledWeak'];old=s['means']['FixedSplit'];ratio=lambda key:a[key]/b[key]
    d={'overall':ratio('RelativeError')<=g['overall_max_ratio'],
       'each_model':all(v['FullCalRefit']['RelativeError']<=g['each_model_max_ratio']*v['DecoupledWeak']['RelativeError'] for v in s['per_model'].values()),
       'weak':ratio('WeakEventError')<=g['weak_max_ratio'],'weak_vs_original':a['WeakEventError']<=g['weak_vs_original_max_ratio']*old['WeakEventError'],
       'avo':ratio('AVOGradientError')<=g['avo_max_ratio'],'noise':ratio('NoiseWindowRMS')<=g['noise_max_ratio'],
       'fault':ratio('FaultNeighborhoodError')<=g['fault_max_ratio'],'pinchout':ratio('PinchoutTailError')<=g['pinchout_max_ratio'],
       'solver':vcheck['optimization_flags']==0 and vcheck['max_normal']<1e-8 and vcheck['baseline_flags']==0}
    d['all_pass']=all(d.values());return d
main['gates']=gates(main)
stress['gates']={'overall':stress['means']['FullCalRefit']['RelativeError']<=stress['means']['DecoupledWeak']['RelativeError'],'weak':stress['means']['FullCalRefit']['WeakEventError']<=stress['means']['DecoupledWeak']['WeakEventError'],'solver':vcheck['optimization_flags']==0}
stress['gates']['all_pass']=all(stress['gates'].values())
for s,kind in [(main,'normal'),(stress,'colored')]:
    rr=[r for r in vrows if r['Suite']==kind];base={r['CaseID']:r for r in rr if r['Method']=='DecoupledWeak'}
    s['case_wins']=sum(float(r['RelativeError'])<float(base[r['CaseID']]['RelativeError']) for r in rr if r['Method']=='FullCalRefit')
controlmean={k:finite_mean(cr,k) for k in keys};dev['means']['FullCalLinear']=controlmean
frozen=json.loads((ROOT/'core_p1/inputs/frozen_refit_profile.json').read_text(encoding='utf-8'))
for run in [Path(dev['run']),validation,api]:
    manifest=json.loads((run/'SOURCE_MANIFEST.json').read_text(encoding='utf-8'));source=run/'source/core_p1'
    assert all(hashlib.sha256((source/p).read_bytes()).hexdigest()==value for p,value in manifest.items())
    assert all(hashlib.sha256((source/p).read_bytes()).hexdigest()==value for p,value in frozen['algorithm_sha256'].items())
assert all(hashlib.sha256((ROOT/'core_p1'/p).read_bytes()).hexdigest()==value for p,value in frozen['algorithm_sha256'].items())
pilot=readrun('pilot_run.json');pc=json.loads((pilot/'matlab/checks.json').read_text(encoding='utf-8'));assert pc['complete'] and pc['cases']==9
pr=readrows(pilot);ids={r['CaseID'] for r in pr};plan=json.loads((ROOT/'core_p1/inputs/adaptive_cv_plan.json').read_text(encoding='utf-8'))
oldrows=list(csv.DictReader((Path(plan['cache'])/'metrics.csv').open(encoding='utf-8-sig')))
pilotreport=summarize(pr+[r for r in oldrows if r['CaseID'] in ids]);pilotreport['scope']='9-case development feasibility screen; not a 60-case or held-out validation result'
ready=dev['gates']['FullCalRefit']['all_pass'] and main['gates']['all_pass'] and stress['gates']['all_pass']
result={'method':'FullCalRefit','ready_for_research_use':ready,'development':dev,'validation':main,'colored_stress':stress,'cv_pilot':pilotreport,'api':api_check,'validation_checks':vcheck,'frozen_algorithm_unchanged':True,'validation_run':str(validation),'control_run':str(control),'pilot_run':str(pilot),'scope':'Five synthetic geological models in the same forward simulator; zero-phase Ricker wavelet, independent known-geometry calibration, 11 traces. No field-data generalization claim.'}
(ROOT/'analysis/round1_summary.json').write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf-8')
status={'method':'FullCalRefit','ready_for_research_use':ready,'research_only':True,'development_cases':60,'new_geology_cases':40,'colored_stress_cases':10,'gates':{'development':dev['gates']['FullCalRefit'],'validation':main['gates'],'colored_stress':stress['gates']},'api_replay_relative_error':api_check['api_replay_relative_error'],'scope':result['scope']}
(ROOT/'core_p1/inputs/round1_validation_status.json').write_text(json.dumps(status,ensure_ascii=False,indent=2),encoding='utf-8')
lines=['# 统一优化结果','','最终冻结方法：完整校准＋非平稳反射系数估计＋振幅重估（FullCalRefit）。同一入口、同一阈值与求解规则适用于所有模型；事件位置来自观测，不使用真实事件数、位置、幅度或真实Q。','', '| 资料 | 方法 | 整体误差 | 弱反射误差 | AVO梯度误差 | 噪声窗RMS |','|---|---|---:|---:|---:|---:|']
for label,s in [('开发60组',dev),('新增模型40组',main),('有色噪声10组',stress)]:
    for m in ['FixedSplit','DecoupledWeak','FullCalRefit']:
        lines.append('| '+label+' | '+m+' | '+' | '.join(f'{s["means"][m][k]:.6g}' for k in ['RelativeError','WeakEventError','AVOGradientError','NoiseWindowRMS'])+' |')
lines+=['','开发/新增模型/有色噪声预设门槛：'+str([dev['gates']['FullCalRefit']['all_pass'],main['gates']['all_pass'],stress['gates']['all_pass']])+'。','',f'整体误差低于现版的案例数：开发 {dev["case_wins"]["FullCalRefit"]}/60；新增模型 {main["case_wins"]}/40；有色噪声 {stress["case_wins"]}/10。','', '| 开发集控制 | 整体误差 | 弱反射误差 |','|---|---:|---:|']
for m in ['DecoupledWeak','FullCalLinear','SameCalRefit','FullCalRefit']:lines.append(f'| {m} | {dev["means"][m]["RelativeError"]:.6f} | {dev["means"][m]["WeakEventError"]:.6f} |')
lines+=['','FullCalLinear仅改变校准规则，SameCalRefit仅改反射系数重估；FullCalRefit同时采用两者。FullCalLinear为冻结方案后补充的归因对照，未用于选择冻结方案。与现版共用校准观测，但Section域的Q/频率处理不同，不能称为完全相同的反演算子。','', '| CV可行性筛选（9组） | 整体误差 | 弱反射误差 |','|---|---:|---:|']
for m in ['DecoupledWeak','AdaptiveJointCV','RobustJointCV']:lines.append(f'| {m} | {pilotreport["means"][m]["RelativeError"]:.6f} | {pilotreport["means"][m]["WeakEventError"]:.6f} |')
lines+=['','CV支线仅做9组开发可行性筛选，不把它当作60组结果。六组保存后进程被外部中断，续跑复用原结果，只计算剩余三组。冻结使用反射系数域方法；普通/鲁棒CV保留源码和全量筛选记录。','',f'公共API输出相对复现误差 {api_check["api_replay_relative_error"]:.3g}；真实Q输入拒绝测试通过。冻结源码哈希与开发、验证、API运行均一致。',f'新增验证最大振幅重估正规方程残差 {vcheck["max_normal"]:.3g}，稀疏阶段未收敛标志 {vcheck["optimization_flags"]}。','', '适用范围：同一仿真器内5个地质模型、给定采集几何及独立校准资料、11道、零相位Ricker子波假设。跨模型收益已检验；尚无实测资料或任意子波形状的验证，不将这些合成试验写成普遍适用的结论。','', '绘图不做逐道归一化，振幅图和误差图分别使用明确的公共色标。所有观测、输出、评分、源码快照及失败运行保留在本目录。']
(ROOT/'ROUND1_RESULTS_ZH.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
print(json.dumps({'ready_for_research_use':ready,'main':main,'colored':stress,'control':controlmean},ensure_ascii=False,indent=2))
