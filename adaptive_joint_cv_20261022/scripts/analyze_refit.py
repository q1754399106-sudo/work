from pathlib import Path
import csv,json,math,statistics,hashlib
ROOT=Path(__file__).resolve().parents[1]
run=Path(json.loads((ROOT/'refit_run.json').read_text(encoding='utf-8'))['run'])
checks=json.loads((run/'matlab/checks.json').read_text(encoding='utf-8'))
assert checks['complete'] and checks['cases']==60 and checks['rows']==120
rows=list(csv.DictReader((run/'matlab/metrics.csv').open(encoding='utf-8-sig')))
plan=json.loads((run/'source/core_p1/inputs/adaptive_cv_plan.json').read_text(encoding='utf-8'))
control=list(csv.DictReader((Path(plan['cache'])/'metrics.csv').open(encoding='utf-8-sig')))
methods=['FixedSplit','DecoupledWeak','SameCalRefit','FullCalRefit']
keys=['RelativeError','EventAmplitudeError','WeakEventError','AVOGradientError','NoiseWindowRMS','FaultNeighborhoodError','PinchoutTailError']
allrows=control+rows
def mean(method,key,model=None,subset=None):
    values=[float(r[key]) for r in allrows if r['Method']==method and (model is None or r['Model']==model) and (subset is None or r['CaseID'] in subset) and math.isfinite(float(r[key]))]
    return statistics.mean(values) if values else None
models=sorted({r['Model'] for r in rows});summary={m:{k:mean(m,k) for k in keys} for m in methods};g=plan['gates'];gates={}
for m in methods[2:]:
    ratio=lambda key:summary[m][key]/summary['DecoupledWeak'][key]
    d={'overall':ratio('RelativeError')<=g['overall_max_ratio'],
       'each_model':all(mean(m,'RelativeError',model)<=g['each_model_max_ratio']*mean('DecoupledWeak','RelativeError',model) for model in models),
       'weak':ratio('WeakEventError')<=g['weak_max_ratio'],
       'weak_vs_original':summary[m]['WeakEventError']<=g['weak_vs_original_max_ratio']*summary['FixedSplit']['WeakEventError'],
       'avo':ratio('AVOGradientError')<=g['avo_max_ratio'],
       'noise':ratio('NoiseWindowRMS')<=g['noise_max_ratio'],
       'fault':ratio('FaultNeighborhoodError')<=g['fault_max_ratio'],
       'pinchout':ratio('PinchoutTailError')<=g['pinchout_max_ratio'],
       'solver':checks['optimization_flags']==0 and checks['max_normal']<1e-8}
    d['all_pass']=all(d.values());gates[m]=d
base={r['CaseID']:r for r in control if r['Method']=='DecoupledWeak'}
result={'run':str(run),'cases':60,'means':summary,'gates':gates,'per_model':{model:{m:{k:mean(m,k,model) for k in keys} for m in methods} for model in models},
        'case_wins':{m:sum(float(r['RelativeError'])<float(base[r['CaseID']]['RelativeError']) for r in rows if r['Method']==m) for m in methods[2:]},'checks':checks}
out=ROOT/'analysis';out.mkdir(exist_ok=True)
(out/'refit_development_summary.json').write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf-8')
lines=['# 统一优化开发集结果','','同一60组观测、同一份校准资料。FullCalRefit在所有采集域使用完整校准估计；SameCalRefit沿用现版的分域估计。因此FullCalRefit与现版在Section域的算子并非相同，不能把全部收益归于振幅重估。所有候选规则均跨模型固定，未使用真实事件位置或振幅。','', '| 方法 | 整体误差 | 弱反射误差 | AVO梯度误差 | 噪声窗RMS |','|---|---:|---:|---:|---:|']
for m in methods:lines.append('| '+m+' | '+' | '.join(f'{summary[m][k]:.6g}' for k in ['RelativeError','WeakEventError','AVOGradientError','NoiseWindowRMS'])+' |')
lines+=['','| 候选 | 全部门槛 | 整体误差更低的案例数 |','|---|---|---:|']
for m in methods[2:]:lines.append(f'| {m} | {"通过" if gates[m]["all_pass"] else "未通过"} | {result["case_wins"][m]}/60 |')
lines+=['','| 模型 | 现版整体误差 | 同校准重估 | 完整校准重估 |','|---|---:|---:|---:|']
for model in models:lines.append('| '+model+' | '+' | '.join(f'{mean(m,"RelativeError",model):.6f}' for m in methods[1:])+' |')
lines+=['','本报告仅为开发集证据。通过后须冻结规则，完成新增模型验证与公共API复现，才可作为本研究适用范围内的整体替换方案。零相位Ricker子波假设、独立校准资料、已知采集几何及11道输入约束必须随版本说明保留。','',f'运行目录：`{run}`。']
(ROOT/'DEVELOPMENT_RESULTS_ZH.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
print(json.dumps({'cases':60,'means':summary,'gates':gates,'case_wins':result['case_wins']},ensure_ascii=False,indent=2))
