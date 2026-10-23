from pathlib import Path
import csv,json,math,statistics,hashlib
ROOT=Path(__file__).resolve().parents[1]
run=Path(json.loads((ROOT/'latest_run.json').read_text(encoding='utf-8'))['run'])
checks=json.loads((run/'matlab/checks.json').read_text(encoding='utf-8'))
assert checks['complete'] and checks['cases']==60 and checks['rows']==180
rows=list(csv.DictReader((run/'matlab/metrics.csv').open(encoding='utf-8-sig')))
consensus=ROOT/'consensus_run.json'
if consensus.exists():
    crun=Path(json.loads(consensus.read_text(encoding='utf-8'))['run'])
    cchecks=json.loads((crun/'matlab/checks.json').read_text(encoding='utf-8'))
    assert cchecks['complete'] and cchecks['cases']==60
    rows+=list(csv.DictReader((crun/'matlab/metrics.csv').open(encoding='utf-8-sig')))
metrics=['RelativeError','EventAmplitudeError','WeakEventError','AVOGradientError','NoiseWindowRMS','FaultNeighborhoodError','PinchoutTailError']
def mean(method,key,model=None):
    v=[float(r[key]) for r in rows if r['Method']==method and (model is None or r['Model']==model) and math.isfinite(float(r[key]))]
    return statistics.mean(v) if v else None
models=sorted({r['Model'] for r in rows});methods=['Baseline','SourceMixture','SourceOnly']
if consensus.exists():methods.append('Consensus')
summary={m:{key:mean(m,key) for key in metrics} for m in methods}
plan=json.loads((run/'source/core_p1/inputs/universal_plan.json').read_text(encoding='utf-8'))
control=list(csv.DictReader((Path(plan['cache'])/'metrics.csv').open(encoding='utf-8-sig')))
control=[r for r in control if r['Method']=='FixedSplit']
assert len(control)==60
summary['OriginalFixedSplit']={key:statistics.mean(float(r[key]) for r in control if math.isfinite(float(r[key]))) for key in metrics}
g=plan['gates'];gates={}
for m in methods[1:]:
    ratio=lambda key:summary[m][key]/summary['Baseline'][key]
    gates[m]={
        'overall':ratio('RelativeError')<=g['overall_relative_error_max_ratio'],
        'each_model':all(mean(m,'RelativeError',model)<=g['each_model_relative_error_max_ratio']*mean('Baseline','RelativeError',model) for model in models),
        'weak':ratio('WeakEventError')<=g['weak_error_max_ratio'],
        'avo':ratio('AVOGradientError')<=g['avo_gradient_max_ratio'],
        'noise':ratio('NoiseWindowRMS')<=g['noise_window_rms_max_ratio'],
    }
    gates[m]['all_pass']=all(gates[m].values())
out=ROOT/'analysis';out.mkdir(exist_ok=True)
result={'run':str(run),'cases':60,'means':summary,'gates':gates,'models':{model:{m:mean(m,'RelativeError',model) for m in methods} for model in models}}
result['per_model_metrics']={model:{m:{key:mean(m,key,model) for key in metrics} for m in methods} for model in models}
base_by_case={r['CaseID']:r for r in rows if r['Method']=='Baseline'}
result['waveform_case_wins']={m:sum(float(r['RelativeError'])<float(base_by_case[r['CaseID']]['RelativeError']) for r in rows if r['Method']==m) for m in methods[1:]}
(out/'summary.json').write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf-8')
lines=['# 统一方法开发集筛选结果','','全部60组资料使用同一组固定规则。以下数值越低越好。','', '| 方法 | 整体误差 | 事件振幅误差 | 弱反射误差 | AVO梯度误差 | 噪声窗RMS |','|---|---:|---:|---:|---:|---:|']
for m in ['OriginalFixedSplit']+methods:
    lines.append('| '+m+' | '+' | '.join(f'{summary[m][key]:.6g}' for key in metrics[:5])+' |')
lines+=['','| 候选 | 整体 | 各模型 | 弱反射 | AVO | 噪声 | 全部通过 |','|---|---|---|---|---|---|---|']
for m,d in gates.items():lines.append('| '+m+' | '+' | '.join('通过' if d[k] else '未通过' for k in ['overall','each_model','weak','avo','noise','all_pass'])+' |')
lines+=['','| 模型 | '+' | '.join(methods)+' |','|---|'+'---:|'*len(methods)]
for model in models:lines.append('| '+model+' | '+' | '.join(f'{mean(m,"RelativeError",model):.6f}' for m in methods)+' |')
if not any(d['all_pass'] for d in gates.values()):
    lines+=['','本轮候选均未通过预设门槛，不能宣布统一优化成功，也不将其设为默认入口；不开展以晋级为目的的新增模型验证。保留全部失败结果，后续应改变机制，再进行新的开发与冻结验证。']
if consensus.exists():lines+=['','Consensus为固定50%混合频谱与50%纯频谱输出的平均。该补充候选是在第一轮部分开发结果出现取舍后登记的，属于第二轮开发，不是第一轮预登记候选。没有搜索融合权重，没有按模型或场景选择分支。分量的求解残差可检查；融合结果本身不是某一个正规方程的解。']
lines+=['','研究候选 `universal_candidate.m` 只读取观测构建的上下文和校准配置，不读取模型名称、场景或评分真值。实验驱动程序读取真值的用途仅为重建评分模板。','', '最终目标仍是一套统一入口和规则，参数仅随观测与校准资料自动估计。本轮未产生通过验收的整体替换版。新增验证模型保持未运行状态，避免将开发集上的选择当作泛化证据。','',f'原始运行：`{run}`。源码及预设门槛随运行快照保存；本报告由 metrics.csv 自动生成。']
(ROOT/'RESULTS_ZH.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
print(json.dumps(result,ensure_ascii=False,indent=2))
