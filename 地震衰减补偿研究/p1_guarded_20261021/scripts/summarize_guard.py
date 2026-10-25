from pathlib import Path
import argparse, collections, csv, hashlib, json, math, statistics, zipfile

ROOT=Path(__file__).resolve().parents[1]
p=argparse.ArgumentParser(); p.add_argument('--run',type=Path,required=True); a=p.parse_args(); run=a.run.resolve()
with (run/'matlab/metrics.csv').open(encoding='utf-8-sig',newline='') as f: rows=list(csv.DictReader(f))
rr=[r for r in rows if r['Route']=='Split']
execution=json.loads((run/'execution.json').read_text(encoding='utf-8'))
drift=[rel for rel,sha in execution['source_hashes'].items() if hashlib.sha256((ROOT/rel).read_bytes()).hexdigest()!=sha]
groups=collections.defaultdict(list)
for r in rr: groups[r['Suite'],r['Seed'],r['Method']].append(r)
metrics=('RelativeError','EventAmplitudeError','WeakEventError','FaultNeighborhoodError','AVOInterceptError','AVOGradientError')
records=[]
for (suite,seed,method),rs in sorted(groups.items()):
    rec={'suite':suite,'seed':seed,'method':method,'observations':len(rs),
         'fallbacks':sum(float(r['FallbackUsed'])==1 for r in rs),'compensated':sum(float(r['CompensationApplied'])==1 for r in rs)}
    for m in metrics:
        v=[float(r[m]) for r in rs if math.isfinite(float(r[m]))]
        rec[m]=statistics.mean(v) if v else None
    records.append(rec)
guard=[r for r in rr if r['Method']=='GuardedNoise']
summary={'rows':len(rows),'observations':len(guard),'compensated_observations':sum(float(r['CompensationApplied'])==1 for r in guard),
         'fallback_observations':sum(float(r['FallbackUsed'])==1 for r in guard),'source_drift':drift,'groups':records,
         'max_normal_residual':max(float(r['NormalResidual']) for r in rows),
         'solver_failures':sum(float(r['NormalResidual'])>2e-7 or float(r['FinalFlag'])!=0 for r in rows),
         'execution_status':execution['status']}
summary['acceptance']={}
for path in sorted(run.glob('acceptance_*/acceptance.json')):
    ac=json.loads(path.read_text(encoding='utf-8'))
    summary['acceptance'][path.parent.name]={'status':ac['status'],'failed_checks':[c for c in ac['checks'] if not c['pass_']]}
(run/'summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding='utf-8')
lookup={(r['suite'],r['seed'],r['method']):r for r in records}
lines=['# 保守策略试验结果','', '**结论：不提升为推荐默认版。** +3% 速度条件下可以避免严重失控，但 −3% 条件下过度拒绝了仍有用的补偿结果。',
       '',f"状态 {execution['status']}；耗时 {execution['elapsed_seconds']:.1f} 秒；{len(rows)} 行指标。",
       '',f"GuardedNoise 实际补偿 {summary['compensated_observations']}/{len(guard)} 组观测，回退 {summary['fallback_observations']}/{len(guard)} 组。",
       '',f"所有分支残差或标志失败 {summary['solver_failures']} 行；执行后核心源码哈希变化 {len(drift)} 个文件。",'',
       '## 主场景平均误差','', '| 观测种子 | Frozen | 仅噪声 | 联合 | 保守策略 | 未补偿 | 保守策略相对 Frozen 改善 |',
       '|---|---:|---:|---:|---:|---:|---:|']
methods=('Frozen','NoiseUncertaintyOnly','CombinedUncertainty','GuardedNoise','NoCompensation')
for seed in sorted({r['seed'] for r in records if r['suite']=='Main'}):
    v=[lookup['Main',seed,m]['RelativeError'] for m in methods]
    lines.append('| '+seed+' | '+' | '.join(f'{x:.6f}' for x in v)+f' | {100*(1-v[3]/v[0]):+.2f}% |')
lines+=['','## 速度压力场景','', '| 条件 | 种子 | Frozen | 仅噪声 | 联合 | 保守策略 | 回退 |',
        '|---|---|---:|---:|---:|---:|---:|']
for suite,seed in sorted({(r['suite'],r['seed']) for r in records if r['suite']!='Main'}):
    vals=[lookup[suite,seed,m]['RelativeError'] for m in methods[:4]]
    lines.append('| '+suite+' | '+seed+' | '+' | '.join(f'{x:.6f}' for x in vals)+' | '+str(lookup[suite,seed,'GuardedNoise']['fallbacks'])+' |')
lines+=['','## 局部指标与校准诊断','']
for r in records:
    if r['suite']=='Main' and r['method'] in ('Frozen','NoiseUncertaintyOnly','CombinedUncertainty','GuardedNoise'):
        lines.append(f"- 种子 {r['seed']} / {r['method']}：弱反射误差 {r['WeakEventError']:.6f}，断层邻域误差 {r['FaultNeighborhoodError']:.6f}，AVO 梯度误差 {r['AVOGradientError']:.6f}。")
for path in sorted((run/'matlab').glob('velocity_stress_*.json')):
    j=json.loads(path.read_text(encoding='utf-8')); est=j['calibration']
    lines.append(f"- {path.name}：速度系数 {j['vp_multiplier']}，Q={est['effective_Q']:.6f}；可靠性：{json.dumps(est['reliability'],ensure_ascii=False)}。")
lines+=['','## 验收与补偿覆盖','']
for name,ac in summary['acceptance'].items():
    lines.append(f"- {name}：{ac['status']}。")
    for ch in ac['failed_checks']: lines.append('  - '+ch['name']+'：'+json.dumps(ch['detail'],ensure_ascii=False))
lines+=['','## 解释边界','',
        '两个相反速度条件的校准归一化损失都约为 8，但直接联合反演一个严重失控、另一个仍有用。因此该损失不能单独充当重建质量判据；不能通过调一个门槛就宣称解决问题。',
        '',
        '异常校准的回退只是原始观测透传，没有恢复衰减损失。必须把回退率与补偿误差一起报告。',
        '', 'Section 数据保留 Frozen，局部指标相同来自保守分支选择，不是获得新的弱反射恢复能力。',
        '', '只测一个模型、一次校准、两个观测种子与两个速度条件；不提供普遍鲁棒性或统计显著性结论。',
        '', '校准门槛是参考此前失败制定的探索性工程规则，尚需在独立模型与新校准上验证误拒率、漏拒率。']
(ROOT/'RESULTS_ZH.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
print(json.dumps({k:v for k,v in summary.items() if k!='groups'},ensure_ascii=False,indent=2))
