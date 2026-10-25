from pathlib import Path
import csv, json, statistics

HERE = Path(__file__).resolve().parent
outdir = HERE / 'matlab_smoke_corrected'
status = json.loads((outdir / 'status.json').read_text(encoding='utf-8'))
with (outdir / 'metrics.csv').open(encoding='utf-8-sig', newline='') as f:
    rows = list(csv.DictReader(f))
split = {(r['Scenario'],r['Method']):r for r in rows if r['Route']=='Split'}
scenes = sorted({r['Scenario'] for r in rows})
scene_stats = []
for scene in scenes:
    b = float(split[scene,'Frozen']['RelativeError'])
    c = float(split[scene,'CombinedUncertainty']['RelativeError'])
    scene_stats.append({'scenario':scene,'frozen':b,'combined':c,'improvement_percent':100*(1-c/b)})
bmean = statistics.mean(r['frozen'] for r in scene_stats)
cmean = statistics.mean(r['combined'] for r in scene_stats)
violations = [{k:r[k] for k in ('Scenario','Method','Route','NormalResidual','FinalFlag')} for r in rows
              if float(r['NormalResidual'])>2e-7 or float(r['FinalFlag'])!=0]
primary = [r for r in rows if r['Method']=='CombinedUncertainty']
summary = {'execution':status,'scene_stats':scene_stats,'frozen_mean':bmean,'combined_mean':cmean,
           'improvement_percent':100*(1-cmean/bmean),'solver_violations':violations,
           'combined_max_normal_residual':max(float(r['NormalResidual']) for r in primary),
           'notes':'One model, one seed: descriptive smoke test only, not statistical evidence of generalization.'}
(HERE / 'smoke_summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding='utf-8')
lines=['','## 本机试跑实测结果','',
       f"状态：{status['status']}；{status['rows']} 行指标；MATLAB 内计算耗时 {status['elapsed_seconds']:.1f} 秒。",'',
       '| 场景 | Frozen Split 相对误差 | Combined Split 相对误差 | 改善 |',
       '|---|---:|---:|---:|']
for r in scene_stats:
    lines.append(f"| {r['scenario']} | {r['frozen']:.6f} | {r['combined']:.6f} | {r['improvement_percent']:+.2f}% |")
lines += ['', f'五场景平均误差 {bmean:.6f} → {cmean:.6f}，相对改善 {summary["improvement_percent"]:.2f}%。仅是单模型单种子的描述性结果。',
          '', f'联合分支两路线最大正规方程相对残差 {summary["combined_max_normal_residual"]:.3e}。超过 2e-7 门槛或求解标志非零的记录共 {len(violations)} 行：','']
for r in violations:
    lines.append(f"- {r['Scenario']} / {r['Method']} / {r['Route']}：NormalResidual={float(r['NormalResidual']):.3e}，FinalFlag={r['FinalFlag']}。")
lines += ['', '无结构分支的实际残差失败与源码缺陷一致；不能因为进程 exit code=0 或 FinalFlag=0 就判全部求解通过。详见 `smoke_summary.json`、`matlab_smoke_corrected/metrics.csv` 与 `matlab_run2_stdout.log`。']
with (HERE / 'REVIEW_ZH.md').open('a',encoding='utf-8') as f:
    f.write('\n'.join(lines)+'\n')
print(json.dumps(summary,ensure_ascii=False,indent=2))
