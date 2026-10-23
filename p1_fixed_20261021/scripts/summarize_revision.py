from pathlib import Path
import argparse, csv, hashlib, json, math, statistics, zipfile

ROOT=Path(__file__).resolve().parents[1]
parser=argparse.ArgumentParser()
parser.add_argument('--run',type=Path,required=True)
parser.add_argument('--previous',type=Path,required=True)
a=parser.parse_args()
run=a.run.resolve()
def read(path):
    with path.open(encoding='utf-8-sig',newline='') as f: return list(csv.DictReader(f))
rows=read(run/'matlab/metrics.csv'); old=read(a.previous)
execution=json.loads((run/'execution.json').read_text(encoding='utf-8'))
quality=json.loads((run/'matlab/checks.json').read_text(encoding='utf-8'))
reg=json.loads((run/'validation/regression_checks.json').read_text(encoding='utf-8'))
source_drift=[rel for rel,sha in execution['source_hashes'].items()
              if hashlib.sha256((ROOT/rel).read_bytes()).hexdigest()!=sha]
summary={'execution_status':execution['status'],'rows':len(rows),'source_changed_after_start':source_drift,
         'quality':quality,'regressions':reg,'methods':{},'scenes':[]}
methods=sorted({r['Method'] for r in rows})
for method in methods:
    rr=[r for r in rows if r['Method']==method and r['Route']=='Split' and r['Suite']=='Main']
    oo=[r for r in old if r['Method']==method and r['Route']=='Split' and r['Suite']=='Main']
    summary['methods'][method]={'main_mean':statistics.mean(float(r['RelativeError']) for r in rr),
                                'previous_main_mean':statistics.mean(float(r['RelativeError']) for r in oo) if oo else None}
for scene in sorted({r['Scenario'] for r in rows if r['Suite']=='Main'}):
    r={r['Method']:r for r in rows if r['Suite']=='Main' and r['Scenario']==scene and r['Route']=='Split'}
    b=float(r['Frozen']['RelativeError']); c=float(r['CombinedUncertainty']['RelativeError'])
    summary['scenes'].append({'scene':scene,'frozen':b,'combined':c,'improvement_percent':100*(1-c/b)})
summary['main_improvement_percent']=100*(1-summary['methods']['CombinedUncertainty']['main_mean']/summary['methods']['Frozen']['main_mean'])
summary['solver_violations']=[{k:r[k] for k in ('Suite','Scenario','Method','Route','NormalResidual','FinalFlag')}
    for r in rows if not math.isfinite(float(r['NormalResidual'])) or float(r['NormalResidual'])>2e-7 or float(r['FinalFlag'])!=0]
summary['acceptance']={}
for path in sorted(run.glob('acceptance_*/acceptance.json')):
    ac=json.loads(path.read_text(encoding='utf-8'))
    summary['acceptance'][path.parent.name]={'status':ac['status'],'checks':len(ac['checks']),
        'failed_checks':[c for c in ac['checks'] if not c['pass_']]}
summary['velocity_stress']=[]
for path in sorted((run/'matlab').glob('velocity_stress_*.json')):
    ac=json.loads(path.read_text(encoding='utf-8')); est=ac['calibration']
    rr={r['Method']:float(r['RelativeError']) for r in rows if r['Suite']=='StressVelocityError' and r['Route']=='Split'}
    summary['velocity_stress'].append({'vp_multiplier':ac['vp_multiplier'],'operator_change':ac['operator_frobenius_change'],
        'effective_Q':est['effective_Q'],'Q_boundary':est['q_boundary'],'optimizer_exitflag':est['optimizer_exitflag'],
        'errors':rr})
(run/'summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding='utf-8')
lines=['# MATLAB 修改后重跑结果','',f"进程状态：{execution['status']}，总耗时 {execution['elapsed_seconds']:.1f} 秒；共 {len(rows)} 行指标。",
       '',f"全部记录的最大正规方程残差：{quality['max_normal_residual']:.3e}；残差或标志失败记录：{len(summary['solver_violations'])}。",
       '',f"执行开始后源码哈希变化：{len(source_drift)} 个文件。",'',
       '| 场景 | Frozen 误差 | 修改后 Combined 误差 | 改善 |','|---|---:|---:|---:|']
for r in summary['scenes']: lines.append(f"| {r['scene']} | {r['frozen']:.6f} | {r['combined']:.6f} | {r['improvement_percent']:+.2f}% |")
lines+=['',f"五场景平均相对改善 {summary['main_improvement_percent']:.2f}%；单模型、单校准、单观测种子，不能据此证明普遍改善。",'',
        '| 方法 | 修改前五场景均值 | 修改后五场景均值 |','|---|---:|---:|']
for name,r in summary['methods'].items():
    prev='—' if r['previous_main_mean'] is None else f"{r['previous_main_mean']:.6f}"
    lines.append(f"| {name} | {prev} | {r['main_mean']:.6f} |")
lines+=['','## 修正后的压力测试与验收','']
for ac in summary['velocity_stress']:
    b=ac['errors']['Frozen']; c=ac['errors']['CombinedUncertainty']
    lines += [f"实际速度扰动 +{100*(ac['vp_multiplier']-1):.0f}%，厚度保持不变；重建后算子差异 Frobenius 范数 {ac['operator_change']:.6f}。",
              '',f"Frozen 相对误差 {b:.6f}，Combined 相对误差 {c:.6f}，候选较 Frozen 退化 {100*(c/b-1):.2f}%。两者误差都大于 1。",
              '',f"该条件下 Q={ac['effective_Q']:.6f}，Q 触边={ac['Q_boundary']}，校准优化器 exitflag={ac['optimizer_exitflag']}。求解器收敛不等于校准可靠或波形恢复成功。",
              '', '这是速度与相应垂直走时/校准失配的联合压力条件，不是保持全部事件走时不变的纯偏移距动校正实验。']
for name,ac in summary['acceptance'].items():
    lines += ['',f"{name}：**{ac['status']}**，{ac['checks']-len(ac['failed_checks'])}/{ac['checks']} 项通过。"]
    for ch in ac['failed_checks']:
        lines.append('- '+ch['name']+'：'+json.dumps(ch['detail'],ensure_ascii=False))
lines+=['','回归检查全部通过：'+str(reg['pass'])+'。',
        '', '新增目标匹配对照用于检验收益来源，不能只凭 StrongerDamping 对照排除所有阻尼解释。',
        '', '方法为不确定性平均的经验模态正则化；未宣称实现完整数据误差协方差。校准协方差包含经验下限。',
        '', '完整九模型、多校准、多种子网格未在本次执行；历史扩展集的两个局部失败不能因本次小范围回归通过而宣布消失。',
        '', '修订版采用独立 MATLAB 实现；未对原包缺失依赖的 Python 数值双生补造验证。']
(ROOT/'RESULTS_ZH.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
manifest=[]
for base in ('core_p1','scripts'):
    for p in sorted((ROOT/base).rglob('*')):
        if p.is_file() and '__pycache__' not in p.parts:
            manifest.append({'path':p.relative_to(ROOT).as_posix(),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()})
(ROOT/'SOURCE_MANIFEST.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
print(json.dumps({k:summary[k] for k in ('execution_status','rows','source_changed_after_start','quality','main_improvement_percent','solver_violations')},ensure_ascii=False,indent=2))
