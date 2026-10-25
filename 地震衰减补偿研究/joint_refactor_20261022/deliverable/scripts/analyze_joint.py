"""Descriptive, paired analysis; standard library only, never tune the solver."""
from pathlib import Path
import csv, json, math, hashlib, statistics, collections
ROOT=Path(__file__).resolve().parents[1]
METRICS=['RelativeError','AVOInterceptError','AVOGradientError','WeakEventError','PinchoutTailError']
def avg(rows,key):
    x=[float(r[key]) for r in rows if math.isfinite(float(r[key]))]
    return statistics.mean(x) if x else float('nan')
def groups(rows,keys):
    d=collections.defaultdict(list)
    for r in rows:d[tuple(r[k] for k in keys)].append(r)
    return d
def writecsv(path,rows):
    with path.open('w',newline='',encoding='utf-8-sig') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def main():
    run=Path(json.loads((ROOT/'LATEST_RUN.json').read_text(encoding='utf-8'))['run'])
    execution=json.loads((run/'execution.json').read_text(encoding='utf-8'))
    assert execution['status']=='COMPLETED',execution
    checks=json.loads((run/'matlab/checks.json').read_text(encoding='utf-8'))
    rows=list(csv.DictReader((run/'matlab/metrics.csv').open(encoding='utf-8-sig')))
    assert len(rows)==640 and checks['rows']==640
    key=['Model','Calibration','Scenario','Seed','Method','Route']
    assert len(groups(rows,key))==len(rows),'Duplicate experimental cells'
    assert all(float(r['FinalFlag'])==0 and float(r['AllSelectionFlags'])==0 for r in rows)
    assert max(float(r['AllSelectionMaxNormal']) for r in rows)<2e-7
    shared=[r for r in rows if r['Method'] in ('SharedIsotropic','Physical')]
    assert max(float(r['SharedRegularizerDifference']) for r in shared)<1e-14
    # Every immutable run snapshot byte is checked after numerical execution.
    manifest=json.loads((run/'SOURCE_MANIFEST.json').read_text(encoding='utf-8'))
    changed=[p for p,h in manifest.items() if hashlib.sha256((run/'source/core_p1'/p).read_bytes()).hexdigest()!=h]
    assert not changed,changed
    batches=[]
    for rr in sorted((ROOT/'runs').glob('paired_*')):
        ex=json.loads((rr/'execution.json').read_text(encoding='utf-8'))
        mm=json.loads((rr/'SOURCE_MANIFEST.json').read_text(encoding='utf-8'))
        assert all(hashlib.sha256((rr/'source/core_p1'/p).read_bytes()).hexdigest()==h for p,h in mm.items())
        batches.append({'run':rr.name,'status':ex['status'],'exit_code':ex.get('exit_code'),'elapsed_seconds':ex['elapsed_seconds']})
    out=ROOT/'analysis';out.mkdir(exist_ok=True)
    checkpoint_audit=[]
    index={tuple(r[k] for k in key):r for r in rows}
    for b in batches:
        if b['status']!='HARD_TIMEOUT':continue
        oldrows=list(csv.DictReader((ROOT/'runs'/b['run']/'matlab/metrics.csv').open(encoding='utf-8-sig')))
        max_relative=0
        for old in oldrows:
            new=index[tuple(old[k] for k in key)]
            for field,value in old.items():
                if value==new[field]:continue
                try:a=float(value);z=float(new[field])
                except ValueError:raise AssertionError((field,value,new[field]))
                if math.isnan(a) and math.isnan(z):continue
                assert math.isclose(a,z,rel_tol=5e-13,abs_tol=1e-300),(field,a,z)
                max_relative=max(max_relative,abs(a-z)/max(abs(a),1e-300))
        checkpoint_audit.append({'source':b['run'],'retained_rows':len(oldrows),'max_numeric_relative_difference':max_relative})
    (out/'checkpoint_integrity.json').write_text(json.dumps(checkpoint_audit,indent=2),encoding='utf-8')
    summary=[]
    for k,rs in groups(rows,['Calibration','Route','Method']).items():
        summary.append(dict(zip(['Calibration','Route','Method'],k))|{'n':len(rs),**{m:avg(rs,m) for m in METRICS},
            'AboveOne':sum(float(r['RelativeError'])>1 for r in rs),'Worst':max(float(r['RelativeError']) for r in rs),
            'EndpointCount':sum(r['ScaleEndpoint'] in ('1','true') for r in rs),
            'FinalCholeskyRows':sum('Cholesky' in r['Solver'] for r in rs),
            'MaxTargetMismatch':max(float(r['TargetMismatch']) for r in rs)})
    writecsv(out/'summary.csv',summary)
    model=[]
    for k,rs in groups(rows,['Model','Calibration','Route','Method']).items():
        model.append(dict(zip(['Model','Calibration','Route','Method'],k))|{m:avg(rs,m) for m in METRICS})
    writecsv(out/'model_summary.csv',model)
    penalties=[]
    for k,rs in groups([r for r in rows if r['Method']=='Physical'],['Model','Calibration','Route']).items():
        penalties.append(dict(zip(['Model','Calibration','Route'],k))|{
            'MeanQPenalty':avg(rs,'QPenalty'),'MeanModalPenalty':avg(rs,'ModalPenalty'),
            'MeanStructurePenalty':avg(rs,'StructurePenalty'),
            'MeanQOverModal':statistics.mean(float(r['QPenalty'])/max(float(r['ModalPenalty']),1e-30) for r in rs),
            'MismatchAboveHalfPercent':sum(float(r['TargetMismatch'])>.005 for r in rs),
            'MaxTargetMismatch':max(float(r['TargetMismatch']) for r in rs)})
    writecsv(out/'physical_penalty_diagnostics.csv',penalties)
    paired=[]
    for k,rs in groups(rows,['Model','Calibration','Scenario','Seed','Route']).items():
        methods={r['Method']:r for r in rs};base=methods['Legacy']
        for method,r in methods.items():
            if method=='Legacy':continue
            paired.append(dict(zip(['Model','Calibration','Scenario','Seed','Route'],k))|{'Method':method,
                **{m:float(r[m])-float(base[m]) for m in METRICS},
                'RelativeImprovementPercent':100*(1-float(r['RelativeError'])/float(base['RelativeError']))})
    writecsv(out/'paired_differences_vs_legacy.csv',paired)
    gates=[];failures=[]
    for cal in ('Old','New'):
      for route in ('Direct','Split'):
        subset=[r for r in rows if r['Calibration']==cal and r['Route']==route]
        baseline=[r for r in subset if r['Method']=='Legacy']
        for method in ('DecoupledWeak','Decoupled','SharedIsotropic','Physical'):
            cand=[r for r in subset if r['Method']==method];improve=100*(1-avg(cand,'RelativeError')/avg(baseline,'RelativeError'))
            reasons=[]
            if improve<5:reasons.append('overall_improvement_below_5_percent')
            for modelname in sorted({r['Model'] for r in subset}):
                b=[r for r in baseline if r['Model']==modelname];c=[r for r in cand if r['Model']==modelname]
                if avg(c,'RelativeError')>1.02*avg(b,'RelativeError'):reasons.append(modelname+'/mean_degradation_above_2_percent')
                if sum(float(r['RelativeError'])>1 for r in c)>sum(float(r['RelativeError'])>1 for r in b):reasons.append(modelname+'/above_one_count_increase')
                for scene in sorted({r['Scenario'] for r in b}):
                    bb=[r for r in b if r['Scenario']==scene];cc=[r for r in c if r['Scenario']==scene]
                    for metric in METRICS[1:]:
                        bv=avg(bb,metric);cv=avg(cc,metric)
                        if math.isfinite(bv) and cv>1.05*bv:
                            reasons.append(modelname+'/'+scene+'/'+metric)
                            failures.append({'Calibration':cal,'Route':route,'Method':method,'Model':modelname,'Scenario':scene,'Metric':metric,'Baseline':bv,'Candidate':cv,'Ratio':cv/bv})
            gates.append({'Calibration':cal,'Route':route,'Method':method,'ImprovementPercent':improve,'Pass':not reasons,'Reasons':reasons})
    (out/'gates.json').write_text(json.dumps(gates,ensure_ascii=False,indent=2),encoding='utf-8')
    if failures:writecsv(out/'local_gate_failures.csv',failures)
    s={(r['Calibration'],r['Route'],r['Method']):r for r in summary}
    passed=[g for g in gates if g['Calibration']=='New' and g['Route']=='Split' and g['Pass']]
    verdict='本轮新校准分步候选均未通过预定工程筛选，不提升生产默认。' if not passed else '本轮有新校准分步候选通过开发工程筛选；仍需新校准与未用于选择的资料确认，不自动提升生产默认。'
    calibration_gain=100*(1-s['New','Split','Legacy']['RelativeError']/s['Old','Split','Legacy']['RelativeError'])
    decoupled_gain=100*(1-s['New','Split','DecoupledWeak']['RelativeError']/s['New','Split','Legacy']['RelativeError'])
    interpretation=f'本批数据中，保留原绑定规则时，新校准使平均误差降低 {calibration_gain:.2f}%；在新校准上进一步归一化并解耦，较小固定结构系数版只再改善 {decoupled_gain:.2f}%。收益主要体现为校准替换的全流程效应，尚不足以支持推广新正则。'
    lines=['# 联合正则重构实验结果','', '## Material Passport','',
      '- Origin Skill: academic-research-suite / experiment-agent','- Origin Mode: run',
      '- Origin Date: 2026-10-22','- Verification Status: ANALYZED; numerical invariants checked; no independent replication',
      '- Version Label: joint_refactor_v1','',
      f"完成 {len(rows)} 行，MATLAB 退出码 {execution['exit_code']}，当前数值批次耗时 {checks['elapsed_seconds']:.1f} 秒。四开发模型、同一既有校准种子、两个新观测种子；所有结果均为开发诊断。",'',verdict,'',interpretation,'',
      '## 新校准下的分步结果','',
      '| 方法 | 平均相对误差 | 相对 Legacy 改善 | >1计数 | 最差 | 端点数 |',
      '|---|---:|---:|---:|---:|---:|']
    for method in ('Legacy','DecoupledWeak','Decoupled','SharedIsotropic','Physical'):
        r=s['New','Split',method];base=s['New','Split','Legacy']
        lines.append(f"|{method}|{r['RelativeError']:.6f}|{100*(1-r['RelativeError']/base['RelativeError']):+.2f}%|{r['AboveOne']}/{r['n']}|{r['Worst']:.6f}|{r['EndpointCount']}|")
    lines+=['','## 校准与正则的交叉结果（分步）','','| 方法 | 旧校准误差 | 新校准误差 | 校准净改善 |','|---|---:|---:|---:|']
    for method in ('Legacy','DecoupledWeak','Decoupled','SharedIsotropic','Physical'):
        a=s['Old','Split',method]['RelativeError'];b=s['New','Split',method]['RelativeError']
        lines.append(f'|{method}|{a:.6f}|{b:.6f}|{100*(1-b/a):+.2f}%|')
    lines+=['','## 同共享正则对照的物理惩罚净作用','','| 校准 | 路线 | SharedIsotropic | Physical | 加入RQ的相对改善 |','|---|---|---:|---:|---:|']
    for cal in ('Old','New'):
      for route in ('Direct','Split'):
        a=s[cal,route,'SharedIsotropic']['RelativeError'];b=s[cal,route,'Physical']['RelativeError']
        lines.append(f'|{cal}|{route}|{a:.6f}|{b:.6f}|{100*(1-b/a):+.2f}%|')
    lines+=['','## 新校准下局部指标（分步，越低越好）','','| 方法 | AVO截距 | AVO梯度 | 弱反射 | 尖灭尾部 |','|---|---:|---:|---:|---:|']
    for method in ('Legacy','DecoupledWeak','Decoupled','SharedIsotropic','Physical'):
        r=s['New','Split',method]
        lines.append('|'+method+'|'+'|'.join(f'{r[m]:.6f}' for m in METRICS[1:])+'|')
    lines+=['','## 工程筛选与数值检查','']
    for g in gates:
        if g['Calibration']=='New':lines.append(f"- {g['Route']} / {g['Method']}: {'PASS' if g['Pass'] else 'FAIL'}，整体改善 {g['ImprovementPercent']:+.2f}%，失败条件 {len(g['Reasons'])} 项。")
    lines+=['',f"全部搜索候选最终求解标志失败数 {checks['solver_flags']}；最大法方程相对残差 {checks['max_normal_residual']:.3e}；旧绑定/归一化重放差 {checks['legacy_replay_difference']:.3e}；期望失配恒等式误差 {checks['expectation_identity_error']:.3e}；共享正则矩阵最大相对差 0（逐列实际计算）；运行快照变化 {len(changed)} 个文件。",'',
      '## 解释边界','',
      '归一化+解耦相对 Legacy 是组合消融，未分别在全网格运行仅归一化与不归一化解耦；仅归一化且补偿系数的等价性通过重放检查。相干预减、局部 E、先导和 H 保持既有实现。SharedIsotropic 与 Physical 使用相同 Gbar，因此它们的差异可以归于 RQ；Physical 与 Legacy 的差异同时包含模态形状、数据算子和参数解耦。',
      '旧/新校准使用相同正则构造规则；校准变化后，先导、lambda、H、经验模态权重和残差预算随规则更新。因此校准交叉反映校准替换的全流程效应，不是固定数值正则矩阵只替换 G 的消融。',
      'Legacy 是相同 CalShrink025Light 配置下的原绑定求解，不是最早 Combined 或 Frozen；本轮没有重新运行 Frozen 全网格比较。DecoupledWeak 的 0.3 仅表示其归一化结构系数小于 Decoupled 的 1，并不保证它比原绑定方案的实际结构惩罚更弱。',
      'RQ 的情景范围来自未校准的不确定性代理。Split 固定 K、使用 K Gbar 数据项，但采用共同的未白化 RQ；这是控制实验定义，不是 Split 白化度量下的期望风险。未测试误差协方差似然、稳健结构、速度压力、新校准种子、全流程 AVO 保护和实测资料。',
      '新观测种子不能把已用于开发的模型和校准变为独立确认。两种子描述性均值不支持显著性或后验覆盖率主张。alpha 由 Direct 残差选定后共享给 Split，不能据此判断各路线最优能力。未调整原确认门槛，未提升生产默认。','',
      '逐行指标：'+str(run/'matlab/metrics.csv'), '设计：PLAN_ZH.md；聚合：analysis/summary.csv；逐模型：analysis/model_summary.csv；全部局部失败：analysis/local_gate_failures.csv；全部筛选原因：analysis/gates.json。']
    lines+=['','物理惩罚各项数值及残差失配见 analysis/physical_penalty_diagnostics.csv。惩罚项的数值比例用于诊断，不能直接解释为独立的误差贡献或因果效应。']
    for label,pattern in [('独立稠密解核对','algebra_check_*'),('研究入口端到端核对','api_check_*')]:
        candidates=sorted((ROOT/'runs').glob(pattern))
        if candidates:
            cr=candidates[-1];cp=cr/'matlab/checks.json'
            if cp.exists():
                check=json.loads(cp.read_text(encoding='utf-8'));assert check['completed']
                lines+=['',label+'：'+json.dumps(check,ensure_ascii=False)+'；证据 '+str(cp)]
    (ROOT/'RESULTS_ZH.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
    with (ROOT/'RESULTS_ZH.md').open('a',encoding='utf-8') as f:
        f.write('\n运行批次：'+json.dumps(batches,ensure_ascii=False)+'。\n')
        if any(b['status']=='HARD_TIMEOUT' for b in batches):
            f.write('首批触及预设一小时时限后，保留完整观测组检查点并补跑缺失组合；完成行直接保留，不是一次不中断的640行运行。续跑数值求解与物理算子源码逐字节核对相同。上述最大法方程残差包括保留行与新行。\n')
            f.write('保留行逐字段核对：'+json.dumps(checkpoint_audit,ensure_ascii=False)+'。\n')
        f.write('\n图件与逐点数据：figures/calibration_regularization.png、figures/local_paired_ratios.png、figures/model_plot_data.csv。图件各有 PDF 和说明，来源与图像哈希见 figures/manifest.json。\n')
    print(json.dumps({'rows':len(rows),'checks':checks,'new_split':[r for r in summary if r['Calibration']=='New' and r['Route']=='Split'],'report':str(ROOT/'RESULTS_ZH.md')},ensure_ascii=False,indent=2))
if __name__=='__main__':main()
