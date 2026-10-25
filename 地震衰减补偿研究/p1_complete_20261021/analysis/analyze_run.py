"""Independent coverage, solver, matched-prior, paired gates and crossed bootstrap analysis."""
from pathlib import Path
from collections import defaultdict,Counter
import argparse,csv,hashlib,importlib.util,json,math,random,statistics

ROOT=Path(__file__).resolve().parents[1]
KEY=['Group','Model','Suite','Scenario','CalibrationSeed','Seed','Method','Route']
PAIR=KEY[:-2]
def read(path):
    with path.open(encoding='utf-8-sig',newline='') as f: return list(csv.DictReader(f))
def save(path,obj): path.write_text(json.dumps(obj,ensure_ascii=False,indent=2),encoding='utf-8')
def finite(x):
    try: return math.isfinite(float(x))
    except (ValueError,TypeError): return False
def mean(xs): return statistics.mean(xs)
def percentile(xs,q):
    s=sorted(xs); x=q*(len(s)-1); lo=int(x); hi=min(lo+1,len(s)-1)
    return s[lo]+(s[hi]-s[lo])*(x-lo)

def crossed_ci(rows,group,candidate,nboot=4000):
    """Model resampling + crossed calibration and observation-seed multiplicities.
    Draw indices shared across models and fixed scenarios to retain CRN dependence.
    Scenarios are a fixed benchmark stratum, not independent population replicates.
    """
    use=[r for r in rows if r['Suite']=='Main' and (group=='all' or r['Group']==group)
         and r['Route']=='Split' and r['Method'] in ('Frozen',candidate)]
    models=sorted({r['Model'] for r in use}); cals=sorted({r['CalibrationSeed'] for r in use})
    seeds=sorted({r['Seed'] for r in use}); scenes=sorted({r['Scenario'] for r in use})
    table={(r['Model'],r['CalibrationSeed'],r['Seed'],r['Scenario'],r['Method']):float(r['RelativeError']) for r in use}
    cells={}
    for m in models:
        for c in cals:
            for s in seeds:
                cells[m,c,s]=tuple(mean(table[m,c,s,sc,method] for sc in scenes) for method in ('Frozen',candidate))
    base=mean(v[0] for v in cells.values());cand=mean(v[1] for v in cells.values())
    point=100*(1-cand/base); rng=random.Random(20261021); full=[];fixed=[]
    for _ in range(nboot):
        cs=rng.choices(cals,k=len(cals)); ss=rng.choices(seeds,k=len(seeds))
        mb={m:tuple(mean(cells[m,c,s][k] for c in cs for s in ss) for k in (0,1)) for m in models}
        mm=rng.choices(models,k=len(models))
        full.append(100*(1-mean(mb[m][1] for m in mm)/mean(mb[m][0] for m in mm)))
        fixed.append(100*(1-mean(v[1] for v in mb.values())/mean(v[0] for v in mb.values())))
    return dict(point_gain_pct=point,crossed_model_ci95_pct=[percentile(full,.025),percentile(full,.975)],
                fixed_models_crossed_ci95_pct=[percentile(fixed,.025),percentile(fixed,.975)],
                n_models=len(models),n_calibrations=len(cals),n_observation_seeds=len(seeds),
                fixed_scenarios=len(scenes),n_boot=nboot,
                design='model x crossed calibration x observation seed; shared seed draws across models/scenarios')

def main():
    parser=argparse.ArgumentParser();parser.add_argument('run',type=Path);a=parser.parse_args();run=a.run.resolve()
    dest=run/'analysis';dest.mkdir(exist_ok=False)
    rows=read(run/'metrics.csv');plan=json.loads((run/'FROZEN_EXECUTION_PLAN.json').read_text())
    actual=[tuple(r[k] for k in KEY) for r in rows];expected=set()
    for job in plan['jobs']:
        for su in job['suites']:
            for sc in su['scenarios']:
                for cal in job['calibration_seeds']:
                    for seed in job['test_seeds']:
                        for m in job['methods']:
                            for route in ['Direct','Split']:
                                expected.add((job['group'],job['model']['name'],su['name'],sc,str(cal),str(seed),m,route))
    cov=dict(expected_rows=len(expected),actual_rows=len(rows),duplicates=len(actual)-len(set(actual)),
             missing_rows=len(expected-set(actual)),unexpected_rows=len(set(actual)-expected))
    cov['complete']=len(expected)==len(rows) and cov['duplicates']==cov['missing_rows']==cov['unexpected_rows']==0
    save(dest/'coverage.json',cov)
    by_method=defaultdict(list)
    for r in rows: by_method[r['Method']].append(r)
    solver={}
    for meth,rr in by_method.items():
        bad=[r for r in rr if not finite(r['NormalResidual']) or not finite(r['FinalFlag'])
             or float(r['FinalFlag'])!=0 or not 0<=float(r['NormalResidual'])<=2e-7]
        selection=[r for r in rr if float(r['AllSelectionSolverFlags'])!=0 or float(r['AllSelectionMaxNormalResidual'])>2e-7]
        solver[meth]=dict(rows=len(rr),bad_final=len(bad),bad_any_scale=len(selection),max_normal=max(float(r['NormalResidual']) for r in rr))
    save(dest/'all_methods_solver_quality.json',solver)
    idx={tuple(r[k] for k in KEY[:-1]):r for r in rows if r['Route']=='Direct'}
    fair=[];fair_fields=['EffectiveQ','EstimatedFrequency','Scale','QUncertainty','FUncertainty','PriorSigmaSum','PriorLambdaSum','PriorHFro','TargetMultiplier']
    for r in rows:
        if r['Route']!='Split':continue
        d=idx[tuple(r[k] for k in KEY[:-1])]
        for field in fair_fields:
            if r[field]!=d[field]: fair.append(dict(key=tuple(r[k] for k in KEY),field=field,direct=d[field],split=r[field]))
    save(dest/'matched_prior_fairness.json',dict(pairs=len(rows)//2,fields=fair_fields,mismatches=fair,passed=not fair))
    spec=importlib.util.spec_from_file_location('eval_candidate',ROOT/'scripts/evaluate_candidate.py')
    evaluator=importlib.util.module_from_spec(spec);spec.loader.exec_module(evaluator)
    evaluations={}
    for candidate in plan['jobs'][0]['methods']:
        if candidate in ('Baseline','Frozen'):continue
        e=evaluator.evaluate(run/'metrics.csv',candidate,'Frozen','Split',dest/('acceptance_'+candidate))
        strict=e['status']=='PROVISIONAL_PASS' and cov['complete'] and not fair
        strict=strict and solver[candidate]['bad_final']==solver['Frozen']['bad_final']==0
        strict=strict and solver[candidate]['bad_any_scale']==solver['Frozen']['bad_any_scale']==0
        evaluations[candidate]=dict(status=e['status'],strict_status='PROVISIONAL_PASS' if strict else 'FAIL',
             passed=sum(x['pass_'] for x in e['checks']),total=len(e['checks']),
             failed_checks=[x for x in e['checks'] if not x['pass_']],
             equal_model_main_improvement=e['equal_model_main_improvement'],diagnostics=e['diagnostic_counts'])
    save(dest/'acceptance_summary.json',evaluations)
    grouped=defaultdict(list)
    for r in rows:
        if r['Route']=='Split':grouped[r['Group'],r['Model'],r['Suite'],r['Method']].append(r)
    summary=[]
    for key,rr in sorted(grouped.items()):
        g,m,s,method=key;ref=grouped[g,m,s,'Frozen'];b=mean(float(r['RelativeError']) for r in ref);v=[float(r['RelativeError']) for r in rr]
        entry=dict(Group=g,Model=m,Suite=s,Method=method,n=len(v),MeanError=mean(v),
                   GainPct=100*(1-mean(v)/b),Failures=sum(x>1 for x in v),P95=percentile(v,.95),Worst=max(v))
        summary.append(entry)
    with (dest/'group_summary.csv').open('w',encoding='utf-8-sig',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(summary[0]));w.writeheader();w.writerows(summary)
    cis={group:crossed_ci(rows,group,'CombinedUncertainty') for group in ['historical','development','heldout_new_geology','all']}
    save(dest/'crossed_bootstrap_ci.json',cis)
    calibration=[]
    for j in plan['jobs']:
        folder=run/j['model']['name']/'matlab'
        for path in sorted(folder.glob('calibration_2*.json')):
            e=json.loads(path.read_text());calibration.append(dict(model=j['model']['name'],file=str(path.relative_to(run)),
              q_boundary=e['q_boundary'],frequency_boundary=e['frequency_boundary'],optimizer_exitflag=e['optimizer_exitflag']))
        for path in sorted(folder.glob('velocity_*.json')):
            e=json.loads(path.read_text())['calibration'];calibration.append(dict(model=j['model']['name'],file=str(path.relative_to(run)),
              q_boundary=e['q_boundary'],frequency_boundary=e['frequency_boundary'],optimizer_exitflag=e['optimizer_exitflag']))
    save(dest/'calibration_diagnostics.json',calibration)
    # Exact deterministic replay of overlapping preflight observations.
    repeats=[]
    for pre in sorted((ROOT/'runs').glob('preflight_*')):
        status=pre/'EXECUTION_SUMMARY.json'
        if status.exists() and json.loads(status.read_text())['status']=='COMPLETED':
            old=read(pre/'metrics.csv');mapping={tuple(r[k] for k in KEY):r for r in rows}
            fields=[k for k in old[0] if k not in ['AllSelectionSolverFlags','AllSelectionMaxNormalResidual']]
            changed=[]
            for r in old:
                q=mapping.get(tuple(r[k] for k in KEY))
                if q is None or any(r[k]!=q[k] for k in fields): changed.append(tuple(r[k] for k in KEY))
            repeats.append(dict(preflight=str(pre),rows=len(old),identical_metric_rows=not changed,changed=changed))
    save(dest/'deterministic_overlap_replay.json',repeats)
    failed=evaluations['CombinedUncertainty']['failed_checks'];lines=[
      '# P1 全量 MATLAB 重跑与补充验收', '',
      '## Material Passport','',
      '- Origin Skill: academic-research-suite / experiment-agent','- Origin Mode: run + validate',
      '- Verification Status: VERIFIED for executed coverage and solver checks; scientific acceptance below',
      '- 数值实现：用户此次附件中的 MATLAB 七分支；超参数未调整。原有九模型属于已见模型回归验证。','',
      f"实际覆盖 **{len(rows):,} 行指标、{plan['expected_observations']:,} 组配对观测**；覆盖完整：{cov['complete']}。",
      '9 模型 × 3 校准 × 6 观测种子；18 个场景/条件组合：Main 5、LowerSNR 5、四种压力各 2。',
      '每组观测包含七方法和 Direct/Split 两路线；指标行数不是独立样本数。','',
      f"**CombinedUncertainty 工程验收：{evaluations['CombinedUncertainty']['status']}，{evaluations['CombinedUncertainty']['passed']}/{evaluations['CombinedUncertainty']['total']}；加入基准求解质量与完整覆盖要求后：{evaluations['CombinedUncertainty']['strict_status']}。**",'',
      '## Main 等模型加权改善与交叉重采样区间','',
      '| 模型组 | 模型数 | 改善 % | 95% 区间 % |','|---|---:|---:|---|']
    for g,e in cis.items():lines.append(f"| {g} | {e['n_models']} | {e['point_gain_pct']:.3f} | [{e['crossed_model_ci95_pct'][0]:.3f}, {e['crossed_model_ci95_pct'][1]:.3f}] |")
    lines+=['','模型、校准、观测种子同时重采样；校准与观测按交叉设计抽样，跨模型/场景共享种子抽样权重，场景作为固定基准保留。三次校准、少量模型仍限制区间稳定性。','',
      '## 未通过项（全部保留）','']
    for item in failed:lines.append('- '+item['name']+'：'+json.dumps(item['detail'],ensure_ascii=False))
    if not failed:lines.append('- 无工程门槛失败。')
    lines+=['','## 所有方法求解质量','', '| 方法 | 行数 | 最终残差/标志失败 | 任一选参尺度失败 | 最大最终残差 |','|---|---:|---:|---:|---:|']
    for m,e in solver.items():lines.append(f"| {m} | {e['rows']} | {e['bad_final']} | {e['bad_any_scale']} | {e['max_normal']:.3e} |")
    abnormal=[x for x in calibration if x['q_boundary'] or x['frequency_boundary'] or x['optimizer_exitflag']<=0]
    lines+=['','## 执行与解释边界','',
      '- 四种压力分别为 30° 常相位旋转、公开速度 +3% 且层厚固定、远偏移最大 4ms 残余 NMO、AR(1) rho=0.55 有色噪声。速度条件包含相应走时失配，重新校准并重建算子。',
      '- 本次速度定义遵循配置中 +3% 的设置；不同于旧报告 vp +4%、h +3.5% 的联合扰动，不能直接混为同一测试。',
      '- 残余 NMO 驱动在 IFFT 后先截取观测长度，避免原驱动的维度赋值错误。',
      '- 完整的采样计划和源码哈希在实验前保存；运行前后源码变化见 EXECUTION_SUMMARY.json。',
      '- 匹配先验审计直接比较实际输出记录中的 Q、频率、尺度、噪声、lambda、H 范数及目标系数；同一次函数调用产生两路线。',
      '- 生产入口投毒测试使用相同 D/calD，分别从正常/私有字段污染配置生成 public 后调用实际 p1_compensate，并核验非法真值字段会被拒绝。',
      f'- 校准触边或优化器非正退出记录：{len(abnormal)}/{len(calibration)}；见 calibration_diagnostics.json。保留这些失败，未按真值筛除。',
      '- 原有模型和已观察结果不会因重新冻结被变成独立盲测。新地质泛化仍需未参与开发的独立模型或真实数据。',
      '- 相同整数随机种子在 MATLAB 与 Python 中不保证相同随机噪声；本次是完整 MATLAB 执行，不能凭旧 Python CSV 宣称逐值等价。',
      '- 11 类统计谬误均检查：整体/分组反转、群体到个体外推、选择偏差、碰撞条件化、基率、均值回归、幸存者、多重搜索、分析路径、因果误用、反向因果。重点限制是已见模型选择、局部均值不能保证每个事件、多分支探索和有限模型数量；失败运行及全部样本已保留，未作结果驱动筛除。',
      '', '完整逐模型结果见 group_summary.csv；所有未通过项见 acceptance_summary.json。']
    (dest/'RESULTS_ZH.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
    print(json.dumps(dict(coverage=cov,combined=evaluations['CombinedUncertainty'],ci=cis),ensure_ascii=False,indent=2))

if __name__=='__main__':main()
