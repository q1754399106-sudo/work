"""Paired engineering acceptance checks; not a proof of blind testing or significance."""
from pathlib import Path
import argparse, collections, csv, html, json, math, statistics

ROOT=Path(__file__).resolve().parents[1]
METRICS=['RelativeError','EventAmplitudeError','WeakEventError','FaultNeighborhoodError',
         'PinchoutTailError','AVOInterceptError','AVOGradientError']

def finite(value):
    try: return math.isfinite(float(value))
    except (ValueError,TypeError): return False

def mean(values): return statistics.mean(values)

def percentile(values,p):
    a=sorted(values); x=(len(a)-1)*p; lo=int(x); hi=min(lo+1,len(a)-1)
    return a[lo]+(a[hi]-a[lo])*(x-lo)

def evaluate(path,candidate,baseline,route,output,policy=None):
    output=Path(output)
    if output.exists(): raise ValueError('Use fresh output directory')
    output.mkdir(parents=True)
    with Path(path).open(encoding='utf-8-sig',newline='') as f: rows=list(csv.DictReader(f))
    required={'Model','Suite','Scenario','Seed','Method','Route','NormalResidual','FinalFlag',*METRICS}
    if not rows or not required.issubset(rows[0]): raise ValueError('Missing required columns or empty data')
    if policy is None:
        policy=json.loads((ROOT/'NEXT_EXPERIMENT.json').read_text(encoding='utf-8'))['acceptance']
    keys=['Model','Suite','Scenario','Seed']
    if 'CalibrationSeed' in rows[0]: keys.append('CalibrationSeed')
    def index(method,whichroute):
        out={}
        for r in rows:
            if r['Method']!=method or r['Route']!=whichroute: continue
            k=tuple(r[f] for f in keys)
            if k in out: raise ValueError(f'Duplicate case: {method}/{whichroute}/{k}')
            if any(not x for x in k): raise ValueError('Empty pairing key')
            out[k]=r
        return out
    b=index(baseline,route); c=index(candidate,route)
    bd=index(baseline,'Direct'); cd=index(candidate,'Direct')
    if not b or set(b)!=set(c): raise ValueError('Missing/unpaired baseline or candidate observations')
    checks=[]; groups=[]; paired=[]
    def check(name,passed,detail):
        checks.append(dict(name=name,pass_=bool(passed),detail=detail))
    check('matched_prior_direct_present',set(cd)==set(c) and set(bd)==set(b),
          'CSV presence only; identical prior/Q and tuning-budget provenance must be audited separately.')
    for k in b:
        for r in [b[k],c[k]]:
            for field in ['RelativeError','NormalResidual','FinalFlag','EventAmplitudeError']:
                if not finite(r[field]): raise ValueError(f'Nonfinite {field} at {k}')
            if float(r['RelativeError'])<0: raise ValueError('Negative error')
        expected=['EventAmplitudeError']
        if b[k]['Scenario'].endswith('AVO'): expected+=['AVOInterceptError','AVOGradientError']
        if b[k]['Scenario']=='WeakSection': expected+=['WeakEventError']
        if b[k]['Scenario']=='PinchoutSection': expected+=['PinchoutTailError']
        if b[k]['Scenario']=='FaultSection': expected+=['FaultNeighborhoodError']
        for metric in expected:
            if not finite(b[k][metric]) or not finite(c[k][metric]):
                raise ValueError(f'Missing applicable local metric {metric} at {k}')
            if float(b[k][metric])<0 or float(c[k][metric])<0:
                raise ValueError(f'Negative error metric {metric} at {k}')
        p={f:b[k][f] for f in keys}
        for metric in METRICS:
            bv=b[k][metric]; cv=c[k][metric]
            if finite(bv)!=finite(cv): raise ValueError(f'Asymmetric missing metric {metric} at {k}')
            p[metric+'_baseline']=bv; p[metric+'_candidate']=cv
            p[metric+'_difference']=float(cv)-float(bv) if finite(cv) else ''
        paired.append(p)
    scenarios={r['Scenario'] for r in c.values()}
    check('local_scenario_coverage',{'WeakSection','PinchoutSection','FaultSection'}.issubset(scenarios)
          and any(x.endswith('AVO') for x in scenarios),sorted(scenarios))
    modelkeys=sorted({(r['Suite'],r['Model']) for r in c.values()})
    main_b=[]; main_c=[]
    for suite,model in modelkeys:
        kk=[k for k in c if c[k]['Suite']==suite and c[k]['Model']==model]
        vb=[float(b[k]['RelativeError']) for k in kk]; vc=[float(c[k]['RelativeError']) for k in kk]
        mb,mc=mean(vb),mean(vc)
        group=dict(Suite=suite,Model=model,n=len(kk),baseline=mb,candidate=mc,
                   improvement_percent=100*(1-mc/mb) if mb else None,
                   baseline_median=statistics.median(vb),candidate_median=statistics.median(vc),
                   baseline_p90=percentile(vb,.9),candidate_p90=percentile(vc,.9),
                   baseline_p95=percentile(vb,.95),candidate_p95=percentile(vc,.95),
                   baseline_worst=max(vb),candidate_worst=max(vc),
                   baseline_failures=sum(x>1 for x in vb),candidate_failures=sum(x>1 for x in vc))
        groups.append(group)
        check(f'{suite}/{model}/mean_noninferiority',mc<=mb*(1+policy['any_suite_model_mean_degradation_max']),group)
        check(f'{suite}/{model}/failure_count',group['candidate_failures']<=group['baseline_failures'],
              {'baseline':group['baseline_failures'],'candidate':group['candidate_failures'],'n':len(kk)})
        if suite=='Main': main_b.append(mb); main_c.append(mc)
        for metric in METRICS[1:]:
            local=[k for k in kk if finite(b[k][metric])]
            if not local: continue
            x=mean([float(b[k][metric]) for k in local]); y=mean([float(c[k][metric]) for k in local])
            rel=policy['pinchout_degradation_relative_max'] if metric=='PinchoutTailError' else policy['event_and_weak_and_AVO_error_degradation_max']
            allowance=max(abs(x)*rel,policy['local_near_zero_absolute_allowance'])
            check(f'{suite}/{model}/{metric}',y<=x+allowance,
                  {'baseline':x,'candidate':y,'allowance':allowance,'n':len(local)})
    if not main_b: raise ValueError('Main suite is required')
    gain=1-mean(main_c)/mean(main_b)
    check('main_equal_model_improvement',gain>=policy['main_equal_model_relative_error_improvement_min'],
          {'baseline':mean(main_b),'candidate':mean(main_c),'relative_improvement':gain})
    solverrows=list(c.values())+list(cd.values())
    solverok=all(finite(r['NormalResidual']) and finite(r['FinalFlag']) and float(r['FinalFlag'])==0
                 and 0<=float(r['NormalResidual'])<=policy['normal_residual_max'] for r in solverrows)
    check('candidate_solver_quality',solverok,'Direct and requested route included')
    status='PROVISIONAL_PASS' if all(x['pass_'] for x in checks) else 'FAIL'
    flagged=[]
    for r in c.values():
        if finite(r.get('ScaleEndpoint')) and float(r['ScaleEndpoint'])!=0: flagged.append('scale endpoint')
        if finite(r.get('TargetMismatch')) and float(r['TargetMismatch'])>.005: flagged.append('residual target mismatch')
    result=dict(status=status,candidate=candidate,baseline=baseline,route=route,pairs=len(c),
                equal_model_main_improvement=gain,groups=groups,checks=checks,diagnostic_counts=dict(collections.Counter(flagged)),
                scope='Engineering gates on supplied paired metrics. No automatic proof of independent geology, absence of leakage or statistical significance.',
                unverified_requirements=['new geological holdout','repeated independent calibration','identical information/priors across routes',
                                         'source/config hashes and actual execution','hierarchical uncertainty intervals'])
    (output/'acceptance.json').write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf-8')
    with (output/'paired_differences.csv').open('w',encoding='utf-8-sig',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(paired[0])); writer.writeheader(); writer.writerows(paired)
    lines=[f'# {candidate} 对比 {baseline}：效果评价', '', f'**结论：{status}。** 配对观测 {len(c)} 组，路线 {route}。',
           f'主试验等模型平均波形误差相对变化：改善 {100*gain:.2f}%（负数表示退化）。',
           '', '这是工程门槛检查；尚需核对新地质独立性、校准重复、信息公平性与统计区间。通过不自动替换主版。',
           '', '| Suite / 组别 | Model / 模型 | Baseline / 原版 | Candidate / 候选 | Improvement / 改善% | Worst / 候选最差 | Failures / 候选误差>1 |',
           '|---|---|---:|---:|---:|---:|---:|']
    for g in groups:
        lines.append(f"|{g['Suite']}|{g['Model']}|{g['baseline']:.6f}|{g['candidate']:.6f}|{g['improvement_percent']:.2f}|{g['candidate_worst']:.6f}|{g['candidate_failures']}/{g['n']}|")
    lines+=['','未通过项目：']
    lines += [f"- {x['name']}：{json.dumps(x['detail'],ensure_ascii=False)}" for x in checks if not x['pass_']] or ['- 无数值门槛失败；仍需完成独立性审计。']
    lines+=['',f'触边/目标失配诊断：{result["diagnostic_counts"]}。','完整检查见 acceptance.json；所有逐观测配对差保存在 paired_differences.csv。']
    (output/'EFFECTS_ZH.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
    body='<h1>'+html.escape(candidate+' vs '+baseline)+'</h1><p><b>'+status+'</b>：主试验等模型改善 '+f'{100*gain:.2f}%</p>'
    body+='<p>负改善表示退化。工程门槛与独立性审计分开；未通过项不得省略。</p><table><tr><th>Suite / 组别</th><th>Model / 模型</th><th>Baseline / 原版</th><th>Candidate / 候选</th><th>Gain / 改善%</th><th>Worst / 最差</th></tr>'
    for g in groups:
        body+='<tr>'+''.join('<td>'+html.escape(str(x))+'</td>' for x in [g['Suite'],g['Model'],round(g['baseline'],6),round(g['candidate'],6),round(g['improvement_percent'],2),round(g['candidate_worst'],6)])+'</tr>'
    body+='</table><h2>完整检查 / All checks</h2><ul>'
    for x in checks:
        body+='<li>'+('PASS' if x['pass_'] else '<b>FAIL</b>')+' '+html.escape(x['name'])+'<pre>'+html.escape(json.dumps(x['detail'],ensure_ascii=False))+'</pre></li>'
    body+='</ul><p>源码/真值隔离/地质独立性/分层区间须另行审计。</p>'
    (output/'effects.html').write_text('<!doctype html><meta charset="utf-8"><title>补偿效果评价</title><style>body{max-width:1100px;margin:32px auto;font:16px/1.6 system-ui;padding:16px}table{border-collapse:collapse;width:100%}td,th{padding:9px;border:1px solid #ccd}pre{white-space:pre-wrap;overflow-wrap:anywhere}b{color:#963d22}</style>'+body,encoding='utf-8')
    return result

if __name__=='__main__':
    p=argparse.ArgumentParser(); p.add_argument('--metrics',type=Path,required=True)
    p.add_argument('--candidate',required=True); p.add_argument('--baseline',default='Frozen')
    p.add_argument('--route',choices=['Direct','Split'],default='Split'); p.add_argument('--output',type=Path,required=True)
    a=p.parse_args()
    try:
        result=evaluate(a.metrics,a.candidate,a.baseline,a.route,a.output)
        print(json.dumps({k:result[k] for k in ['status','pairs','equal_model_main_improvement']}))
    except (ValueError,KeyError) as exc:
        raise SystemExit('INVALID_INPUT: '+str(exc))
