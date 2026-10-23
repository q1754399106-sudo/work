"""Report terminal models without changing frozen methods or final analysis."""
from pathlib import Path
import argparse, csv, datetime, hashlib, importlib.util, json, statistics, sys
sys.dont_write_bytecode = True

def read(p): return json.loads(p.read_text(encoding='utf-8-sig'))
def save(p, value): p.write_text(json.dumps(value, ensure_ascii=False, indent=2)+'\n', encoding='utf-8')
def rows(p):
    with p.open(encoding='utf-8-sig', newline='') as f: return list(csv.DictReader(f))
def module(name, p):
    spec=importlib.util.spec_from_file_location(name,p)
    out=importlib.util.module_from_spec(spec); spec.loader.exec_module(out); return out
def main():
    p=argparse.ArgumentParser(); p.add_argument('run',type=Path); a=p.parse_args(); run=a.run.resolve()
    plan=read(run/'FROZEN_EXECUTION_PLAN.json'); manifest=read(run/'SOURCE_MANIFEST.json')
    def unchanged():
        current={p.relative_to(run/'source').as_posix():hashlib.sha256(p.read_bytes()).hexdigest()
                 for p in (run/'source').rglob('*') if p.is_file() and '__pycache__' not in p.parts}
        return current==manifest
    assert unchanged(), 'Frozen source changed'
    audit=module('frozen_analysis',run/'source/scripts/analyze_validation.py')
    gate=module('frozen_gate',run/'source/scripts/evaluate_candidate.py')
    destination=run/'model_reports'; destination.mkdir(exist_ok=True); reports=[]
    for job in plan['jobs']:
        name=job['model']['name']; model=run/name; ex=model/'execution.json'
        if not ex.is_file() or read(ex)['status']!='COMPLETED': continue
        out=destination/name
        if (out/'SUMMARY.json').is_file():
            reports.append(read(out/'SUMMARY.json')); continue
        assert not out.exists(), 'Incomplete previous report: inspect before resuming'
        rr=rows(model/'matlab/metrics.csv')
        for r in rr: r['RowProvenance']='new_MATLAB_execution'
        replay_exact=None
        if job['reuse_original']:
            original=run.parents[2]/'p1_complete_20261021/runs/full_20261021T145453Z/metrics.csv'
            assert hashlib.sha256(original.read_bytes()).hexdigest()==plan['original_metrics_sha256']
            prior=[r for r in rows(original) if r['Model']==name and r['Method'] in ['Frozen','CombinedUncertainty']]
            index={tuple(r[k] for k in audit.KEY):r for r in prior}
            replay=rows(model/'matlab/replay_metrics.csv')
            replay_exact=bool(replay) and all(tuple(r[k] for k in audit.KEY) in index and
                all(r.get(k)==v for k,v in index[tuple(r[k] for k in audit.KEY)].items()) for r in replay)
            assert replay_exact, 'Original comparator replay differs'
            rr.extend(dict(r,CompensationApplied='1',ComputationAliasOf='',RowProvenance='completed_original_MATLAB_run') for r in prior)
        one=dict(plan,jobs=[job],expected_rows=3600)
        cov=audit.coverage(rr,one); assert cov['complete'], 'Model coverage incomplete'
        bad=[{k:r[k] for k in audit.KEY} for r in rr if audit.bad_solver(r,plan['acceptance']['normal_residual_max'])]
        out.mkdir(); fields=list(dict.fromkeys(k for r in rr for k in r))
        with (out/'metrics.csv').open('w',encoding='utf-8-sig',newline='') as f:
            w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(rr)
        result=gate.evaluate(out/'metrics.csv','CombinedUncertaintyV2','Frozen','Split',out/'vs_frozen',plan['acceptance'])
        suites=[]
        for suite in job['suites']:
            selected=[r for r in rr if r['Suite']==suite['name'] and r['Route']=='Split']
            means={m:statistics.mean(float(r['RelativeError']) for r in selected if r['Method']==m)
                   for m in ['Frozen','CombinedUncertainty','CombinedUncertaintyV2']}
            suites.append(dict(suite=suite['name'],means=means,
                gain_vs_original=1-means['CombinedUncertaintyV2']/means['CombinedUncertainty'],
                gain_vs_frozen=1-means['CombinedUncertaintyV2']/means['Frozen']))
        summary=dict(model=name,created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
            coverage=cov,all_method_solver_failures=bad,source_unchanged=unchanged(),
            original_replay_exact=replay_exact,suites=suites,
            single_model_frozen_gate_status=result['status'],failed_gates=[x for x in result['checks'] if not x['pass_']],
            scope='Descriptive completed-model report. Frozen gates vs Frozen applied to this model only; 5% Main gate is also a final cohort requirement. No final cohort acceptance or cross-model generalization claimed. No tuning on confirmation outcomes.')
        assert summary['source_unchanged']
        save(out/'SUMMARY.json',summary)
        lines=[f'# {name} 单模型结果','', '仅评价已完成模型；最终分组验收和跨模型统计仍待全部完成。', '',
            '| 条件 | 原版误差 | 新版误差 | 相对原版改善 | 相对 Frozen 改善 |', '|---|---:|---:|---:|---:|']
        for s in suites:
            lines.append(f"|{s['suite']}|{s['means']['CombinedUncertainty']:.6f}|{s['means']['CombinedUncertaintyV2']:.6f}|{100*s['gain_vs_original']:.2f}%|{100*s['gain_vs_frozen']:.2f}%|")
        lines+=['',f"单模型门槛：{result['status']}；所有方法求解失败行数：{len(bad)}。负改善表示变差。", '', '单模型未通过项目：']
        lines += [f"- {x['name']}：{json.dumps(x['detail'],ensure_ascii=False)}" for x in summary['failed_gates']] or ['- 无。']
        (out/'RESULTS_ZH.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
        reports.append(summary)
    save(destination/'INDEX.json',dict(reports=reports,final_acceptance=False))
    print(json.dumps([dict(model=r['model'],main=r['suites'][0],status=r['single_model_frozen_gate_status'],failed=r['failed_gates']) for r in reports],ensure_ascii=True,indent=2))
if __name__=='__main__': main()
