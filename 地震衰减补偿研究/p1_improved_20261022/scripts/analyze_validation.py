"""Audit frozen confirmation coverage, numerical quality and crossed uncertainty.

No candidate selection or tuning is performed here. Failed gates are retained.
"""
from pathlib import Path
import argparse,collections,csv,hashlib,importlib.util,itertools,json,math,random,statistics

KEY=['Model','Suite','Scenario','CalibrationSeed','Seed','Method','Route']
def read(p):return json.loads(p.read_text(encoding='utf-8'))
def save(p,v):p.write_text(json.dumps(v,ensure_ascii=False,indent=2),encoding='utf-8')
def csvread(p):
 with p.open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))
def csvwrite(p,rows):
 with p.open('w',encoding='utf-8-sig',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def coverage(rows,plan):
 expected=set()
 for j in plan['jobs']:
  for s in j['suites']:
   for sc,c,o,m,r in itertools.product(s['scenarios'],j['calibration_seeds'],j['test_seeds'],plan['methods'],['Direct','Split']):
    expected.add((j['model']['name'],s['name'],sc,str(c),str(o),m,r))
 actual=collections.Counter(tuple(r[k] for k in KEY) for r in rows)
 missing=expected-set(actual);extra=set(actual)-expected;dups=[k for k,n in actual.items() if n!=1]
 return dict(complete=not(missing or extra or dups) and len(rows)==plan['expected_rows'],
  expected=len(expected),actual=len(rows),missing=sorted(missing),extra=sorted(extra),duplicates=sorted(dups))
def bad_solver(r,limit):
 try:
  values=[float(r[k]) for k in ['FinalFlag','NormalResidual','AllSelectionSolverFlags','AllSelectionMaxNormalResidual']]
  return not(all(math.isfinite(x) for x in values) and values[0]==values[2]==0 and 0<=values[1]<=limit and 0<=values[3]<=limit)
 except (ValueError,KeyError):return True
def quantile(a,p):
 a=sorted(a);x=(len(a)-1)*p;i=int(x);return a[i]+(a[min(i+1,len(a)-1)]-a[i])*(x-i)
def crossed_bootstrap(rows,jobs,baseline,draws=3000):
 """Resample model, calibration and observation axes; keep methods paired.

 Calibration and observation draws are shared within each cohort, matching
 common random numbers. Separate calls give independent cohort intervals.
 Scenes are a fixed stratum, not a fourth source of independent replication.
 """
 idx={tuple(r[k] for k in KEY):float(r['RelativeError']) for r in rows if r['Suite']=='Main' and r['Route']=='Split'}
 models=[j['model']['name'] for j in jobs];cals=[str(x) for x in jobs[0]['calibration_seeds']];seeds=[str(x) for x in jobs[0]['test_seeds']]
 assert all([str(x) for x in j['calibration_seeds']]==cals and [str(x) for x in j['test_seeds']]==seeds for j in jobs)
 scenes=next(s['scenarios'] for s in jobs[0]['suites'] if s['name']=='Main');means={}
 for m,c,o,method in itertools.product(models,cals,seeds,[baseline,'CombinedUncertaintyV2']):
  means[m,c,o,method]=statistics.mean(idx[m,'Main',s,c,o,method,'Split'] for s in scenes)
 def gain(mm,cc,oo):
  b=sum(means[m,c,o,baseline] for m,c,o in itertools.product(mm,cc,oo));v=sum(means[m,c,o,'CombinedUncertaintyV2'] for m,c,o in itertools.product(mm,cc,oo))
  return 1-v/b
 rng=random.Random(20261022);values=[]
 for _ in range(draws):values.append(gain(rng.choices(models,k=len(models)),rng.choices(cals,k=len(cals)),rng.choices(seeds,k=len(seeds))))
 return dict(baseline=baseline,candidate='CombinedUncertaintyV2',relative_improvement=gain(models,cals,seeds),
  crossed_percentile_95_interval=[quantile(values,.025),quantile(values,.975)],draws=draws,models=len(models),calibrations=len(cals),observation_seeds=len(seeds),
  scope='Model/calibration/observation crossed bootstrap; paired common draws within cohort; fixed Main scenarios; synthetic simulator only')
def main():
 p=argparse.ArgumentParser();p.add_argument('run',type=Path);p.add_argument('--draws',type=int,default=3000);a=p.parse_args();run=a.run.resolve()
 plan=read(run/'FROZEN_EXECUTION_PLAN.json');ex=read(run/'EXECUTION_SUMMARY.json');assert ex['status']=='COMPLETED','Numerical execution incomplete'
 dest=run/'analysis';dest.mkdir(exist_ok=False);rows=csvread(run/'metrics.csv');cov=coverage(rows,plan);save(dest/'coverage.json',cov)
 if not cov['complete']:raise SystemExit('Incomplete/duplicate coverage; see analysis/coverage.json')
 manifest=read(run/'SOURCE_MANIFEST.json');changed=[k for k,h in manifest.items() if not (run/'source'/k).is_file() or hashlib.sha256((run/'source'/k).read_bytes()).hexdigest()!=h]
 extra=[p.relative_to(run/'source').as_posix() for p in (run/'source').rglob('*') if p.is_file() and '__pycache__' not in p.parts and p.relative_to(run/'source').as_posix() not in manifest]
 replay=read(run/'ORIGINAL_REPLAY_AUDIT.json');quality={}
 for method in plan['methods']:
  rr=[r for r in rows if r['Method']==method];bad=[{k:r[k] for k in KEY} for r in rr if bad_solver(r,plan['acceptance']['normal_residual_max'])]
  quality[method]=dict(rows=len(rr),passed=not bad,failures=bad)
 integrity=dict(source_changed=changed,source_extra=extra,original_baseline_replay=replay,all_methods_solver_quality=quality,
  compensation_applied_all=all(float(r.get('CompensationApplied','0'))==1 for r in rows))
 integrity['passed']=not(changed or extra) and replay['exact_match'] and all(x['passed'] for x in quality.values()) and integrity['compensation_applied_all'];save(dest/'integrity.json',integrity)
 spec=importlib.util.spec_from_file_location('gate',run/'source/scripts/evaluate_candidate.py');gate=importlib.util.module_from_spec(spec);spec.loader.exec_module(gate)
 cohorts={'prospective_synthetic':[j for j in plan['jobs'] if not j['reuse_original']], 'original_nine_regression':[j for j in plan['jobs'] if j['reuse_original']]};results={};intervals={}
 for cohort,jj in cohorts.items():
  names={j['model']['name'] for j in jj};rr=[r for r in rows if r['Model'] in names];path=dest/(cohort+'.csv');csvwrite(path,rr);results[cohort]={}
  for method in ['CombinedUncertainty','FrozenJointCalibration','CombinedJointCalibration','CombinedUncertaintyV2']:
   g=gate.evaluate(path,method,'Frozen','Split',dest/cohort/method,plan['acceptance']);results[cohort][method]=dict(status=g['status'],gain_vs_frozen=g['equal_model_main_improvement'],failed=[x for x in g['checks'] if not x['pass_']])
  intervals[cohort]={b:crossed_bootstrap(rr,jj,b,a.draws) for b in ['Frozen','CombinedUncertainty','FrozenJointCalibration','CombinedJointCalibration']}
 save(dest/'crossed_bootstrap.json',intervals)
 status='PASS' if integrity['passed'] and all(results[c]['CombinedUncertaintyV2']['status']=='PROVISIONAL_PASS' and intervals[c]['CombinedUncertainty']['crossed_percentile_95_interval'][0]>0 for c in cohorts) else 'FAIL'
 summary=dict(status=status,coverage=cov,integrity_passed=integrity['passed'],cohorts=results,bootstrap=intervals,
  scope='Original nine models are regression. Four prospective synthetic models were frozen before selection and simulated after policy freeze; same simulator, not independent field validation.',
  unverified='Layerwise Q uncertainty coverage and arbitrary layerwise velocity/geometry recovery are not established.')
 save(dest/'acceptance.json',summary)
 lines=['# CombinedUncertainty 改进版确认结果','',f'综合判定：**{status}**。覆盖 {len(rows)} 行；完整覆盖：{cov["complete"]}；源码、重放及所有方法求解质量：{integrity["passed"]}。','',
  '原九模型属于回归试验。四个提前封存的新合成模型属于前瞻确认，仍使用相同模拟器，不能称为外部实测验证。','',
  '| 数据组 | 对照 | 主试验相对改善 | 交叉 bootstrap 95% 区间 |','|---|---|---:|---:|']
 for c,bs in intervals.items():
  for b,z in bs.items():
   lo,hi=z['crossed_percentile_95_interval'];lines.append(f'|{c}|{b}|{100*z["relative_improvement"]:.2f}%|[{100*lo:.2f}%, {100*hi:.2f}%]|')
 lines+=['','工程门槛失败项（保留局部误差、各模型及压力条件的失败）：']
 for c,z in results.items():
  for x in z['CombinedUncertaintyV2']['failed']:lines.append(f'- {c} / {x["name"]}：{json.dumps(x["detail"],ensure_ascii=False)}')
 lines+=['','完整验收与对照消融见 acceptance.json，各模型局部指标见分组目录。分层 Q 区间覆盖率、任意分层速度和几何反演尚未得到验证。']
 (dest/'RESULTS_ZH.md').write_text('\n'.join(lines)+'\n',encoding='utf-8');print(json.dumps(dict(status=status,rows=len(rows),integrity=integrity['passed'])))
 return 0
if __name__=='__main__':raise SystemExit(main())
