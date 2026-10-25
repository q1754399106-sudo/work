"""Assess the fixed development candidates; never imports holdout results."""
from pathlib import Path
import argparse,csv,hashlib,importlib.util,json,statistics
ROOT=Path(__file__).resolve().parents[1]
def save(path,obj):path.write_text(json.dumps(obj,ensure_ascii=False,indent=2),encoding='utf-8')
def main():
 p=argparse.ArgumentParser();p.add_argument('run',type=Path);a=p.parse_args();run=a.run.resolve()
 execution=json.loads((run/'execution.json').read_text(encoding='utf-8'))
 if execution['status']!='COMPLETED':raise SystemExit('Numerical experiment is not complete')
 plan=json.loads((run/'matlab/executed_plan.json').read_text(encoding='utf-8'))
 dest=run/'analysis';dest.mkdir(exist_ok=False)
 spec=importlib.util.spec_from_file_location('gate',ROOT/'scripts/evaluate_candidate.py');gate=importlib.util.module_from_spec(spec);spec.loader.exec_module(gate)
 with (run/'matlab/metrics.csv').open(encoding='utf-8-sig',newline='') as f:rows=list(csv.DictReader(f))
 methods=['CombinedUncertainty','CombinedJointCalibration']+[x['name'] for x in plan['profiles']]
 result={}
 for method in methods:
  data=gate.evaluate(run/'matlab/metrics.csv',method,'Frozen','Split',dest/method)
  re=[r for r in rows if r['Method']==method]
  bad=[r for r in re if float(r['FinalFlag'])!=0 or float(r['NormalResidual'])>2e-7
       or float(r['AllSelectionSolverFlags'])!=0 or float(r['AllSelectionMaxNormalResidual'])>2e-7]
  result[method]=dict(status=data['status'],passed=sum(x['pass_'] for x in data['checks']),total=len(data['checks']),
   failed=[x for x in data['checks'] if not x['pass_']],gain_vs_frozen=data['equal_model_main_improvement'],solver_failures=len(bad))
 old=result['CombinedUncertainty']['gain_vs_frozen'];feasible=[]
 candidates=list(plan['profiles'])
 if plan.get('include_original_prior_candidate'):
  candidates.append(dict(name='CombinedJointCalibration',confidence_floor=1,section_modal_blend=1,section_target=-1,section_gamma=10))
 for profile in candidates:
  x=result[profile['name']];gain=(1-(1-x['gain_vs_frozen'])/(1-old));x['gain_vs_original_combined']=gain
  # Do not trade a primary error regression for a new label or a pass-through.
  x['feasible_for_validation']=x['status']=='PROVISIONAL_PASS' and x['solver_failures']==0 and gain>=0
  if x['feasible_for_validation']:feasible.append((x['gain_vs_frozen'],profile))
 selected=max(feasible,key=lambda x:x[0])[1] if feasible else None
 summary=dict(development_only=True,rows=len(rows),results=result,selected_profile=selected,
              selection_rule='All original engineering gates and solver quality; Main mean must improve vs original Combined; among feasible profiles minimize Main equal-model error',
              source_manifest_sha256=hashlib.sha256((run/'SOURCE_MANIFEST.json').read_bytes()).hexdigest())
 save(dest/'development_selection.json',summary)
 print(json.dumps(summary,ensure_ascii=True,indent=2))
 return 0
if __name__=='__main__':raise SystemExit(main())
