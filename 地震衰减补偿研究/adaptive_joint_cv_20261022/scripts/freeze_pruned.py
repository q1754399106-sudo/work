from pathlib import Path
import csv,json,hashlib,datetime
from pruning_gates import summarize,gates
ROOT=Path(__file__).resolve().parents[1];run=Path(json.loads((ROOT/'pruning_run.json').read_text(encoding='utf-8'))['run'])
checks=json.loads((run/'matlab/checks.json').read_text(encoding='utf-8'));assert checks['complete'] and checks['cases']==110
rows=list(csv.DictReader((run/'matlab/metrics.csv').open(encoding='utf-8-sig')));main=summarize([r for r in rows if r['Suite']=='normal']);stress=summarize([r for r in rows if r['Suite']=='colored'])
assert main['cases']==100 and stress['cases']==10
g=json.loads((ROOT/'core_p1/inputs/adaptive_cv_plan.json').read_text(encoding='utf-8'))['gates'];main['gates']=gates(main,g,checks)
stress['passes']=all(stress['means']['PrunedFullCalRefit'][k]<=stress['means']['DecoupledWeak'][k] for k in ['RelativeError','WeakEventError']) and checks['optimization_flags']==0
result={'main':main,'colored':stress,'checks':checks,'run':str(run),'scope':'Expanded development after round-one validation failure; none of these 110 cases is now an unseen test.'}
(ROOT/'analysis/pruned_development_summary.json').write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf-8')
assert main['gates']['all_pass'] and stress['passes'],'Pruned candidate did not pass development gates'
files=['pruned_compensate.m','matlab_src/prune_dictionary.m','matlab_src/reflectivity_prune_prepared.m','matlab_src/refit_dictionary.m','matlab_src/refit_lasso.m','matlab_src/reflectivity_refit_prepared.m','matlab_src/p1_prepare_candidate.m','matlab_src/p1_noise_uncertainty.m','inputs/candidate_policy.json']
hashes={p:hashlib.sha256((ROOT/'core_p1'/p).read_bytes()).hexdigest() for p in files}
for p,v in hashes.items():
    if p=='pruned_compensate.m':continue
    assert hashlib.sha256((run/'source/core_p1'/p).read_bytes()).hexdigest()==v
frozen={'method':'PrunedFullCalRefit','status':'FROZEN_BEFORE_SECOND_ROUND_VALIDATION','frozen_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'development_run':str(run),'algorithm_sha256':hashes,'development_cases':110,'all_development_gates_pass':True,'validation_plan_sha256':hashlib.sha256((ROOT/'core_p1/inputs/pruned_validation_plan.json').read_bytes()).hexdigest()}
path=ROOT/'core_p1/inputs/frozen_pruned_profile.json';assert not path.exists();path.write_text(json.dumps(frozen,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps({'main':main,'colored':stress,'frozen':True},ensure_ascii=False,indent=2))
