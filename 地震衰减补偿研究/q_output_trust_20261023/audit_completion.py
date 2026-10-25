from pathlib import Path
import json,hashlib,csv,datetime
import numpy as np
ROOT=Path(__file__).resolve().parent;digest=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
spec=json.loads((ROOT/'SPEC.json').read_text());selected=json.loads((ROOT/'SELECTED.json').read_text());frozen=json.loads((ROOT/'FROZEN_CANDIDATE.json').read_text());old=json.loads((ROOT/'BASELINE_HASHES.json').read_text());checks={}
checks['all_frozen_sources_inputs_targets_unchanged']=all(digest(ROOT/p)==h for p,h in frozen['hashes'].items())
checks['122_prior_MATLAB_dependencies_unchanged']=all(digest(ROOT/'core'/p)==h for p,h in old.items())
checks['standalone_controls']=json.loads((ROOT/'analysis/standalone_controls.json').read_text())['passed']
checks['independent_conditional_covariance_check']=json.loads((ROOT/'analysis/covariance_dense_symmetric.json').read_text())['passed']
checks['prepared_cache_replays_selected_estimator']=json.loads((ROOT/'analysis/cache_replay_check.json').read_text())['passed']
stages=dict(development=['Development'],validation=[v['name'] for v in spec['variants']],sensitivity=['QPoint080','QPoint125']);total=0;waveforms=0;records={};guardvalues=[];caches=[]
for stage,cases in stages.items():
 run=Path(json.loads((ROOT/f'{stage}_run.json').read_text())['run']);ex=json.loads((run/'execution.json').read_text());mf=json.loads((run/'SOURCE_MANIFEST.json').read_text());records[stage]=dict(run=str(run),elapsed_seconds=ex['elapsed_seconds'],cases={})
 checks[f'{stage}_completed_exit_zero']=ex['status']=='COMPLETED' and ex['exit_code']==0
 checks[f'{stage}_executed_source_and_input_hashes']=all(digest(run/'source'/p)==h for p,h in mf.items())
 checks[f'{stage}_same_core_estimator']=digest(run/'source/core/qtrust_refit.m')==digest(ROOT/'core/qtrust_refit.m')
 if stage=='validation':checks['frozen_before_new_model_execution']=datetime.datetime.fromisoformat(frozen['frozen_utc'])<datetime.datetime.fromisoformat(ex['started_utc'])
 methods=['FullCalRefit','PrunedFullCalRefit','QAwareRidge','QAwareCluster','Trust05','Trust10','Trust20'] if stage=='development' else spec['methods']
 for case in cases:
  ck=json.loads((run/case/'checks.json').read_text());records[stage]['cases'][case]=ck
  checks[f'{case}_numerical_and_constraints']=ck['complete'] and ck['optimization_flags']==0 and ck['max_normal']<1e-8 and ck['max_guard_ratio']<=1+1e-6 and ck['controls']['passed']
  inp=json.loads((run/'source'/f'{case}.json').read_text());checks[f'{case}_truth_fields_excluded']=not set(inp['public']).intersection({'Q','fpeak','truth','rho','vs','coeff','event_t'})
  total+=ck['traces']
  for cfg in inp['datasets']:
   data=json.loads((run/case/(cfg['name']+'.json')).read_text())
   cache=run/case/(cfg['name']+'_prepared.mat');assert cache.is_file()
   caches.append(dict(stage=stage,case=case,condition=cfg['name'],path=cache.relative_to(ROOT).as_posix(),bytes=cache.stat().st_size,sha256=digest(cache),source_manifest_sha256=digest(run/'SOURCE_MANIFEST.json'),scope='Prepared dictionary, calibration uncertainty and observed-data support; no target arrays'))
   for m in methods:
    a=np.array(data[m]);assert a.shape==(inp['public']['nt'],22) and np.isfinite(a).all();waveforms+=a.shape[1]
   for m,di in data['diagnostics'].items():
    for item in di:
     path=np.atleast_1d(item['objective_path']);assert np.all(np.diff(path)<=1e-9*max(path[0],1e-30))
     if 'kappa' in item:assert item['max_protected_change_sd']<=item['kappa']+1e-7;guardvalues.append(item['max_protected_change_sd']/item['kappa'])
checks['all_220_trace_conditions']=total==220;checks['all_1364_method_waveforms']=waveforms==1364
checks['all_cumulative_output_guards_satisfied']=max(guardvalues)<=1+1e-6
for name,n in [('development_metrics.csv',16),('validation_metrics.csv',42),('sensitivity_metrics.csv',18),('previous_round_output_diagnosis.csv',8)]:
 with (ROOT/'analysis'/name).open(encoding='utf-8-sig',newline='') as f:r=list(csv.DictReader(f))
 checks[f'{name}_complete']=len(r)==n
base=json.loads((ROOT/'inputs'/f'{spec["variants"][0]["name"]}.json').read_text())
for case,factor in [('QPoint080',.8),('QPoint125',1.25)]:
 p=json.loads((ROOT/'inputs'/f'{case}.json').read_text());checks[f'{case}_same_observation_registered_scale']=p['q_point_scale']==factor and len(p['datasets'])==1 and np.array_equal(np.array(p['datasets'][0]['D']),np.array(base['datasets'][0]['D']))
seed=lambda s:{z for v in s['variants'] for z in [v['seed'],v['calibration_seed']]+v['observation_seeds']}
prior_path=ROOT/'vendor/round1_SPEC.json';prior=json.loads(prior_path.read_text());prior_provenance=json.loads((ROOT/'vendor/round1_SPEC_PROVENANCE.json').read_text())
checks['round1_seed_source_matches_its_frozen_hash']=digest(prior_path)==prior_provenance['expected_sha256']
checks['new_seeds_disjoint_from_round1']=not seed(spec).intersection(seed(prior))
for name in ['RESULTS_ZH.md','README_ZH.md','index.html','analysis/validation_summary.json','analysis/sensitivity_summary.json','figures/manifest.json']:
 checks[f'artifact_{name}']=(ROOT/name).is_file() and (ROOT/name).stat().st_size>0
fm=json.loads((ROOT/'figures/manifest.json').read_text());checks['figure_hashes_match']=all(digest(ROOT/p['path'])==p['sha256'] for p in fm['plots'])
result=dict(completed=all(checks.values()),checked_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),checks=checks,trace_conditions=total,method_waveforms=waveforms,selected_kappa=selected['kappa'],execution=records,note='Experimental completion is distinct from scientific novelty and performance promotion.')
(ROOT/'CACHE_MANIFEST.json').write_text(json.dumps(dict(caches=caches,note='Reusable prepared-data caches for future development. Reuse requires matching case, observation, point-Q setting and solver source; these are not independent data.'),indent=2))
(ROOT/'COMPLETION_AUDIT.json').write_text(json.dumps(result,indent=2));print(json.dumps(dict(completed=result['completed'],trace_conditions=total,method_waveforms=waveforms,failed=[k for k,v in checks.items() if not v]),indent=2));assert result['completed']
