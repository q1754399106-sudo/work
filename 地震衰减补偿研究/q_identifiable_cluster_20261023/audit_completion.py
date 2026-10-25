from pathlib import Path
import json,hashlib,datetime,csv
import numpy as np

ROOT=Path(__file__).resolve().parent
digest=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
spec=json.loads((ROOT/'SPEC.json').read_text(encoding='utf-8'))
frozen=json.loads((ROOT/'FROZEN_CANDIDATE.json').read_text(encoding='utf-8'))
checks={}
checks['frozen_sources_inputs_and_truth_unchanged']=all(digest(ROOT/p)==h for p,h in frozen['hashes'].items())
old=json.loads((ROOT/'FROZEN_BASELINE_HASHES.json').read_text(encoding='utf-8'))
checks['original_baseline_MATLAB_unchanged']=all(digest(ROOT/'core'/p)==h for p,h in old.items())
records={};expected=dict(development=['Development'],validation=[v['name'] for v in spec['variants']],sensitivity=['QPoint080','QPoint125'])
trace_conditions=0;method_traces=0
for stage,cases in expected.items():
    run=Path(json.loads((ROOT/f'{stage}_run.json').read_text(encoding='utf-8'))['run'])
    ex=json.loads((run/'execution.json').read_text(encoding='utf-8'))
    checks[f'{stage}_exit_zero']=ex['status']=='COMPLETED' and ex['exit_code']==0
    if stage=='validation':checks['frozen_before_validation_run']=datetime.datetime.fromisoformat(frozen['frozen_utc'])<datetime.datetime.fromisoformat(ex['started_utc'])
    mf=json.loads((run/'SOURCE_MANIFEST.json').read_text(encoding='utf-8'))
    checks[f'{stage}_exact_executed_source_snapshot']=all(digest(run/'source'/p)==h for p,h in mf.items())
    records[stage]=dict(run=str(run),elapsed_seconds=ex['elapsed_seconds'],cases={})
    for case in cases:
        ck=json.loads((run/case/'checks.json').read_text(encoding='utf-8'))
        checks[f'{case}_numerical']=ck['complete'] and ck['optimization_flags']==0 and ck['max_normal']<1e-8 and ck['zero_Q_covariance_limit']<1e-10 and ck['zero_input_exact']
        if stage!='development':checks[f'{case}_limiting_controls']=ck['zero_Q_ridge_refit_limit']<1e-6 and ck['synthetic_controls']['duplicate_rank_after']==1
        inp=json.loads((ROOT/'inputs'/f'{case}.json').read_text(encoding='utf-8'))
        checks[f'{case}_public_input_truth_isolation']=not set(inp['public']).intersection({'Q','fpeak','truth','rho','vs','coeff','event_t'})
        trace_conditions+=ck['traces']
        for cfg in inp['datasets']:
            r=json.loads((run/case/(cfg['name']+'.json')).read_text(encoding='utf-8'))
            assert np.array(r['Observed']).shape==(inp['public']['nt'],22)
            for method in spec['methods']:
                a=np.array(r[method]);assert a.shape==(inp['public']['nt'],22) and np.isfinite(a).all()
                method_traces+=a.shape[1]
        records[stage]['cases'][case]=ck
nominal_input=json.loads((ROOT/'inputs/FreshFoldLowQ.json').read_text(encoding='utf-8'))
for case,factor in [('QPoint080',.8),('QPoint125',1.25)]:
    pert=json.loads((ROOT/'inputs'/f'{case}.json').read_text(encoding='utf-8'))
    checks[f'{case}_same_observation_as_nominal']=np.array_equal(np.array(pert['datasets'][0]['D']),np.array(nominal_input['datasets'][0]['D'])) and len(pert['datasets'])==1
    checks[f'{case}_registered_scale']=pert['q_point_scale']==factor
seed_values=[v[k] for v in spec['variants'] for k in ['seed','calibration_seed']]+[s for v in spec['variants'] for s in v['observation_seeds']]
checks['new_model_and_noise_seeds_distinct']=len(seed_values)==len(set(seed_values))
checks['all_220_trace_conditions']=trace_conditions==220
checks['all_1320_method_waveforms']=method_traces==1320
for name,n in [('development_metrics.csv',14),('validation_metrics.csv',42),('sensitivity_metrics.csv',18)]:
    with (ROOT/'analysis'/name).open(encoding='utf-8-sig',newline='') as f:rows=list(csv.DictReader(f))
    checks[f'{name}_complete']=len(rows)==n
for name in ['RESULTS_ZH.md','METHOD_ZH.md','README_ZH.md','index.html','analysis/validation_summary.json','analysis/sensitivity_summary.json','figures/manifest.json']:
    checks[f'artifact_{name}']=(ROOT/name).is_file() and (ROOT/name).stat().st_size>0
fm=json.loads((ROOT/'figures/manifest.json').read_text(encoding='utf-8'))
checks['all_rendered_figure_hashes']=all(digest(ROOT/p['path'])==p['sha256'] for p in fm['plots'])
result=dict(completed=all(checks.values()),checked_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),checks=checks,trace_conditions=trace_conditions,method_waveforms=method_traces,execution=records,
    interpretation='Experiment completion is distinct from candidate promotion; inspect validation gates for performance. Hash checks prove code/data consistency, not scientific novelty or physical validity.')
(ROOT/'COMPLETION_AUDIT.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
print(json.dumps(dict(completed=result['completed'],trace_conditions=trace_conditions,method_waveforms=method_traces,failed=[k for k,v in checks.items() if not v]),indent=2))
assert result['completed']
