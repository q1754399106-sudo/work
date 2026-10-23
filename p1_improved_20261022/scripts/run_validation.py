"""Freeze the development-selected candidate and run reserved + regression jobs.

Numerical baselines from the completed original nine-model run are reused only
for exactly matching cases and independently replayed on every scene/condition.
New conditions and new models execute both original controls again.
"""
from pathlib import Path
import argparse,concurrent.futures,csv,datetime,hashlib,json,platform,shutil,subprocess,sys,time
ROOT=Path(__file__).resolve().parents[1]
ORIGINAL=ROOT.parent/'p1_complete_20261021/runs/full_20261021T145453Z'
MATLAB='D:/Program Files (x86)/MATLAB2024B/bin/matlab.exe'
METHODS=['Frozen','CombinedUncertainty','FrozenJointCalibration','CombinedJointCalibration','CombinedUncertaintyV2']
KEY=['Model','Suite','Scenario','CalibrationSeed','Seed','Method','Route']
def read(p):return json.loads(p.read_text(encoding='utf-8'))
def save(p,x):
 t=p.with_suffix(p.suffix+'.tmp');t.write_text(json.dumps(x,ensure_ascii=False,indent=2),encoding='utf-8');t.replace(p)
def hashes(folder):return {p.relative_to(folder).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(folder.rglob('*')) if p.is_file() and '__pycache__' not in p.parts}
def csvread(p):
 with p.open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))
def csvwrite(p,rows):
 fields=list(dict.fromkeys(k for r in rows for k in r))
 with p.open('w',encoding='utf-8-sig',newline='') as f:
  w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(rows)
def jobs():
 reserved=read(ROOT/'RESERVED_VALIDATION_DESIGN.json');expected=(ROOT/'RESERVED_VALIDATION_DESIGN.sha256').read_text(encoding='ascii').strip()
 assert hashlib.sha256((ROOT/'RESERVED_VALIDATION_DESIGN.json').read_bytes()).hexdigest()==expected,'Reserved models changed'
 original=read(ROOT/'core_p1/inputs/p1_experiment_spec.json');main=original['scenarios'];avo=original['stress_scenarios']
 common=dict(vp_multiplier=1.,h_multiplier=1.,phase_deg=0.,max_nmo_shift_s=0.,ar1_rho=0.)
 def suite(name,kind,scenes,snr=18,**kw):return dict(common,name=name,kind=kind,scenarios=scenes,snr_db=snr,**kw)
 suites=[suite('Main','normal',main),suite('LowerSNR','normal',main,12),suite('StressNonRickerPhase','phase',avo,phase_deg=30.),
   suite('StressVelocityError','velocity',avo,vp_multiplier=1.03),suite('StressResidualNMO','nmo',avo,max_nmo_shift_s=.004),
   suite('StressColoredNoise','colored',avo,ar1_rho=.55),suite('StressVelocityMinus3','velocity',avo,vp_multiplier=.97)]
 result=[dict(model=m,group='prospective_synthetic',calibration_seeds=reserved['calibration_seeds'],test_seeds=reserved['observation_seeds'],suites=suites,reuse_original=False) for m in reserved['models']]
 for group,models in original['models'].items():
  for model in models:result.append(dict(model=model,group=group,calibration_seeds=original['calibration_seeds'],test_seeds=original['test_seeds'],suites=suites,reuse_original=True))
 return result
def execute(job,run,timeout):
 d=run/job['model']['name'];d.mkdir();save(d/'job.json',job);q=lambda p:p.as_posix().replace("'","''")
 cmd=[MATLAB,'-wait','-nojvm','-singleCompThread','-batch',f"addpath('{q(run/'source/core_p1')}');run_validation_job('{q(d/'job.json')}','{q(d/'matlab')}');"]
 status=dict(status='RUNNING',model=job['model']['name'],command=cmd,started_utc=datetime.datetime.now(datetime.timezone.utc).isoformat());start=time.monotonic()
 try:
  with (d/'stdout.log').open('wb') as out,(d/'stderr.log').open('wb') as err:
   p=subprocess.Popen(cmd,stdout=out,stderr=err,creationflags=getattr(subprocess,'CREATE_NO_WINDOW',0));status['pid']=p.pid;save(d/'execution.json',status)
   print(f"STARTED {job['model']['name']} pid={p.pid}",flush=True)
   try:
    status['exit_code']=p.wait(timeout=timeout);status['status']='COMPLETED' if p.returncode==0 else 'FAILED'
   except subprocess.TimeoutExpired:
    print(f'Hard timeout: terminating only owned batch tree {p.pid}',flush=True);subprocess.run(['taskkill','/PID',str(p.pid),'/T','/F'],capture_output=True);status['status']='HARD_TIMEOUT'
 except Exception as exc:status.update(status='FAILED',error=repr(exc))
 status['elapsed_seconds']=time.monotonic()-start;save(d/'execution.json',status)
 print(f"FINISHED {job['model']['name']} {status['status']} {status['elapsed_seconds']:.1f}s",flush=True);return status
def merge(run,plan):
 assert hashlib.sha256((ORIGINAL/'metrics.csv').read_bytes()).hexdigest()==plan['original_metrics_sha256'],'Original comparator data changed after freeze'
 prior=[r for r in csvread(ORIGINAL/'metrics.csv') if r['Method'] in ['Frozen','CombinedUncertainty']]
 index={tuple(r[k] for k in KEY):r for r in prior};rows=[];replay_differences=[];replay_count=0;fresh_count=0
 for job in plan['jobs']:
  d=run/job['model']['name'];new=csvread(d/'matlab/metrics.csv');rows.extend(new);fresh_count+=len(new)
  if job['reuse_original']:
   imported=[dict(r,CompensationApplied='1',ComputationAliasOf='',RowProvenance='completed_original_MATLAB_run') for r in prior if r['Model']==job['model']['name']]
   rows.extend(imported)
   for r in csvread(d/'matlab/replay_metrics.csv'):
    old=index.get(tuple(r[k] for k in KEY));replay_count+=1
    changed=list(r) if old is None else [k for k in old if r.get(k)!=old[k]]
    if changed:replay_differences.append(dict(key=tuple(r[k] for k in KEY),changed_fields=changed,values={k:dict(old=old.get(k) if old else None,replay=r.get(k)) for k in changed}))
  for r in new:r['RowProvenance']='new_MATLAB_execution'
 csvwrite(run/'metrics.csv',rows)
 save(run/'ORIGINAL_REPLAY_AUDIT.json',dict(rows=replay_count,exact_match=not replay_differences,differences=replay_differences,original_run=str(ORIGINAL),original_csv_sha256=hashlib.sha256((ORIGINAL/'metrics.csv').read_bytes()).hexdigest()))
 return dict(rows=len(rows),new_numerical_rows=fresh_count,imported_baseline_rows=len(rows)-fresh_count,replay_rows=replay_count,baseline_replay_exact=not replay_differences)
def main():
 p=argparse.ArgumentParser();p.add_argument('--workers',type=int,default=2);p.add_argument('--timeout',type=float,default=43200);p.add_argument('--run-pointer',type=Path);a=p.parse_args()
 policy=read(ROOT/'core_p1/inputs/candidate_policy.json');assert policy['promotion_status']=='PENDING_INDEPENDENT_CONFIRMATION'
 readiness=read(ROOT/'CONFIRMATION_READINESS.json');assert readiness['status']=='READY_FOR_FROZEN_CONFIRMATION','Production and harness verification must pass first'
 assert readiness['policy_sha256']==hashlib.sha256((ROOT/'core_p1/inputs/candidate_policy.json').read_bytes()).hexdigest(),'Policy changed after verification'
 verified=Path(readiness['stages']['validation_harness_smoke']['run']);verified_manifest=read(verified/'SOURCE_MANIFEST.json')
 assert hashes(ROOT/'core_p1')==verified_manifest,'Candidate code changed after numerical harness verification'
 stamp=datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%SZ');run=ROOT/'runs'/('validation_'+stamp);run.mkdir(parents=True)
 for folder in ['core_p1','scripts']:shutil.copytree(ROOT/folder,run/'source'/folder,ignore=shutil.ignore_patterns('__pycache__'))
 for filename in ['NEXT_EXPERIMENT.json','RESERVED_VALIDATION_DESIGN.json','RESERVED_VALIDATION_DESIGN.sha256']:
  shutil.copy2(ROOT/filename,run/'source'/filename)
 manifest=hashes(run/'source');save(run/'SOURCE_MANIFEST.json',manifest)
 jj=jobs();plan=dict(stage='frozen_candidate_confirmation_and_regression',policy=policy,jobs=jj,methods=METHODS,
   created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),acceptance=read(ROOT/'NEXT_EXPERIMENT.json')['acceptance'],reserved_design_sha256=(ROOT/'RESERVED_VALIDATION_DESIGN.sha256').read_text(encoding='ascii').strip(),
   original_metrics_sha256=hashlib.sha256((ORIGINAL/'metrics.csv').read_bytes()).hexdigest(),production_and_harness_readiness=readiness,
   expected_observations=sum(len(j['calibration_seeds'])*len(j['test_seeds'])*sum(len(s['scenarios']) for s in j['suites']) for j in jj),
   original_nine_models_are_regression=True,prospective_models_are_same_simulator_not_external_field_data=True,
   original_comparators_reused_only_at_identical_inputs_with_exact_scene_replay=True)
 plan['expected_rows']=plan['expected_observations']*len(METHODS)*2;save(run/'FROZEN_EXECUTION_PLAN.json',plan)
 save(run/'RUN_METADATA.json',dict(python=sys.version,platform=platform.platform(),matlab_executable=MATLAB,
   workers=a.workers,per_model_timeout_seconds=a.timeout,created_utc=plan['created_utc'],source_manifest_sha256=hashlib.sha256((run/'SOURCE_MANIFEST.json').read_bytes()).hexdigest(),
   policy_sha256=hashlib.sha256((run/'source/core_p1/inputs/candidate_policy.json').read_bytes()).hexdigest(),
   original_run=str(ORIGINAL),reserved_models_frozen_before_selection=True,original_models_are_regression=True,
   actual_matlab_versions_recorded_in_per_model_checks=True))
 shutil.copy2(ROOT/'RESERVED_VALIDATION_DESIGN.json',run/'RESERVED_VALIDATION_DESIGN.json')
 if a.run_pointer:save(a.run_pointer,dict(run=str(run)))
 print(f'RUN_DIRECTORY {run}',flush=True);print(f"EXPECTED observations={plan['expected_observations']} rows={plan['expected_rows']}",flush=True)
 with concurrent.futures.ThreadPoolExecutor(max_workers=a.workers) as pool:
  fs=[pool.submit(execute,j,run,a.timeout) for j in jj];statuses=[f.result() for f in concurrent.futures.as_completed(fs)]
 after=hashes(run/'source');changes=[k for k in set(manifest)|set(after) if manifest.get(k)!=after.get(k)]
 good=all(s['status']=='COMPLETED' for s in statuses);summary=dict(status='COMPLETED' if good else 'FAILED',tasks=statuses,source_changed_during_run=changes)
 if good:summary.update(merge(run,plan));summary['expected_rows']=plan['expected_rows'];summary['coverage_complete']=summary['rows']==plan['expected_rows']
 save(run/'EXECUTION_SUMMARY.json',summary);print(json.dumps({k:v for k,v in summary.items() if k!='tasks'},ensure_ascii=True),flush=True)
 return 0 if good and not changes else 1
if __name__=='__main__':raise SystemExit(main())
