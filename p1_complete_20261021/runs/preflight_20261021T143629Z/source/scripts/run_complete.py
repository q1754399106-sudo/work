"""Local MATLAB experiment runner. Only standard library; immutable per-run source snapshot."""
from pathlib import Path
import argparse, concurrent.futures, csv, datetime, hashlib, json, os, shutil, subprocess, time

ROOT=Path(__file__).resolve().parents[1]
MATLAB=Path('D:/Program Files (x86)/MATLAB2024B/bin/matlab.exe')

def write_json(path,obj):
    tmp=path.with_suffix(path.suffix+'.tmp')
    tmp.write_text(json.dumps(obj,ensure_ascii=False,indent=2),encoding='utf-8')
    tmp.replace(path)

def hashes():
    return {p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest()
            for folder in ['core_p1','scripts'] for p in sorted((ROOT/folder).rglob('*'))
            if p.is_file() and '__pycache__' not in p.parts}

def make_jobs(preflight):
    spec=json.loads((ROOT/'core_p1/inputs/p1_experiment_spec.json').read_text())
    scenarios=spec['scenarios']; avo=spec['stress_scenarios']
    common=dict(vp_multiplier=1.,h_multiplier=1.,phase_deg=0.,max_nmo_shift_s=0.,ar1_rho=0.)
    def suite(name,kind,scenes,snr=18,**kw):
        return dict(common,name=name,kind=kind,scenarios=scenes,snr_db=snr,**kw)
    suites=[suite('Main','normal',scenarios),suite('LowerSNR','normal',scenarios,12),
            suite('StressNonRickerPhase','phase',avo,phase_deg=30.),
            suite('StressVelocityError','velocity',avo,vp_multiplier=1.03),
            suite('StressResidualNMO','nmo',avo,max_nmo_shift_s=.004),
            suite('StressColoredNoise','colored',avo,ar1_rho=.55)]
    methods=[m['name'] for m in spec['methods']]
    jobs=[]
    for group,models in spec['models'].items():
        for model in models:
            jobs.append(dict(model=model,group=group,calibration_seeds=spec['calibration_seeds'],
                             test_seeds=spec['test_seeds'],methods=methods,suites=suites,preflight=preflight))
    if preflight:
        jobs=jobs[:1]
        jobs[0].update(calibration_seeds=spec['calibration_seeds'][:1],test_seeds=spec['test_seeds'][:1],
                       suites=[suites[0],suites[3],suites[4]])
    return jobs

def execute(job,run,timeout):
    name=job['model']['name']; task=run/name; task.mkdir()
    jobfile=task/'job.json'; write_json(jobfile,job)
    quote=lambda p:p.as_posix().replace("'","''")
    command=[str(MATLAB),'-wait','-singleCompThread','-batch',
             f"addpath('{quote(ROOT/'core_p1')}'); run_complete_job('{quote(jobfile)}','{quote(task/'matlab')}');"]
    status=dict(status='RUNNING',model=name,command=command,started_utc=datetime.datetime.now(datetime.timezone.utc).isoformat())
    start=time.monotonic()
    try:
        with (task/'stdout.log').open('wb') as out,(task/'stderr.log').open('wb') as err:
            proc=subprocess.Popen(command,stdout=out,stderr=err,creationflags=getattr(subprocess,'CREATE_NO_WINDOW',0))
            status['pid']=proc.pid; write_json(task/'execution.json',status)
            print(f'STARTED {name} pid={proc.pid}',flush=True)
            try:
                status['exit_code']=proc.wait(timeout=timeout)
                status['status']='COMPLETED' if proc.returncode==0 else 'FAILED'
            except subprocess.TimeoutExpired:
                print(f'HARD TIMEOUT {name}; terminating only owned process tree {proc.pid}',flush=True)
                subprocess.run(['taskkill','/PID',str(proc.pid),'/T','/F'],capture_output=True)
                status['status']='HARD_TIMEOUT'
    except Exception as exc:
        status.update(status='FAILED',error=repr(exc))
    status['elapsed_seconds']=time.monotonic()-start
    write_json(task/'execution.json',status)
    print(f"FINISHED {name} {status['status']} {status['elapsed_seconds']:.1f}s",flush=True)
    return status

def main():
    parser=argparse.ArgumentParser(); parser.add_argument('--stage',choices=['preflight','full'],required=True)
    parser.add_argument('--workers',type=int,default=2); parser.add_argument('--timeout',type=float,default=43200)
    a=parser.parse_args(); jobs=make_jobs(a.stage=='preflight')
    stamp=datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%SZ')
    run=ROOT/'runs'/(a.stage+'_'+stamp); run.mkdir(parents=True)
    before=hashes(); write_json(run/'SOURCE_MANIFEST.json',before)
    plan=dict(stage=a.stage,jobs=jobs,source_hashes=before,
              created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
              prior_results_seen=True,hyperparameters_changed=False,
              original_models_are_regression_not_blind_holdout=True,
              velocity_stress='inverse vp * 1.03, thickness fixed, recalibrate Q/f and rebuild all operators',
              lower_snr='12dB all five scenes; full calibration x observation cross product',
              engine='new isolated hidden batch processes; original interactive MATLAB untouched',
              expected_observations=sum(len(j['calibration_seeds'])*len(j['test_seeds'])*sum(len(s['scenarios']) for s in j['suites']) for j in jobs))
    plan['expected_rows']=plan['expected_observations']*len(jobs[0]['methods'])*2
    write_json(run/'FROZEN_EXECUTION_PLAN.json',plan)
    # Persist the actual source used, independent of later working-directory edits.
    shutil.copytree(ROOT/'core_p1',run/'source'/'core_p1')
    shutil.copytree(ROOT/'scripts',run/'source'/'scripts',ignore=shutil.ignore_patterns('__pycache__'))
    print(f'RUN_DIRECTORY {run}',flush=True)
    print(f"EXPECTED observations={plan['expected_observations']} rows={plan['expected_rows']}",flush=True)
    with concurrent.futures.ThreadPoolExecutor(max_workers=a.workers) as pool:
        futures=[pool.submit(execute,j,run,a.timeout) for j in jobs]
        statuses=[f.result() for f in concurrent.futures.as_completed(futures)]
    after=hashes(); changed=[k for k in set(before)|set(after) if before.get(k)!=after.get(k)]
    good=all(s['status']=='COMPLETED' for s in statuses)
    result=dict(status='COMPLETED' if good else 'FAILED',tasks=statuses,source_changed_during_run=changed)
    if good:
        files=[run/j['model']['name']/'matlab/metrics.csv' for j in jobs]
        rows=[]
        for path in files:
            with path.open(encoding='utf-8-sig',newline='') as f: rows.extend(csv.DictReader(f))
        with (run/'metrics.csv').open('w',encoding='utf-8-sig',newline='') as f:
            writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
        result.update(rows=len(rows),expected_rows=plan['expected_rows'],coverage_complete=len(rows)==plan['expected_rows'])
    write_json(run/'EXECUTION_SUMMARY.json',result)
    print(json.dumps({k:v for k,v in result.items() if k!='tasks'},ensure_ascii=False),flush=True)
    return 0 if good and not changed else 1

if __name__=='__main__': raise SystemExit(main())
