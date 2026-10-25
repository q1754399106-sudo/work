from pathlib import Path
import sys, json, hashlib, shutil, subprocess, time, datetime

ROOT=Path(__file__).resolve().parent

def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    stage=sys.argv[1];spec=json.loads((ROOT/'SPEC.json').read_text(encoding='utf-8'))
    assert stage in ['development','validation','sensitivity']
    cases=['Development'] if stage=='development' else (['QPoint080','QPoint125'] if stage=='sensitivity' else [v['name'] for v in spec['variants']])
    if stage in ['validation','sensitivity']:
        frozen=json.loads((ROOT/'FROZEN_CANDIDATE.json').read_text(encoding='utf-8'))
        assert all(digest(ROOT/p)==h for p,h in frozen['hashes'].items()),'Candidate changed after freeze'
    stamp=datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%SZ')
    run=ROOT/'runs'/f'{stage}_{stamp}';run.mkdir(parents=True)
    shutil.copytree(ROOT/'core',run/'source/core')
    (run/'source/SPEC.json').write_bytes((ROOT/'SPEC.json').read_bytes())
    (ROOT/f'{stage}_run.json').write_text(json.dumps({'run':str(run)},indent=2),encoding='utf-8')
    original=json.loads((ROOT/'FROZEN_BASELINE_HASHES.json').read_text(encoding='utf-8'))
    assert all(digest(run/'source/core'/p)==h for p,h in original.items())
    manifest={p.relative_to(run/'source').as_posix():digest(p) for p in (run/'source').rglob('*') if p.is_file()}
    (run/'SOURCE_MANIFEST.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
    q=lambda p:str(p).replace('\\','/').replace("'","''")
    statement=f"addpath('{q(run/'source/core')}');"+''.join(f"run_qcluster_case('{q(run/case)}','{q(ROOT/'inputs'/f'{case}.json')}');" for case in cases)
    cmd=['D:/Program Files (x86)/MATLAB2024B/bin/matlab.exe','-wait','-nojvm','-singleCompThread','-batch',statement]
    start=time.monotonic();status=dict(stage=stage,status='STARTING',cases=cases,started_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),command=cmd,timeout_seconds=5400)
    print('RUN',str(run),flush=True)
    with (run/'stdout.log').open('wb') as out,(run/'stderr.log').open('wb') as err:
        proc=subprocess.Popen(cmd,stdout=out,stderr=err,creationflags=getattr(subprocess,'CREATE_NO_WINDOW',0))
        status.update(status='RUNNING',pid=proc.pid)
        (run/'execution.json').write_text(json.dumps(status,indent=2),encoding='utf-8')
        try:
            code=proc.wait(timeout=5400);status.update(exit_code=code,status='COMPLETED' if code==0 else 'FAILED')
        except subprocess.TimeoutExpired:
            print('Hard timeout: stopping owned process tree',proc.pid,flush=True)
            subprocess.run(['taskkill','/PID',str(proc.pid),'/T','/F'],capture_output=True)
            status.update(status='HARD_TIMEOUT',exit_code=None)
    status['elapsed_seconds']=time.monotonic()-start
    (run/'execution.json').write_text(json.dumps(status,indent=2),encoding='utf-8')
    print(json.dumps(status,indent=2),flush=True)
    return 0 if status['status']=='COMPLETED' else 1
if __name__=='__main__':raise SystemExit(main())
