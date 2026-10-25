from pathlib import Path
import sys,json,hashlib,shutil,subprocess,datetime,time
ROOT=Path(__file__).resolve().parent
digest=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 stage=sys.argv[1];assert stage in ['development','validation','sensitivity'];spec=json.loads((ROOT/'SPEC.json').read_text(encoding='utf-8'))
 kappa=1.;cases=['Development']
 if stage in ['validation','sensitivity']:
  frozen=json.loads((ROOT/'FROZEN_CANDIDATE.json').read_text(encoding='utf-8'));assert all(digest(ROOT/p)==h for p,h in frozen['hashes'].items())
  kappa=json.loads((ROOT/'SELECTED.json').read_text(encoding='utf-8'))['kappa'];cases=['QPoint080','QPoint125'] if stage=='sensitivity' else [v['name'] for v in spec['variants']]
 stamp=datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%SZ');run=ROOT/'runs'/f'{stage}_{stamp}';run.mkdir(parents=True)
 shutil.copytree(ROOT/'core',run/'source/core');shutil.copy2(ROOT/'SPEC.json',run/'source/SPEC.json')
 for case in cases:shutil.copy2(ROOT/'inputs'/f'{case}.json',run/'source'/f'{case}.json')
 original=json.loads((ROOT/'BASELINE_HASHES.json').read_text(encoding='utf-8'));assert all(digest(run/'source/core'/p)==h for p,h in original.items())
 manifest={p.relative_to(run/'source').as_posix():digest(p) for p in (run/'source').rglob('*') if p.is_file()};(run/'SOURCE_MANIFEST.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
 (ROOT/f'{stage}_run.json').write_text(json.dumps(dict(run=str(run)),indent=2),encoding='utf-8')
 q=lambda p:str(p).replace('\\','/').replace("'","''")
 expr=f"addpath('{q(run/'source/core')}');"+''.join(f"run_output_trust_case('{q(run/case)}','{q(run/'source'/f'{case}.json')}','{stage}',{kappa});" for case in cases)
 cmd=['D:/Program Files (x86)/MATLAB2024B/bin/matlab.exe','-wait','-nojvm','-singleCompThread','-batch',expr]
 start=time.monotonic();status=dict(stage=stage,cases=cases,kappa=kappa,status='STARTING',started_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),command=cmd,timeout_seconds=5400)
 print('RUN',str(run),flush=True)
 with (run/'stdout.log').open('wb') as out,(run/'stderr.log').open('wb') as err:
  proc=subprocess.Popen(cmd,stdout=out,stderr=err,creationflags=getattr(subprocess,'CREATE_NO_WINDOW',0));status.update(pid=proc.pid,status='RUNNING');(run/'execution.json').write_text(json.dumps(status,indent=2),encoding='utf-8')
  try:
   code=proc.wait(timeout=5400);status.update(exit_code=code,status='COMPLETED' if code==0 else 'FAILED')
  except subprocess.TimeoutExpired:
   print('Hard timeout: stopping owned process tree',proc.pid,flush=True);subprocess.run(['taskkill','/PID',str(proc.pid),'/T','/F'],capture_output=True);status.update(exit_code=None,status='HARD_TIMEOUT')
 status['elapsed_seconds']=time.monotonic()-start;(run/'execution.json').write_text(json.dumps(status,indent=2),encoding='utf-8');print(json.dumps(status,indent=2),flush=True)
 return 0 if status['status']=='COMPLETED' else 1
if __name__=='__main__':raise SystemExit(main())
