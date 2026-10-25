"""Replay the predetermined prospective Section case through frozen API.

Uses the original saved calibration observation; no candidate tuning or retry.
"""
from pathlib import Path
import argparse,csv,datetime,hashlib,json,shutil,subprocess,time
DRIVER_DIR=Path(__file__).resolve().parent
KEY=['Model','Suite','Scenario','CalibrationSeed','Seed','Method','Route']
def hashes(folder):return {p.relative_to(folder).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(folder.rglob('*')) if p.is_file() and '__pycache__' not in p.parts}
def read(p):return json.loads(p.read_text(encoding='utf-8'))
def save(p,x):p.write_text(json.dumps(x,ensure_ascii=False,indent=2),encoding='utf-8')
def csvread(p):
 with p.open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))
def main():
 p=argparse.ArgumentParser();p.add_argument('run',type=Path);a=p.parse_args();run=a.run.resolve()
 assert read(run/'EXECUTION_SUMMARY.json')['status']=='COMPLETED','Wait for completed numerical confirmation'
 manifest=read(run/'SOURCE_MANIFEST.json');assert hashes(run/'source')==manifest,'Frozen source differs before replay'
 dest=run/'deterministic_replay';dest.mkdir();source=dest/'driver_source';source.mkdir()
 for name in ['replay_confirmation.py','replay_confirmation_case.m']:shutil.copy2(DRIVER_DIR/name,source/name)
 save(dest/'DRIVER_MANIFEST.json',hashes(source));q=lambda p:p.as_posix().replace("'","''")
 cmd=[read(run/'RUN_METADATA.json')['matlab_executable'],'-wait','-nojvm','-singleCompThread','-batch',f"addpath('{q(source)}');replay_confirmation_case('{q(run)}','{q(dest/'matlab')}');"]
 status=dict(status='RUNNING',command=cmd,started_utc=datetime.datetime.now(datetime.timezone.utc).isoformat());start=time.monotonic()
 with (dest/'stdout.log').open('wb') as out,(dest/'stderr.log').open('wb') as err:
  proc=subprocess.Popen(cmd,stdout=out,stderr=err,creationflags=getattr(subprocess,'CREATE_NO_WINDOW',0));status['pid']=proc.pid;save(dest/'execution.json',status)
  try:code=proc.wait(timeout=1800)
  except subprocess.TimeoutExpired:
   subprocess.run(['taskkill','/PID',str(proc.pid),'/T','/F'],capture_output=True);status.update(status='HARD_TIMEOUT',elapsed_seconds=time.monotonic()-start);save(dest/'execution.json',status);return 1
 status.update(status='COMPLETED' if code==0 else 'FAILED',exit_code=code,elapsed_seconds=time.monotonic()-start);save(dest/'execution.json',status)
 if code:return code
 original={tuple(r[k] for k in KEY):r for r in csvread(run/'metrics.csv')};replayed=csvread(dest/'matlab/metrics.csv');differences=[]
 for r in replayed:
  key=tuple(r[k] for k in KEY);old=original.get(key);changed=list(r) if old is None else [k for k in old if old[k]!=r.get(k)]
  if changed:differences.append(dict(key=key,fields=changed,values={k:dict(original=old.get(k) if old else None,replay=r.get(k)) for k in changed}))
 unchanged=hashes(run/'source')==manifest;passed=len(replayed)==2 and not differences and unchanged and read(dest/'matlab/checks.json')['passed']
 result=dict(passed=passed,rows=len(replayed),all_csv_fields_exact=not differences,differences=differences,frozen_source_unchanged=unchanged,
  scope='Predetermined first prospective model / first calibration and observation / Main WeakSection, Direct and Split production API. No claim that every confirmation case was replayed.')
 save(dest/'REPLAY_AUDIT.json',result);print(json.dumps(result,ensure_ascii=True,indent=2));return 0 if passed else 1
if __name__=='__main__':raise SystemExit(main())
