"""Wait for selected dev policy, then verify API and development-only harness.

No reserved simulations or retries; retained logs diagnose any failed stage.
"""
from pathlib import Path
import argparse,ctypes,datetime,hashlib,json,subprocess,sys,time
from run_validation import csvread,KEY,ORIGINAL
ROOT=Path(__file__).resolve().parents[1]
def read(p):return json.loads(p.read_text(encoding='utf-8'))
def write(p,v):
 t=p.with_suffix(p.suffix+'.tmp');t.write_text(json.dumps(v,ensure_ascii=False,indent=2),encoding='utf-8');t.replace(p)
def alive(pid):
 k=ctypes.windll.kernel32;k.OpenProcess.restype=ctypes.c_void_p;h=k.OpenProcess(0x1000,False,pid)
 if not h:return False
 code=ctypes.c_ulong();ok=k.GetExitCodeProcess(ctypes.c_void_p(h),ctypes.byref(code));k.CloseHandle(ctypes.c_void_p(h));return bool(ok and code.value==259)
def main():
 p=argparse.ArgumentParser();p.add_argument('--policy-ready',action='store_true');a=p.parse_args()
 state=ROOT/'CONFIRMATION_READINESS.json'
 if state.exists() and a.policy_ready:
  prior=read(state);assert prior['status']=='DEVELOPMENT_REVIEW_REQUIRED' and not prior.get('stages'),'Review existing numerical verification before retrying'
  archive=ROOT/'driver_attempts'/('readiness_before_policy_'+datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%SZ')+'.json');archive.parent.mkdir(exist_ok=True);write(archive,prior)
 write(state,dict(status='CHECKING_SELECTED_POLICY' if a.policy_ready else 'WAITING_FOR_DEVELOPMENT',goal_complete=False));start=time.monotonic()
 while not a.policy_ready:
  dev=read(ROOT/'DEVELOPMENT_PIPELINE.json');status=dev['status']
  if status=='CANDIDATE_SELECTED_AWAITING_CONFIRMATION':break
  if status!='LAYERED_DEVELOPMENT_RUNNING' or not alive(38956):
   write(state,dict(status='DEVELOPMENT_REVIEW_REQUIRED',development=dev,automatic_retry=False));return 1
  if time.monotonic()-start>10800:
   write(state,dict(status='OBSERVATION_TIMEOUT_REVIEW_REQUIRED',processes_not_stopped=True));return 1
  time.sleep(20)
 policy=ROOT/'core_p1/inputs/candidate_policy.json';policydata=read(policy);assert policydata['promotion_status']=='PENDING_INDEPENDENT_CONFIRMATION'
 selected=read(Path(policydata['selected_on_development_run'])/'analysis/development_selection.json')
 assert selected['selected_profile']==policydata['selected_profile'] and selected['results'][selected['selected_profile']['name']]['feasible_for_validation'],'Policy lacks feasible completed development evidence'
 policyhash=hashlib.sha256(policy.read_bytes()).hexdigest();stages={}
 for entry,label in [('verify_v2','production_verification'),('run_validation_smoke','validation_harness_smoke')]:
  pointer=ROOT/(label.upper()+'_RUN.json')
  if pointer.exists():
   write(state,dict(status='EXISTING_VERIFICATION_REQUIRES_REVIEW',pointer=str(pointer),automatic_retry=False));return 1
  write(state,dict(status=label.upper()+'_RUNNING',stages=stages,policy_sha256=policyhash))
  with (ROOT/(label+'_runner_stdout.log')).open('wb') as out,(ROOT/(label+'_runner_stderr.log')).open('wb') as err:
   done=subprocess.run([sys.executable,str(ROOT/'scripts/run_local.py'),'--entry',entry,'--label',label,'--timeout','3600','--run-pointer',str(pointer)],stdout=out,stderr=err)
  run=Path(read(pointer)['run']) if pointer.exists() else None
  if done.returncode or run is None:
   write(state,dict(status='VERIFICATION_FAILED',stage=label,run=str(run),exit_code=done.returncode,automatic_retry=False));return 1
  checks=read(run/'matlab/checks.json');assert checks['passed'];stages[label]=dict(run=str(run),checks=checks)
  manifest=read(run/'SOURCE_MANIFEST.json');assert all(hashlib.sha256((run/'source/core_p1'/k).read_bytes()).hexdigest()==h for k,h in manifest.items())
 smoke=Path(stages['validation_harness_smoke']['run']);old={tuple(r[k] for k in KEY):r for r in csvread(ORIGINAL/'metrics.csv')};differences=[]
 replay=csvread(smoke/'matlab/job/replay_metrics.csv')
 for r in replay:
  key=tuple(r[k] for k in KEY);b=old.get(key);changed=list(r) if b is None else [k for k in b if r.get(k)!=b[k]]
  if changed:differences.append(dict(key=key,changed_fields=changed,values={k:dict(old=b.get(k) if b else None,replay=r.get(k)) for k in changed}))
 audit=dict(rows=len(replay),exact_match=not differences,differences=differences);write(smoke/'ORIGINAL_REPLAY_AUDIT.json',audit)
 unchanged=policyhash==hashlib.sha256(policy.read_bytes()).hexdigest()
 passed=not differences and len(replay)==12 and unchanged
 write(state,dict(status='READY_FOR_FROZEN_CONFIRMATION' if passed else 'REPLAY_OR_POLICY_REVIEW_REQUIRED',stages=stages,replay=audit,policy_sha256=policyhash,policy_unchanged=unchanged,full_confirmation_started=False,goal_complete=False))
 return 0 if passed else 1
if __name__=='__main__':
 try:code=main()
 except Exception as exc:
  state=ROOT/'CONFIRMATION_READINESS.json';prior=read(state) if state.exists() else {}
  write(state,dict(status='READINESS_DRIVER_FAILED',error=repr(exc),previous=prior,automatic_retry=False))
  raise
 raise SystemExit(code)
