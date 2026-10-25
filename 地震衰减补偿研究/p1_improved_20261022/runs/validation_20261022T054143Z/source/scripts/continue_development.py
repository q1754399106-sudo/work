"""Serial local experiment handoff; never repeats failed MATLAB runs.

Keeps the goal moving without periodic AI invocations. Full validation is still
required after the candidate is selected; this script cannot mark the goal done.
"""
from pathlib import Path
import argparse,ctypes,datetime,hashlib,json,subprocess,sys,time
ROOT=Path(__file__).resolve().parents[1]
def read(p):return json.loads(p.read_text(encoding='utf-8'))
def write(p,v):
 t=p.with_suffix(p.suffix+'.tmp');t.write_text(json.dumps(v,ensure_ascii=False,indent=2),encoding='utf-8');t.replace(p)
def alive(pid):
 k=ctypes.windll.kernel32;k.OpenProcess.restype=ctypes.c_void_p
 h=k.OpenProcess(0x1000,False,int(pid))
 if not h:return False
 code=ctypes.c_ulong();ok=k.GetExitCodeProcess(ctypes.c_void_p(h),ctypes.byref(code));k.CloseHandle(ctypes.c_void_p(h))
 return bool(ok and code.value==259)
def main():
 p=argparse.ArgumentParser();p.add_argument('--prior-run',type=Path,required=True);p.add_argument('--layered-study',type=Path,required=True);a=p.parse_args()
 statefile=ROOT/'DEVELOPMENT_PIPELINE.json';start=time.monotonic();waits=[a.prior_run.resolve(),a.layered_study.resolve()]
 write(statefile,dict(status='WAITING_FOR_CURRENT_DEVELOPMENT',runs=[str(x) for x in waits],started_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()))
 while True:
  statuses=[read(x/'execution.json') for x in waits]
  if all(x['status']=='COMPLETED' for x in statuses):break
  failures=[dict(run=str(r),execution=s) for r,s in zip(waits,statuses) if s['status'] not in ['RUNNING','COMPLETED'] or (s['status']=='RUNNING' and not alive(s['pid']))]
  if failures:
   write(statefile,dict(status='DEVELOPMENT_FAILURE_REVIEW_REQUIRED',failures=failures,automatic_retry=False));return 1
  if time.monotonic()-start>10800:
   write(statefile,dict(status='WAIT_TIMEOUT_REVIEW_REQUIRED',workers_not_terminated=True));return 1
  time.sleep(20)
 # Apply the declared gates to the complete first screen and retain its failures.
 if not (a.prior_run/'analysis/development_selection.json').exists():
  with (a.prior_run/'analysis_stdout.log').open('wb') as out,(a.prior_run/'analysis_stderr.log').open('wb') as err:
   result=subprocess.run([sys.executable,str(ROOT/'scripts/analyze_development.py'),str(a.prior_run)],stdout=out,stderr=err)
  if result.returncode:
   write(statefile,dict(status='FIRST_SCREEN_ANALYSIS_FAILED',exit_code=result.returncode));return 1
 # Change the development model family, seeds, candidate configurations and gates
 # only as stated: here only the calibration model is extended; all are dev data.
 plan=read(ROOT/'core_p1/inputs/development_trial.json');plan['calibration_mode']='layered';plan['include_original_prior_candidate']=True
 write(ROOT/'core_p1/inputs/development_trial.json',plan)
 pointer=ROOT/'LAYERED_DEVELOPMENT_RUN.json'
 write(statefile,dict(status='LAYERED_DEVELOPMENT_RUNNING',earlier_screen=str(a.prior_run),layered_calibration_study=str(a.layered_study)))
 with (ROOT/'layered_development_runner_stdout.log').open('wb') as out,(ROOT/'layered_development_runner_stderr.log').open('wb') as err:
  result=subprocess.run([sys.executable,str(ROOT/'scripts/run_local.py'),'--entry','run_development_trial','--label','layered_development_screen','--timeout','9000','--run-pointer',str(pointer)],stdout=out,stderr=err)
 run=Path(read(pointer)['run']) if pointer.exists() else None
 if result.returncode or run is None:
  write(statefile,dict(status='LAYERED_DEVELOPMENT_FAILED',run=str(run),exit_code=result.returncode,automatic_retry=False));return 1
 with (run/'analysis_stdout.log').open('wb') as out,(run/'analysis_stderr.log').open('wb') as err:
  result=subprocess.run([sys.executable,str(ROOT/'scripts/analyze_development.py'),str(run)],stdout=out,stderr=err)
 if result.returncode:
  write(statefile,dict(status='LAYERED_SCREEN_ANALYSIS_FAILED',run=str(run),exit_code=result.returncode));return 1
 selection=read(run/'analysis/development_selection.json');chosen=selection['selected_profile']
 if chosen is None:
  write(statefile,dict(status='NO_FEASIBLE_CANDIDATE_CONTINUE_DEVELOPMENT',run=str(run),selection=str(run/'analysis/development_selection.json'),goal_complete=False));return 0
 policy=dict(version='CombinedUncertaintyV2-development-selected',calibration_mode='layered',selected_profile=chosen,
             selected_on_development_run=str(run),selection_rule=selection['selection_rule'],
             development_manifest_sha256=selection['source_manifest_sha256'],
             created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
             promotion_status='PENDING_INDEPENDENT_CONFIRMATION',frozen_before_reserved_validation=True)
 write(ROOT/'core_p1/inputs/candidate_policy.json',policy)
 write(statefile,dict(status='CANDIDATE_SELECTED_AWAITING_CONFIRMATION',run=str(run),policy=policy,goal_complete=False))
 return 0
if __name__=='__main__':raise SystemExit(main())
