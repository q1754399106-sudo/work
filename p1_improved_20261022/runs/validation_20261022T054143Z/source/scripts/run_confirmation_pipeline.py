"""Frozen numerical confirmation followed by analysis; no retries or promotion.

Stops for actual figure inspection before packaging/final delivery.
"""
from pathlib import Path
import datetime,json,subprocess,sys
ROOT=Path(__file__).resolve().parents[1]
STATE=ROOT/'CONFIRMATION_PIPELINE.json'
def read(p):return json.loads(p.read_text(encoding='utf-8'))
def save(v):
 v['updated_utc']=datetime.datetime.now(datetime.timezone.utc).isoformat();temp=STATE.with_suffix('.tmp');temp.write_text(json.dumps(v,ensure_ascii=False,indent=2),encoding='utf-8');temp.replace(STATE)
def call(cmd,stdout,stderr):
 with stdout.open('wb') as out,stderr.open('wb') as err:return subprocess.run(cmd,cwd=ROOT,stdout=out,stderr=err).returncode
def main():
 assert not STATE.exists(),'Review existing confirmation rather than repeating it'
 pointer=ROOT/'FULL_CONFIRMATION_RUN.json';assert not pointer.exists()
 assert read(ROOT/'CONFIRMATION_READINESS.json')['status']=='READY_FOR_FROZEN_CONFIRMATION'
 save(dict(status='NUMERICAL_CONFIRMATION_RUNNING',workers=2,automatic_retry=False,goal_complete=False))
 code=call([sys.executable,str(ROOT/'scripts/run_validation.py'),'--workers','2','--run-pointer',str(pointer)],ROOT/'full_confirmation_runner_stdout.log',ROOT/'full_confirmation_runner_stderr.log')
 run=Path(read(pointer)['run']) if pointer.exists() else None
 if code or run is None:
  save(dict(status='NUMERICAL_DRIVER_REVIEW_REQUIRED',run=str(run),exit_code=code,automatic_retry=False));return 1
 summary=read(run/'EXECUTION_SUMMARY.json')
 if not(summary['coverage_complete'] and summary['baseline_replay_exact'] and not summary['source_changed_during_run']):
  save(dict(status='NUMERICAL_INTEGRITY_REVIEW_REQUIRED',run=str(run),summary=summary,automatic_retry=False));return 1
 save(dict(status='ANALYSIS_RUNNING',run=str(run),goal_complete=False))
 code=call([sys.executable,str(run/'source/scripts/analyze_validation.py'),str(run)],run/'analysis_stdout.log',run/'analysis_stderr.log')
 if code:
  save(dict(status='ANALYSIS_REVIEW_REQUIRED',run=str(run),exit_code=code,automatic_retry=False));return 1
 save(dict(status='FIGURES_RUNNING',run=str(run),method_acceptance=read(run/'analysis/acceptance.json')['status'],goal_complete=False))
 code=call([str(ROOT/'.figure_env/Scripts/python.exe'),str(run/'source/scripts/plot_validation.py'),str(run)],run/'figures_stdout.log',run/'figures_stderr.log')
 if code:
  save(dict(status='FIGURES_REVIEW_REQUIRED',run=str(run),exit_code=code,automatic_retry=False));return 1
 save(dict(status='READY_FOR_VISUAL_AND_COMPLETION_REVIEW',run=str(run),method_acceptance=read(run/'analysis/acceptance.json')['status'],figures_manual_review='PENDING',package_created=False,goal_complete=False))
 return 0
if __name__=='__main__':
 try:code=main()
 except Exception as exc:
  previous=read(STATE) if STATE.exists() else {};save(dict(status='PIPELINE_DRIVER_REVIEW_REQUIRED',error=repr(exc),previous=previous,automatic_retry=False));raise
 raise SystemExit(code)
