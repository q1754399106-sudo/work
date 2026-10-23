"""One new dev screen; no calibration refitting or automatic experiment retry."""
from pathlib import Path
import argparse,csv,hashlib,json,shutil,subprocess,sys
ROOT=Path(__file__).resolve().parents[1]
def read(p):return json.loads(p.read_text(encoding='utf-8'))
def save(p,v):p.write_text(json.dumps(v,ensure_ascii=False,indent=2),encoding='utf-8')
def rows(p):
 with p.open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))
def main():
 parser=argparse.ArgumentParser();parser.add_argument('--resume-prelaunch',action='store_true');args=parser.parse_args()
 state=ROOT/'SHRINKAGE_DEVELOPMENT_PIPELINE.json'
 if state.exists():
  assert args.resume_prelaunch and read(state)['status']=='DRIVER_FAILED_REVIEW_REQUIRED' and not (ROOT/'SHRINKAGE_DEVELOPMENT_RUN.json').exists(),'Inspect previous numerical run; no automatic retry'
  archive=ROOT/'driver_attempts/shrinkage_prelaunch_1';archive.mkdir(parents=True)
  for name in ['SHRINKAGE_DEVELOPMENT_PIPELINE.json','shrinkage_controller_stdout.log','shrinkage_controller_stderr.log']:
   if (ROOT/name).exists():shutil.copy2(ROOT/name,archive/name)
  oldplan=ROOT/'core_p1/inputs/cached_shrinkage/trial_plan.json'
  if oldplan.exists():shutil.copy2(oldplan,archive/'trial_plan.json')
 source=Path(read(ROOT/'LAYERED_DEVELOPMENT_RUN.json')['run']);assert read(source/'execution.json')['status']=='COMPLETED'
 inputs=ROOT/'core_p1/inputs/cached_shrinkage';inputs.mkdir(exist_ok=args.resume_prelaunch);plan=read(source/'matlab/executed_plan.json');plan['profiles']=[]
 for key in ['calibration_seeds','test_seeds']:
  values=plan[key] if isinstance(plan[key],list) else [plan[key]];plan[key]=[int(x) for x in values]
 for weight in [.25,.5]:
  for light in [False,True]:
   plan['profiles'].append(dict(name=f'CalShrink{round(weight*100):03d}'+('Light' if light else 'Original'),confidence_floor=.2 if light else 1,
    section_modal_blend=.25 if light else 1,section_target=1 if light else -1,section_gamma=3 if light else 10,section_calibration_weight=weight))
 plan.update(calibration_mode='layered_section_shrinkage',include_original_prior_candidate=False,development_only=True,
  cached_calibrations_original_run=str(source),offset_candidates_alias='CombinedJointCalibration; candidate changes Section only')
 save(inputs/'trial_plan.json',plan);shutil.copy2(source/'matlab/metrics.csv',inputs/'original_development_metrics.csv')
 names=['Dev_ModerateQ','Dev_DeepContrast']+[m['name'] for m in plan['additional_models']]
 for name in names:shutil.copy2(source/'matlab'/f'calibration_{name}_Main_{plan["calibration_seeds"][0]}.json',inputs/f'calibration_{name}.json')
 manifest={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs.iterdir() if p.is_file()};save(inputs/'INPUT_PROVENANCE.json',dict(original_run=str(source),files=manifest))
 pointer=ROOT/'SHRINKAGE_DEVELOPMENT_RUN.json';save(state,dict(status='RUNNING',calibration_refitting=False,expected_new_section_rows=192,expected_total_rows=672))
 with (ROOT/'shrinkage_screen_runner_stdout.log').open('wb') as out,(ROOT/'shrinkage_screen_runner_stderr.log').open('wb') as err:
  done=subprocess.run([sys.executable,str(ROOT/'scripts/run_local.py'),'--entry','run_cached_shrinkage_trial','--label','shrinkage_development_screen','--timeout','3600','--run-pointer',str(pointer)],stdout=out,stderr=err)
 run=Path(read(pointer)['run']) if pointer.exists() else None
 if done.returncode or run is None:save(state,dict(status='FAILED_REVIEW_REQUIRED',run=str(run),automatic_retry=False));return 1
 frozen=run/'source/core_p1/inputs/cached_shrinkage';frozen_plan=read(frozen/'trial_plan.json');old=rows(frozen/'original_development_metrics.csv');fresh=rows(run/'matlab/metrics_section_only.csv');assert len(fresh)==192
 merged=[dict(r,RowProvenance='completed_development_run',ComputationAliasOf='',CalibrationBlend='') for r in old if r['Method'] in ['Frozen','CombinedUncertainty','CombinedJointCalibration']]+fresh
 for profile in frozen_plan['profiles']:
  for r in old:
   if r['Method']=='CombinedJointCalibration' and not r['Scenario'].endswith('Section'):
    merged.append(dict(r,Method=profile['name'],Preset=profile['name'],RowProvenance='exact_Offset_alias_Section_only_change',ComputationAliasOf='CombinedJointCalibration',CalibrationBlend='1'))
 key=['Model','Suite','Scenario','CalibrationSeed','Seed','Method','Route'];assert len(merged)==672 and len({tuple(r[k] for k in key) for r in merged})==672
 with (run/'matlab/metrics.csv').open('w',encoding='utf-8-sig',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(dict.fromkeys(k for r in merged for k in r)));w.writeheader();w.writerows(merged)
 with (run/'analysis_stdout.log').open('wb') as out,(run/'analysis_stderr.log').open('wb') as err:
  done=subprocess.run([sys.executable,str(ROOT/'scripts/analyze_development.py'),str(run)],stdout=out,stderr=err)
 if done.returncode:save(state,dict(status='ANALYSIS_FAILED_REVIEW_REQUIRED',run=str(run),automatic_retry=False));return 1
 selection=read(run/'analysis/development_selection.json')
 save(state,dict(status='DEVELOPMENT_FEASIBLE_REVIEW_API' if selection['selected_profile'] else 'NO_FEASIBLE_CANDIDATE_CONTINUE_DEVELOPMENT',run=str(run),selection=str(run/'analysis/development_selection.json'),selected_profile=selection['selected_profile'],goal_complete=False))
 return 0
if __name__=='__main__':
 try:code=main()
 except Exception as exc:
  save(ROOT/'SHRINKAGE_DEVELOPMENT_PIPELINE.json',dict(status='DRIVER_FAILED_REVIEW_REQUIRED',error=repr(exc),automatic_retry=False));raise
 raise SystemExit(code)
