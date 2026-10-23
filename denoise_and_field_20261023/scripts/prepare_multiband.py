from pathlib import Path
import json,shutil,datetime
import numpy as np
from multiband_filter import suppress_multiband,MULTIBAND_PARAMETERS
ROOT=Path(__file__).resolve().parents[1];OLD=ROOT.parent/'complex_fieldlike_20261023';archive=ROOT/'attempt_v1_broadband';archive.mkdir(exist_ok=True)
for rel in ['PLAN.json','prepare.log','analysis/preprocess_diagnostics.json','analysis/clean_input_controls.json','analysis/Moderate_preprocess.npz','analysis/Strong_preprocess.npz','analysis/Moderate_clean_control.npz','analysis/Strong_clean_control.npz','core_p1/inputs/denoised_input.json']:
 p=ROOT/rel;target=archive/rel
 assert p.resolve().is_relative_to(ROOT.resolve()) and target.resolve().is_relative_to(ROOT.resolve())
 if p.exists() and not target.exists():target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,target)
shutil.copy2(ROOT/'scripts/coherent_filter.py',archive/'coherent_filter.py')
plan=json.loads((archive/'PLAN.json').read_text(encoding='utf-8'));plan.update(attempt='v2 band-limited detection, preconditioned, same base thresholds',registered_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),parameters=MULTIBAND_PARAMETERS,development_note='Same synthetic observations are development; broadband version detected nothing. Both attempts retained. No held-out generalization claim; field records processed after freezing this rule.')
(ROOT/'PLAN.json').write_text(json.dumps(plan,ensure_ascii=False,indent=2),encoding='utf-8')
payload=json.loads((OLD/'core_p1/inputs/complex_input.json').read_text(encoding='utf-8'));diagnostics={};x=np.arange(99)*30.
for data in payload['datasets']:
 name=data['name'];D=np.asarray(data['D']);out,removed,diag=suppress_multiband(D,payload['public']['dt'],x);data['D']=out.tolist();diagnostics[name]=diag;np.savez_compressed(ROOT/'analysis'/(name+'_preprocess.npz'),original=D,denoised=out,removed=removed);print(name,'removed fraction',diag['removed_norm_fraction'],'components',[len(b['diagnostics'].get('components',[])) for b in diag['bands']],flush=True)
(ROOT/'core_p1/inputs/denoised_input.json').write_text(json.dumps(payload,separators=(',',':')),encoding='utf-8');(ROOT/'analysis/preprocess_diagnostics.json').write_text(json.dumps(diagnostics,ensure_ascii=False,indent=2),encoding='utf-8')
truth=np.load(OLD/'truth/model_truth.npz');controls={}
for name in ['Moderate','Strong']:
 target=truth[name+'_target'];filtered,removed,diag=suppress_multiband(target,payload['public']['dt'],x);controls[name]=dict(relative_change=float(np.linalg.norm(filtered-target)/np.linalg.norm(target)),diagnostics=diag);np.savez_compressed(ROOT/'analysis'/(name+'_clean_control.npz'),truth=target,filtered=filtered,removed=removed);print('Clean control',name,controls[name]['relative_change'],flush=True)
(ROOT/'analysis/clean_input_controls.json').write_text(json.dumps(controls,ensure_ascii=False,indent=2),encoding='utf-8')
