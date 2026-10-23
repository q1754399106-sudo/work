from pathlib import Path
import json,datetime,hashlib,shutil
import numpy as np
from coherent_filter import suppress_coherent,PARAMETERS
ROOT=Path(__file__).resolve().parents[1];OLD=ROOT.parent/'complex_fieldlike_20261023'
for p in [ROOT/'analysis',ROOT/'inputs',ROOT/'figures']:p.mkdir(exist_ok=True)
if not (ROOT/'PLAN.json').exists():
 plan=dict(registered_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),method='Data-only long-aperture linear coherence packet estimation and conservative subtraction',parameters=PARAMETERS,source_model='Previous frozen complex field-like 99 trace model',conditions=['Moderate','Strong'],truth_use='Scoring only, not input to filter or compensated solvers',gates=dict(clean_input_relative_change_max=.02,per_condition_method_waveform_ratio_max=1.02,low_energy_error_ratio_max=1.05),promotion='Research attempt; retain failures, do not change default if any gate fails',real_data='USGS Alaska 31-81 records 111,120,130 first 96 traces; Yilmaz wz25 entire gather. No reference truth or trusted Q calibration; denoising-only field experiment, no invented Q or recovery-error claims.')
 (ROOT/'PLAN.json').write_text(json.dumps(plan,ensure_ascii=False,indent=2),encoding='utf-8')
shutil.copytree(OLD/'core_p1',ROOT/'core_p1',dirs_exist_ok=True);shutil.copy2(OLD/'scripts/run_local.py',ROOT/'scripts/run_local.py')
hashes=json.loads((OLD/'FROZEN_SOLVER_HASHES.json').read_text(encoding='utf-8'))
assert all(hashlib.sha256((ROOT/'core_p1'/p).read_bytes()).hexdigest()==h for p,h in hashes.items())
(ROOT/'FROZEN_SOLVER_HASHES.json').write_text(json.dumps(hashes,indent=2),encoding='utf-8')
run=Path(json.loads((OLD/'run.json').read_text(encoding='utf-8'))['run']);shutil.copy2(run/'matlab/calibration.json',ROOT/'core_p1/inputs/shared_calibration.json')
payload=json.loads((OLD/'core_p1/inputs/complex_input.json').read_text(encoding='utf-8'));originals={d['name']:np.asarray(d['D']) for d in payload['datasets']};diagnostics={};x=np.arange(99)*30.
for data in payload['datasets']:
 name=data['name'];D=originals[name];out,removed,diag=suppress_coherent(D,payload['public']['dt'],x);data['D']=out.tolist();diagnostics[name]=diag
 np.savez_compressed(ROOT/'analysis'/(name+'_preprocess.npz'),original=D,denoised=out,removed=removed)
 print(name,json.dumps(diag,ensure_ascii=True),flush=True)
(ROOT/'core_p1/inputs/denoised_input.json').write_text(json.dumps(payload,separators=(',',':')),encoding='utf-8')
(ROOT/'analysis/preprocess_diagnostics.json').write_text(json.dumps(diagnostics,ensure_ascii=False,indent=2),encoding='utf-8')
# Clean-input control is evaluated only after all noisy outputs have been frozen.
truth=np.load(OLD/'truth/model_truth.npz');controls={}
for name in originals:
 target=truth[name+'_target'];filtered,removed,diag=suppress_coherent(target,payload['public']['dt'],x)
 controls[name]=dict(relative_change=float(np.linalg.norm(filtered-target)/np.linalg.norm(target)),diagnostics=diag);np.savez_compressed(ROOT/'analysis'/(name+'_clean_control.npz'),truth=target,filtered=filtered,removed=removed)
(ROOT/'analysis/clean_input_controls.json').write_text(json.dumps(controls,ensure_ascii=False,indent=2),encoding='utf-8')
print('Frozen preprocessing done; clean-control results saved separately.',flush=True)
