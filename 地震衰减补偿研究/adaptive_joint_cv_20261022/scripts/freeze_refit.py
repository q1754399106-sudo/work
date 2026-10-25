from pathlib import Path
import json,hashlib,datetime
ROOT=Path(__file__).resolve().parents[1]
summary=json.loads((ROOT/'analysis/refit_development_summary.json').read_text(encoding='utf-8'))
assert summary['gates']['FullCalRefit']['all_pass']
core=ROOT/'core_p1';run=Path(summary['run']);old=run/'source/core_p1'
files=['reflectivity_compensate.m','matlab_src/refit_dictionary.m','matlab_src/refit_lasso.m','matlab_src/reflectivity_refit_prepared.m','matlab_src/p1_prepare_candidate.m','matlab_src/p1_joint_calibrate.m','matlab_src/p1_layered_calibrate.m','matlab_src/p1_noise_uncertainty.m','inputs/candidate_policy.json']
hashes={p:hashlib.sha256((core/p).read_bytes()).hexdigest() for p in files}
assert all((old/p).exists() and hashlib.sha256((old/p).read_bytes()).hexdigest()==v for p,v in hashes.items())
result={'method':'FullCalRefit','status':'FROZEN_BEFORE_VALIDATION','frozen_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'development_run':str(run),'development_cases':60,'all_development_gates_pass':True,'development_summary_sha256':hashlib.sha256((ROOT/'analysis/refit_development_summary.json').read_bytes()).hexdigest(),'algorithm_sha256':hashes,'selection':'One profile for every model/domain. Lowest mean waveform error among the two predefined refit profiles passing every gate. No validation-based retuning.'}
target=core/'inputs/frozen_refit_profile.json'
assert not target.exists()
target.write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps(result,ensure_ascii=False,indent=2))
