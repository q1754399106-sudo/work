from pathlib import Path
import json,csv,math,hashlib,datetime
ROOT=Path(__file__).resolve().parents[1];old=ROOT/'runs/cv_pilot_20261022T155015Z';new=Path(json.loads((ROOT/'pilot_run.json').read_text(encoding='utf-8'))['run'])
prior=list(csv.DictReader((old/'matlab/metrics.csv').open(encoding='utf-8-sig')));resumed=list(csv.DictReader((new/'matlab/metrics.csv').open(encoding='utf-8-sig')));assert len(prior)==12 and len(resumed)==18
for a,b in zip(prior,resumed):
    assert a.keys()==b.keys()
    for k in a:
        if a[k]==b[k]:continue
        x=float(a[k]);y=float(b[k]);assert (math.isnan(x) and math.isnan(y)) or math.isclose(x,y,rel_tol=1e-13,abs_tol=1e-15)
files=list((old/'matlab/cases').glob('*'));assert all(hashlib.sha256(p.read_bytes()).hexdigest()==hashlib.sha256((new/'matlab/cases'/p.name).read_bytes()).hexdigest() for p in files)
algorithms=['p1_solve_branch.m','adaptive_cv_prepared.m','cv_masks.m','cv_training_models.m','cv_normalize_context.m','cv_solve.m']
assert all(hashlib.sha256((old/'source/core_p1/matlab_src'/p).read_bytes()).hexdigest()==hashlib.sha256((new/'source/core_p1/matlab_src'/p).read_bytes()).hexdigest() for p in algorithms)
meta=old/'execution.json';data=json.loads(meta.read_text(encoding='utf-8'))
if data['status']=='RUNNING':
    (old/'execution_stale_original.json').write_bytes(meta.read_bytes());data['status']='INTERRUPTED';data['recovery_note']='MATLAB reported 0x40010004 after six saved cases; launcher was no longer present. Six cases preserved bit-for-bit; remaining three solved in the resume run.';data['recovered_into']=str(new);data['detected_utc']=datetime.datetime.now(datetime.timezone.utc).isoformat();meta.write_text(json.dumps(data,ensure_ascii=False,indent=2),encoding='utf-8')
result={'passed':True,'checkpoint_cases':6,'checkpoint_scoring_rows':12,'resumed_total_cases':9,'cached_output_files_unchanged':len(files),'algorithm_hashes_unchanged':algorithms}
(ROOT/'analysis/checkpoint_audit.json').write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf-8');print(json.dumps(result,indent=2))
