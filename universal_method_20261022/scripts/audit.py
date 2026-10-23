from pathlib import Path
import json,hashlib,math,csv
ROOT=Path(__file__).resolve().parents[1]
run=Path(json.loads((ROOT/'latest_run.json').read_text(encoding='utf-8'))['run'])
source=run/'source/core_p1'
manifest=json.loads((run/'SOURCE_MANIFEST.json').read_text(encoding='utf-8'))
assert all(hashlib.sha256((source/k).read_bytes()).hexdigest()==v for k,v in manifest.items())
plan=json.loads((source/'inputs/universal_plan.json').read_text(encoding='utf-8'));cache=Path(plan['cache'])
policy=json.loads((source/'inputs/candidate_policy.json').read_text(encoding='utf-8'))
qerr=ferr=0.;count=0
for p in (cache/'cases').glob('*.json'):
    meta=json.loads(p.read_text(encoding='utf-8'))['meta']
    a=json.loads((cache/f"calibration_{meta['model']}.json").read_text(encoding='utf-8'))
    fixed=a['diagnostics']['fixed_estimate'];base=a['estimate']
    q=base.get('layer_Q',[base['effective_Q']]*4);f=base['estimated_frequency']
    if meta['domain']=='Section':
        w=policy['selected_profile']['section_calibration_weight']
        q=[math.exp((1-w)*math.log(fixed['effective_Q'])+w*math.log(v)) for v in q]
        f=math.exp((1-w)*math.log(fixed['estimated_frequency'])+w*math.log(f))
    qerr=max(qerr,max(abs(x-y) for x,y in zip(q,meta['Q_used'])))
    ferr=max(ferr,abs(f-meta['estimated_frequency']));count+=1
assert count==60 and qerr<1e-10 and ferr<1e-10
rows=list(csv.DictReader((run/'matlab/metrics.csv').open(encoding='utf-8-sig')))
assert len(rows)==180 and len({r['CaseID'] for r in rows})==60
for case in {r['CaseID'] for r in rows}:
    assert {r['Method'] for r in rows if r['CaseID']==case}=={'Baseline','SourceMixture','SourceOnly'}
checks=json.loads((run/'matlab/checks.json').read_text(encoding='utf-8'))
assert checks['complete'] and checks['max_normal']<2e-7 and checks['baseline_score_replay']<1e-12
report={'source_hashes_unchanged':True,'paired_cases':count,'same_calibration_Q_max_difference':qerr,'same_calibration_frequency_max_difference':ferr,'baseline_score_replay':checks['baseline_score_replay'],'max_normal_residual':checks['max_normal'],'algorithm_rule':'Same fixed profile across all models and scenarios; truth/scene used by scoring driver only'}
cp=ROOT/'consensus_run.json'
if cp.exists():
    crun=Path(json.loads(cp.read_text(encoding='utf-8'))['run'])
    cm=json.loads((crun/'SOURCE_MANIFEST.json').read_text(encoding='utf-8'))
    assert all(hashlib.sha256((crun/'source/core_p1'/k).read_bytes()).hexdigest()==v for k,v in cm.items())
    cc=json.loads((crun/'matlab/checks.json').read_text(encoding='utf-8'))
    assert cc['complete'] and cc['cases']==60 and cc['weights']==[.5,.5]
    report['consensus_source_hashes_unchanged']=True
    report['consensus_cases']=60
(ROOT/'analysis/fairness_audit.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
print(json.dumps(report,indent=2))
