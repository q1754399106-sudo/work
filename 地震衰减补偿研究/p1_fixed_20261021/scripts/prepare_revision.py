from pathlib import Path
import json

root=Path(__file__).resolve().parents[1]
p=root/'core_p1/run_p1_research.m'
p.write_text(p.read_text(encoding='utf-8').split('function run_p1_grid(')[0],encoding='utf-8')
p=root/'core_p1/inputs/p1_experiment_spec.json'
spec=json.loads(p.read_text(encoding='utf-8'))
spec['version']='2026-10-21-P1-fixed-v1'
spec['frozen_before_test']=False
spec['revision_note']='Bug fixes and uncertainty-averaged empirical weights; revised after observing historical tests. Requires fresh future holdouts.'
for name,role in [('TargetMatched','Frozen prior with candidate discrepancy target'),('CombinedFixedTarget','Candidate weights with Frozen discrepancy target')]:
    if not any(m['name']==name for m in spec['methods']):
        spec['methods'].append({'name':name,'role':role})
for collection in ('methods','stress_suites'):
    fields=sorted(set().union(*(r.keys() for r in spec[collection])))
    spec[collection]=[{k:r.get(k) for k in fields} for r in spec[collection]]
spec['p1_hyperparameters'].pop('offset_avo_ortho_weight_qonly',None)
spec['p1_hyperparameters'].pop('q_perturbation_factors',None)
spec['p1_hyperparameters'].pop('q_perturbation_weights',None)
p.write_text(json.dumps(spec,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
