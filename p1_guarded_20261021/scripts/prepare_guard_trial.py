from pathlib import Path
import json
root=Path(__file__).resolve().parents[1]
p=root/'core_p1/inputs/p1_experiment_spec.json'
spec=json.loads(p.read_text(encoding='utf-8'))
spec['version']='2026-10-21-guard-trial-v1'
spec['frozen_before_test']=False
spec['revision_note']='Exploratory guarded routing based on previously observed failures; no claim of blind testing.'
for name in ('GuardedNoise','NoCompensation'):
    if not any(r['name']==name for r in spec['methods']): spec['methods'].append({'name':name,'role':'guard_policy_trial'})
fields=sorted(set().union(*(r.keys() for r in spec['methods'])))
spec['methods']=[{k:r.get(k) for k in fields} for r in spec['methods']]
p.write_text(json.dumps(spec,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
