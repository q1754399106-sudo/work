from pathlib import Path
import json,csv,hashlib
import numpy as np
from analyze import score,masks,writecsv
ROOT=Path(__file__).resolve().parent
main=Path(json.loads((ROOT/'validation_run.json').read_text(encoding='utf-8'))['run'])
run=Path(json.loads((ROOT/'sensitivity_run.json').read_text(encoding='utf-8'))['run'])
status=json.loads((run/'execution.json').read_text(encoding='utf-8'));assert status['status']=='COMPLETED' and status['exit_code']==0
spec=json.loads((ROOT/'SPEC.json').read_text(encoding='utf-8'));truth=np.load(ROOT/'truth/FreshFoldLowQ.npz')
t=truth['t'];live=truth['live'];target=truth['Moderate_target'];low=masks(target,t,live)
nominal=json.loads((main/'FreshFoldLowQ/Moderate.json').read_text(encoding='utf-8'))
rows=[];checks=[]
for case,factor in [('Nominal',1.),('QPoint080',.8),('QPoint125',1.25)]:
    if case=='Nominal':result=nominal
    else:
        checks.append(json.loads((run/case/'checks.json').read_text(encoding='utf-8')))
        result=json.loads((run/case/'Moderate.json').read_text(encoding='utf-8'))
    for method in spec['methods']:
        s=np.asarray(result[method]);base=np.asarray(nominal[method])
        rows.append(dict(Case=case,QPointScale=factor,Method=method,**score(s,target,t,live,low),
            OutputChangeFromNominal=float(np.linalg.norm(s-base)/max(np.linalg.norm(base),1e-30))))
writecsv(ROOT/'analysis/sensitivity_metrics.csv',rows)
summary={m:dict(worst_error=max(r['LiveTraceError'] for r in rows if r['Method']==m),
    mean_error=float(np.mean([r['LiveTraceError'] for r in rows if r['Method']==m])),
    maximum_output_change=max(r['OutputChangeFromNominal'] for r in rows if r['Method']==m)) for m in spec['methods']}
(ROOT/'analysis/sensitivity_summary.json').write_text(json.dumps(dict(scope='First fresh geology, same Moderate observation, nominal and artificial point-Q scales 0.8/1.25; not independent validation cases or true-Q perturbations.',summary=summary,checks=checks),indent=2),encoding='utf-8')
print(json.dumps(summary,indent=2))
