from pathlib import Path
import json,numpy as np
from analyze_round2 import score,masks,writecsv
ROOT=Path(__file__).resolve().parent
spec=json.loads((ROOT/'SPEC.json').read_text());case=spec['variants'][0]['name'];vrun=Path(json.loads((ROOT/'validation_run.json').read_text())['run']);run=Path(json.loads((ROOT/'sensitivity_run.json').read_text())['run']);ex=json.loads((run/'execution.json').read_text());assert ex['status']=='COMPLETED' and ex['exit_code']==0
truth=np.load(ROOT/'truth'/f'{case}.npz');target=truth['Moderate_target'];t=truth['t'];live=truth['live'];low=masks(target,t,live);nominal=json.loads((vrun/case/'Moderate.json').read_text());rows=[];checks=[]
for name,factor in [('Nominal',1.),('QPoint080',.8),('QPoint125',1.25)]:
 data=nominal if name=='Nominal' else json.loads((run/name/'Moderate.json').read_text())
 if name!='Nominal':checks.append(json.loads((run/name/'checks.json').read_text()))
 for method in spec['methods']:
  s=np.array(data[method]);ref=np.array(nominal[method]);rows.append(dict(Case=name,QPointScale=factor,Method=method,**score(s,target,t,live,low),OutputChange=float(np.linalg.norm(s-ref)/max(np.linalg.norm(ref),1e-30))))
writecsv(ROOT/'analysis/sensitivity_metrics.csv',rows)
summary={m:dict(worst_error=max(r['LiveTraceError'] for r in rows if r['Method']==m),mean_error=float(np.mean([r['LiveTraceError'] for r in rows if r['Method']==m])),max_output_change=max(r['OutputChange'] for r in rows if r['Method']==m)) for m in spec['methods']}
(ROOT/'analysis/sensitivity_summary.json').write_text(json.dumps(dict(summary=summary,checks=checks,scope='Same Moderate observation, estimated Q scaled 0.8/1/1.25; not independent observations or physical true-Q tests.'),indent=2),encoding='utf-8');print(json.dumps(summary,indent=2))
