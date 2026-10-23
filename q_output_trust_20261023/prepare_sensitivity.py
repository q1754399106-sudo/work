from pathlib import Path
import json
ROOT=Path(__file__).resolve().parent
spec=json.loads((ROOT/'SPEC.json').read_text());case=spec['variants'][0]['name'];run=Path(json.loads((ROOT/'validation_run.json').read_text())['run']);ex=json.loads((run/'execution.json').read_text());assert ex['status']=='COMPLETED'
base=json.loads((ROOT/'inputs'/f'{case}.json').read_text());base['cached_calibration']=json.loads((run/case/'calibration.json').read_text());base['datasets']=[base['datasets'][0]]
for name,factor in [('QPoint080',.8),('QPoint125',1.25)]:
 d=dict(base);d['case_name']=name;d['q_point_scale']=factor;p=ROOT/'inputs'/f'{name}.json';assert not p.exists();p.write_text(json.dumps(d,separators=(',',':')),encoding='utf-8')
print('Q-point sensitivity inputs ready')
