from pathlib import Path
import json
ROOT=Path(__file__).resolve().parent
mainrun=Path(json.loads((ROOT/'validation_run.json').read_text(encoding='utf-8'))['run'])
status=json.loads((mainrun/'execution.json').read_text(encoding='utf-8'));assert status['status']=='COMPLETED'
base=json.loads((ROOT/'inputs/FreshFoldLowQ.json').read_text(encoding='utf-8'))
cal=json.loads((mainrun/'FreshFoldLowQ/calibration.json').read_text(encoding='utf-8'))
base['cached_calibration']=cal;base['datasets']=[base['datasets'][0]]
for scale,label in [(.8,'QPoint080'),(1.25,'QPoint125')]:
    d=dict(base);d['case_name']=label;d['q_point_scale']=scale
    p=ROOT/'inputs'/f'{label}.json';assert not p.exists();p.write_text(json.dumps(d,separators=(',',':')),encoding='utf-8')
print('Registered Q-point sensitivity inputs ready')
