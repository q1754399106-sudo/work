from pathlib import Path
import json,hashlib
ROOT=Path(__file__).resolve().parents[1]
original=json.loads((ROOT/'refit_run.json').read_text(encoding='utf-8'))['run'].replace('\\','/')
oldroot=original.split('/runs/')[0]
plan=json.loads((ROOT/'core_p1/inputs/adaptive_cv_plan.json').read_text(encoding='utf-8'));oldcache=plan['cache'].replace('\\','/')
def relocate(value):
    if isinstance(value,dict):return {k:relocate(v) for k,v in value.items()}
    if isinstance(value,list):return [relocate(v) for v in value]
    if isinstance(value,str):
        s=value.replace('\\','/')
        if s.startswith(oldroot+'/'):return ROOT.as_posix()+s[len(oldroot):]
        if s==oldcache or s.startswith(oldcache+'/'):return (ROOT/'reproduction_inputs').as_posix()+s[len(oldcache):]
    return value
backups=ROOT/'relocation_originals';log=[]
for path in list(ROOT.glob('*_run.json'))+list((ROOT/'core_p1/inputs').glob('*.json')):
    before=path.read_bytes();data=json.loads(before.decode('utf-8-sig'));after=relocate(data)
    if after==data:continue
    relative=path.relative_to(ROOT);target=backups/relative;target.parent.mkdir(parents=True,exist_ok=True)
    if not target.exists():target.write_bytes(before)
    path.write_text(json.dumps(after,ensure_ascii=False,indent=2),encoding='utf-8')
    log.append({'path':relative.as_posix(),'before_sha256':hashlib.sha256(before).hexdigest(),'after_sha256':hashlib.sha256(path.read_bytes()).hexdigest()})
(ROOT/'relocation_log.json').write_text(json.dumps(log,ensure_ascii=False,indent=2),encoding='utf-8')
print(f'Updated {len(log)} current input/pointer files; historical runs unchanged.')
