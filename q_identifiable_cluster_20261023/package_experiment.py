from pathlib import Path
import zipfile,json,hashlib
ROOT=Path(__file__).resolve().parent
audit=json.loads((ROOT/'COMPLETION_AUDIT.json').read_text(encoding='utf-8'));assert audit['completed']
dest=ROOT.parent/'q_identifiable_cluster_experiment_20261023.zip'
assert not dest.exists(),'Refusing to overwrite prior package'
paths=[]
for p in ROOT.rglob('*'):
    if not p.is_file() or '__pycache__' in p.parts:continue
    rel=p.relative_to(ROOT)
    if rel.parts[0]=='generated' and p.name not in ['VARIANT.json','PLAN.json','generated_source.py']:continue
    paths.append(p)
with zipfile.ZipFile(dest,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
    for p in paths:z.write(p,p.relative_to(ROOT.parent).as_posix())
with zipfile.ZipFile(dest) as z:
    assert z.testzip() is None
out=dict(path=str(dest),bytes=dest.stat().st_size,sha256=hashlib.sha256(dest.read_bytes()).hexdigest(),files=len(paths),
    contents='All actual inversion/scoring inputs and 22-trace targets, exact run outputs/snapshots/logs, code, frozen manifests, CSV and figures. Full 99-trace generator scratch arrays remain locally in generated/ and are excluded; generator source and variant plans retained.')
(ROOT/'PACKAGE_MANIFEST.json').write_text(json.dumps(out,indent=2),encoding='utf-8')
print(json.dumps(out,indent=2))
