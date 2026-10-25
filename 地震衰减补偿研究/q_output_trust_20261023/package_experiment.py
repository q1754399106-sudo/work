from pathlib import Path
import json,zipfile,hashlib
ROOT=Path(__file__).resolve().parent
assert json.loads((ROOT/'COMPLETION_AUDIT.json').read_text())['completed']
dest=ROOT.parent/'q_output_trust_experiment_20261023.zip';assert not dest.exists();paths=[]
for p in ROOT.rglob('*'):
 if not p.is_file() or '__pycache__' in p.parts or p.name.endswith('_prepared.mat'):continue
 rel=p.relative_to(ROOT)
 if rel.parts[0]=='generated' and p.name not in ['VARIANT.json','PLAN.json','generated_source.py']:continue
 paths.append(p)
with zipfile.ZipFile(dest,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
 for p in paths:z.write(p,p.relative_to(ROOT.parent).as_posix())
with zipfile.ZipFile(dest) as z:assert z.testzip() is None
out=dict(path=str(dest),bytes=dest.stat().st_size,sha256=hashlib.sha256(dest.read_bytes()).hexdigest(),files=len(paths),contents='All inversion/scoring inputs and targets, original JSON outputs, source/input snapshots, logs, CSV, figures and code. Rebuildable prepared MAT caches and full-99-trace generator scratch arrays remain local; generator code and variant plans retained.')
(ROOT/'PACKAGE_MANIFEST.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2))
