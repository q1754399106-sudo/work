from pathlib import Path
import argparse, difflib, hashlib, json, zipfile

ROOT=Path(__file__).resolve().parents[1]
p=argparse.ArgumentParser(); p.add_argument('--run',type=Path,required=True); p.add_argument('--original',type=Path,required=True)
a=p.parse_args(); original=a.original.resolve(); run=a.run.resolve()
changes=[]
for path in sorted((ROOT/'core_p1').rglob('*')):
    if not path.is_file(): continue
    relative=path.relative_to(ROOT); before=original/relative
    old=before.read_text(encoding='utf-8').splitlines(keepends=True) if before.exists() else []
    new=path.read_text(encoding='utf-8').splitlines(keepends=True)
    changes.extend(difflib.unified_diff(old,new,fromfile='original/'+relative.as_posix(),tofile='fixed/'+relative.as_posix()))
(ROOT/'CHANGES.patch').write_text(''.join(changes),encoding='utf-8')
files=[]
for base in (ROOT/'core_p1',ROOT/'scripts',run):
    files += [q for q in base.rglob('*') if q.is_file() and '__pycache__' not in q.parts]
files += [ROOT/q for q in ('README_ZH.md','RESULTS_ZH.md','SOURCE_MANIFEST.json','NEXT_EXPERIMENT.json','CHANGES.patch')]
files=sorted(set(files))
manifest=[{'path':f.relative_to(ROOT).as_posix(),'bytes':f.stat().st_size,'sha256':hashlib.sha256(f.read_bytes()).hexdigest()} for f in files]
(ROOT/'PACKAGE_MANIFEST.json').write_text(json.dumps({'file_count':len(manifest),'files':manifest},indent=2),encoding='utf-8')
files.append(ROOT/'PACKAGE_MANIFEST.json')
target=ROOT.parent/'p1_fixed_20261021.zip'
if target.exists(): raise SystemExit('Refusing to overwrite existing delivery ZIP')
with zipfile.ZipFile(target,'w',compression=zipfile.ZIP_DEFLATED) as z:
    for f in files: z.write(f,arcname=ROOT.name+'/'+f.relative_to(ROOT).as_posix())
with zipfile.ZipFile(target) as z: assert z.testzip() is None
print(json.dumps({'zip':str(target),'files':len(files),'bytes':target.stat().st_size},ensure_ascii=False))
