from pathlib import Path
import argparse, difflib, hashlib, json, zipfile

ROOT=Path(__file__).resolve().parents[1]
p=argparse.ArgumentParser(); p.add_argument('--run',type=Path,required=True); a=p.parse_args(); run=a.run.resolve()
old=ROOT.parent/'p1_fixed_20261021'
source=[]; diff=[]
for file in sorted((ROOT/'core_p1').rglob('*')):
    if not file.is_file(): continue
    relative=file.relative_to(ROOT)
    source.append({'path':relative.as_posix(),'sha256':hashlib.sha256(file.read_bytes()).hexdigest()})
    before=old/relative
    previous=before.read_text(encoding='utf-8').splitlines(keepends=True) if before.exists() else []
    diff.extend(difflib.unified_diff(previous,file.read_text(encoding='utf-8').splitlines(keepends=True),
                                  fromfile='fixed/'+relative.as_posix(),tofile='guarded/'+relative.as_posix()))
(ROOT/'SOURCE_MANIFEST.json').write_text(json.dumps(source,indent=2),encoding='utf-8')
(ROOT/'CHANGES.patch').write_text(''.join(diff),encoding='utf-8')
files=[]
for base in (ROOT/'core_p1',ROOT/'scripts',run):
    files.extend(q for q in base.rglob('*') if q.is_file() and '__pycache__' not in q.parts)
files.extend(ROOT/q for q in ('README_ZH.md','RESULTS_ZH.md','SOURCE_MANIFEST.json','CHANGES.patch','NEXT_EXPERIMENT.json'))
files=sorted(set(files))
manifest=[{'path':q.relative_to(ROOT).as_posix(),'bytes':q.stat().st_size,'sha256':hashlib.sha256(q.read_bytes()).hexdigest()} for q in files]
(ROOT/'PACKAGE_MANIFEST.json').write_text(json.dumps({'file_count':len(manifest),'files':manifest},indent=2),encoding='utf-8')
files.append(ROOT/'PACKAGE_MANIFEST.json')
target=ROOT.parent/'p1_guarded_20261021.zip'
if target.exists(): raise SystemExit('Refusing to overwrite delivery ZIP')
with zipfile.ZipFile(target,'w',compression=zipfile.ZIP_DEFLATED) as z:
    for q in files: z.write(q,ROOT.name+'/'+q.relative_to(ROOT).as_posix())
with zipfile.ZipFile(target) as z:
    assert z.testzip() is None
    assert all(hashlib.sha256(z.read(ROOT.name+'/'+q['path'])).hexdigest()==q['sha256'] for q in manifest)
print(json.dumps({'zip':str(target),'verified_files':len(manifest),'bytes':target.stat().st_size},ensure_ascii=False))
