from pathlib import Path
import json,hashlib,zipfile
ROOT=Path(__file__).resolve().parents[1]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 run=Path(json.loads((ROOT/'LATEST_RUN.json').read_text(encoding='utf-8'))['run'])
 assert json.loads((run/'execution.json').read_text())['status']=='COMPLETED'
 assert json.loads((run/'matlab/checks.json').read_text())['rows']==120
 figs=json.loads((ROOT/'figures/manifest.json').read_text(encoding='utf-8'))
 assert len(figs['case_scales'])==60 and len(figs['figures'])==368
 assert all(sha(ROOT/'figures'/f['file'])==f['sha256'] for f in figs['figures'])
 source=json.loads((run/'SOURCE_MANIFEST.json').read_text())
 assert all(sha(run/'source/core_p1'/k)==v for k,v in source.items())
 paths=[p for p in sorted(ROOT.rglob('*')) if p.is_file() and '__pycache__' not in p.parts]
 paths=[p for p in paths if p.name!='PACKAGE_MANIFEST.json']
 manifest={p.relative_to(ROOT).as_posix():sha(p) for p in paths}
 (ROOT/'PACKAGE_MANIFEST.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2),encoding='utf-8')
 paths.append(ROOT/'PACKAGE_MANIFEST.json')
 archive=ROOT.parent/'fixed_vs_new_gathers_comparison_20261022.zip';temporary=archive.with_suffix('.building.zip');assert not temporary.exists()
 with zipfile.ZipFile(temporary,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
  for p in paths:z.write(p,p.relative_to(ROOT).as_posix())
 with zipfile.ZipFile(temporary) as z:
  assert z.testzip() is None
  assert all(hashlib.sha256(z.read(k)).hexdigest()==v for k,v in manifest.items())
 temporary.replace(archive)
 print(json.dumps({'archive':str(archive),'bytes':archive.stat().st_size,'files':len(paths),'sha256':sha(archive)},ensure_ascii=False,indent=2))
if __name__=='__main__':main()
