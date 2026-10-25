from pathlib import Path
import json,zipfile,hashlib,sys
ROOT=Path(__file__).resolve().parents[1]
status=json.loads((ROOT/'core_p1/inputs/validation_status.json').read_text(encoding='utf-8'));assert status['ready_for_research_use'] and status['method']=='PrunedFullCalRefit'
def make(path,files):
    assert not path.exists();manifest={}
    with zipfile.ZipFile(path,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
        for p in sorted(files):
            name=p.relative_to(ROOT).as_posix();data=p.read_bytes();z.writestr(name,data);manifest[name]=hashlib.sha256(data).hexdigest()
        z.writestr('PACKAGE_SHA256.json',json.dumps(manifest,ensure_ascii=False,indent=2))
    with zipfile.ZipFile(path) as z:
        assert z.testzip() is None
        assert all(hashlib.sha256(z.read(name)).hexdigest()==sha for name,sha in manifest.items())
    return {'path':str(path),'files':len(manifest),'bytes':path.stat().st_size,'sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'verified':True}
core=[p for p in (ROOT/'core_p1').rglob('*') if p.is_file() and '__pycache__' not in p.parts]
docs=[ROOT/p for p in ['START_HERE.md','METHOD_ZH.md','RESULTS_ZH.md','ROUND1_RESULTS_ZH.md']]
if '--results-only' in sys.argv:
    code=json.loads((ROOT/'analysis/code_package.json').read_text(encoding='utf-8'))
    assert json.loads((ROOT/'analysis/delivery_check.json').read_text(encoding='utf-8'))['passed']
else:
    code=make(ROOT.parent/'optimized_seismic_method_20261023.zip',core+docs+[ROOT/'figures/overall_comparison.png'])
    (ROOT/'analysis/code_package.json').write_text(json.dumps(code,ensure_ascii=False,indent=2),encoding='utf-8')
if '--code-only' in sys.argv:print(json.dumps(code,ensure_ascii=False,indent=2));raise SystemExit(0)
allfiles=[p for p in ROOT.rglob('*') if p.is_file() and '__pycache__' not in p.parts and p.name!='PACKAGE_RESULTS.json']
full=make(ROOT.parent/'optimized_seismic_experiments_20261023.zip',allfiles)
(ROOT/'PACKAGE_RESULTS.json').write_text(json.dumps({'code':code,'experiments':full},ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps({'code':code,'experiments':full},ensure_ascii=False,indent=2))
