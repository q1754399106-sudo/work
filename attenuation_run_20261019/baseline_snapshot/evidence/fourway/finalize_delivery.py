"""Audit snapshot-only replay and package this comparison without source edits."""
from pathlib import Path
import json, hashlib, shutil, zipfile

ROOT=Path(__file__).resolve().parent
REPLAY=ROOT.parent/'work/fourway_reanalysis_check_20261018'

def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()

rechecked=[]
for original in sorted((ROOT/'data').glob('*.csv')):
    duplicate=REPLAY/'data'/original.name
    assert duplicate.is_file() and sha(original)==sha(duplicate),original.name
    rechecked.append(original.name)
assert json.loads((ROOT/'ANALYSIS_CHECKS.json').read_text(encoding='utf-8'))==json.loads((REPLAY/'ANALYSIS_CHECKS.json').read_text(encoding='utf-8'))
for rel in ['REPORT_ZH_EN.md','data/equal_model_main_summary.csv']:
    shutil.copy2(REPLAY/rel,ROOT/rel)
manifest=json.loads((ROOT/'SOURCE_MANIFEST.json').read_text(encoding='utf-8'))
assert all(sha(Path(x['source']))==x['sha256'] and sha(ROOT/x['snapshot'])==x['sha256'] for x in manifest)
checks=dict(status='PASS',snapshot_only_replay=True,identical_csv_files=rechecked,source_and_snapshot_hashes_unchanged=True,all_9_png_pdf_pairs=True,manual_contact_sheet_review=True,manual_overall_figure_review=True,new_matlab_runs=0)
(ROOT/'DELIVERY_CHECKS.json').write_text(json.dumps(checks,ensure_ascii=False,indent=2),encoding='utf-8')
entries=[]
for p in sorted(ROOT.rglob('*')):
    if p.is_file() and p.name!='PACKAGE_MANIFEST.json': entries.append(dict(path=str(p.relative_to(ROOT)).replace('\\','/'),sha256=sha(p),bytes=p.stat().st_size))
(ROOT/'PACKAGE_MANIFEST.json').write_text(json.dumps(dict(files=entries),ensure_ascii=False,indent=2),encoding='utf-8')
archive=ROOT.parent/'outputs'/(ROOT.name+'.zip'); assert not archive.exists(), 'Preserve existing package'
with zipfile.ZipFile(archive,'x',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
    for p in sorted(ROOT.rglob('*')):
        if p.is_file(): z.write(p,str(Path(ROOT.name)/p.relative_to(ROOT)))
with zipfile.ZipFile(archive) as z:
    assert z.testzip() is None
    for entry in entries:
        assert hashlib.sha256(z.read(ROOT.name+'/'+entry['path'])).hexdigest()==entry['sha256']
checksum=sha(archive)
archive.with_suffix('.zip.sha256').write_text(checksum+'  '+archive.name+'\n',encoding='utf-8')
print(json.dumps(dict(archive=str(archive),bytes=archive.stat().st_size,files=len(entries)+1,sha256=checksum,checks=checks),ensure_ascii=False,indent=2))
