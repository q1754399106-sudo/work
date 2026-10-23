from pathlib import Path
import hashlib,json,zipfile
from PIL import Image

R=Path(__file__).resolve().parent
plots=json.loads((R/'PLOT_MANIFEST.json').read_text(encoding='utf-8'))
assert len(plots)==38
views=json.loads((R/'GALLERY_VIEWS.json').read_text(encoding='utf-8'));assert len(views)==36
for p in plots:
    path=R/'figures'/p['file']
    assert hashlib.sha256(path.read_bytes()).hexdigest()==p['sha256']
    with Image.open(path) as im:im.verify()
    with Image.open(path) as im:
        assert im.size==(p['width'],p['height']) and abs(im.info['dpi'][0]-400)<1
    if p['view']=='detail':assert (p['width'],p['height'])==(7600,2960)
for v in views.values():assert (R/'figures'/v['file']).is_file()
assert json.loads((R/'checks.json').read_text())['optimization_flags']==0
assert json.loads((R/'API_CHECK.json').read_text())['maximum_relative_difference']==0
for file,h in json.loads((R/'TRIAL_CODE_HASHES.json').read_text()).items():assert hashlib.sha256((R/file).read_bytes()).hexdigest()==h
for file,h in json.loads((R/'PUBLIC_API_HASHES.json').read_text()).items():assert hashlib.sha256((R/file).read_bytes()).hexdigest()==h

def flatten(value):
    if isinstance(value,list):
        for z in value:yield from flatten(z)
    elif isinstance(value,(int,float)):yield value

gains=[]
for p in (R/'results').glob('*/Q_*.json'):
    for record in json.loads(p.read_text(encoding='utf-8'))['records']:
        for d in record['diagnostics']:gains.extend(flatten(d['group_gain']))
assert gains and all(0<=g<=1 for g in gains)
for p in R.glob('*.py'):compile(p.read_text(encoding='utf-8'),str(p),'exec')
qa=dict(figures_checked=38,gallery_views=36,detail_pixels=[7600,2960],dpi=400,
        gain_groups=len(gains),minimum_gain=min(gains),maximum_gain=max(gains),
        zero_gains=sum(g==0 for g in gains),source_and_api_hashes_match=True)
(R/'QA_CHECK.json').write_text(json.dumps(qa,indent=2),encoding='utf-8')
(R/'RUN_LEDGER.json').write_text(json.dumps(dict(status='complete',completed_blocks=180,expected_blocks=180,
    gallery=str(R/'index.html'),figures=38,storage=str(R),inference_rule='one frozen candidate',api_verified=True),indent=2),encoding='utf-8')

excluded={'PACKAGE_MANIFEST.json','DELIVERY_CHECK.json'}
files=[p for p in R.rglob('*') if p.is_file() and p.name not in excluded and '__pycache__' not in p.parts]
baseline=R.parent/'q_dense_gathers_20261023'
files+=list((baseline/'results').rglob('*.json'))
records=[dict(path=p.relative_to(R.parent).as_posix(),bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in files]
(R/'PACKAGE_MANIFEST.json').write_text(json.dumps(records,ensure_ascii=False,indent=2),encoding='utf-8')
archive=R.parent/'soft_refit_complete_20261023.zip'
with zipfile.ZipFile(archive,'w',zipfile.ZIP_DEFLATED,compresslevel=5) as z:
    for p in files:z.write(p,p.relative_to(R.parent).as_posix())
    z.write(R/'PACKAGE_MANIFEST.json','soft_refit_20261023/PACKAGE_MANIFEST.json')
with zipfile.ZipFile(archive) as z:
    assert z.testzip() is None
    for item in records:assert hashlib.sha256(z.read(item['path'])).hexdigest()==item['sha256']
delivery=dict(complete=True,files=len(files)+1,bytes=archive.stat().st_size,archive=str(archive),figures=38,all_hashes_verified=True)
(R/'DELIVERY_CHECK.json').write_text(json.dumps(delivery,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps(delivery,ensure_ascii=False))
