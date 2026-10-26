from pathlib import Path
import sys,json,hashlib,math
from PIL import Image,ImageOps,ImageDraw
root=Path(__file__).resolve().parent
sys.path.insert(0,str(root.parent/'novelty_evidence_20261025/cache/deps'))
import numpy as np
from scipy.io import loadmat

cases=json.loads((root/'CASE_INDEX.json').read_text(encoding='utf-8'))
figures=json.loads((root/'FIGURE_INDEX.json').read_text(encoding='utf-8'))
assert len(cases)==35 and len({r['model'] for r in cases})==29
assert len(figures)==38 and len({r['stem'] for r in figures})==38
manifest=json.loads((root/'SOURCE_MANIFEST.json').read_text(encoding='utf-8'))
for r in manifest['files']:
    p=root/r['snapshot'];assert hashlib.sha256(p.read_bytes()).hexdigest()==r['sha256']
for r in cases:
    s=loadmat(root/'results'/f'{r["case"]}.mat');D,T,X=[s[k] for k in ['D','target','compensated']]
    assert D.shape==T.shape==X.shape and np.isfinite(D).all() and np.isfinite(T).all() and np.isfinite(X).all()
    assert abs(np.linalg.norm(X-T)/np.linalg.norm(T)-r['output_error'])<1e-12
    assert abs(np.linalg.norm(D-T)/np.linalg.norm(T)-r['input_error'])<1e-12
    assert r['solver']['converged'] and r['solver']['kkt_ratio']<=.02
    for a in r['sources']:
        assert hashlib.sha256((root.parent/a['path']).read_bytes()).hexdigest()==a['sha256']
    if r['case'].endswith('_Qminus30'):
        ref=loadmat(root/'results'/f'{r["model"]}.mat');assert np.array_equal(D,ref['D']) and np.array_equal(T,ref['target'])
metadata=[]
for r in figures:
    stem=root/'figures'/r['stem']
    for ext in ['.png','.pdf','.fig']:assert stem.with_suffix(ext).stat().st_size>1000
    assert stem.with_suffix('.pdf').read_bytes().startswith(b'%PDF-')
    with Image.open(stem.with_suffix('.png')) as im:
        im.load();assert im.width>1800 and im.height>1400
        assert all(abs(x-180)<.2 for x in im.info['dpi'])
        metadata.append(dict(stem=r['stem'],pixels=im.size,dpi=im.info['dpi']))
qa=root/'cache/qa';qa.mkdir(parents=True,exist_ok=True)
for page in range(math.ceil(len(figures)/6)):
    sheet=Image.new('RGB',(1440,1560),'#e9edf1');draw=ImageDraw.Draw(sheet)
    for i,r in enumerate(figures[page*6:page*6+6]):
        with Image.open(root/'figures'/f'{r["stem"]}.png') as im:thumb=ImageOps.contain(im.convert('RGB'),(700,475))
        x=(i%2)*720+10;y=(i//2)*520+30;sheet.paste(thumb,(x,y));draw.text((x,y-20),r['stem'],fill='black')
    sheet.save(qa/f'contact_{page+1}.jpg',quality=93)
report=dict(cases=len(cases),models=29,figures=len(figures),source_files=len(manifest['files']),
            all_converged=True,max_kkt=max(r['solver']['kkt_ratio'] for r in cases),
            input_output_metrics_verified=True,H_same_observation_verified=True,numpy=np.__version__,exports=metadata)
(root/'OUTPUT_VALIDATION.json').write_text(json.dumps(report,indent=2,ensure_ascii=False),encoding='utf-8')
print('PASS: 35 cases / 29 models / 38 triplets; metrics, source hashes, convergence and paired H inputs verified.')
