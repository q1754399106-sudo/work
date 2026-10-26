"""Verify source integrity, exported triplets and image metadata; make QA contact sheets."""
from pathlib import Path
import hashlib, json, math
from PIL import Image, ImageOps, ImageDraw

root=Path(__file__).resolve().parent
manifest=json.loads((root/'SOURCE_MANIFEST.json').read_text(encoding='utf-8'))
for r in manifest['files']:
    data=(root/r['snapshot']).read_bytes()
    assert len(data)==r['bytes'] and hashlib.sha256(data).hexdigest()==r['sha256'], r['source']
figs=json.loads((root/'FIGURE_INDEX.json').read_text(encoding='utf-8'))
checks=json.loads((root/'VALIDATION.json').read_text(encoding='utf-8'))
assert checks['complex_V3_max_error']==0
assert all(checks['H_base_match'])
assert len(figs)==26
qa=root/'cache/qa';qa.mkdir(parents=True,exist_ok=True)
metadata=[]
for r in figs:
    stem=root/'figures'/r['name']
    for ext in ['.png','.pdf','.fig']:
        p=stem.with_suffix(ext);assert p.is_file() and p.stat().st_size>1000
    assert stem.with_suffix('.pdf').read_bytes().startswith(b'%PDF-')
    with Image.open(stem.with_suffix('.png')) as im:
        im.load();assert im.width>=1800 and im.height>=1000
        assert all(199<x<201 for x in im.info['dpi'])
        metadata.append(dict(name=r['name'],width=im.width,height=im.height,dpi=im.info['dpi'],mode=im.mode))
for page in range(math.ceil(len(figs)/6)):
    sheet=Image.new('RGB',(1440,1560),'#e9edf1');draw=ImageDraw.Draw(sheet)
    for cell,r in enumerate(figs[page*6:page*6+6]):
        with Image.open(root/r['png']) as im:
            thumb=ImageOps.contain(im.convert('RGB'),(700,475))
        xx=(cell%2)*720+10;yy=(cell//2)*520+30
        sheet.paste(thumb,(xx,yy));draw.text((xx,yy-20),r['name'],fill='black')
    sheet.save(qa/f'contact_{page+1}.jpg',quality=92)
report=dict(figures=len(figs),model_configurations=29,source_files_verified=len(manifest['files']),
            exports=metadata,checks=checks)
(root/'OUTPUT_VALIDATION.json').write_text(json.dumps(report,indent=2,ensure_ascii=False),encoding='utf-8')
print(f'PASS: {len(figs)} PNG/PDF/FIG triplets; {len(manifest["files"])} source hashes; all model checks.')
