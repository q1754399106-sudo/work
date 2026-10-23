from pathlib import Path
import json, re
from PIL import Image
from pypdf import PdfReader
import pypdfium2 as pdfium

root = Path(__file__).resolve().parent
manifest = json.loads((root / 'PLOT_MANIFEST.json').read_text(encoding='utf-8'))
assert len(manifest) == 32
pdfs = []
for entry in manifest:
    with Image.open(root / 'figures' / entry['file']) as im:
        assert im.size == (entry['width_pixels'], entry['height_pixels'])
        assert abs(im.info['dpi'][0] - 400) < 1
    if entry['view'] == 'detail':
        assert (entry['width_pixels'], entry['height_pixels']) == (6400, 2720)
    if entry['pdf_file']:
        path = root / 'figures' / entry['pdf_file']
        reader = PdfReader(path)
        assert len(reader.pages) == 1
        assert len(reader.pages[0].extract_text()) > 20
        pdfs.append(path.name)
assert len(pdfs) == 24
html = (root / 'index.html').read_text(encoding='utf-8')
views = json.loads(re.search(r'const views=(.*?);const \$', html).group(1))
assert len(views) == 32
assert all((root / 'figures' / v['file']).is_file() for v in views.values())
qa = root / 'qa'
qa.mkdir(exist_ok=True)
for name in ['complex_Strong_Q+0_detail.pdf', 'primary_Strong_wiggle_0.pdf']:
    doc = pdfium.PdfDocument(str(root / 'figures' / name))
    bitmap = doc[0].render(scale=1.5)
    bitmap.to_pil().save(qa / (Path(name).stem + '_pdf_preview.png'))
    bitmap.close()
    doc.close()
result = dict(pngs_checked=32, pdfs_checked=24, gallery_views_checked=32,
              detail_pixels=[6400, 2720], dpi=400, pdf_previews=2)
(qa / 'CHECKS.json').write_text(json.dumps(result, indent=2), encoding='utf-8')
print(result)
