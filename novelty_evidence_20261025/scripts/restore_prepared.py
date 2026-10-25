"""Restore immutable archived states to their ignored runtime-cache locations."""
from common import *
manifest=json.loads((ROOT/'inputs/prepared_manifest.json').read_text(encoding='utf-8'))
for name,meta in manifest.items():
    data=(ROOT/'inputs/prepared'/name).read_bytes()
    assert hashlib.sha256(data).hexdigest()==meta['sha256'] and len(data)==meta['bytes']
    dst=REPO/meta['source']
    if dst.exists():assert hashlib.sha256(dst.read_bytes()).hexdigest()==meta['sha256'],f'Existing different state retained: {dst}'
    else:dst.parent.mkdir(parents=True,exist_ok=True);dst.write_bytes(data)
print('Verified/restored',len(manifest),'prepared states')
