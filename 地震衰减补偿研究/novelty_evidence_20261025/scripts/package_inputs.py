"""Archive only the precomputed states required to reproduce this comparison."""
from common import *
import shutil,platform,scipy
manifest={}
for model in ['F1','F2','F3','F4','V4','V5']:
    for noise in ['Moderate','Strong']:
        name=f'{model}_Q+0_{noise}_nominal.mat'
        src=REPO/'round3_20261024/cache'/name
        dst=ROOT/'inputs/prepared'/name;dst.parent.mkdir(parents=True,exist_ok=True)
        data=src.read_bytes();digest=hashlib.sha256(data).hexdigest()
        if dst.exists():assert hashlib.sha256(dst.read_bytes()).hexdigest()==digest
        else:dst.write_bytes(data)
        manifest[name]=dict(bytes=len(data),sha256=digest,source=str(src.relative_to(REPO)),kind='Archived MATLAB refit state, not a new solve')
writejson(ROOT/'inputs/prepared_manifest.json',manifest)
original_inputs={}
for p in (ROOT/'cache/inputs').rglob('*.mat'):
    rel=str(p.relative_to(ROOT/'cache/inputs')).replace('\\','/')
    data=p.read_bytes();assert data==gitbytes(rel),f'Changed Git input: {rel}'
    original_inputs[rel]=dict(bytes=len(data),sha256=hashlib.sha256(data).hexdigest())
writejson(ROOT/'inputs/git_inputs_manifest.json',dict(commit='198f04b76d54c5ca77cd7b59df3a5653d6001a3e',files=original_inputs))
writejson(ROOT/'analysis/environment.json',dict(python=sys.version,numpy=np.__version__,scipy=scipy.__version__,platform=platform.platform(),blas_threads=1))
print('Prepared states:',len(manifest),'bytes:',sum(v['bytes'] for v in manifest.values()))
