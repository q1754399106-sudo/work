from pathlib import Path
import subprocess,hashlib,json

ROOT=Path(__file__).resolve().parent
REPO=ROOT.parents[1]
REV=subprocess.check_output(['git','rev-parse','HEAD'],cwd=REPO).decode().strip()
PREFIX='地震衰减补偿研究/'
records=[]

def extract(rel):
    data=subprocess.check_output(['git','show',f'{REV}:{PREFIX}{rel}'],cwd=REPO)
    if data.startswith(b'version https://git-lfs.github.com/spec/v1'):
        oid=next(x.split(':',1)[1] for x in data.decode().splitlines() if x.startswith('oid sha256:'))
        obj=REPO/'.git/lfs/objects'/oid[:2]/oid[2:4]/oid
        data=obj.read_bytes() if obj.exists() else subprocess.check_output(['git','lfs','smudge'],input=data,cwd=REPO)
        assert hashlib.sha256(data).hexdigest()==oid
    return data

def save(rel,target=None):
    data=extract(rel);target=target or ('cache/inputs/'+rel)
    out=ROOT/target;out.parent.mkdir(parents=True,exist_ok=True);out.write_bytes(data)
    records.append(dict(source=PREFIX+rel,snapshot=target,sha256=hashlib.sha256(data).hexdigest(),bytes=len(data)))

paths=subprocess.check_output(['git','ls-tree','-r','--name-only','-z',REV,'--',PREFIX+'compensation_package_20261024'],cwd=REPO).decode().split('\0')
for path in filter(None,paths):
    rel=path[len(PREFIX):]
    if rel.endswith(('.m','.json')):save(rel,'source/'+rel)
for fname in ['benchmark.py','common.py']:
    save('novelty_evidence_20261025/scripts/'+fname,'source/sc/scripts/'+fname)
save('novelty_evidence_20261025/analysis/selection.json','source/sc/analysis/selection.json')
for name in ['F1','F2','F3','F4','V4','V5']:
    save(f'round3_20261024/dev/{name}.mat');save(f'round3_20261024/dev/{name}_truth.mat')
    save(f'round3_20261024/results/{name}_Q+0_Moderate.mat')
for name in ['V1','V2','V3']:
    save(f'coherence_soft_20261023/inputs/validation/{name}.mat');save(f'coherence_soft_20261023/truth/{name}_truth.mat')
    save(f'coherence_soft_20261023/results/{name}_Q+0_Moderate.mat')
for k in range(1,5):save(f'b0_world_truth_20261024/inputs/U{k}.mat');save(f'b0_world_truth_20261024/truth/U{k}_truth.mat')
for k in range(1,7):save(f'fd2d_20261025/models/H{k}.mat');save(f'fd2d_20261025/truth/H{k}_truth.mat')
(ROOT/'SOURCE_MANIFEST.json').write_text(json.dumps(dict(revision=REV,files=records),ensure_ascii=False,indent=2),encoding='utf-8')
print('Snapshots',len(records),'bytes',sum(x['bytes'] for x in records))
