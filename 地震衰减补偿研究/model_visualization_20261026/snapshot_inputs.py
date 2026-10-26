"""Snapshot exact archived inputs without changing shared sparse checkout."""
from pathlib import Path
import subprocess, hashlib, json

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
PREFIX = '地震衰减补偿研究/'

def git(*args, **kw):
    return subprocess.check_output(['git', *args], cwd=REPO, **kw)

revision = git('rev-parse', 'HEAD').decode().strip()
paths = [
    'calq_confirm_20261025/inputs/p1_experiment_spec.json',
    'calq_confirm_20261025/source/config_model.m',
    'calq_confirm_20261025/core/gen_complex_model.m',
    'fd2d_20261025/scripts/fd_setup.m',
    'fd2d_20261025/scripts/fd_geology.m',
    'round3_20261024/make_F34.m',
    'bae_q_20261024/make_independent.m',
    'b0_world_truth_20261024/scripts/make_unseen.m',
]
for n in ['F1','F2','F3','F4','V4','V5']:
    paths += [f'round3_20261024/dev/{n}.mat', f'round3_20261024/dev/{n}_truth.mat']
for n in ['V1','V2','V3']:
    paths += [f'coherence_soft_20261023/inputs/validation/{n}.mat',
              f'coherence_soft_20261023/truth/{n}_truth.mat']
for i in range(1,5):
    paths += [f'b0_world_truth_20261024/inputs/U{i}.mat',f'b0_world_truth_20261024/truth/U{i}_truth.mat']
for i in range(1,7):
    paths += [f'fd2d_20261025/truth/H{i}_truth.mat']
records = []
for rel in paths:
    src = PREFIX + rel
    data = git('show', f'{revision}:{src}')
    if data.startswith(b'version https://git-lfs.github.com/spec/v1'):
        oid = next(s.split(':',1)[1] for s in data.decode().splitlines() if s.startswith('oid sha256:'))
        obj = REPO / '.git/lfs/objects' / oid[:2] / oid[2:4] / oid
        data = obj.read_bytes() if obj.exists() else git('lfs','smudge', input=data)
        assert hashlib.sha256(data).hexdigest() == oid, src
    out = HERE / 'inputs' / rel
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_bytes(data)
    records.append(dict(source=src, snapshot=str(out.relative_to(HERE)).replace('\\','/'),
                        sha256=hashlib.sha256(data).hexdigest(), bytes=len(data)))
(HERE/'SOURCE_MANIFEST.json').write_text(json.dumps(dict(revision=revision,files=records),indent=2,ensure_ascii=False),encoding='utf-8')
print(f'Snapshotted {len(records)} inputs at {revision[:12]} ({sum(r["bytes"] for r in records)/1e6:.1f} MB)')
