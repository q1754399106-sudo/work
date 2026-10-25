"""Check actual F1 cross-group noise covariance using direct noise maps."""
from pathlib import Path
import sys, json, hashlib
ROOT = Path(__file__).resolve().parent
TASK = ROOT.parent/'joint_event_shrinkage_20261025'
sys.path.insert(0, str(TASK/'scripts'))
import joint
from common import load_scene, dictionary, dense, vec
import numpy as np

cm, obs, sigma, st, _ = load_scene('F1', 'Moderate')
prepared = joint.prepare(cm, obs, sigma, st)
B, W, P, scale = dictionary(cm)
R = dense(st['R0'])
eligible = [j for j, t in enumerate(prepared['traces']) if t is not None]
indices = [eligible[0], eligible[len(eligible)//2], eligible[-1]]
rows = []
for j in indices:
    s = vec(st['supports'][j], int)-1
    gap = max(1, int(np.floor(.5/(cm['fpeak']*cm['dt'])+.5)))
    labels = np.cumsum(np.r_[0, np.diff(s)>gap])
    C = B[:, s]
    singular = np.linalg.svd(C, compute_uv=False)
    top, bottom = singular[0]**2, singular[-1]**2
    ridge = max(1e-8*top, (top-1e4*bottom)/(1e4-1))
    # Independent dense solve rather than the submitted SVD covariance product.
    K = np.linalg.solve(C.T@C + ridge*np.eye(len(s)), C.T)
    maps = []
    waves = []
    for k in range(int(labels.max())+1):
        mask = labels == k
        maps.append((sigma[j]*(P[:, s[mask]]@K[mask, :])).ravel())
        waves.append(W[:, s[mask]]@R[s[mask], j])
    maps = np.asarray(maps)
    waves = np.asarray(waves)
    N = maps@maps.T
    G = waves@waves.T
    submitted = prepared['traces'][j]
    ne = float(np.linalg.norm(N-submitted['N'])/max(np.linalg.norm(N), 1e-30))
    ge = float(np.linalg.norm(G-submitted['G'])/max(np.linalg.norm(G), 1e-30))
    assert ne < 1e-9 and ge < 1e-9, (j, ne, ge)
    rows.append(dict(trace=j, groups=len(waves), noise_relative_difference=ne,
                     waveform_gram_relative_difference=ge))
result = dict(passed=True, scene='F1_Moderate',
              reference='dense ridge solve and explicit output noise maps',
              joint_sha256=hashlib.sha256((TASK/'scripts/joint.py').read_bytes()).hexdigest(),
              traces=rows)
(ROOT/'covariance_review.json').write_text(json.dumps(result, indent=2), encoding='utf-8')
print(json.dumps(result, indent=2))
