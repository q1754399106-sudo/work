"""Independent small-dimensional active-set enumeration of the submitted QP.

This check does not use the candidate's matrix assembly or optimizer as its
reference. It exhausts all lower/free/upper active sets for a convex box QP.
"""
import os
for name in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS'):
    os.environ[name] = '1'
from pathlib import Path
import sys, itertools, json, hashlib
ROOT = Path(__file__).resolve().parent
TASK = ROOT.parent / 'joint_event_shrinkage_20261025'
sys.path.insert(0, str(TASK / 'scripts'))
import joint
import numpy as np

def enumerate_box(H, b):
    n = len(b)
    best = None
    for state in itertools.product((-1, 0, 1), repeat=n):
        state = np.array(state)
        free = state == 0
        x = (state == 1).astype(float)
        if free.any():
            x[free] = np.linalg.lstsq(H[np.ix_(free, free)],
                                    (b - H @ x)[free], rcond=None)[0]
        if np.any(x < -1e-9) or np.any(x > 1 + 1e-9):
            continue
        value = float(x @ H @ x - 2 * b @ x)
        if best is None or value < best[0]:
            best = value, x.copy()
    return best

def run():
    rng = np.random.default_rng(2026102517)
    records = []
    for case in range(18):
        n = 4
        Y = rng.normal(size=(7, n))
        L = rng.normal(size=(6, n))
        if case % 3 == 0:
            Y[:, 1] = -.95 * Y[:, 0] + .01 * Y[:, 1]
        G = Y.T @ Y
        N = .2 * L.T @ L
        lam = np.exp(rng.uniform(-.69, 2.77, n))
        tr = dict(G=G, N=N, E=np.diag(G), V=np.diag(N), lam=lam)
        for mode in ('Joint', 'GOnly', 'NOnly'):
            eta = (.25, .5, 1)[case % 3]
            GG = ((1-eta)*np.diag(np.diag(G)) + eta*G
                  if mode != 'NOnly' else np.diag(np.diag(G)))
            NN = ((1-eta)*np.diag(np.diag(N)) + eta*N
                  if mode != 'GOnly' else np.diag(np.diag(N)))
            D = np.diag(np.sqrt(lam))
            H = GG + D @ NN @ D
            b = GG @ np.ones(n)
            exact, x = enumerate_box(H, b)
            got, info = joint.solve_trace(tr, eta, mode)
            difference = abs(float(got @ H @ got - 2*b @ got) - exact)
            scaled = difference / max(1, abs(exact))
            assert scaled < 1e-9, (case, mode, scaled)
            assert np.max(np.abs(got-x)) < 1e-5
            records.append(dict(case=case, mode=mode, eta=eta,
                                objective_relative_difference=scaled,
                                gain_max_difference=float(np.max(np.abs(got-x)))))
    result = dict(reference='exhaustive active-set enumeration', cases=len(records),
                  passed=True, submitted_joint_sha256=hashlib.sha256(
                      (TASK/'scripts/joint.py').read_bytes()).hexdigest(),
                  max_objective_relative_difference=max(r['objective_relative_difference'] for r in records),
                  max_gain_difference=max(r['gain_max_difference'] for r in records),
                  records=records)
    (ROOT/'solver_review.json').write_text(json.dumps(result, indent=2), encoding='utf-8')
    print(json.dumps({k:v for k,v in result.items() if k!='records'}, indent=2))

if __name__ == '__main__':
    run()
