"""Check mlrandn.randn against MATLAB-generated noise archived in bae_q_20261024 (F1/F2, two SNR
levels each). make_independent.m: rng(seed,'twister'); u=randn(nt,nx); AR(1) 0.45 in time;
u=.9u+.1 mean(u,2); N=u*scale. The script inverts mixing and AR filter, fits the one scale per
gather by least squares and reports the largest absolute deviation from mlrandn.

Usage: python3 validate_randn.py <dir with bae_q_20261024 inputs/independent and truth/independent>
"""
import sys
from pathlib import Path
import numpy as np, scipy.io as sio
sys.path.insert(0, str(Path(__file__).resolve().parent))
from mlrandn import randn

base = Path(sys.argv[1])
worst = 0.0
for model, seed in [("F1", 2026102451), ("F2", 2026102461)]:
    s = sio.loadmat(base / "inputs" / "independent" / f"{model}.mat", squeeze_me=True, struct_as_record=False)
    t = sio.loadmat(base / "truth" / "independent" / f"{model}_truth.mat", squeeze_me=True, struct_as_record=False)
    D0 = t["clean"]
    for k, o in enumerate(s["obs"]):
        N = o.D - D0
        m = N.mean(1, keepdims=True)
        u = (N - .1 * m) / .9
        u0 = u.copy(); u0[1:] = u[1:] - .45 * u[:-1]
        v = u0.flatten(order="F")
        z, _ = randn(seed + 1 + k, v.size)
        scale = v @ z / (z @ z)
        err = np.abs(v / scale - z).max()
        worst = max(worst, err)
        print(f"{model} {o.name}: seed {seed + 1 + k}, {v.size} normals, max |recovered - mlrandn| = {err:.2e}")
print(f"worst {worst:.2e}")
assert worst < 1e-12
