"""Frequency-domain Q-compensated reverse time migration (exploratory, synthetic): cross-check of the Kirchhoff
pilot with wave-equation operators.

For the migration model (smooth c0, smooth Q times m) and amplitude factor lam:
  lam = +1  both wavefields amplitude-boosted (conventional Q-RTM);  lam = 0  dispersion only (no boost).
Per frequency:  I(x, f) = sum_shots Re[conj(S) R],  S = A^-1 src W,  R = conj(A^-1 conj(d)).
Images per frequency are stored, so the low-pass stabilisation (cosine taper from fc to fc + 8 Hz) is applied
afterwards. Signal and a unit-variance noise record are migrated separately (the operator is linear).
Reference: the attenuation-free data migrated with the acoustic smooth operator.
usage: rtm.py            -> CACHE/rtm.npz  (ops, nf, 151, NZ) float32 for signal and noise
"""
import time
from multiprocessing import Pool
from common import *

M = model()
OPS = [("ref", 1.0, 0.0, 0.0)] + [("m%g_lam%d" % (m, lam), m, 1.0 / m, -lam / m) for m in MS for lam in (1, 0)]
_d = np.load(CACHE / "data.npz"); DV, DR = _d["visco"], _d["ref"]
NOISE = (lambda r: r.standard_normal(DV.shape) + 1j * r.standard_normal(DV.shape))(np.random.default_rng([SEED, 99])) / np.sqrt(2)


def one(j):
    f = FREQ[j]; lap = laplacian(f); n = lap.shape[0]; ns = len(IXS); rec = gidx(IXR, ZSR)
    b = np.zeros((n, 3 * ns), complex); b[gidx(IXS, ZSR), np.arange(ns)] = -ricker(f) / H ** 2
    out = np.zeros((len(OPS), 2, len(IXR), NZ), np.float32)
    for o, (name, m, d, kap) in enumerate(OPS):
        dat = DR[j] if name == "ref" else DV[j]
        b[:, ns:] = 0
        for s in range(ns):
            b[rec, ns + s] = np.conj(dat[s]); b[rec, 2 * ns + s] = np.conj(NOISE[j, s])
        u = interior(factor(lap, M["cs"], M["qis"], f, d, kap).solve(b))[::2]         # (151, NZ, 3 ns)
        S = u[..., :ns]
        out[o, 0] = (np.conj(S) * np.conj(u[..., ns:2 * ns])).real.sum(-1)
        out[o, 1] = (np.conj(S) * np.conj(u[..., 2 * ns:])).real.sum(-1)
    return out


if __name__ == "__main__":
    t0 = time.time()
    with Pool(4) as p:
        out = np.array(p.map(one, range(len(FREQ)), chunksize=1))                     # (nf, ops, 2, nx, nz)
    np.savez(CACHE / "rtm.npz", img=out.transpose(1, 2, 0, 3, 4), ops=[o[0] for o in OPS])
    print("done %.0f s" % (time.time() - t0))
