"""Frequency-domain viscoacoustic shot records (exploratory, synthetic).

D_visco = W(f) [G(true c, true Q) - G(smooth c, smooth Q)]      reflections with attenuation and dispersion
D_ref   = W(f) [G(true c, no Q)   - G(smooth c, no Q)]          the same reflections without attenuation
Arrays (nf, nshot, nrec) complex128 in CACHE/data.npz.
"""
import time
from multiprocessing import Pool
from common import *

M = model()


def one(j):
    f = FREQ[j]; lap = laplacian(f); n = lap.shape[0]
    b = np.zeros((n, len(IXS)), complex); b[gidx(IXS, ZSR), np.arange(len(IXS))] = -1 / H ** 2
    rec = gidx(IXR, ZSR); g = []
    for (c, qi, d, kap) in ((M["c"], M["qi"], 1, 1), (M["cs"], M["qis"], 1, 1), (M["c"], M["qi"], 0, 0), (M["cs"], M["qis"], 0, 0)):
        g.append(factor(lap, c, qi, f, d, kap).solve(b)[rec].T)          # (nshot, nrec)
    w = ricker(f)
    return w * (g[0] - g[1]), w * (g[2] - g[3])


if __name__ == "__main__":
    t0 = time.time()
    with Pool(4) as p:
        out = p.map(one, range(len(FREQ)), chunksize=1)
    np.savez(CACHE / "data.npz", visco=np.array([o[0] for o in out]), ref=np.array([o[1] for o in out]))
    v = np.array([o[0] for o in out]); r = np.array([o[1] for o in out])
    print("done %.0f s; rms visco/ref by frequency:" % (time.time() - t0))
    for j in range(0, len(FREQ), 7):
        print("  f=%5.2f  %.3f" % (FREQ[j], np.sqrt((abs(v[j]) ** 2).mean() / (abs(r[j]) ** 2).mean())))
