"""Ceiling check on used models: sat and split at every lambda on the grid, oracle per cell, against the band-limited inverse Q (k = 3). Descriptive.
usage: dev_ceiling.py M -> results/dev_ceiling_<M>.json"""
import json, sys, time
import numpy as np
import step2b as s
from step2b import NEW, columns, metrics, spec, LX, GAMMA
import gate1 as g
from gate1 import n, NX, I, a0_ax, atten
import theory as th, stable, satdamp as sd, bandlim as bl

NEW.update(M1=np.array([40.0, 50, 70, 100]), M10=np.array([45.0, 70, 55, 85]), M11=np.array([30.0, 30, 60, 100])); s.SEEDS = 5; FP = 60.0; g.RIC = spec(FP); th.S2 = g.RIC ** 2
GRID = np.logspace(-6, 0, 13)[::-1]; KEYS = ("ga", "eg", "egl", "mg", "cg", "wave"); CELLS = [(sc, qf) for sc in ("events_aligned", "thin_aligned", "dense_aligned", "dense_nmo") for qf in (1.0, 0.8, 1.2)]


def main(M):
    rows = []; p = s.OUT / f"dev_ceiling_{M}.json"; t0 = time.monotonic()
    for sc, qf in CELLS:
        Q = NEW[M] * qf; D, Cl, Tr, sig, meta, keep = columns(M, sc, FP); nc = D.shape[2]; a0, ax = a0_ax(Q); G0 = atten(a0); Gx = [atten(ax[:, i]) for i in range(NX)]; Gd = [atten(a0 + ax[:, i]) for i in range(NX)]; GtG = [G.T @ G for G in Gd]; GtD = [Gd[i].T @ D[:, i, :] for i in range(NX)]; Y = np.stack([np.linalg.solve(G.T @ G + LX * I, G.T @ D[:, i, :]) for i, G in enumerate(Gx)], axis=1); target = np.sum(sig ** 2, axis=0) * n; lab = dict(M=M, scene=sc, qf=qf)
        def rec(m, lam, O, extra=None):
            mt = metrics(O, O, Tr, keep, FP); rows.append(dict(lab, m=m, lam=float(lam), **{k: [float(x) for x in mt[k]] for k in KEYS}, **(extra or {})))
        for lam in GRID:
            wn = stable.bound(a0, ax, lam, GAMMA); O = np.zeros_like(D)
            for i in range(NX):
                L = sd.Lmat(ax[:, i], sd.fb_floor(np.full(n, sd.FF[-1]), ax[:, i], wn[:, i])); O[:, i, :] = np.linalg.solve(GtG[i] + lam * (L.T @ L), GtD[i])
            res = sum(np.sum(((Gd[i] @ O[:, i, :] - D[:, i, :]) * keep[:, i, None]) ** 2, axis=0) for i in range(NX)); rec("sat", lam, O, dict(res=(res / target).tolist()))
            F0 = np.linalg.solve(G0.T @ G0 + lam * I, G0.T); rec("split", lam, np.stack([F0 @ Y[:, i, :] for i in range(NX)], axis=1))
        rec("band3", 0.0, np.stack([bl.apply(D[:, :, c], Q, sig[:, c], "snr", 3.0) for c in range(nc)], 2))
        print(M, sc, qf, f"{time.monotonic() - t0:.0f}s", flush=True); p.write_text(json.dumps(rows), encoding="utf-8")


if __name__ == "__main__":
    main(sys.argv[1])
