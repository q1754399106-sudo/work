"""Two-layer anchoring test (PLAN_LAYER_ZH.md): QScan with separate Q factors above / below the model boundary.

usage: layered.py [--smoke]    worker; jobs claimed with results6/<job>.claim, resumable.
Datasets: experiment-14 recipe (gen.signal on structure line 10227), c 0.7 / 1.0 / 1.4, new seeds 9961-9964,
noise x 0.5 (true Q proportional: u1 = u2 = u). Input wavelet: lib.wavelets (window 1.2-1.6 s, absorption at 1.4 s)
at the fixed anchor c_a = 0.7 / 1.0 / 1.4. Scanned model: Q x [c1, c1, c2, c2], i.e. c1 above the boundary at
2.041 s (layers 1-2) and c2 below (layers 3-4). Lasso objective of lib.qscan at lambda 1/2 on a 2-D grid of
v = 1/c (the QGRID ratio, extended to c 0.30-4.0); 3 x 3 pattern search from (1, 1) until the centre is the
minimum, then the vertex of a quadratic fitted to the final 3 x 3 block (per-axis parabolas as fallback).
-> results6/<dataset>__a<c_a>.json
"""
import os, sys, json, time
from env15 import np, lib, semi, sc11, bm, gen, ROOT

STRUCT, NOISE = 10227, [10149, 10273, 10355]
MULT, LAM, MAXITER = 0.5, 0.5, 20000
ANCHORS = (0.7, 1.0, 1.4); SEEDS = (9961, 9962, 9963, 9964)
OUT = ROOT / "results6"; bm.MAXITER = MAXITER
R = float(lib.QGRID[1] / lib.QGRID[0])
CG = 1.0 * R ** np.arange(-10, 11); CG = CG[(CG >= 0.299) & (CG <= 4.01)]   # c grid, includes 1
VG = 1.0 / CG


def jobs():
    return [(f"c{c}_s{s}", c, s, a) for s in SEEDS for c in (0.7, 1.0, 1.4) for a in ANCHORS]


def search(D, cm, Q, pp, wv):
    sigma, K, live = pp["sigma"], pp["K"], pp["live"]; n = D.shape[0]; s2 = np.where(live, sigma ** 2, np.inf)
    J, X, conv = {}, {}, {}

    def ev(i, j):
        if (i, j) in J:
            return J[(i, j)]
        cc = dict(cm); cc["Q"] = Q * np.array([CG[i], CG[i], CG[j], CG[j]])
        B, _, scale = lib.dictionary(cc, wv["w_in"], wv["w_out"])
        pen1 = np.sqrt(2 * np.log(n)) * sigma[None, :] * np.median(scale) / scale[:, None]
        near = min(X, key=lambda k: abs(k[0] - i) + abs(k[1] - j)) if X else None
        Xs, info = bm.solve(B, np.eye(n), D, sigma, scale, K, LAM, 0.0, live, start=X[near] if near else None)
        res = D - B @ Xs; X[(i, j)] = Xs; conv[(i, j)] = bool(info["converged"])
        J[(i, j)] = float(np.sum((0.5 * np.sum(res ** 2, axis=0) + LAM * np.sum(pen1 * np.abs(Xs), axis=0))[live]
                                 / s2[live]))
        return J[(i, j)]

    m = len(CG); c = (int(np.argmin(np.abs(CG - 1))),) * 2; moves = 0
    while True:
        blk = [(i, j) for i in range(c[0] - 1, c[0] + 2) for j in range(c[1] - 1, c[1] + 2) if 0 <= i < m and 0 <= j < m]
        best = min(blk, key=lambda p: ev(*p))
        if best == c or moves > 40:
            break
        c = best; moves += 1
    edge = c[0] in (0, m - 1) or c[1] in (0, m - 1)
    v = [float(VG[c[0]]), float(VG[c[1]])]; how = "grid"
    if not edge:
        pts = [(i, j) for i in range(c[0] - 1, c[0] + 2) for j in range(c[1] - 1, c[1] + 2)]
        x = np.array([VG[i] - VG[c[0]] for i, _ in pts]); y = np.array([VG[j] - VG[c[1]] for _, j in pts])
        z = np.array([J[q] for q in pts]); M = np.c_[np.ones(9), x, y, x * x, x * y, y * y]
        p = np.linalg.lstsq(M, z, rcond=None)[0]; H = np.array([[2 * p[3], p[4]], [p[4], 2 * p[5]]])
        lo = [VG[c[0] + 1] - VG[c[0]], VG[c[1] + 1] - VG[c[1]]]; hi = [VG[c[0] - 1] - VG[c[0]], VG[c[1] - 1] - VG[c[1]]]
        if np.all(np.linalg.eigvalsh(H) > 0):
            d = -np.linalg.solve(H, p[1:3])
            if all(lo[k] <= d[k] <= hi[k] for k in range(2)):
                v = [float(VG[c[0]] + d[0]), float(VG[c[1]] + d[1])]; how = "quad2d"
        if how == "grid":
            for k in range(2):
                idx = [c[k] + 1, c[k], c[k] - 1]
                zz = [J[(i, c[1])] if k == 0 else J[(c[0], i)] for i in idx]
                a = np.polyfit(VG[idx], zz, 2)
                if a[0] > 0:
                    v[k] = float(np.clip(-a[1] / (2 * a[0]), VG[idx[0]], VG[idx[2]]))
            how = "axis"
    trust = all(conv.get((i, j), True) for i in range(c[0] - 1, c[0] + 2) for j in range(c[1] - 1, c[1] + 2))
    pts = [dict(v1=float(VG[i]), v2=float(VG[j]), J=J[(i, j)], converged=conv[(i, j)]) for (i, j) in sorted(J)]
    return dict(v1=v[0], v2=v[1], vertex=how, edge=bool(edge), moves=moves, neighbourhood_converged=bool(trust),
                n_evals=len(J), points=pts)


def run(did, c, s, ca, ctx):
    cm = ctx["cm"]; Q = np.asarray(cm["Q"], float)
    g = gen.signal(ctx, c, s); D = g["x"] + MULT * g["level"] * ctx["noise"][0]; pp = sc11.prepare_sc(D, cm)
    wv = lib.wavelets(D, semi.TW, dict(cm, Q=Q * ca)); t0 = time.monotonic(); res = search(D, cm, Q, pp, wv)
    return dict(id=did, c_true=c, u=1.0 / c, seed=s, c_anchor=ca, v_anchor=1.0 / ca, noise_mult=MULT, lam=LAM,
                level=g["level"], seconds=time.monotonic() - t0, **res)


if __name__ == "__main__":
    OUT.mkdir(exist_ok=True); ctx = gen.prepare(STRUCT, NOISE)
    if "--smoke" in sys.argv:          # development seed 9940, not a test dataset; timing and search check only
        rec = run("c1.0_s9940", 1.0, 9940, 1.0, ctx); rec.pop("points")
        print(json.dumps(rec, indent=1)); sys.exit()
    for did, c, s, ca in jobs():
        job = f"{did}__a{ca}"; path = OUT / f"{job}.json"
        if path.exists():
            continue
        try:
            os.close(os.open(OUT / f"{job}.claim", os.O_CREAT | os.O_EXCL | os.O_WRONLY))
        except FileExistsError:
            continue
        rec = run(did, c, s, ca, ctx)
        path.write_text(json.dumps(rec, indent=1), encoding="utf-8"); (OUT / f"{job}.claim").unlink()
        print(job, f"v1 {rec['v1']:.3f} v2 {rec['v2']:.3f}", rec["vertex"], f"{rec['n_evals']} evals",
              f"{rec['seconds']:.0f}s", flush=True)
