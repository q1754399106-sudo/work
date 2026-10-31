"""Iterated re-anchoring test (PLAN_ITER_ZH.md): feed the selected Q factor back as the wavelet anchor.

usage: iterate.py [--smoke]    worker; jobs claimed with results5/<job>.claim, resumable.
Datasets: experiment-14 recipe (gen.signal on structure line 10227), c 0.7 / 1.0 / 1.4, noise x 0.5, full analysis
window 1.2-3.0 s; seeds 9941-9944 (the datasets of the analysis-window test trunc.py) and new seeds 9951-9954.
Step k: estimate the input wavelet at anchor v_a = v_{k-1} (v_0 = 1, the current pipeline), evaluate the lasso
objective of lib.qscan on grid points around the expected minimum, and take the vertex of the parabola
in v = 1/c through the discrete minimum and its two neighbours -> v_k (continuous).
Arms (lambda 1/2 unless noted; all 24 datasets unless noted):
  "main"   wavelet window 1.2-1.6 s, absorption at 1.4 s (kappa 0.675), 8 steps;
  "mid"    window 1.8-2.2 s, absorption at 2.0 s (kappa 0.96), 8 steps;
  "deep"   window 2.4-2.8 s, absorption at 2.6 s (kappa 1.14), 6 steps;
  "main16" as main at lambda 1/16, 6 steps, seeds 9941-9944 only (secondary, run last).
-> results5/<dataset>__<arm>.json
"""
import os, sys, json, time
from env15 import np, lib, semi, sc11, bm, gen, ROOT
from anchor2 import wavelets_at

STRUCT, NOISE = 10227, [10149, 10273, 10355]
MULT, MAXITER = 0.5, 20000
ARMS = {"main": dict(win=(1.2, 1.6, 1.4), steps=8, lam=0.5), "mid": dict(win=(1.8, 2.2, 2.0), steps=8, lam=0.5),
        "deep": dict(win=(2.4, 2.8, 2.6), steps=6, lam=0.5), "main16": dict(win=(1.2, 1.6, 1.4), steps=6, lam=2 ** -4)}
SEEDS = (9941, 9942, 9943, 9944, 9951, 9952, 9953, 9954)
OUT = ROOT / "results5"; bm.MAXITER = MAXITER
VG = 1.0 / np.asarray(lib.QGRID, float)          # v on the grid (decreasing with index)


def jobs():
    out = [(f"c{c}_s{s}", c, s, arm) for arm in ("main", "mid", "deep") for s in SEEDS for c in (0.7, 1.0, 1.4)]
    return out + [(f"c{c}_s{s}", c, s, "main16") for s in SEEDS[:4] for c in (0.7, 1.0, 1.4)]


def objective(D, cm, Q, pp, wv, j, start, lam):
    sigma, K, live = pp["sigma"], pp["K"], pp["live"]; n = D.shape[0]
    cc = dict(cm); cc["Q"] = Q * lib.QGRID[j]
    B, _, scale = lib.dictionary(cc, wv["w_in"], wv["w_out"])
    pen1 = np.sqrt(2 * np.log(n)) * sigma[None, :] * np.median(scale) / scale[:, None]
    X, info = bm.solve(B, np.eye(n), D, sigma, scale, K, lam, 0.0, live, start=start)
    res = D - B @ X; s2 = np.where(live, sigma ** 2, np.inf)
    J = float(np.sum((0.5 * np.sum(res ** 2, axis=0) + lam * np.sum(pen1 * np.abs(X), axis=0))[live] / s2[live]))
    return J, X, bool(info["converged"]), int(info["iterations"])


def vertex(v3, J3):
    """Vertex of the parabola through three (v, J) points, clipped to their range."""
    a = np.polyfit(v3, J3, 2)
    if a[0] <= 0:
        return float(v3[int(np.argmin(J3))])
    return float(np.clip(-a[1] / (2 * a[0]), min(v3), max(v3)))


def step(D, cm, Q, pp, wv, v_guess, warm, lam):
    """Local scan around v_guess; extend until the discrete minimum is interior (or at the grid edge)."""
    J, conv, its = {}, {}, {}
    j0 = int(np.argmin(np.abs(VG - v_guess))); todo = list(range(max(0, j0 - 2), min(len(VG), j0 + 3)))
    while True:
        for j in todo:
            if j in J:
                continue
            nb = [k for k in (j - 1, j + 1) if k in warm]
            start = warm[j] if j in warm else (warm[nb[0]] if nb else None)
            J[j], warm[j], conv[j], its[j] = objective(D, cm, Q, pp, wv, j, start, lam)
        jm = min(J, key=J.get)
        if jm in (0, len(VG) - 1):
            edge = True; break
        need = [k for k in (jm - 1, jm + 1) if k not in J]
        if not need:
            edge = False; break
        todo = need
    if edge:
        v = float(VG[jm])
    else:
        idx = [jm - 1, jm, jm + 1]; v = vertex(VG[idx], np.array([J[k] for k in idx]))
    local = {int(k): dict(v=float(VG[k]), J=J[k], converged=conv[k], iterations=its[k]) for k in sorted(J)}
    trust = all(conv[k] for k in (jm - 1, jm, jm + 1) if k in conv)
    return v, dict(v_grid=float(VG[jm]), edge=edge, neighbourhood_converged=trust, points=local)


def run(did, c, s, arm, ctx):
    cm = ctx["cm"]; Q = np.asarray(cm["Q"], float); cfg = ARMS[arm]
    g = gen.signal(ctx, c, s); D = g["x"] + MULT * g["level"] * ctx["noise"][0]; pp = sc11.prepare_sc(D, cm)
    v_a, hist, warm = 1.0, [], {}
    for k in range(1, cfg["steps"] + 1):
        t0 = time.monotonic(); wv = wavelets_at(D, semi.TW, dict(cm, Q=Q / v_a), *cfg["win"])
        v, diag = step(D, cm, Q, pp, wv, v_a, warm, cfg["lam"])
        hist.append(dict(k=k, v_anchor=v_a, v_sel=v, seconds=time.monotonic() - t0, **diag))
        print(did, arm, k, f"v_a {v_a:.4f} -> v {v:.4f}", f"{hist[-1]['seconds']:.0f}s", flush=True)
        if diag["edge"]:
            break
        v_a = v
    return dict(id=did, arm=arm, c_true=c, u=1.0 / c, seed=s, noise_mult=MULT, lam=cfg["lam"], window=list(cfg["win"]),
                level=g["level"], history=hist)


if __name__ == "__main__":
    OUT.mkdir(exist_ok=True); ctx = gen.prepare(STRUCT, NOISE)
    if "--smoke" in sys.argv:          # development seed 9940, not a test dataset; one step, timing only
        ARMS["main"]["steps"] = 1; t0 = time.monotonic(); rec = run("c1.0_s9940", 1.0, 9940, "main", ctx)
        print(json.dumps(rec["history"][0], indent=1)[:1500], f"total {time.monotonic() - t0:.0f}s"); sys.exit()
    for did, c, s, arm in jobs():
        job = f"{did}__{arm}"; path = OUT / f"{job}.json"
        if path.exists():
            continue
        try:
            os.close(os.open(OUT / f"{job}.claim", os.O_CREAT | os.O_EXCL | os.O_WRONLY))
        except FileExistsError:
            continue
        t0 = time.monotonic(); rec = run(did, c, s, arm, ctx); rec["seconds"] = time.monotonic() - t0
        path.write_text(json.dumps(rec, indent=1), encoding="utf-8"); (OUT / f"{job}.claim").unlink()
        print(job, "done", f"{rec['seconds']:.0f}s", flush=True)
