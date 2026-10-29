"""Experiment 14 runs (PLAN_ZH.md).

usage: run14.py dev K N     the 20 sections of experiment 12 (development set): QScan-W at all weights
       run14.py test K N    the new held-out sets (groups T, L, N): QScan, QScan-W at all weights, SR
Worker K of N, resumable; results/<stage>/<id>.json.

QScan    lib.qscan unchanged (weight 0.5, one input wavelet estimated with the model Q, i.e. factor 1)
QScan-W  the same objective, solver, sigma, K and grid, but the input wavelet re-estimated with Q x c for every
         scanned factor c (lib.wavelets with the scaled Q); FISTA limit 20000 steps, warm start from the solution
         of the previous factor at the same weight.
SR       sq.estimate (experiment 10), nboot 0.
"""
import sys, json, time
from env14 import np, lib, semi, sq, bm, sc11, gen, synth, wq, ROOT, E12OUT

LAMS = [0.5, 2 ** -2, 2 ** -4, 2 ** -6]
MAXITER = 20000
STRUCT, NOISE = 10227, [10149, 10273, 10355]
C_T = [0.6, 0.8, 1.0, 1.2, 1.45]; SEEDS_T = [9901, 9902, 9903, 9904]
C_S = [0.75, 1.2]; WELLS_L = ["F-12", "F-14"]; SEED_L = 9911; SEEDS_N = [9921, 9922]; NOISE_MULT = [0.5, 2.0]
OUT = ROOT / "results"


def test_sets():
    out = [(f"T_c{c}_s{s}", "T", None, 1.0, c, s) for c in C_T for s in SEEDS_T]
    out += [(f"L_{w}_c{c}_s{SEED_L}", "L", w, 1.0, c, SEED_L) for w in WELLS_L for c in C_S]
    out += [(f"N{m}_c{c}_s{s}", "N", None, m, c, s) for m in NOISE_MULT for c in C_S for s in SEEDS_N]
    return out


def dev_sets():
    return sorted(q for q in E12OUT.glob("T_*.npz") if "__" not in q.name)


def qscan_w(D, cm, Q, pp):
    """QScan-W: objective of lib.qscan per weight and factor, input wavelet re-estimated for every factor."""
    sigma, K, live = pp["sigma"], pp["K"], pp["live"]; n = D.shape[0]
    J = np.zeros((len(LAMS), len(lib.QGRID))); C = np.zeros(J.shape, bool); IT = np.zeros(J.shape, int)
    s2 = np.where(live, sigma ** 2, np.inf); prev = [None] * len(LAMS)
    for j, c in enumerate(lib.QGRID):
        cc = dict(cm); cc["Q"] = np.asarray(Q, float) * c
        wv = lib.wavelets(D, semi.TW, cc); B, _, scale = lib.dictionary(cc, wv["w_in"], wv["w_out"])
        pen1 = np.sqrt(2 * np.log(n)) * sigma[None, :] * np.median(scale) / scale[:, None]
        for i, lam in enumerate(LAMS):
            X, info = bm.solve(B, np.eye(n), D, sigma, scale, K, lam, 0.0, live, start=prev[i]); prev[i] = X
            res = D - B @ X
            J[i, j] = float(np.sum((0.5 * np.sum(res ** 2, axis=0) + lam * np.sum(pen1 * np.abs(X), axis=0))[live] / s2[live]))
            C[i, j], IT[i, j] = bool(info["converged"]), int(info["iterations"])
    sel = {str(lam): float(lib.QGRID[int(np.argmin(J[i]))]) for i, lam in enumerate(LAMS)}
    return sel, dict(J=J.tolist(), converged=C.tolist(), iterations=IT.tolist())


def signal_L(ctx, w, c, seed):
    """run11.signal_L (experiment 11): the long Gardner log of well w at the centre trace."""
    cm = ctx["cm"]; f = np.fft.rfftfreq(8192, semi.FINE); rng = np.random.default_rng(seed)
    R = synth.reflectivity(ctx, w, rng); S = semi.source_spectrum(ctx["w_in"], f, 0.0, cm, c)
    x = semi.model(R, semi.tstar_true(cm, c), S, f, float(cm["fref"]))
    b = ctx["band"]; k = float(np.sqrt(ctx["Creal"][b].sum() / semi.band_power(x, 1.2, 1.6)[b].sum())); x = x * k
    N0, bq = ctx["noise"][0], ctx["bq"]
    err = [float(np.sqrt(np.mean((np.log(np.maximum(lib.coherence_snr(x + L * N0, semi.TW, 2.0, 3.0)[1][bq], 1e-3))
                                  - np.log(np.maximum(ctx["snr_real"][bq], 1e-3))) ** 2))) for L in gen.LEVEL_GRID]
    return dict(x=x, level=float(gen.LEVEL_GRID[int(np.argmin(err))]))


def dev(k, n):
    out = OUT / "dev"; out.mkdir(parents=True, exist_ok=True); cm = semi.curtain_cm(); Q = np.asarray(cm["Q"], float)
    for i, p in enumerate(dev_sets()):
        did = p.stem; path = out / f"{did}.json"
        if i % n != k or path.exists():
            continue
        t0 = time.monotonic(); D = np.load(p)["Input"].astype(float); pp = sc11.prepare_sc(D, cm)
        sel, diag = qscan_w(D, cm, Q, pp)
        rec = dict(id=did, c_true=float(did.split("_")[1][1:]), lams=LAMS, qscan_w=sel, diag=diag, seconds=time.monotonic() - t0)
        path.write_text(json.dumps(rec, indent=1), encoding="utf-8"); print(did, f"{rec['seconds']:.0f}s", flush=True)


def test(k, n):
    out = OUT / "test"; out.mkdir(parents=True, exist_ok=True); ctx = gen.prepare(STRUCT, NOISE)
    cm = ctx["cm"]; Q = np.asarray(cm["Q"], float); f = np.fft.rfftfreq(sq.NFFT, sq.DT); wins = sq.windows()
    T = sq.mw_table(cm, Q, f, wins); P = wq.picks(); W = {w: wq.well(w, P) for w in WELLS_L}
    for i, (did, grp, well, mult, c, seed) in enumerate(test_sets()):
        path = out / f"{did}.json"
        if i % n != k or path.exists():
            continue
        t0 = time.monotonic()
        g = signal_L(ctx, W[well], c, seed) if grp == "L" else gen.signal(ctx, c, seed)
        D = g["x"] + mult * g["level"] * ctx["noise"][0]
        sr = sq.estimate(D, cm, Q, table=T, wins=wins, nboot=0, seed=seed)
        pp = sc11.prepare_sc(D, cm); fq, qi = sc11.qscan(D, cm, Q, pp)
        sel, diag = qscan_w(D, cm, Q, pp)
        rec = dict(id=did, group=grp, well=well, noise_mult=mult, c_true=c, seed=seed, level=g["level"], lines=ctx["lines"],
                   c_sr=float(sr["c"]), qscan=fq, qscan_boundary=bool(qi["at_boundary"]), lams=LAMS, qscan_w=sel,
                   diag=diag, seconds=time.monotonic() - t0)
        path.write_text(json.dumps(rec, indent=1), encoding="utf-8"); print(did, f"{rec['seconds']:.0f}s", flush=True)


bm.MAXITER = MAXITER

if __name__ == "__main__":
    (dev if sys.argv[1] == "dev" else test)(int(sys.argv[2]), int(sys.argv[3]))
