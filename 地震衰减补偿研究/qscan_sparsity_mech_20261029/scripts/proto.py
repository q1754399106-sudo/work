"""Prototype (exploratory): does the identifiability of Q by the QScan lasso criterion depend on how sparse the
reflectivity is, while SR does not?

Prediction from the mechanism: with a wrong Q the reflectivity absorbs the attenuation error as a spectral tilt.
For a sparse series the tilt smears every spike and raises the L1 norm, so the criterion has its minimum near the
true Q; for a dense Gaussian series L1 is proportional to L2, measures energy rather than sparsity, and the
minimum is set by the dictionary normalisation, not by the true Q (selection compressed towards a fixed value).

Families, same lateral structure (tau of STRUCT) and same rms as the log reflectivity of experiment 12:
  log    semi.reflectivity, exactly as experiments 8-12
  sparse Bernoulli-Gaussian at 1 ms, density P_SPARSE, laterally identical before the structure (eta = 0)
  gauss  white Gaussian at 1 ms with the eta mix and lateral smoothing of semi.reflectivity
Everything else (source, true t*, noise level calibration, noise line) is gen.signal / run12 of experiment 12.
Methods: QScan objective at weights 0.5 (lib.qscan) and 2^-5, SR (sq.estimate, nboot 0).

usage: proto.py K N    worker K of N; writes proto/<id>.json
"""
import sys, json, time
from env13 import np, lib, semi, sq, bm, sc11, gen, ROOT

STRUCT, NOISE = 10283, [10107, 10235, 10319]          # experiment 12
C = [0.65, 1.05, 1.5]; SEEDS = [9801, 9802]; FAMILIES = ["log", "sparse", "gauss"]
LAMS = [0.5, 2 ** -5]; P_SPARSE = 0.01
OUT = ROOT / "proto"


def datasets():
    return [(f"{fam}_c{c}_s{s}", fam, c, s) for fam in FAMILIES for c in C for s in SEEDS]


def deposit(base, tau):
    """semi.reflectivity's placement of base (len(TF), nx) along the structure tau."""
    nf, nx = base.shape; R = np.zeros((nf, nx)); TF, TW, FINE = semi.TF, semi.TW, semi.FINE
    for j in range(nx):
        pos = (TF + np.interp(TF, TW, tau[:, j]) - TF[0]) / FINE
        i = np.floor(pos).astype(int); fr = pos - i
        for ii, ww in ((i, 1 - fr), (i + 1, fr)):
            ok = (ii >= 0) & (ii < nf); np.add.at(R[:, j], ii[ok], base[ok, j] * ww[ok])
    return R


def reflectivity(fam, tau, rng):
    if fam == "log":
        return semi.reflectivity(tau, rng)[0]
    nf, nx = len(semi.TF), tau.shape[1]; rms = float(np.sqrt(np.mean(semi.log_reflectivity()[1] ** 2)))
    if fam == "sparse":
        r0 = rng.normal(size=nf) * (rng.random(nf) < P_SPARSE); base = np.repeat(r0[:, None], nx, axis=1)
    else:
        r0 = rng.normal(size=nf); ind = rng.normal(size=(nf, nx))
        k = np.arange(-15, 16); g = np.exp(-0.5 * (k / 5.0) ** 2); g /= g.sum()
        ind = np.stack([np.convolve(np.pad(row, 15, mode="reflect"), g, mode="valid") for row in ind])
        ind *= np.sqrt(np.mean(r0 ** 2) / np.mean(ind ** 2))
        base = np.sqrt(1 - semi.ETA ** 2) * r0[:, None] + semi.ETA * ind
    base *= rms / np.sqrt(np.mean(base ** 2))
    return deposit(base, tau)


def signal(ctx, fam, c, seed):
    """gen.signal with the reflectivity family swapped."""
    cm = ctx["cm"]; f = np.fft.rfftfreq(8192, semi.FINE); rng = np.random.default_rng(seed)
    R = reflectivity(fam, ctx["tau"], rng); S = semi.source_spectrum(ctx["w_in"], f, 0.0, cm, c)
    x = semi.model(R, semi.tstar_true(cm, c), S, f, float(cm["fref"]))
    b = ctx["band"]; k = float(np.sqrt(ctx["Creal"][b].sum() / semi.band_power(x, 1.2, 1.6)[b].sum())); x = x * k
    N0, bq = ctx["noise"][0], ctx["bq"]
    err = [float(np.sqrt(np.mean((np.log(np.maximum(lib.coherence_snr(x + L * N0, semi.TW, 2.0, 3.0)[1][bq], 1e-3))
                                  - np.log(np.maximum(ctx["snr_real"][bq], 1e-3))) ** 2))) for L in gen.LEVEL_GRID]
    return x, float(gen.LEVEL_GRID[int(np.argmin(err))])


def scan(D, cm, Q, pp):
    wv, sigma, K, live = pp["wv"], pp["sigma"], pp["K"], pp["live"]; n = D.shape[0]
    J = np.zeros((len(LAMS), len(lib.QGRID))); L1 = np.zeros_like(J); R = np.zeros_like(J); C_ = np.zeros(J.shape, bool)
    s2 = np.where(live, sigma ** 2, np.inf); w = np.sqrt(2 * np.log(n)) * sigma[None, :]
    for j, c in enumerate(lib.QGRID):
        cc = dict(cm); cc["Q"] = np.asarray(Q, float) * c
        B, _, scale = lib.dictionary(cc, wv["w_in"], wv["w_out"])
        for i, lam in enumerate(LAMS):
            X, info = bm.solve(B, np.eye(n), D, sigma, scale, K, lam, 0.0, live)
            res = D - B @ X; pw = w * np.median(scale) / scale[:, None]
            R[i, j] = float(np.sum((0.5 * np.sum(res ** 2, axis=0))[live] / s2[live]))
            L1[i, j] = float(np.sum(np.sum(pw * np.abs(X), axis=0)[live] / s2[live]))
            J[i, j], C_[i, j] = R[i, j] + lam * L1[i, j], bool(info["converged"])
    return J, R, L1, C_


def main(k, n):
    OUT.mkdir(parents=True, exist_ok=True); ctx = gen.prepare(STRUCT, NOISE)
    cm = ctx["cm"]; Q = np.asarray(cm["Q"], float); f = np.fft.rfftfreq(sq.NFFT, sq.DT); wins = sq.windows()
    T = sq.mw_table(cm, Q, f, wins)
    for i, (did, fam, c, seed) in enumerate(datasets()):
        path = OUT / f"{did}.json"
        if i % n != k or path.exists():
            continue
        t0 = time.monotonic(); x, level = signal(ctx, fam, c, seed); D = x + level * ctx["noise"][0]
        sr = sq.estimate(D, cm, Q, table=T, wins=wins, nboot=0, seed=seed)
        pp = sc11.prepare_sc(D, cm); J, R, L1, Cv = scan(D, cm, Q, pp)
        sel = {str(lam): float(lib.QGRID[int(np.argmin(J[i]))]) for i, lam in enumerate(LAMS)}
        rec = dict(id=did, family=fam, c_true=c, seed=seed, level=level, c_sr=float(sr["c"]), selected=sel,
                   grid=lib.QGRID.tolist(), lams=LAMS, J=J.tolist(), residual=R.tolist(), l1=L1.tolist(),
                   converged=Cv.tolist(), seconds=time.monotonic() - t0)
        path.write_text(json.dumps(rec, indent=1), encoding="utf-8"); print(did, f"{rec['seconds']:.0f}s", flush=True)


if __name__ == "__main__":
    main(int(sys.argv[1]), int(sys.argv[2]))
