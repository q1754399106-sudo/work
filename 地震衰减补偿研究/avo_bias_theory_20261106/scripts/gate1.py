"""Gate 1 (GATE1_ZH.md): does the split's AVO advantage survive NMO stretch, dense reflectivity, thin beds and an
independent forward model?  usage: gate1.py smoke | check | run M   (M = M1, M2, M3)  -> results/gate1_<M>.json
Forward: reflector-by-reflector synthesis in the frequency domain. Inverse: the original column-varying matrices."""
import json, sys, time
from pathlib import Path
import numpy as np
from core import C, atten, ricker

OUT = Path(__file__).resolve().parents[1] / "results"; OUT.mkdir(exist_ok=True)
n = C.nt; T = C.t; DT = C.dt; OFF = C.offset; NX = len(OFF); I = np.eye(n)
LAMS = np.logspace(-5, 0, 9); LX = 1e-3; FC = (40.0, 60.0, 80.0); GRID = np.logspace(-7, 0, 41)
MODELS = {"M1": np.array([40.0, 50, 70, 100]), "M2": 0.7 * np.array([40.0, 50, 70, 100]), "M3": np.array([100.0, 70, 50, 40])}
SCENES = ["events_aligned", "thin_aligned", "dense_aligned", "events_nmo", "thin_nmo", "dense_nmo"]
SNRS = (18, 12); QF = (0.8, 1.0, 1.2); SEEDS = 10; MUTE = 1.3; BAND = (10.0, 60.0); WIN = (0.25, 1.10); NSYN = 4096
SYN_T = NSYN * DT
FS = np.fft.rfftfreq(NSYN, DT)
RIC = np.abs(np.fft.rfft(np.roll(ricker(np.arange(NSYN) * DT - SYN_T / 2), -NSYN // 2)))   # zero-phase amplitude spectrum
BOUND = np.r_[0, np.cumsum(2 * C.h / C.vp[:3]), np.inf]


def geometry(t0s, x):
    """tau (nt, 4) two-way vertical time per layer, p, sqrt term; exact ray through layers to reflector at t0."""
    t0s = np.atleast_1d(t0s); tau = np.maximum(0, np.minimum(t0s[:, None], BOUND[1:]) - BOUND[:-1]); tau[~np.isfinite(tau)] = 0
    p = np.zeros(len(t0s))
    if x > 0:
        vmax = np.array([C.vp[tau[k] > 0].max() if (tau[k] > 0).any() else C.vp[0] for k in range(len(t0s))])
        lo = np.zeros(len(t0s)); hi = (1 - 1e-12) / vmax
        for _ in range(60):
            pm = 0.5 * (lo + hi); xx = np.sum(tau * C.vp ** 2 * pm[:, None] / np.sqrt(np.maximum(1 - (pm[:, None] * C.vp) ** 2, 1e-12)), axis=1)
            lo = np.where(xx < x, pm, lo); hi = np.where(xx < x, hi, pm)
        p = 0.5 * (lo + hi)
    root = np.sqrt(np.maximum(1 - (p[:, None] * C.vp) ** 2, 1e-12))
    return tau, p, root


GEO = {}
for i, x in enumerate(OFF):
    tau, p, root = geometry(T, x); GEO[i] = dict(tau=tau, p=p, root=root, tx=np.sum(tau / root, axis=1))


def a0_ax(Q):
    a0 = GEO[0]["tau"] @ (1 / Q); ax = np.stack([np.sum(GEO[i]["tau"] / (Q * GEO[i]["root"]), axis=1) - a0 for i in range(NX)], axis=1)
    return a0, ax


VLOC = np.array([C.vp[int(np.searchsorted(BOUND[1:], t, side="right"))] if t < BOUND[-2] else C.vp[3] for t in T])
SIN2 = np.stack([(GEO[i]["p"] * VLOC) ** 2 for i in range(NX)], axis=1)


def rpp(j, sin_i):
    """Exact Zoeppritz P-P reflection at interface j (between layers j and j+1) for incidence sin_i."""
    a1, a2, b1, b2, r1, r2 = C.vp[j], C.vp[j + 1], C.vs[j], C.vs[j + 1], C.rho[j], C.rho[j + 1]
    sin_i = np.atleast_1d(sin_i); out = np.zeros(len(sin_i))
    for k, s in enumerate(sin_i):
        p = s / a1; c1 = np.sqrt(1 - s ** 2 + 0j); c2 = np.sqrt(1 - (p * a2) ** 2 + 0j); d1 = np.sqrt(1 - (p * b1) ** 2 + 0j); d2 = np.sqrt(1 - (p * b2) ** 2 + 0j)
        A = np.array([[-p * a1, -d1, p * a2, d2], [c1, -p * b1, c2, -p * b2],
                      [2 * r1 * b1 ** 2 * p * c1, r1 * b1 * (1 - 2 * b1 ** 2 * p ** 2), 2 * r2 * b2 ** 2 * p * c2, r2 * b2 * (1 - 2 * b2 ** 2 * p ** 2)],
                      [-r1 * a1 * (1 - 2 * b1 ** 2 * p ** 2), 2 * r1 * b1 ** 2 * p * d1, r2 * a2 * (1 - 2 * b2 ** 2 * p ** 2), -2 * r2 * b2 ** 2 * p * d2]], complex)
        b = np.array([p * a1, c1, 2 * r1 * b1 ** 2 * p * c1, r1 * a1 * (1 - 2 * b1 ** 2 * p ** 2)], complex)
        out[k] = np.linalg.solve(A, b)[0].real
    return out


def reflectors(scene, rng):
    """List: times t0 (on the grid), R (n_ref, NX)."""
    kind = scene.split("_")[0]; ts, rs = [], []
    ei = np.round(C.event_t / DT).astype(int)
    for j, i0 in enumerate(ei):
        sin_i = np.array([GEO[i]["p"][i0] * C.vp[j] for i in range(NX)]); r = rpp(j, sin_i)
        if kind == "thin":
            ts += [i0, i0 + 8]; rs += [r, -r]
        else:
            ts.append(i0); rs.append(r)
    if kind == "dense":
        idx = np.where(rng.random(n) < 0.3)[0]; idx = idx[(idx > 10) & (idx < n - 10)]
        A = rng.normal(0, 0.05, len(idx)); B = -0.8 * A + rng.normal(0, 0.03, len(idx))
        ts += list(idx); rs += list(A[:, None] + B[:, None] * SIN2[idx])
    return np.array(ts), np.array(rs)


def synth(idx, R, Q, nmo):
    """Attenuated gather (pre-NMO if nmo) and the clean truth at t0."""
    a0, ax = a0_ax(Q); b = -np.pi * FS + 2j * FS * np.log(np.maximum(FS, 1e-300) / C.fref); b[0] = 0
    d = np.zeros((NSYN, NX)); tru = np.zeros((n, NX)); e0 = np.exp(-2j * np.pi * np.outer(FS, T[idx]))
    for i in range(NX):
        tk = GEO[i]["tx"][idx] if nmo else T[idx]; a = a0[idx] + ax[idx, i]
        X = RIC[:, None] * np.exp(np.outer(b, a)) * np.exp(-2j * np.pi * np.outer(FS, tk)); d[:, i] = np.fft.irfft(X @ R[:, i], NSYN)
        tru[:, i] = np.fft.irfft((RIC[:, None] * e0) @ R[:, i], NSYN)[:n]
    return d, tru


def to_t0(d, nmo):
    """NMO to t0 and stretch mute (or plain crop)."""
    if not nmo:
        return d[:n], np.ones((n, NX), bool)
    M = 4; Nu = M * NSYN; X = np.fft.rfft(d, axis=0); Xu = np.zeros((Nu // 2 + 1, NX), complex); Xu[:X.shape[0]] = X; du = np.fft.irfft(Xu, Nu, axis=0) * M; tu = np.arange(Nu) * DT / M
    out = np.zeros((n, NX)); keep = np.ones((n, NX), bool)
    for i in range(NX):
        tx = GEO[i]["tx"]; out[:, i] = np.interp(tx, tu, du[:, i], right=0.0); keep[:, i] = tx <= MUTE * np.maximum(T, DT)
    return out * keep, keep


def bandpass(X):
    F = np.fft.rfft(X, 2048, axis=0); f = np.fft.rfftfreq(2048, DT)
    g = np.clip((f - BAND[0] + 5) / 5, 0, 1) * np.clip((BAND[1] + 8 - f) / 8, 0, 1)
    shape = (-1,) + (1,) * (X.ndim - 1)
    return np.fft.irfft(F * g.reshape(shape), 2048, axis=0)[:n]


TW = (T >= WIN[0]) & (T <= WIN[1])


def avo(Y, keep):
    """Per-sample fit y = A + B sin^2 over kept offsets; Y (n, NX, ...). A, B and the valid rows."""
    W = keep.astype(float); S = SIN2 * W; n0 = W.sum(axis=1); sx = S.sum(axis=1); sxx = (SIN2 * S).sum(axis=1); det = n0 * sxx - sx ** 2
    ok = (n0 >= 4) & (det > 1e-12) & TW
    ext = (slice(None), slice(None)) + (None,) * (Y.ndim - 2)
    sy = (Y * W[ext]).sum(axis=1); sxy = (Y * S[ext]).sum(axis=1); ds = np.where(ok, det, 1.0)
    e = (slice(None),) + (None,) * (Y.ndim - 2)
    return (sxx[e] * sy - sx[e] * sxy) / ds[e], (n0[e] * sxy - sx[e] * sy) / ds[e], ok


def score_batch(O, Tru, keep, Atrue, Btrue, ok):
    """O (n, NX, ncol), Tru (n, NX, ncol): metrics per column."""
    mask = (keep & TW[:, None])[:, :, None]; wave = np.sqrt(np.sum(((O - Tru) * mask) ** 2, axis=(0, 1)) / np.sum((Tru * mask) ** 2, axis=(0, 1)))
    A, B, _ = avo(bandpass(O), keep); r = ok[:, None]
    ba = np.sum(A * Atrue * r, axis=0) / np.sum(Atrue ** 2 * r, axis=0); bg = np.sum(B * Btrue * r, axis=0) / np.sum(Btrue ** 2 * r, axis=0)
    cg = np.sum(B * Btrue * r, axis=0) / np.sqrt(np.sum(B ** 2 * r, axis=0) * np.sum(Btrue ** 2 * r, axis=0))
    return dict(wave=wave, mg=np.abs(bg / ba - 1), ma=np.abs(ba - 1), mb=np.abs(bg - 1), cg=cg, ba=ba, bg=bg)


def build_data(Mname, scene):
    """All columns of a scene: (snr, seed). Data and truth per column."""
    Q = MODELS[Mname]; nmo = scene.endswith("nmo"); cols = []
    for seed in range(SEEDS):
        rng = np.random.default_rng(10000 * list(MODELS).index(Mname) + 100 * SCENES.index(scene) + seed); idx, R = reflectors(scene, rng); dp, tru = synth(idx, R, Q, nmo)
        rms = np.sqrt(np.sum(dp[:int(2.4 / DT)] ** 2, axis=0) / n)
        for snr in SNRS:
            noise = rms / 10 ** (snr / 20) * rng.standard_normal(dp.shape); y, keep = to_t0(dp + noise, nmo); cols.append((snr, seed, y, tru, keep, rms / 10 ** (snr / 20)))
    return cols


def solve_all(G, w, lam, D):
    return np.linalg.solve(G.T @ G + lam * np.diag(w), G.T @ D)


def disc_batch(G, Y, target):
    U, s, Vt = np.linalg.svd(G, full_matrices=False); b = U.T @ Y; res0 = np.sum(Y ** 2, axis=0) - np.sum(b ** 2, axis=0); lam = np.full(Y.shape[1], GRID[0])
    fl = (GRID[:, None] / (s[None, :] ** 2 + GRID[:, None]))
    for j in range(Y.shape[1]):
        r = np.sum((fl * b[:, j]) ** 2, axis=1) + res0[j]; g = np.where(r <= target[j] ** 2)[0]; lam[j] = GRID[g[-1]] if len(g) else GRID[0]
    return Vt.T @ (s[:, None] / (s[:, None] ** 2 + lam) * b)


def run_scene(Mname, qf, scene, rows, cols):
    Q = MODELS[Mname] * qf; a0, ax = a0_ax(Q); G0 = atten(a0); Gx = [atten(ax[:, i]) for i in range(NX)]; Gd = [atten(a0 + ax[:, i]) for i in range(NX)]
    keep = cols[0][4]; ncol = len(cols); Dat = np.stack([c[2] for c in cols], axis=2); Tru = np.stack([c[3] for c in cols], axis=2); snr = np.array([c[0] for c in cols]); seed = np.array([c[1] for c in cols])
    At, Bt, ok = avo(bandpass(Tru), keep)
    Bx = [np.linalg.solve(g.T @ g + LX * I, g.T) for g in Gx]; Y1 = np.stack([Bx[i] @ Dat[:, i, :] for i in range(NX)], axis=1)

    def rec(m, lam, O):
        sc = score_batch(O, Tru, keep, At, Bt, ok)
        for c in range(ncol):
            rows.append(dict(M=Mname, qf=qf, scene=scene, snr=int(snr[c]), seed=int(seed[c]), m=m, lam=float(lam), **{k: float(v[c]) for k, v in sc.items()}))

    rec("observed", 0.0, Dat)
    for lam in LAMS:
        O = np.zeros_like(Dat)
        F0 = np.linalg.solve(G0.T @ G0 + lam * I, G0.T)
        for i in range(NX):
            O[:, i, :] = F0 @ Y1[:, i, :]
        rec("split", lam, O)
        for name, fc in (("direct", None),) + tuple((f"consistent{f:.0f}", f) for f in FC):
            for i in range(NX):
                w = np.ones(n) if fc is None else np.exp(-2 * np.pi * fc * ax[:, i]); O[:, i, :] = solve_all(Gd[i], w, lam, Dat[:, i, :])
            rec(name, lam, O)
    O1 = np.zeros_like(Dat); O2 = np.zeros_like(Dat)
    for i in range(NX):
        target = np.array([c[5][i] for c in cols]) * np.sqrt(n); O1[:, i, :] = disc_batch(Gd[i], Dat[:, i, :], target)
        L = np.linalg.cholesky(Bx[i] @ Bx[i].T + 1e-6 * I); O2[:, i, :] = disc_batch(np.linalg.solve(L, G0), np.linalg.solve(L, Y1[:, i, :]), target)
    rec("direct_disc", 0.0, O1); rec("covsplit_disc", 0.0, O2)


def check():
    """orig scene (K1b reproduction): the inverse matrices as forward, template G/A of the earlier pilot. Already-seen numbers."""
    from core import reference
    truth = np.stack([reference(x)[0] for x in OFF], axis=1); S2 = np.stack([reference(x)[2] for x in OFF], axis=1); Q = MODELS["M1"]; a0, ax = a0_ax(Q)
    Gd = [atten(a0 + ax[:, i]) for i in range(NX)]; clean = np.stack([Gd[i] @ truth[:, i] for i in range(NX)], axis=1); sig = np.linalg.norm(clean, axis=0) / np.sqrt(n) / 10 ** (18 / 20)

    def ga(D):
        amp = np.zeros((3, NX))
        for j, te in enumerate(C.event_t):
            k = np.abs(T - te) <= 0.07; w = ricker(T[k] - te); amp[j] = (w @ D[k]) / (w @ w)
        return np.array([np.polyfit(S2[j], amp[j], 1)[::-1] for j in range(3)])

    fT = ga(truth); G0 = atten(a0); Gx = [atten(ax[:, i]) for i in range(NX)]; res = {"split": [], "direct": []}
    F0 = np.linalg.solve(G0.T @ G0 + 1e-3 * I, G0.T); Bx = [np.linalg.solve(g.T @ g + LX * I, g.T) for g in Gx]
    for s in range(12):
        d = clean + sig * np.random.default_rng(s).standard_normal(clean.shape)
        D1 = np.stack([F0 @ (Bx[i] @ d[:, i]) for i in range(NX)], axis=1); D2 = np.stack([solve_all(Gd[i], np.ones(n), 1e-4, d[:, [i]])[:, 0] for i in range(NX)], axis=1)
        for k, D in (("split", D1), ("direct", D2)):
            f = ga(D); res[k].append(np.mean(np.abs(f[:, 1] / f[:, 0] - fT[:, 1] / fT[:, 0])))
    print("orig scene, template G/A error: split", round(float(np.mean(res["split"])), 3), "(K1b 0.040); direct", round(float(np.mean(res["direct"])), 3), "(K1b 0.387)", flush=True)
    ok = np.linalg.norm(np.fft.irfft(RIC, NSYN)[0:1]) > 0; print("forward parameter check: Ricker spectrum peak at", float(FS[np.argmax(RIC)]), "Hz", flush=True)


def smoke():
    t0 = time.monotonic(); rows = []
    for scene in ("dense_nmo",):
        cols = build_data("M1", scene)[:4]; print("data built", round(time.monotonic() - t0, 1), "s", flush=True); run_scene("M1", 1.0, scene, rows, cols)
    v = [r["mg"] for r in rows]; print("rows", len(rows), "finite", bool(np.all(np.isfinite(v))), "elapsed", round(time.monotonic() - t0, 1), "s; no scores printed", flush=True)


def run(Mname):
    rows = []; path = OUT / f"gate1_{Mname}.json"
    for scene in SCENES:
        cols = build_data(Mname, scene)
        for qf in QF:
            t0 = time.monotonic(); run_scene(Mname, qf, scene, rows, cols); path.write_text(json.dumps(rows), encoding="utf-8"); print(Mname, scene, qf, f"{time.monotonic() - t0:.0f}s", flush=True)


if __name__ == "__main__":
    m = sys.argv[1]
    if m == "smoke":
        smoke()
    elif m == "check":
        check()
    else:
        run(sys.argv[2])
