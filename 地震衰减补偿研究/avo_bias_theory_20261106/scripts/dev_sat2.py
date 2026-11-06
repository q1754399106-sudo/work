"""Development on used models (M1, M10, M11), all step 2 scenes, descriptive. Saturated damping with one parameter gamma:
(G'G + lam L'L) s = G'd, L symbol exp[-pi min(f, f_n) a_x], exp(-2 pi f_n a_x) = w_n(t, x) (noise bound of stable.py).
Old metrics (E_G, M_G, corr G, wave) plus: per-event band-limited G/A against truth, noisy and noise-free, and E_G with a local amplitude normalisation.
usage: dev_sat2.py M -> results/dev_sat2_<M>.json.gz"""
import gzip, json, sys, time
from pathlib import Path
import numpy as np
import gate1 as g
from gate1 import C, NX, n, I, DT, T, a0_ax, atten, avo, MODELS, reflectors, synth, to_t0, SEEDS, SNRS, TW, MUTE
import step2 as s2
from step2 import GRID
import theory as th, stable, satdamp as sd
from gate1b import EV
from diag_aligned import bp

OUT = Path(__file__).resolve().parents[1] / "results"; LX = 1e-3; GAM = (6.0, 10.0, 15.0, 25.0); BASE = ("method", "split", "direct"); BANDS = {"b": (10.0, 60.0), "f": None}; SIG = int(round(0.05 / DT))


def columns(M, scene):
    Q = MODELS[M]; nmo = scene.endswith("nmo"); D, C0, Tr, Sg, meta = [], [], [], [], []; keep = None
    for seed in range(SEEDS):
        rng = np.random.default_rng(10000 * list(MODELS).index(M) + 100 * g.SCENES.index(scene) + seed); idx, R = reflectors(scene, rng); dp, tru = synth(idx, R, Q, nmo); rms = np.sqrt(np.sum(dp[:int(2.4 / DT)] ** 2, axis=0) / n); cl, keep = to_t0(dp, nmo)
        for snr in SNRS:
            sg = rms / 10 ** (snr / 20); y, _ = to_t0(dp + sg * rng.standard_normal(dp.shape), nmo); D.append(y); C0.append(cl); Tr.append(tru); Sg.append(sg); meta.append((int(snr), int(seed)))
    return np.stack(D, 2), np.stack(C0, 2), np.stack(Tr, 2), np.stack(Sg, 1), meta, keep


def smooth(X):
    k = np.exp(-0.5 * (np.arange(-3 * SIG, 3 * SIG + 1) / SIG) ** 2); F = np.fft.rfft(X, 2 * n, axis=0); K = np.fft.rfft(k, 2 * n)
    return np.roll(np.fft.irfft(F * K.reshape((-1,) + (1,) * (X.ndim - 1)), 2 * n, axis=0), -3 * SIG, axis=0)[:n]


def metrics(O, Oc, Tru, keep):
    """per-column dict. wave on O; everything else after the 10-60 band unless named."""
    out = {}; mask = (keep & TW[:, None])[:, :, None]; out["wave"] = np.sqrt(np.sum(((O - Tru) * mask) ** 2, axis=(0, 1)) / np.sum((Tru * mask) ** 2, axis=(0, 1)))
    for b, band in BANDS.items():
        At, Bt, ok = avo(bp(Tru, band), keep); A, B, _ = avo(bp(O, band), keep); Ac, Bc, _ = avo(bp(Oc, band), keep); r = ok[:, None]; dot = lambda X, Y: np.sum(X * Y * r, axis=0); nb = dot(Bt, Bt)
        ba = dot(A, At) / dot(At, At); bg = dot(B, Bt) / nb
        out["eg" + b] = np.sqrt(dot(B / ba - Bt, B / ba - Bt) / nb); out["mg" + b] = bg / ba - 1; out["cg" + b] = dot(B, Bt) / np.sqrt(dot(B, B) * nb); out["ba" + b] = ba; out["bg" + b] = bg
        eps = 1e-3 * np.max(At ** 2, axis=0); sa = np.clip(smooth(A * At) / (smooth(At ** 2) + eps), 0.2, 2.0); sac = np.clip(smooth(Ac * At) / (smooth(At ** 2) + eps), 0.2, 2.0)
        out["egl" + b] = np.sqrt(dot(B / sa - Bt, B / sa - Bt) / nb); out["mgl" + b] = dot(B / sa, Bt) / nb - 1; out["mglc" + b] = dot(Bc / sac, Bt) / nb - 1
        e = []; ec = []
        for k0 in EV:      # per-event G/A against the truth over a window around the event
            w = (np.abs(T - T[k0]) <= 0.05)[:, None] & r; d = lambda X, Y: np.sum(X * Y * w, axis=0); tg = d(Bt, At) / d(At, At); e.append(d(B, At) / d(A, At) - tg); ec.append(d(Bc, At) / d(Ac, At) - tg)
        e = np.array(e); ec = np.array(ec); out["ga" + b] = np.mean(np.abs(e), axis=0); out["gac" + b] = np.mean(np.abs(ec), axis=0)
        for j in range(3):
            out[f"gac{j}{b}"] = ec[j]
    return out


def apply_base(m, lam, Q, Y):
    a0, ax = a0_ax(Q); out = np.zeros_like(Y)
    if m == "split":
        G0 = atten(a0); F0 = np.linalg.solve(G0.T @ G0 + lam * I, G0.T)
        for i in range(NX):
            G = atten(ax[:, i]); out[:, i, :] = F0 @ np.linalg.solve(G.T @ G + LX * I, G.T @ Y[:, i, :])
        return out
    w = np.ones((n, NX))
    if m == "method":
        fs = np.array([th.fc_star(a0[k], lam) for k in range(n)]); w = np.maximum(np.exp(-2 * np.pi * fs[:, None] * ax), stable.bound(a0, ax, lam, 6.0))
    for i in range(NX):
        G = atten(a0 + ax[:, i]); out[:, i, :] = np.linalg.solve(G.T @ G + lam * np.diag(w[:, i]), G.T @ Y[:, i, :])
    return out


def rec(rows, lab, m, O, Oc, Tru, keep, meta, lams, fb=0):
    s = metrics(O, Oc, Tru, keep)
    for c in range(O.shape[2]):
        rows.append(dict(lab, m=m, snr=meta[c][0], seed=meta[c][1], lam=float(lams[c]), fallback=fb, **{k: float(v[c]) for k, v in s.items()}))


def cell(M, scene, qf, ref, rows):
    Q = MODELS[M] * qf; D, Cl, Tr, sig, meta, keep = columns(M, scene); nc = D.shape[2]; lab = dict(M=M, scene=scene, qf=qf); a0, ax = a0_ax(Q)
    for m in BASE:
        lams = np.array([ref[(scene, qf, *meta[c], m)]["lam"] for c in range(nc)]); O = np.zeros_like(D); Oc = np.zeros_like(D)
        for lam in np.unique(lams):
            k = np.where(lams == lam)[0]; O[:, :, k] = apply_base(m, lam, Q, D[:, :, k]); Oc[:, :, k] = apply_base(m, lam, Q, Cl[:, :, k])
        rec(rows, lab, m, O, Oc, Tr, keep, meta, lams)
    Gd = [atten(a0 + ax[:, i]) for i in range(NX)]; GtG = [G.T @ G for G in Gd]; GtD = [Gd[i].T @ D[:, i, :] for i in range(NX)]; GtC = [Gd[i].T @ Cl[:, i, :] for i in range(NX)]; target = np.sum(sig ** 2, axis=0) * n
    for gam in GAM:
        ch = [None] * nc; out = np.zeros_like(D); outc = np.zeros_like(D); last = None
        for lam in GRID:
            wn = stable.bound(a0, ax, lam, gam); O = np.zeros_like(D); Oc = np.zeros_like(D)
            for i in range(NX):
                fn = sd.fb_floor(np.full(n, sd.FF[-1]), ax[:, i], wn[:, i]); L = sd.Lmat(ax[:, i], fn); K = GtG[i] + lam * (L.T @ L); O[:, i, :] = np.linalg.solve(K, GtD[i]); Oc[:, i, :] = np.linalg.solve(K, GtC[i])
            res = sum(np.sum(((Gd[i] @ O[:, i, :] - D[:, i, :]) * keep[:, i, None]) ** 2, axis=0) for i in range(NX)); last = (O, Oc)
            for c in range(nc):
                if ch[c] is None and res[c] <= target[c]:
                    ch[c] = float(lam); out[:, :, c] = O[:, :, c]; outc[:, :, c] = Oc[:, :, c]
            if all(x is not None for x in ch):
                break
        fb = sum(x is None for x in ch)
        for c in range(nc):
            if ch[c] is None:
                ch[c] = float(GRID[-1]); out[:, :, c] = last[0][:, :, c]; outc[:, :, c] = last[1][:, :, c]
        rec(rows, lab, f"sat{gam:.0f}", out, outc, Tr, keep, meta, ch, fb)


def main(M):
    R0 = json.loads(gzip.open(OUT / f"step2_{M}.json.gz", "rt", encoding="utf-8").read()); ref = {(r["scene"], r["qf"], r["snr"], r["seed"], r["m"]): r for r in R0}; rows = []; t0 = time.monotonic(); p = OUT / f"dev_sat2_{M}.json.gz"
    for scene in g.SCENES:
        for qf in (1.0, 0.8, 1.2):
            cell(M, scene, qf, ref, rows); print(M, scene, qf, f"{time.monotonic() - t0:.0f}s", flush=True); gzip.open(p, "wt", encoding="utf-8").write(json.dumps(rows))


if __name__ == "__main__":
    main(sys.argv[1])
