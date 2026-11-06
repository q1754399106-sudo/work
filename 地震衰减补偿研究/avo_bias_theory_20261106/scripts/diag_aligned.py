"""Post hoc diagnosis of step 2, aligned scenes, Q correct, used models only (descriptive, no rule).
Why is per-sample E_G of the method worse than the split when the template bias is the same?
Uses the lambda each method chose in step 2. usage: diag_aligned.py M -> results/diag_aligned_<M>.json
Exact split of the error: E_G^2 = M_G^2 + |R|^2 / |Bt|^2, R = residual orthogonal to Bt, with R = shape (signal) + noise (+ a small cross term)."""
import gzip, json, sys
from pathlib import Path
import numpy as np
import gate1 as g
from gate1 import C, NX, n, I, DT, a0_ax, atten, avo, MODELS, reflectors, synth, SEEDS, SNRS, TW, FS
import step2 as s2
import theory as th
import stable

OUT = Path(__file__).resolve().parents[1] / "results"; LX = 1e-3; GAMMA = 6.0; METH = ("method", "split", "direct")
BANDS = {"10-60": (10.0, 60.0), "10-40": (10.0, 40.0), "10-100": (10.0, 100.0), "none": None}; FB = np.arange(10, 130, 10.0)
SCENES = ("events_aligned", "thin_aligned", "dense_aligned")


def bp(X, band):
    if band is None:
        return X
    F = np.fft.rfft(X, 2048, axis=0); f = np.fft.rfftfreq(2048, DT); w = np.clip((f - band[0] + 5) / 5, 0, 1) * np.clip((band[1] + 8 - f) / 8, 0, 1)
    return np.fft.irfft(F * w.reshape((-1,) + (1,) * (X.ndim - 1)), 2048, axis=0)[:n]


def columns(M, scene):
    """Reproduces build_data, keeping signal and noise apart."""
    Q = MODELS[M]; D, N, Tr, Sg, meta = [], [], [], [], []
    for seed in range(SEEDS):
        rng = np.random.default_rng(10000 * list(MODELS).index(M) + 100 * g.SCENES.index(scene) + seed); idx, R = reflectors(scene, rng); dp, tru = synth(idx, R, Q, False); rms = np.sqrt(np.sum(dp[:int(2.4 / DT)] ** 2, axis=0) / n)
        for snr in SNRS:
            sg = rms / 10 ** (snr / 20); nz = sg * rng.standard_normal(dp.shape); D.append(dp[:n]); N.append(nz[:n]); Tr.append(tru); Sg.append(sg); meta.append((int(snr), int(seed)))
    return np.stack(D, 2), np.stack(N, 2), np.stack(Tr, 2), np.stack(Sg, 1), meta


def op(m, lam, Q, D):
    a0, ax = a0_ax(Q); out = np.zeros_like(D); nc = D.shape[2]
    if m == "split":
        G0 = atten(a0); F0 = np.linalg.solve(G0.T @ G0 + lam * I, G0.T)
        for i in range(NX):
            G = atten(ax[:, i]); out[:, i, :] = F0 @ np.linalg.solve(G.T @ G + LX * I, G.T @ D[:, i, :])
        return out, {}
    w = np.ones((n, NX)); info = {}
    if m == "method":
        fs = np.array([th.fc_star(a0[k], lam) for k in range(n)]); we = np.exp(-2 * np.pi * fs[:, None] * ax); wn = stable.bound(a0, ax, lam, GAMMA); w = np.maximum(we, wn)
        info = dict(fc=float(np.median(fs[TW])), bound_frac=float(np.mean((wn > we)[TW][:, 1:])), bound_frac_far=float(np.mean((wn > we)[TW][:, -3:])))
    for i in range(NX):
        G = atten(a0 + ax[:, i]); out[:, i, :] = np.linalg.solve(G.T @ G + lam * np.diag(w[:, i]), G.T @ D[:, i, :])
    return out, info


def spectra(Bs, Bt):
    """Transfer of the signal part against truth per 10 Hz band, Re<X, Xt>/<Xt, Xt>; rows = bands."""
    f = np.fft.rfftfreq(2048, DT); M = TW[:, None].astype(float); X = np.fft.rfft(Bs * M, 2048, axis=0); Xt = np.fft.rfft(Bt * M, 2048, axis=0); out = []
    for lo in FB:
        k = (f >= lo) & (f < lo + 10); out.append(np.sum((X[k] * np.conj(Xt[k])).real, axis=0) / np.maximum(np.sum(np.abs(Xt[k]) ** 2, axis=0), 1e-30))
    return np.array(out)


def main(M):
    R0 = json.loads(gzip.open(OUT / f"step2_{M}.json.gz", "rt", encoding="utf-8").read()); ref = {(r["scene"], r["snr"], r["seed"], r["m"]): r for r in R0 if r["qf"] == 1.0 and "aligned" in r["scene"]}
    Q = MODELS[M]; rows = []; keep = np.ones((n, NX), bool)
    for scene in SCENES:
        S, Nz, Tr, sg, meta = columns(M, scene); nc = S.shape[2]
        for m in METH:
            lams = np.array([ref[(scene, *meta[c], m)]["lam"] for c in range(nc)]); Os = np.zeros_like(S); On = np.zeros_like(S); info = {}
            for lam in np.unique(lams):
                k = np.where(lams == lam)[0]; Os[:, :, k], info[float(lam)] = op(m, lam, Q, S[:, :, k]); On[:, :, k], _ = op(m, lam, Q, Nz[:, :, k])
            rec = [dict(M=M, scene=scene, m=m, snr=meta[c][0], seed=meta[c][1], lam=float(lams[c]), **info[float(lams[c])]) for c in range(nc)]
            for name, band in BANDS.items():
                At, Bt, ok = avo(bp(Tr, band), keep); As, Bs, _ = avo(bp(Os, band), keep); An, Bn, _ = avo(bp(On, band), keep); A = As + An; B = Bs + Bn; r = ok[:, None]; dotT = lambda X, Y: np.sum(X * Y * r, axis=0)
                ba = dotT(A, At) / dotT(At, At); bg = dotT(B, Bt) / dotT(Bt, Bt); bas = dotT(As, At) / dotT(At, At); bgs = dotT(Bs, Bt) / dotT(Bt, Bt); bgn = dotT(Bn, Bt) / dotT(Bt, Bt); nBt = dotT(Bt, Bt)
                eg = np.sqrt(dotT(B / ba - Bt, B / ba - Bt) / nBt); shape = np.sqrt(dotT(Bs - bgs * Bt, Bs - bgs * Bt) / nBt) / ba; noise = np.sqrt(dotT(Bn - bgn * Bt, Bn - bgn * Bt) / nBt) / ba
                for c in range(nc):
                    rec[c].update({f"eg_{name}": float(eg[c]), f"mg_{name}": float(abs(bg[c] / ba[c] - 1)), f"sg_{name}": float(bg[c] / ba[c] - 1), f"shape_{name}": float(shape[c]), f"noise_{name}": float(noise[c]), f"mgs_{name}": float(bgs[c] / bas[c] - 1), f"ba_{name}": float(ba[c]), f"bg_{name}": float(bg[c])})
            for c in range(nc):      # reproduction check against the stored score
                rec[c]["stored_eg"] = ref[(scene, *meta[c], m)]["eg"]
            At, Bt, ok = avo(Tr, keep); As, Bs, _ = avo(Os, keep); ga = spectra(As, At); gb = spectra(Bs, Bt)
            for c in range(nc):
                rec[c]["ga"] = ga[:, c].tolist(); rec[c]["gb"] = gb[:, c].tolist()
            rows += rec; print(M, scene, m, flush=True)
    (OUT / f"diag_aligned_{M}.json").write_text(json.dumps(rows), encoding="utf-8")


if __name__ == "__main__":
    main(sys.argv[1])
