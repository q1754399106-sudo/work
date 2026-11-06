"""Step 2 (STEP2_ZH.md): the frozen method in the realistic scenes. usage: step2.py run M | stress -> results/step2_<M>.json.gz, step2_stress.json
No scores printed."""
import gzip, json, sys, time
from pathlib import Path
import numpy as np
import gate1 as g
from gate1 import C, GEO, NX, n, I, DT, a0_ax, atten, avo, bandpass, MODELS, build_data, TW, reflectors, synth, MUTE, T
import theory as th
import stable

OUT = Path(__file__).resolve().parents[1] / "results"; MODELS.update(M4=np.array([60.0, 45, 80, 55]), M5=np.array([30.0, 60, 120, 90]), M10=np.array([45.0, 70, 55, 85]), M11=np.array([30.0, 30, 60, 100]))
LX = 1e-3; GAMMA = 6.0; GRID = np.logspace(-6, 0, 13)[::-1]; QF = (0.8, 1.0, 1.2); SCENES = ("events_aligned", "thin_aligned", "dense_aligned", "events_nmo", "thin_nmo", "dense_nmo"); METH = ("method", "split", "direct")


def score(O, Tru, keep, At, Bt, ok):
    mask = (keep & TW[:, None])[:, :, None]; wave = np.sqrt(np.sum(((O - Tru) * mask) ** 2, axis=(0, 1)) / np.sum((Tru * mask) ** 2, axis=(0, 1)))
    A, B, _ = avo(bandpass(O), keep); r = ok[:, None]; ba = np.sum(A * At * r, axis=0) / np.sum(At ** 2 * r, axis=0); bg = np.sum(B * Bt * r, axis=0) / np.sum(Bt ** 2 * r, axis=0)
    cg = np.sum(B * Bt * r, axis=0) / np.sqrt(np.sum(B ** 2 * r, axis=0) * np.sum(Bt ** 2 * r, axis=0)); eg = np.sqrt(np.sum(((B / ba - Bt) * r) ** 2, axis=0) / np.sum((Bt * r) ** 2, axis=0))
    return dict(wave=wave, eg=eg, mg=np.abs(bg / ba - 1), ma=np.abs(ba - 1), mb=np.abs(bg - 1), cg=cg)


def run_cell(Qm, Dat, Tru, keep, sig, meta, label, rows):
    a0, ax = a0_ax(Qm); G0 = atten(a0); Gx = [atten(ax[:, i]) for i in range(NX)]; Gd = [atten(a0 + ax[:, i]) for i in range(NX)]; GtG = [G.T @ G for G in Gd]; GtD = [Gd[i].T @ Dat[:, i, :] for i in range(NX)]; nc = Dat.shape[2]
    Y = np.stack([np.linalg.solve(G.T @ G + LX * I, G.T @ Dat[:, i, :]) for i, G in enumerate(Gx)], axis=1); At, Bt, ok = avo(bandpass(Tru), keep); target = np.sum(sig ** 2, axis=0) * n; ch = {m: [None] * nc for m in METH}; out = {m: np.zeros_like(Dat) for m in METH}; last = {}
    for lam in GRID:
        fs = np.array([th.fc_star(a0[k], lam) for k in range(n)]); F0 = np.linalg.solve(G0.T @ G0 + lam * I, G0.T); w = np.maximum(np.exp(-2 * np.pi * fs[:, None] * ax), stable.bound(a0, ax, lam, GAMMA))
        O = dict(split=np.stack([F0 @ Y[:, i, :] for i in range(NX)], axis=1), direct=np.stack([np.linalg.solve(GtG[i] + lam * I, GtD[i]) for i in range(NX)], axis=1), method=np.stack([np.linalg.solve(GtG[i] + lam * np.diag(w[:, i]), GtD[i]) for i in range(NX)], axis=1)); last = O
        for m, o in O.items():
            res = sum(np.sum(((Gd[i] @ o[:, i, :] - Dat[:, i, :]) * keep[:, i, None]) ** 2, axis=0) for i in range(NX))
            for c in range(nc):
                if ch[m][c] is None and res[c] <= target[c]:
                    ch[m][c] = float(lam); out[m][:, :, c] = o[:, :, c]
    fb = {}
    for m, o in last.items():
        fb[m] = sum(v is None for v in ch[m])
        for c in range(nc):
            if ch[m][c] is None:
                ch[m][c] = float(GRID[-1]); out[m][:, :, c] = o[:, :, c]
    for m in METH:
        s = score(out[m], Tru, keep, At, Bt, ok)
        for c in range(nc):
            rows.append(dict(label, snr=meta[c][0], seed=meta[c][1], m=m, lam=ch[m][c], fallback=fb[m], **{k: float(v[c]) for k, v in s.items()}))


def run(Mname):
    rows = []; path = OUT / f"step2_{Mname}.json.gz"; t0 = time.monotonic(); Q0 = MODELS[Mname]
    for scene in SCENES:
        g.MODELS = MODELS; cols = build_data(Mname, scene); keep = cols[0][4]; Dat = np.stack([c[2] for c in cols], axis=2); Tru = np.stack([c[3] for c in cols], axis=2); sig = np.stack([c[5] for c in cols], axis=1); meta = [(int(c[0]), int(c[1])) for c in cols]
        for qf in QF:
            run_cell(Q0 * qf, Dat, Tru, keep, sig, meta, dict(M=Mname, scene=scene, qf=qf), rows)
            print(Mname, scene, qf, f"{time.monotonic() - t0:.0f}s", flush=True); gzip.open(path, "wt", encoding="utf-8").write(json.dumps(rows))


def stress():
    from gate1b import nmo_cols, txs, EV
    rows = []; Q0 = MODELS["M1"]; SC = list(g.SCENES)
    for scene in ("events_nmo", "dense_nmo"):
        cols = [c for c in build_data("M1", scene) if c[0] == 18]
        for vs in (0.98, 1.02):
            keep = np.stack([txs(i, vs) <= MUTE * np.maximum(T, DT) for i in range(NX)], axis=1); D, Tr, S, meta = [], [], [], []
            for snr, seed, y, tru, kp, sg in cols:
                rng = np.random.default_rng(100 * SC.index(scene) + seed); idx, R = reflectors(scene, rng); dp, tr = synth(idx, R, Q0, True); d = dp + sg[None, :] * rng.standard_normal(dp.shape)
                D.append(np.stack([nmo_cols(d[:, [i]], txs(i, vs), keep[:, i])[:, 0] for i in range(NX)], axis=1)); Tr.append(tr); S.append(sg); meta.append((snr, seed))
            run_cell(Q0, np.stack(D, axis=2), np.stack(Tr, axis=2), keep, np.stack(S, axis=1), meta, dict(M="M1", scene=scene, qf=1.0, vs=vs), rows); print(scene, vs, flush=True)
    (OUT / "step2_stress.json").write_text(json.dumps(rows), encoding="utf-8")


if __name__ == "__main__":
    if sys.argv[1] == "run":
        run(sys.argv[2])
    else:
        stress()
