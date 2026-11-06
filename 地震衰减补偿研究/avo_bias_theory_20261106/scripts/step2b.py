"""Frozen test of STEP2B_ZH.md: saturated damping against the split, the single-frequency method and one-step, on new models M12-M15, wavelets 60 and 40 Hz.
usage: step2b.py smoke | run M -> results/step2b_<M>.json.gz. No scores printed."""
import gzip, json, sys, time
from pathlib import Path
import numpy as np
import gate1 as g
from gate1 import C, NX, n, I, DT, T, NSYN, a0_ax, atten, avo, reflectors, synth, to_t0, TW, ricker
import theory as th, stable, satdamp as sd

OUT = Path(__file__).resolve().parents[1] / "results"; LX = 1e-3; GAMMA = 6.0; GRID = np.logspace(-6, 0, 13)[::-1]; QF = (1.0, 0.8, 1.2); SNRS = (18, 12); SEEDS = 10; FPS = (60.0, 40.0); METH = ("sat", "split", "method", "direct")
NEW = dict(M12=np.array([55.0, 40, 100, 70]), M13=np.array([28.0, 65, 40, 90]), M14=np.array([70.0, 28, 35, 50]), M15=np.array([40.0, 120, 60, 30])); SCENES = list(g.SCENES); EV = np.round(C.event_t / DT).astype(int); SIG = int(round(0.05 / DT))


def spec(fp):
    return np.abs(np.fft.rfft(np.roll(ricker(np.arange(NSYN) * DT - NSYN * DT / 2, fp), -NSYN // 2)))


def bp(X, band):
    F = np.fft.rfft(X, 2048, axis=0); f = np.fft.rfftfreq(2048, DT); w = np.clip((f - band[0] + 5) / 5, 0, 1) * np.clip((band[1] + 8 - f) / 8, 0, 1)
    return np.fft.irfft(F * w.reshape((-1,) + (1,) * (X.ndim - 1)), 2048, axis=0)[:n]


def columns(Mname, scene, fp):
    Q = NEW[Mname]; nmo = scene.endswith("nmo"); D, C0, Tr, Sg, meta = [], [], [], [], []; keep = None
    for seed in range(SEEDS):
        rng = np.random.default_rng(8_000_000 + int(fp) * 10000 + int(Mname[1:]) * 1000 + SCENES.index(scene) * 100 + seed); idx, R = reflectors(scene, rng); dp, tru = synth(idx, R, Q, nmo); rms = np.sqrt(np.sum(dp[:int(2.4 / DT)] ** 2, axis=0) / n); cl, keep = to_t0(dp, nmo)
        for snr in SNRS:
            sg = rms / 10 ** (snr / 20); y, _ = to_t0(dp + sg * rng.standard_normal(dp.shape), nmo); D.append(y); C0.append(cl); Tr.append(tru); Sg.append(sg); meta.append((int(snr), int(seed)))
    return np.stack(D, 2), np.stack(C0, 2), np.stack(Tr, 2), np.stack(Sg, 1), meta, keep


def smooth(X):
    k = np.exp(-0.5 * (np.arange(-3 * SIG, 3 * SIG + 1) / SIG) ** 2); F = np.fft.rfft(X, 2 * n, axis=0); K = np.fft.rfft(k, 2 * n)
    return np.roll(np.fft.irfft(F * K.reshape((-1,) + (1,) * (X.ndim - 1)), 2 * n, axis=0), -3 * SIG, axis=0)[:n]


def metrics(O, Oc, Tru, keep, fp):
    out = {}; mask = (keep & TW[:, None])[:, :, None]; out["wave"] = np.sqrt(np.sum(((O - Tru) * mask) ** 2, axis=(0, 1)) / np.sum((Tru * mask) ** 2, axis=(0, 1))); band = (10.0, fp); win = 3.0 / fp
    At, Bt, ok = avo(bp(Tru, band), keep); A, B, _ = avo(bp(O, band), keep); Ac, Bc, _ = avo(bp(Oc, band), keep); r = ok[:, None]; dot = lambda X, Y: np.sum(X * Y * r, axis=0); nb = dot(Bt, Bt)
    ba = dot(A, At) / dot(At, At); bg = dot(B, Bt) / nb
    out["eg"] = np.sqrt(dot(B / ba - Bt, B / ba - Bt) / nb); out["mg"] = bg / ba - 1; out["cg"] = dot(B, Bt) / np.sqrt(dot(B, B) * nb); out["ba"] = ba; out["bg"] = bg
    eps = 1e-3 * np.max(At ** 2, axis=0); sa = np.clip(smooth(A * At) / (smooth(At ** 2) + eps), 0.2, 2.0); sac = np.clip(smooth(Ac * At) / (smooth(At ** 2) + eps), 0.2, 2.0)
    out["egl"] = np.sqrt(dot(B / sa - Bt, B / sa - Bt) / nb); out["mglc"] = dot(Bc / sac, Bt) / nb - 1
    e, ec = [], []
    for k0 in EV:
        w = (np.abs(T - T[k0]) <= win)[:, None] & r; d = lambda X, Y: np.sum(X * Y * w, axis=0); tg = d(Bt, At) / d(At, At); e.append(d(B, At) / d(A, At) - tg); ec.append(d(Bc, At) / d(Ac, At) - tg)
    e = np.array(e); ec = np.array(ec); out["ga"] = np.mean(np.abs(e), axis=0); out["gac"] = np.mean(np.abs(ec), axis=0)
    for j in range(3):
        out[f"gac{j}"] = ec[j]
    return out


def cell(Mname, scene, qf, fp, rows):
    Q = NEW[Mname] * qf; D, Cl, Tr, sig, meta, keep = columns(Mname, scene, fp); nc = D.shape[2]; a0, ax = a0_ax(Q); RH = np.concatenate([D, Cl], axis=2); G0 = atten(a0); Gx = [atten(ax[:, i]) for i in range(NX)]; Gd = [atten(a0 + ax[:, i]) for i in range(NX)]; GtG = [G.T @ G for G in Gd]; GtR = [Gd[i].T @ RH[:, i, :] for i in range(NX)]
    Y = np.stack([np.linalg.solve(G.T @ G + LX * I, G.T @ RH[:, i, :]) for i, G in enumerate(Gx)], axis=1); At0 = None; target = np.sum(sig ** 2, axis=0) * n; ch = {m: [None] * nc for m in METH}; out = {m: np.zeros_like(RH) for m in METH}; last = {}
    for lam in GRID:
        todo = [m for m in METH if any(v is None for v in ch[m])]
        if not todo:
            break
        wn = stable.bound(a0, ax, lam, GAMMA) if ("sat" in todo or "method" in todo) else None
        for m in todo:
            if m == "split":
                F0 = np.linalg.solve(G0.T @ G0 + lam * I, G0.T); O = np.stack([F0 @ Y[:, i, :] for i in range(NX)], axis=1)
            elif m == "direct":
                O = np.stack([np.linalg.solve(GtG[i] + lam * I, GtR[i]) for i in range(NX)], axis=1)
            elif m == "method":
                fs = np.array([th.fc_star(a0[k], lam) for k in range(n)]); w = np.maximum(np.exp(-2 * np.pi * fs[:, None] * ax), wn); O = np.stack([np.linalg.solve(GtG[i] + lam * np.diag(w[:, i]), GtR[i]) for i in range(NX)], axis=1)
            else:
                O = np.zeros_like(RH)
                for i in range(NX):
                    L = sd.Lmat(ax[:, i], sd.fb_floor(np.full(n, sd.FF[-1]), ax[:, i], wn[:, i])); O[:, i, :] = np.linalg.solve(GtG[i] + lam * (L.T @ L), GtR[i])
            last[m] = O; res = sum(np.sum(((Gd[i] @ O[:, i, :nc] - D[:, i, :]) * keep[:, i, None]) ** 2, axis=0) for i in range(NX))
            for c in range(nc):
                if ch[m][c] is None and res[c] <= target[c]:
                    ch[m][c] = float(lam); out[m][:, :, c] = O[:, :, c]; out[m][:, :, nc + c] = O[:, :, nc + c]
    for m in METH:
        fb = sum(v is None for v in ch[m])
        for c in range(nc):
            if ch[m][c] is None:
                ch[m][c] = float(GRID[-1]); out[m][:, :, c] = last[m][:, :, c]; out[m][:, :, nc + c] = last[m][:, :, nc + c]
        s = metrics(out[m][:, :, :nc], out[m][:, :, nc:], Tr, keep, fp)
        for c in range(nc):
            rows.append(dict(M=Mname, fp=fp, scene=scene, qf=qf, m=m, snr=meta[c][0], seed=meta[c][1], lam=ch[m][c], fallback=fb, **{k: float(v[c]) for k, v in s.items()}))


def run(Mname, only=None):
    rows = []; t0 = time.monotonic(); p = OUT / f"step2b_{Mname}.json.gz"
    for fp in FPS:
        g.RIC = spec(fp); th.S2 = g.RIC ** 2
        for scene in SCENES:
            for qf in QF:
                cell(Mname, scene, qf, fp, rows); print(Mname, fp, scene, qf, f"{time.monotonic() - t0:.0f}s", flush=True)
                if only:
                    return rows
                gzip.open(p, "wt", encoding="utf-8").write(json.dumps(rows))


if __name__ == "__main__":
    if sys.argv[1] == "smoke":      # finite check only, on a used model name is not possible: uses one new-model cell without storing or printing scores
        r = run("M12", only=True); print("rows", len(r), "finite", bool(np.all(np.isfinite([x["wave"] + x["ga"] + x["gac"] for x in r]))))
    else:
        run(sys.argv[2])
