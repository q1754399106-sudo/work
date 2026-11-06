"""Frozen test of STEP1_ZH.md on new models M6-M9, wavelet 60 and 40 Hz, SNR 18 and 12. -> results/step1.json. No scores printed."""
import json, time
from pathlib import Path
import numpy as np
import gate1 as g
from gate1 import C, GEO, SIN2, NX, n, T, I, DT, NSYN, a0_ax, atten, ricker
import theory as th
import stable

OUT = Path(__file__).resolve().parents[1] / "results"; LX = 1e-3; GAMMA = 6.0; GRID = np.logspace(-6, 0, 13)[::-1]; QF = (0.8, 1.0, 1.2); SEEDS = 20
NEW = dict(M6=np.array([25.0, 40, 90, 60]), M7=np.array([80.0, 35, 45, 120]), M8=np.array([35.0, 35, 35, 35]), M9=np.array([150.0, 90, 30, 60]))
EV = np.round(C.event_t / DT).astype(int)


def spec(fp):
    return np.abs(np.fft.rfft(np.roll(ricker(np.arange(NSYN) * DT - NSYN * DT / 2, fp), -NSYN // 2)))


def template(D, fp):
    amp = np.zeros((3, NX) + D.shape[2:])
    for j, te in enumerate(C.event_t):
        k = np.abs(T - te) <= 0.07; w = ricker(T[k] - te, fp); amp[j] = np.tensordot(w, D[k], axes=(0, 0)) / (w @ w)
    return amp


def ga(a):
    return np.array([(lambda p: p[0] / p[1])(np.polyfit(SIN2[EV[j]], a[j], 1)) for j in range(3)])


rows = []; t0 = time.monotonic(); idx = EV.copy(); R = np.array([g.rpp(j, np.array([GEO[i]["p"][EV[j]] * C.vp[j] for i in range(NX)])) for j in range(3)])
for fp in (60.0, 40.0):
    g.RIC = spec(fp); th.S2 = g.RIC ** 2
    for Mname, Qt in NEW.items():
        dp, tru = g.synth(idx, R, Qt, False); clean = dp[:n]; rms = np.sqrt(np.sum(dp[:int(2.4 / DT)] ** 2, axis=0) / n); gT = ga(template(tru, fp))
        for snr in (18, 12):
            sigma = rms / 10 ** (snr / 20); target2 = float(np.sum(sigma ** 2) * n); rng = np.random.default_rng(99000 + int(fp) * 100 + int(Mname[1]) * 10 + snr); Dat = clean[:, :, None] + sigma[None, :, None] * rng.standard_normal((n, NX, SEEDS))
            for qf in QF:
                Qm = Qt * qf; a0, ax = a0_ax(Qm); G0 = atten(a0); Gx = [atten(ax[:, i]) for i in range(NX)]; Gd = [atten(a0 + ax[:, i]) for i in range(NX)]; GtG = [G.T @ G for G in Gd]; GtD = [Gd[i].T @ Dat[:, i, :] for i in range(NX)]
                Y = np.stack([np.linalg.solve(G.T @ G + LX * I, G.T @ Dat[:, i, :]) for i, G in enumerate(Gx)], axis=1); ch = {m: [None] * SEEDS for m in ("method", "split", "direct")}; out = {m: np.zeros_like(Dat) for m in ch}; last = {}
                for lam in GRID:
                    fs = np.array([th.fc_star(a0[k], lam) for k in range(n)]); F0 = np.linalg.solve(G0.T @ G0 + lam * I, G0.T); w = np.maximum(np.exp(-2 * np.pi * fs[:, None] * ax), stable.bound(a0, ax, lam, GAMMA))
                    O = dict(split=np.stack([F0 @ Y[:, i, :] for i in range(NX)], axis=1), direct=np.stack([np.linalg.solve(GtG[i] + lam * I, GtD[i]) for i in range(NX)], axis=1), method=np.stack([np.linalg.solve(GtG[i] + lam * np.diag(w[:, i]), GtD[i]) for i in range(NX)], axis=1)); last = O
                    for m, o in O.items():
                        res = sum(np.sum((Gd[i] @ o[:, i, :] - Dat[:, i, :]) ** 2, axis=0) for i in range(NX))
                        for s in range(SEEDS):
                            if ch[m][s] is None and res[s] <= target2:
                                ch[m][s] = float(lam); out[m][:, :, s] = o[:, :, s]
                for m, o in last.items():
                    for s in range(SEEDS):
                        if ch[m][s] is None:
                            ch[m][s] = float(GRID[-1]); out[m][:, :, s] = o[:, :, s]
                for m, O in out.items():
                    meas = ga(template(O, fp).mean(axis=2)); wave = float(np.mean(np.sqrt(np.sum((O - tru[:, :, None]) ** 2, axis=(0, 1)) / np.sum(tru ** 2))))
                    rows.append(dict(M=Mname, fp=fp, snr=snr, qf=qf, m=m, bias=float(np.mean(np.abs(meas - gT))), bias_ev=[float(v) for v in meas - gT], wave=wave, lam=float(np.median(ch[m]))))
                print(Mname, fp, snr, qf, f"{time.monotonic() - t0:.0f}s", flush=True); (OUT / "step1.json").write_text(json.dumps(rows), encoding="utf-8")
