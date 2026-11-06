"""Measured template G/A on isolated events with the independent forward model and the real column-varying solvers,
next to the closed-form prediction. usage: measure.py dev   (models M1-M3, already seen) | test (frozen, new models and seeds)"""
import json, sys, time
from pathlib import Path
import numpy as np
import gate1 as g
from gate1 import C, GEO, SIN2, NX, n, T, I, DT, a0_ax, atten, ricker, synth, MODELS
import theory as th

OUT = Path(__file__).resolve().parents[1] / "results"; EV = np.round(C.event_t / DT).astype(int); LX = 1e-3
MODELS.update(M4=np.array([60.0, 45, 80, 55]), M5=np.array([30.0, 60, 120, 90]))
LAMS = (1e-5, 1e-4, 1e-3, 1e-2); QF = (0.8, 1.0, 1.2); FCS = (20.0, 40.0, 60.0, 80.0, 100.0)


def reflector_series():
    idx = EV.copy(); R = np.array([g.rpp(j, np.array([GEO[i]["p"][EV[j]] * C.vp[j] for i in range(NX)])) for j in range(3)]); return idx, R


def template(D):
    amp = np.zeros((3, NX) + D.shape[2:])
    for j, te in enumerate(C.event_t):
        k = np.abs(T - te) <= 0.07; w = ricker(T[k] - te); amp[j] = np.tensordot(w, D[k], axes=(0, 0)) / (w @ w)
    return amp


def ga(amp):
    """Signed G/A per event from amplitudes (3, NX) by LS against sin^2."""
    out = np.zeros(3)
    for j in range(3):
        p = np.polyfit(SIN2[EV[j]], amp[j], 1); out[j] = p[0] / p[1]
    return out


def ridge(G, lam, w, D):
    return np.linalg.solve(G.T @ G + lam * np.diag(w), G.T @ D)


def run(Mname, snr, seeds, base):
    Qt = MODELS[Mname]; idx, R = reflector_series(); dp, tru = synth(idx, R, Qt, False); clean = dp[:n]; rms = np.sqrt(np.sum(dp[:int(2.4 / DT)] ** 2, axis=0) / n); rows = []
    rng = np.random.default_rng(base + 1000 * list(MODELS).index(Mname) + snr); Dat = clean[:, :, None] + (rms / 10 ** (snr / 20))[None, :, None] * rng.standard_normal((n, NX, seeds))
    aT = template(tru); gT = ga(aT)
    for qf in QF:
        Qm = Qt * qf; a0, ax = a0_ax(Qm); G0 = atten(a0); Gx = [atten(ax[:, i]) for i in range(NX)]; Gd = [atten(a0 + ax[:, i]) for i in range(NX)]
        Y = np.stack([np.linalg.solve(G.T @ G + LX * I, G.T @ Dat[:, i, :]) for i, G in enumerate(Gx)], axis=1)
        for lam in LAMS:
            F0 = np.linalg.solve(G0.T @ G0 + lam * I, G0.T); fs_t = np.array([th.fc_star(a0[k], lam) for k in range(n)])
            outs = {"split": np.stack([F0 @ Y[:, i, :] for i in range(NX)], axis=1), "direct": np.stack([ridge(Gd[i], lam, np.ones(n), Dat[:, i, :]) for i in range(NX)], axis=1)}
            outs["fc_star"] = np.stack([ridge(Gd[i], lam, np.exp(-2 * np.pi * fs_t * ax[:, i]), Dat[:, i, :]) for i in range(NX)], axis=1)
            for fc in FCS:
                outs[f"fc{fc:.0f}"] = np.stack([ridge(Gd[i], lam, np.exp(-2 * np.pi * fc * ax[:, i]), Dat[:, i, :]) for i in range(NX)], axis=1)
            for m, O in outs.items():
                A = template(O); wave = float(np.mean(np.sqrt(np.sum((O - tru[:, :, None]) ** 2, axis=(0, 1)) / np.sum(tru ** 2))))
                per = np.array([ga(A[:, :, s]) for s in range(seeds)])  # seeds, events
                meas = ga(A.mean(axis=2))    # G/A of the mean amplitudes: the bias, with little noise influence
                for j in range(3):
                    r = aT[j]; rec = dict(M=Mname, snr=snr, qf=qf, lam=lam, m=m, j=j, truth=float(gT[j]), meas=float(meas[j]), sd=float(per[:, j].std(ddof=1) / np.sqrt(seeds)), mae=float(np.mean(np.abs(per[:, j] - gT[j]))), wave=wave, fc_star=float(fs_t[EV[j]]))
                    if m in ("split", "direct"):
                        rec["pred"] = th.predict_ga(m, j, Qt, Qm, lam, r)[0]
                    elif m == "fc_star":
                        rec["pred"] = th.predict_ga("cons", j, Qt, Qm, lam, r, fc=fs_t[EV[j]])[0]
                    else:
                        rec["pred"] = th.predict_ga("cons", j, Qt, Qm, lam, r, fc=float(m[2:]))[0]
                    if m == "direct" and qf == 1.0:
                        rec["first"] = th.first_order(j, Qm, lam)
                    rows.append(rec)
    return rows


def summary(rows, title):
    rows = [r for r in rows]; print(title)
    for m in ("direct", "split", "fc_star"):
        for lab, sel in (("Q correct", lambda r: r["qf"] == 1.0), ("Q wrong", lambda r: r["qf"] != 1.0)):
            S = [r for r in rows if r["m"] == m and sel(r)]; d = np.array([r["meas"] - r["truth"] for r in S]); p = np.array([r["pred"] - r["truth"] for r in S]); c = np.polyfit(p, d, 1) if np.std(p) > 0 else [np.nan, np.nan]; r2 = np.corrcoef(p, d)[0, 1] ** 2 if np.std(p) > 0 else np.nan
            print(f"  {m:8s} {lab:9s} n {len(S):3d}: measured vs predicted bias: slope {c[0]:6.3f} R2 {r2:6.3f} rmse {np.sqrt(np.mean((d - p) ** 2)):.3f}; mean |measured| {np.mean(np.abs(d)):.3f} mean |pred| {np.mean(np.abs(p)):.3f}")
    F = [r for r in rows if r["m"] == "direct" and r["qf"] == 1.0 and "first" in r]; d = np.array([r["meas"] - r["truth"] for r in F]); a = np.array([r["first"] for r in F]); print(f"  first-order formula, direct, Q correct: slope {np.polyfit(a, d, 1)[0]:.3f} R2 {np.corrcoef(a, d)[0, 1] ** 2:.3f} rmse {np.sqrt(np.mean((d - a) ** 2)):.3f}")
    # optimal fc against fc*: for each case interpolate the zero crossing of the measured bias across the fc grid
    keys = sorted({(r["M"], r["snr"], r["qf"], r["lam"], r["j"]) for r in rows}); ok = []
    for k in keys:
        C_ = {r["m"]: r for r in rows if (r["M"], r["snr"], r["qf"], r["lam"], r["j"]) == k}; xs = np.array(FCS); ys = np.array([C_[f"fc{f:.0f}"]["meas"] - C_[f"fc{f:.0f}"]["truth"] for f in FCS]); cross = None
        for a, b, u, v in zip(xs[:-1], xs[1:], ys[:-1], ys[1:]):
            if u * v <= 0 and u != v:
                cross = a + (b - a) * u / (u - v); break
        ok.append((k, cross, C_["direct"]["fc_star"]))
    for lab, sel in (("Q correct", lambda k: k[2] == 1.0), ("Q wrong", lambda k: k[2] != 1.0)):
        A = [(c, s) for k, c, s in ok if sel(k) and c is not None]; n_all = sum(1 for k, c, s in ok if sel(k))
        if A:
            A = np.array(A); print(f"  optimal fc vs fc* ({lab}): crossing found {len(A)}/{n_all}; corr {np.corrcoef(A[:, 0], A[:, 1])[0, 1]:.3f}; mean |diff| {np.mean(np.abs(A[:, 0] - A[:, 1])):.1f} Hz; mean crossing {A[:, 0].mean():.1f}, mean fc* {A[:, 1].mean():.1f}")
    print("  mean bias magnitude and waveform, by method at each lambda (Q correct):")
    for lam in LAMS:
        line = f"   lam {lam:.0e}"
        for m in ("direct", "split", "fc_star"):
            S = [r for r in rows if r["m"] == m and r["qf"] == 1.0 and r["lam"] == lam]; line += f" | {m} |bias| {np.mean([abs(r['meas'] - r['truth']) for r in S]):.3f} wave {np.mean([r['wave'] for r in S]):.3f}"
        print(line)


def ricker_spec2(fp):
    return np.abs(np.fft.rfft(np.roll(ricker(np.arange(g.NSYN) * DT - g.NSYN * DT / 2, fp), -g.NSYN // 2))) ** 2


def run_disc(Mname, snr, seeds, base):
    """No-truth lambda: one lambda per gather, the largest on the grid whose total residual <= sigma sqrt(n NX). Sigma known."""
    Qt = MODELS[Mname]; idx, R = reflector_series(); dp, tru = synth(idx, R, Qt, False); clean = dp[:n]; rms = np.sqrt(np.sum(dp[:int(2.4 / DT)] ** 2, axis=0) / n); rows = []
    sigma = rms / 10 ** (snr / 20); target2 = float(np.sum(sigma ** 2) * n); rng = np.random.default_rng(base + 500 + 1000 * list(MODELS).index(Mname) + snr); Dat = clean[:, :, None] + sigma[None, :, None] * rng.standard_normal((n, NX, seeds)); aT = template(tru); gT = ga(aT)
    grid = np.logspace(-6, 0, 25)[::-1]
    for qf in QF:
        Qm = Qt * qf; a0, ax = a0_ax(Qm); G0 = atten(a0); Gx = [atten(ax[:, i]) for i in range(NX)]; Gd = [atten(a0 + ax[:, i]) for i in range(NX)]; Y = np.stack([np.linalg.solve(G.T @ G + LX * I, G.T @ Dat[:, i, :]) for i, G in enumerate(Gx)], axis=1)
        chosen = {m: [None] * seeds for m in ("split", "direct", "fc_star")}; out = {m: np.zeros_like(Dat) for m in chosen}
        for lam in grid:
            F0 = np.linalg.solve(G0.T @ G0 + lam * I, G0.T); fs_t = np.array([th.fc_star(a0[k], lam) for k in range(n)])
            O = dict(split=np.stack([F0 @ Y[:, i, :] for i in range(NX)], axis=1), direct=np.stack([ridge(Gd[i], lam, np.ones(n), Dat[:, i, :]) for i in range(NX)], axis=1),
                     fc_star=np.stack([ridge(Gd[i], lam, np.exp(-2 * np.pi * fs_t * ax[:, i]), Dat[:, i, :]) for i in range(NX)], axis=1))
            for m, o in O.items():
                res = sum(np.sum((Gd[i] @ o[:, i, :] - Dat[:, i, :]) ** 2, axis=0) for i in range(NX))
                for s in range(seeds):
                    if chosen[m][s] is None and res[s] <= target2:
                        chosen[m][s] = float(lam); out[m][:, :, s] = o[:, :, s]
        for m, o in O.items():         # no lambda met the target: smallest lambda
            for s in range(seeds):
                if chosen[m][s] is None:
                    chosen[m][s] = float(grid[-1]); out[m][:, :, s] = o[:, :, s]
        for m in chosen:
            A = template(out[m]); meas = ga(A.mean(axis=2)); wave = float(np.mean(np.sqrt(np.sum((out[m] - tru[:, :, None]) ** 2, axis=(0, 1)) / np.sum(tru ** 2))))
            for j in range(3):
                rows.append(dict(M=Mname, snr=snr, qf=qf, m=m, j=j, truth=float(gT[j]), meas=float(meas[j]), wave=wave, lam_med=float(np.median(chosen[m]))))
    return rows


def run_specs(Mname, seeds, base):
    """Descriptive: fc* computed with a wrong wavelet peak frequency (x0.8, x1.2). Q correct, SNR 18."""
    Qt = MODELS[Mname]; idx, R = reflector_series(); dp, tru = synth(idx, R, Qt, False); clean = dp[:n]; rms = np.sqrt(np.sum(dp[:int(2.4 / DT)] ** 2, axis=0) / n); sigma = rms / 10 ** (18 / 20); rng = np.random.default_rng(base + 900 + 1000 * list(MODELS).index(Mname)); Dat = clean[:, :, None] + sigma[None, :, None] * rng.standard_normal((n, NX, seeds)); gT = ga(template(tru)); a0, ax = a0_ax(Qt); Gd = [atten(a0 + ax[:, i]) for i in range(NX)]; rows = []
    for fp in (0.8, 1.0, 1.2):
        spec = ricker_spec2(fp * C.fpeak)
        for lam in (1e-3, 1e-2):
            fs_t = np.array([th.fc_star(a0[k], lam, spec) for k in range(n)]); O = np.stack([ridge(Gd[i], lam, np.exp(-2 * np.pi * fs_t * ax[:, i]), Dat[:, i, :]) for i in range(NX)], axis=1); meas = ga(template(O).mean(axis=2))
            rows += [dict(M=Mname, fp=fp, lam=lam, j=j, bias=float(meas[j] - gT[j]), fc=float(fs_t[EV[j]])) for j in range(3)]
    return rows


if __name__ == "__main__":
    mode = sys.argv[1]; names, snrs, base = (("M1", "M2", "M3"), (60, 18), 31000) if mode == "dev" else (("M4", "M5"), (60, 18), 77000)
    rows = []; t0 = time.monotonic()
    if mode == "test":
        D = []; P = []
        for mname in names:
            D += run_disc(mname, 18, 20, base); P += run_specs(mname, 20, base); print(mname, "disc and spectra", f"{time.monotonic() - t0:.0f}s", flush=True)
        (OUT / "test_disc.json").write_text(json.dumps(D), encoding="utf-8"); (OUT / "test_spec.json").write_text(json.dumps(P), encoding="utf-8")
    for mname in names:
        for snr in snrs:
            rows += run(mname, snr, 20 if snr < 60 else 3, base); print(mname, snr, f"{time.monotonic() - t0:.0f}s", flush=True)
    (OUT / f"{mode}.json").write_text(json.dumps(rows), encoding="utf-8")
    for snr in snrs:
        summary([r for r in rows if r["snr"] == snr], f"\n{mode}, SNR {snr}")
    print("\nExample: M, qf, lam, event | direct truth/meas/pred | fc* meas (fc*)")
    for r in rows:
        if r["m"] == "direct" and r["snr"] == 60 and r["lam"] in (1e-4, 1e-3) and r["M"] == names[0]:
            f = next(x for x in rows if x["m"] == "fc_star" and all(x[k] == r[k] for k in ("M", "snr", "qf", "lam", "j")))
            print(f"{r['M']} {r['qf']} {r['lam']:.0e} {r['j']} | {r['truth']:+.3f} {r['meas']:+.3f} {r['pred']:+.3f} first {r.get('first', float('nan')):+.3f} | {f['meas']:+.3f} ({r['fc_star']:.1f} Hz)")
