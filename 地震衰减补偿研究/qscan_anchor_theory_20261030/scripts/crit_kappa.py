"""Other selection criteria (PLAN_CRIT_ZH.md): does the anchoring coefficient depend on the criterion?

usage: crit_kappa.py    worker; jobs claimed with results17/<job>.claim, resumable.
For every candidate c (QGRID points with 0.52 <= c <= 2.5, i.e. v <= 1.92) one lasso solve as lib.qscan
(lambda 1/2, wavelet at the fixed anchor c_a, default window), warm-started along the grid from small v.
From the same solutions, five selection criteria (each: grid minimum, then the parabola vertex in v):
  J       the QScan lasso objective (reference);
  lp1     ||r||_1 / ||r||_2, r = X / scale (physical reflectivity) -- Ma et al. 2018-type lp sparsity, p = 1,
          made scale-invariant (the raw sum |r| always picks the weakest compensation, a pure scale effect);
  lp05    (sum |r|^0.5)^2 / ||r||_1 (scale-invariant p = 1/2);
  yan     Yan et al. 2023: (1 - phi) eps* + phi sig*, eps = sum residual^2, sig = (sum |r|^0.5)^2,
          min-max over the candidates, phi = 0.5;
  chengD  Cheng 2025-type shallow/deep consistency of the compensated output O = W_out * r: |slope| over
          10-60 Hz of ln(A_deep / A_shallow), A = live-trace mean amplitude spectrum of 1.2-1.8 s / 2.4-3.0 s.
Data: semi-synthetic (experiment-14 recipe, structure 10227, noise line 10149 x 1) c 0.7 / 1.0 / 1.4, seeds
10201-10202, anchors 0.7 / 1.0 / 1.4; field inlines 10141, 10193, 10241, 10305, anchors 0.7 / 1.4.
-> results17/<dataset>__a<c_a>.json
"""
import os, sys, json, time
from env15 import np, lib, semi, sc11, bm, gen, ROOT

sys.path.insert(0, str(ROOT.parent / "field_baselines_volve_20261029" / "scripts"))
import lib13 as L13  # noqa: E402

OUT = ROOT / "results17"; bm.MAXITER = 20000
QG = np.asarray(lib.QGRID, float); IDX = [j for j in range(len(QG)) if QG[j] >= 0.519]; VG = 1.0 / QG
CRITS = ("J", "lp1", "lp05", "yan", "chengD")


def jobs():
    semi_j = [(f"S_c{c}_s{s}", ("S", c, s), a) for s in (10201, 10202) for c in (0.7, 1.0, 1.4) for a in (0.7, 1.0, 1.4)]
    field_j = [(f"F_{ln}", ("F", ln, None), a) for ln in (10141, 10193, 10241, 10305) for a in (0.7, 1.4)]
    return semi_j + field_j


def spec(O, live, t0, t1, nfft=256):
    i = semi.win(semi.TW, t0, t1); X = O[i][:, live] * np.hanning(i.sum())[:, None]
    return np.fft.rfftfreq(nfft, semi.DT), np.sqrt(np.mean(np.abs(np.fft.rfft(X, nfft, axis=0)) ** 2, axis=1))


def scan_all(D, cm, Q, pp, wv):
    sigma, K, live = pp["sigma"], pp["K"], pp["live"]; n = D.shape[0]; s2 = np.where(live, sigma ** 2, np.inf)
    rec = {k: {} for k in ("J", "lp1", "lp05", "eps", "sig", "chengD")}; conv = {}; prev = None
    for j in sorted(IDX, reverse=True):                        # small v first
        cc = dict(cm); cc["Q"] = Q * QG[j]; B, P, scale = lib.dictionary(cc, wv["w_in"], wv["w_out"])
        pen1 = np.sqrt(2 * np.log(n)) * sigma[None, :] * np.median(scale) / scale[:, None]
        X, info = bm.solve(B, np.eye(n), D, sigma, scale, K, 0.5, 0.0, live, start=prev); prev = X
        res = D - B @ X; r = X / scale[:, None]; O = P @ X
        rec["J"][j] = float(np.sum((0.5 * np.sum(res ** 2, axis=0) + 0.5 * np.sum(pen1 * np.abs(X), axis=0))[live] / s2[live]))
        l1 = float(np.sum(np.abs(r[:, live]))); l2 = float(np.sqrt(np.sum(r[:, live] ** 2))); h = float(np.sum(np.sqrt(np.abs(r[:, live]))))
        rec["lp1"][j] = l1 / l2; rec["lp05"][j] = h ** 2 / l1; rec["sig"][j] = h ** 2; rec["eps"][j] = float(np.sum(res[:, live] ** 2))
        f, As = spec(O, live, 1.2, 1.8); _, Ad = spec(O, live, 2.4, 3.0); b = (f >= 10) & (f <= 60)
        rec["chengD"][j] = float(abs(np.polyfit(f[b], np.log(np.maximum(Ad[b], 1e-30) / np.maximum(As[b], 1e-30)), 1)[0]))
        conv[j] = bool(info["converged"])
    e = np.array([rec["eps"][j] for j in IDX]); s = np.array([rec["sig"][j] for j in IDX])
    mm = lambda x: (x - x.min()) / (x.max() - x.min()) if x.max() > x.min() else 0 * x
    yan = 0.5 * mm(e) + 0.5 * mm(s); rec["yan"] = {j: float(v) for j, v in zip(IDX, yan)}
    out = {}
    for k in CRITS:
        vals = np.array([rec[k][j] for j in IDX]); jm = int(np.argmin(vals)); edge = jm in (0, len(IDX) - 1)
        v = float(VG[IDX[jm]])
        if not edge:
            ii = [IDX[jm - 1], IDX[jm], IDX[jm + 1]]; a = np.polyfit(VG[ii], [rec[k][x] for x in ii], 2)
            if a[0] > 0:
                v = float(np.clip(-a[1] / (2 * a[0]), min(VG[ii]), max(VG[ii])))
        out[k] = dict(v_sel=v, v_grid=float(VG[IDX[jm]]), edge=bool(edge),
                      ok=bool(all(conv[x] for x in IDX[max(0, jm - 1):jm + 2])))
    return out, {k: {str(j): v for j, v in d.items()} for k, d in rec.items()}


if __name__ == "__main__":
    OUT.mkdir(exist_ok=True); ctx = gen.prepare(10227, [10149, 10273, 10355]); cm = ctx["cm"]; Q = np.asarray(cm["Q"], float)
    for did, spec_, a in jobs():
        job = f"{did}__a{a}"; path = OUT / f"{job}.json"
        if path.exists():
            continue
        try:
            os.close(os.open(OUT / f"{job}.claim", os.O_CREAT | os.O_EXCL | os.O_WRONLY))
        except FileExistsError:
            continue
        t0 = time.monotonic(); kind, x, s = spec_
        if kind == "S":
            g = gen.signal(ctx, x, s); D = g["x"] + g["level"] * ctx["noise"][0]; u = 1.0 / x
        else:
            D = L13.read_inline(x); u = None
        pp = sc11.prepare_sc(D, cm); wv = lib.wavelets(D, semi.TW, dict(cm, Q=Q * a))
        sel, raw = scan_all(D, cm, Q, pp, wv)
        rec = dict(dataset=did, kind=kind, u=u, c_anchor=a, v_anchor=1 / a, sel=sel, raw=raw, seconds=time.monotonic() - t0)
        path.write_text(json.dumps(rec, indent=1), encoding="utf-8"); (OUT / f"{job}.claim").unlink()
        print(job, {k: round(v["v_sel"], 3) for k, v in sel.items()}, f"{rec['seconds']:.0f}s", flush=True)
