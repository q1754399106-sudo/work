"""Differential field test (BETA_ZH.md step 4): add a known attenuation to real inlines and measure each estimator's response.

usage: beta_diff.py    worker; jobs claimed with results14/<job>.claim, resumable.
Lines 10141, 10193, 10241, 10305. The real section (1.2-3.0 s) is passed through the forward attenuation of
Delta u x curtain t* (operator exp(b(f) Delta_u a(t)) per output time, b as lib.dictionary), Delta u = 0, 0.3, 0.6.
Estimators on each section: SR (sq.estimate, curtain model Q); fixed-anchor QScan with the default wavelet window at
c_a = 1 and with the deep window (2.4-2.8 s) at c_a = 0.7 and 1.4 (lambda 1/2, qprofile.scan).
-> results14/IL_<line>_d<Delta_u>.json
"""
import os, sys, json, time
from env15 import np, lib, semi, sc11, bm, ROOT, sq
from anchor2 import wavelets_at
from qprofile import scan

sys.path.insert(0, str(ROOT.parent / "field_baselines_volve_20261029" / "scripts"))
import lib13 as L13  # noqa: E402

OUT = ROOT / "results14"; bm.MAXITER = 20000
LINES = (10141, 10193, 10241, 10305); DUS = (0.0, 0.3, 0.6)


def attenuate(D, cm, du, nfft=1024):
    n = len(semi.TW); f = np.fft.rfftfreq(nfft, semi.DT); lag = np.arange(n); a = lib.absorption(cm, semi.TW)
    b = -np.pi * f + 2j * f * np.log(np.maximum(f, np.finfo(float).tiny) / float(cm["fref"])); b[0] = 0
    M = np.empty((n, n))
    for i in range(n):
        H = np.exp(b * du * a[i]); H[-1] = H[-1].real; M[i] = np.fft.irfft(H, nfft)[(i - lag) % nfft]
    return M @ D


if __name__ == "__main__":
    OUT.mkdir(exist_ok=True); cm = semi.curtain_cm(); Q = np.asarray(cm["Q"], float)
    for ln in LINES:
        for du in DUS:
            job = f"IL_{ln}_d{du}"; path = OUT / f"{job}.json"
            if path.exists():
                continue
            try:
                os.close(os.open(OUT / f"{job}.claim", os.O_CREAT | os.O_EXCL | os.O_WRONLY))
            except FileExistsError:
                continue
            t0 = time.monotonic(); D0 = L13.read_inline(ln); live = np.abs(D0).sum(axis=0) > 0
            D = attenuate(D0, cm, du) * live[None, :] if du else D0
            sr = sq.estimate(D, cm, Q)["c"]; pp = sc11.prepare_sc(D, cm); rec = dict(inline=ln, du=du, c_sr=sr)
            for tag, ca, win in (("w1.4_a1.0", 1.0, None), ("w2.6_a0.7", 0.7, (2.4, 2.8, 2.6)), ("w2.6_a1.4", 1.4, (2.4, 2.8, 2.6))):
                acm = dict(cm, Q=Q * ca); wv = lib.wavelets(D, semi.TW, acm) if win is None else wavelets_at(D, semi.TW, acm, *win)
                r = scan(D, cm, Q, pp, lambda cc: wv); rec[tag] = dict(v_sel=r["v_sel"], edge=r["edge"], ok=r["neighbourhood_converged"])
            rec["seconds"] = time.monotonic() - t0
            path.write_text(json.dumps(rec, indent=1), encoding="utf-8"); (OUT / f"{job}.claim").unlink()
            print(job, f"c_sr {sr:.3f}", {k: round(v["v_sel"], 3) for k, v in rec.items() if isinstance(v, dict)}, f"{rec['seconds']:.0f}s", flush=True)
