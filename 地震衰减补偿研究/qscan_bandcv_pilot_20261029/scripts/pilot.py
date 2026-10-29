"""Pilot (exploratory, descriptive): select the Q factor by band cross-validation instead of the full-band lasso.

For each factor c in QGRID: solve the lasso of lib.qscan on the low band only (data and dictionary passed through the
same zero-phase low-pass F_L), predict the high band F_H B_c X_c, and score the high-band misfit. The low band carries
little of the attenuation contrast, so the reflectivity it yields cannot absorb a wrong Q in the high band.

Bands (fixed before any run): low = 3-40 Hz with a cosine taper 37.5-42.5 Hz; high = its complement up to the
target-band limit fmax_target of lib.wavelets (5 Hz taper below it). Low-band noise sigma = full-band sigma times the
white-noise energy fraction of the low band. Weights LAMS include 0.5 (the QScan weight). Scores per c and weight:
  E_raw  = sum over live traces of ||F_H D - F_H B X||^2 / sigma^2
  E_gain = the same after one least-squares gain g >= 0 on the prediction (removes lasso shrinkage of the level)
Data: the 20 held-out semi-synthetic sections of experiment 12 (true c 0.65-1.5), same files as qscan_reg_pilot.

usage: pilot.py K N [i,j,...]    worker K of N, or only the listed dataset indices; writes results/<id>.json
(run with 10 workers; workers 4-9 were stopped after their first dataset to free memory, and indices 14-19
were then run by one process)
"""
import sys, json, time
from env13 import np, lib, semi, sc11, bm, ROOT, datasets

LAMS = [0.5, 0.125, 2 ** -5]
F_SPLIT, TAPER, NFFT = 40.0, 5.0, 1024
OUT = ROOT / "results"


def gains(dt, fd):
    f = np.fft.rfftfreq(NFFT, dt)
    gl = np.clip((F_SPLIT + TAPER / 2 - f) / TAPER, 0, 1); gl = 0.5 * (1 - np.cos(np.pi * gl)); gl[f < 3] = 0
    top = np.clip((fd - f) / TAPER, 0, 1); top = 0.5 * (1 - np.cos(np.pi * top))
    gh = (1 - gl) * top; gh[f < 3] = 0
    return f, gl, gh


def band_matrix(g, n):
    """n x n matrix of the zero-phase filter with gain g (zero-padded FFT, no wrap-around)."""
    I = np.eye(n); return np.fft.irfft(np.fft.rfft(I, NFFT, axis=0) * g[:, None], NFFT, axis=0)[:n]


def scan(D, cm, Q, pp):
    wv, sigma, K, live = pp["wv"], pp["sigma"], pp["K"], pp["live"]; n = D.shape[0]
    f, gl, gh = gains(float(cm["dt"]), wv["fmax_target"])
    FL, FH = band_matrix(gl, n), band_matrix(gh, n)
    sig_l = sigma * np.sqrt(np.mean(gl ** 2)); sig_h = sigma * np.sqrt(np.mean(gh ** 2))
    DL, DH = FL @ D, FH @ D; s2 = np.where(live, sig_h ** 2, np.inf)
    shape = (len(LAMS), len(lib.QGRID))
    Er, Eg, G, C, ITS = np.zeros(shape), np.zeros(shape), np.zeros(shape), np.zeros(shape, bool), np.zeros(shape, int)
    for j, c in enumerate(lib.QGRID):
        cc = dict(cm); cc["Q"] = np.asarray(Q, float) * c
        B, _, scale = lib.dictionary(cc, wv["w_in"], wv["w_out"])
        BL, BH = FL @ B, FH @ B
        for i, lam in enumerate(LAMS):
            X, info = bm.solve(BL, np.eye(n), DL, sig_l, scale, K, lam, 0.0, live)
            P = BH @ X; r = DH - P
            Er[i, j] = float(np.sum(np.sum(r ** 2, axis=0)[live] / s2[live]))
            num = float(np.sum((np.sum(DH * P, axis=0))[live] / s2[live])); den = float(np.sum(np.sum(P * P, axis=0)[live] / s2[live]))
            g = max(num / den, 0.0) if den > 0 else 0.0
            Eg[i, j] = float(np.sum(np.sum((DH - g * P) ** 2, axis=0)[live] / s2[live]))
            G[i, j], C[i, j], ITS[i, j] = g, bool(info["converged"]), int(info["iterations"])
    return dict(E_raw=Er, E_gain=Eg, gain=G, converged=C, iterations=ITS,
                bands=dict(split=F_SPLIT, taper=TAPER, fmax_target=float(wv["fmax_target"]), f_peak=float(wv["f_peak"])))


def main(k, n, only=None):
    OUT.mkdir(parents=True, exist_ok=True); cm = semi.curtain_cm(); Q = np.asarray(cm["Q"], float)
    for i, p in enumerate(datasets()):
        did = p.stem; path = OUT / f"{did}.json"
        if (i % n != k if only is None else i not in only) or path.exists():
            continue
        t0 = time.monotonic(); D = np.load(p)["Input"].astype(float)
        pp = sc11.prepare_sc(D, cm); s = scan(D, cm, Q, pp)
        sel = {m: {str(lam): float(lib.QGRID[int(np.argmin(s[m][i]))]) for i, lam in enumerate(LAMS)} for m in ("E_raw", "E_gain")}
        rec = dict(id=did, c_true=float(did.split("_")[1][1:]), grid=lib.QGRID.tolist(), lams=LAMS, selected=sel,
                   bands=s["bands"], seconds=time.monotonic() - t0,
                   **{m: s[m].tolist() for m in ("E_raw", "E_gain", "gain", "converged", "iterations")})
        path.write_text(json.dumps(rec, indent=1), encoding="utf-8"); print(did, f"{rec['seconds']:.0f}s", flush=True)


if __name__ == "__main__":
    main(int(sys.argv[1]), int(sys.argv[2]), [int(x) for x in sys.argv[3].split(",")] if len(sys.argv) > 3 else None)
