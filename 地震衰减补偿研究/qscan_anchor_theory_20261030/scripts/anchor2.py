"""Anchor-window experiment (PLAN_WINDOW_ZH.md): QScan with the input wavelet estimated in a deeper window.

usage: anchor2.py    worker; (dataset, window, anchor) jobs, claimed with results2/<job>.claim, resumable.
wavelets_at is lib.wavelets with the shallow window [t1, t2] and the absorption time tc as parameters.
-> results2/<id>__w<tc>__a<c_a>.json
"""
import os, json, time
from env15 import np, lib, semi, sc11, bm, ROOT, E12OUT

LAM = 0.5; ANCHORS = [0.7, 1.4]; WINDOWS = [(1.8, 2.2, 2.0), (2.4, 2.8, 2.6)]; MAXITER = 20000
IDS = [f"T_c{c}_s{s}" for c in (0.65, 0.85, 1.05, 1.25, 1.5) for s in (9701, 9702)]
OUT = ROOT / "results2"; bm.MAXITER = MAXITER


def wavelets_at(D, tw, cm, t1, t2, tc, nfft=512):
    """lib.wavelets with the wavelet window [t1, t2] and absorption at tc instead of 1.20-1.60 s and 1.40 s."""
    f, snr_s, P_s = lib.coherence_snr(D, tw, t1, t2, nfft)
    _, snr_d, _ = lib.coherence_snr(D, tw, 2.00, 3.00, nfft)
    amp = np.sqrt(P_s) / np.exp(-np.pi * f * lib.absorption(cm, [tc])[0])
    fpk = float(f[np.argmax(np.sqrt(P_s) * (f >= 5))])
    fs = lib.snr_band(f, snr_s, fpk); fd = lib.snr_band(f, snr_d, fpk)

    def cut(fc):
        g = np.ones_like(f); r = (f > fc - 5) & (f < fc)
        g[r] = 0.5 * (1 + np.cos(np.pi * (f[r] - (fc - 5)) / 5)); g[f >= fc] = 0; g[f < 3] = 0
        return g

    a_in = amp * cut(fs); a_out = a_in * cut(min(fd, fs))

    def zp(a):
        w = np.fft.irfft(a, nfft); w = np.roll(w, 100)[:201]
        return w / np.max(np.abs(w))

    w_in, w_out = zp(a_in), zp(a_out)
    w_out *= np.max(np.abs(np.fft.rfft(w_in, nfft))) / np.max(np.abs(np.fft.rfft(w_out, nfft)))
    return dict(w_in=w_in, w_out=w_out, fmax_shallow=fs, fmax_target=fd)


def scan(D, cm, Q, ca, win, pp):
    sigma, K, live = pp["sigma"], pp["K"], pp["live"]; n = D.shape[0]
    ca_cm = dict(cm); ca_cm["Q"] = Q * ca; wv = wavelets_at(D, semi.TW, ca_cm, *win)
    J = np.zeros(len(lib.QGRID)); C = np.zeros(J.shape, bool); IT = np.zeros(J.shape, int)
    s2 = np.where(live, sigma ** 2, np.inf); prev = None
    for j, c in enumerate(lib.QGRID):
        cc = dict(cm); cc["Q"] = Q * c
        B, _, scale = lib.dictionary(cc, wv["w_in"], wv["w_out"])
        pen1 = np.sqrt(2 * np.log(n)) * sigma[None, :] * np.median(scale) / scale[:, None]
        X, info = bm.solve(B, np.eye(n), D, sigma, scale, K, LAM, 0.0, live, start=prev); prev = X
        res = D - B @ X
        J[j] = float(np.sum((0.5 * np.sum(res ** 2, axis=0) + LAM * np.sum(pen1 * np.abs(X), axis=0))[live] / s2[live]))
        C[j], IT[j] = bool(info["converged"]), int(info["iterations"])
    return float(lib.QGRID[int(np.argmin(J))]), dict(J=J.tolist(), converged=C.tolist(), iterations=IT.tolist(),
                                                     fmax_shallow=wv["fmax_shallow"])


if __name__ == "__main__":
    OUT.mkdir(exist_ok=True); cm = semi.curtain_cm(); Q = np.asarray(cm["Q"], float)
    for did in IDS:
        for win in WINDOWS:
            for ca in ANCHORS:
                job = f"{did}__w{win[2]}__a{ca}"; path = OUT / f"{job}.json"
                if path.exists():
                    continue
                try:
                    os.close(os.open(OUT / f"{job}.claim", os.O_CREAT | os.O_EXCL | os.O_WRONLY))
                except FileExistsError:
                    continue
                t0 = time.monotonic(); D = np.load(E12OUT / f"{did}.npz")["Input"].astype(float); pp = sc11.prepare_sc(D, cm)
                sel, diag = scan(D, cm, Q, ca, win, pp)
                rec = dict(id=did, c_true=float(did.split("_")[1][1:]), window=list(win), c_anchor=ca, lam=LAM, c_sel=sel,
                           diag=diag, seconds=time.monotonic() - t0)
                path.write_text(json.dumps(rec, indent=1), encoding="utf-8"); (OUT / f"{job}.claim").unlink()
                print(job, sel, f"{rec['seconds']:.0f}s", flush=True)
