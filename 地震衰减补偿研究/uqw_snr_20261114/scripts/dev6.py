"""Step 1b (development data, exploratory): why does shaping after UQ-W fail when noise dominates? Oracle decomposition.

4x noise (kind res or gau), variant O. UQ-W is linear for a fixed 1/snr, so its signal and noise parts Xs, Xn are known.
  UQW        uqw with the coherence snr              UQWo   uqw with the true snr(t, f)
  H          uqwh.half (lambda 0.5, total spectrum, target |W_out| estimated from the noisy data)
  H_tc       as H, target wavelet estimated from the 1x-noise data (clean target)
  H_s        as H, measured spectrum = the signal part Xs only (oracle signal spectrum)
  H_s_tc     both
  WS         oracle Wiener shaping, lambda 1: G = |W_out| As / (As^2 + An^2), As, An the spectra of Xs, Xn (clean target)
  WS50       the same filter to the power 0.5
  VDB_e      VDB with eps e (1e-6 is the frozen value)
usage: dev6.py KIND ID [ID ...] -> results/dev6/<KIND>_<ID>.json
"""
import sys, json, time
from env import np, lib, semi, score, sc11, TW, DT, snr_rows, CW, load, t2, lib12, SEL, MS, H, DEV, ROOT, F10, WIN, noise_of, spectra, cross, starts, fd

EPS = (1e-6, 1e-4, 1e-3, 1e-2, 3e-2, 0.1, 0.3)


def true_inv(S, N):
    cen, inv = [], []
    for a in starts(TW):
        Ps, Pn = cross(spectra(S, TW, a, a + WIN), 0), cross(spectra(N, TW, a, a + WIN), 0); cen.append(a + WIN / 2)
        inv.append(np.interp(F10, np.fft.rfftfreq(512, DT), Pn / np.maximum(Ps, 1e-30 * Ps.max())))
    inv = np.minimum(np.array(inv), 1e3); cen = np.array(cen)
    return np.stack([np.array([np.interp(t, cen, inv[:, k]) for k in range(len(F10))]) for t in TW])


def shape(X, w_out, lam, Xa=None, Xn=None):
    """uqwh.half generalised: spectrum measured on Xa (default X); with Xn the Wiener form Wo As / (As^2 + An^2)."""
    Xa = X if Xa is None else Xa
    n = X.shape[0]; L = int(round(H.WIN / DT)); step = int(round(H.HOP / DT)); h = np.hanning(L); f = np.fft.rfftfreq(H.NF, DT)
    Wo = np.abs(np.fft.rfft(w_out, H.NF)); Wo = Wo / Wo.max(); band = (f >= H.LO) & (f <= H.HI); out = np.zeros_like(X); ws = np.zeros(n)
    sp = lambda Y, i0, i1, hw: np.sqrt(lib.smooth(np.mean(np.abs(np.fft.rfft(Y[i0:i1] * hw[:, None], H.NF, axis=0)) ** 2, axis=1), H.NSM))
    for s0 in range(-L // 2, n - L // 2, step):
        i0, i1 = max(s0, 0), min(s0 + L, n); hw = h[i0 - s0:i1 - s0]; Fx = np.fft.rfft(X[i0:i1] * hw[:, None], H.NF, axis=0)
        A = sp(Xa, i0, i1, hw); nrm = A[band].max(); A = A / nrm; An = 0 if Xn is None else sp(Xn, i0, i1, hw) / nrm
        G0 = Wo * A / (A ** 2 + An ** 2 + SEL["B3"]["eps"] * np.max(A ** 2)); G = np.where(band, G0 ** lam, G0)
        out[i0:i1] += np.fft.irfft(Fx * G[:, None], H.NF, axis=0)[:i1 - i0]; ws[i0:i1] += hw
    return out / np.maximum(ws, 1e-3 * ws.max())[:, None]


if __name__ == "__main__":
    kind = sys.argv[1]; cm = semi.curtain_cm(); Q = np.asarray(cm["Q"], float); fp = float(cm["fpeak"]); out = ROOT / "results" / "dev6"; out.mkdir(parents=True, exist_ok=True)
    for did in sys.argv[2:]:
        t0 = time.monotonic(); arr, meta = load(did); S = arr["signal"]; N = noise_of(arr, kind, 4.0, 20261114000 + DEV.index(did)); D = S + N
        w_out = sc11.prepare_sc(D, cm)["wv"]["w_out"]; w_cl = sc11.prepare_sc(S + noise_of(arr, kind, 1.0, 20261114000 + DEV.index(did)), cm)["wv"]["w_out"]
        inv = snr_rows(D, TW, F10); inv_o = true_inv(S, N); rec = []
        for m in MS:
            c = m * meta["c_eff_cum"]; u = lambda Y, iv: CW.uqw(Y, cm, Q, c, H.SIG, TW, 1.4, H.GAMMA, inv_snr=iv)
            Xs, Xn = u(S, inv), u(N, inv); X = Xs + Xn; ph = fd.iq_t(D, cm, Q, None, c, 0, TW)
            o = dict(UQW=X, UQWo=u(D, inv_o), H=shape(X, w_out, 0.5), H_tc=shape(X, w_cl, 0.5), H_s=shape(X, w_out, 0.5, Xa=Xs), H_s_tc=shape(X, w_cl, 0.5, Xa=Xs),
                     WS=shape(X, w_cl, 1.0, Xa=Xs, Xn=Xn), WS50=shape(X, w_cl, 0.5, Xa=Xs, Xn=Xn))
            for e in EPS:
                o["VDB_%g" % e] = lib12.tvw(ph, w_out, e, win=SEL["B3"]["win"])
            rec += [dict(method=k, m=m, S=score(Y, arr["ideal"])["S"], T2=t2(Y, arr["R"], fp)) for k, Y in o.items()]
        (out / f"{kind}_{did}.json").write_text(json.dumps(rec), encoding="utf-8"); print(kind, did, "%.0fs" % (time.monotonic() - t0), flush=True)
