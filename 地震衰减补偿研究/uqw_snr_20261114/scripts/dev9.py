"""Step 4 (development data, exploratory): fix the shaping target instead of the operator.

dev6.py: under 4x noise the target wavelet w_out of lib.wavelets collapses (deep edge fd about 42 Hz; shallow edge fs stops at
the noise notch), and shaping toward it removes what UQ-W recovered; with the 1x target, shaping after UQ-W helps. Targets
here use only the noisy input:
  w_in_r   input wavelet of sca_ceiling libc.wavelets(robust=True): gap-tolerant shallow edge fs, NO deep cut
  w_out_r  output wavelet of the same (deep cut min(fd, fs) kept)
Methods (lambda 0.5 = uqwh.half unless noted):
  UQW, H (frozen target, as dev4/5), Hr (w_out_r), Hin (w_in_r), HWin (target |w_in_r| sqrt(w), w = snr / (1 + snr) averaged
  over the window: the expected output amplitude of a Wiener filter), VDB (frozen), VDBr (w_out_r), VDBin (w_in_r).
States: semi O / F x res1, res4, gau4; FD H1-H6 at 30 and 14 dB (target Ricker: Hr = Hin = H there, only HWin differs).
usage: dev9.py semi ID [ID ...] | dev9.py fd  -> results/dev9/
"""
import sys, json, time
from env import np, lib, semi, score, sc11, TW, DT, snr_rows, CW, load, t2, lib12, SEL, tvfilt, CURVE, MS, H, DEV, ROOT, F10, noise_of, fd, BASE
sys.path.insert(0, str(BASE / "fault_q_confirm_20261030" / "scripts")); sys.path.insert(0, str(BASE / "sca_ceiling_20261030" / "scripts"))
import libc  # noqa: E402

OUT = ROOT / "results" / "dev9"; OUT.mkdir(parents=True, exist_ok=True)


def half_t(X, w_t, inv, lam=0.5, sqrtw=True):
    """uqwh.half toward |W_t| (w)^(1/2) per window (sqrtw) or |W_t| (not)."""
    n = X.shape[0]; L = int(round(H.WIN / DT)); step = int(round(H.HOP / DT)); h = np.hanning(L); f = np.fft.rfftfreq(H.NF, DT)
    live = np.abs(X).sum(axis=0) > 0; Wo = np.abs(np.fft.rfft(w_t, H.NF)); Wo = Wo / Wo.max(); band = (f >= H.LO) & (f <= H.HI)
    out = np.zeros_like(X); ws = np.zeros(n); W = 1.0 / (1.0 + inv)
    for s0 in range(-L // 2, n - L // 2, step):
        i0, i1 = max(s0, 0), min(s0 + L, n); hw = h[i0 - s0:i1 - s0]; Fx = np.fft.rfft(X[i0:i1] * hw[:, None], H.NF, axis=0)
        A = np.sqrt(lib.smooth(np.mean(np.abs(Fx[:, live]) ** 2, axis=1), H.NSM)); A = A / A[band].max()
        T = Wo * (np.sqrt(np.interp(f, F10, (W[i0:i1] * hw[:, None]).sum(0) / hw.sum())) if sqrtw else 1.0); T = T / T[band].max()
        G0 = T * A / (A ** 2 + SEL["B3"]["eps"] * np.max(A ** 2)); G = np.where(band, G0 ** lam, G0)
        out[i0:i1] += np.fft.irfft(Fx * G[:, None], H.NF, axis=0)[:i1 - i0]; ws[i0:i1] += hw
    return out / np.maximum(ws, 1e-3 * ws.max())[:, None]


def methods(D, cm, Q, c, tw, anchor, w_out, w_in_r, w_out_r, inv):
    X = CW.uqw(D, cm, Q, c, H.SIG, tw, anchor, H.GAMMA, inv_snr=inv); ph = fd.iq_t(D, cm, Q, None, c, 0, tw); v = lambda w: lib12.tvw(ph, w, SEL["B3"]["eps"], win=SEL["B3"]["win"])
    return dict(UQW=X, H=H.half(X, w_out), Hr=H.half(X, w_out_r), Hin=H.half(X, w_in_r), HWin=half_t(X, w_in_r, inv),
                VDB=v(w_out), VDBr=v(w_out_r), VDBin=v(w_in_r))


if sys.argv[1] == "semi":
    cm = semi.curtain_cm(); Q = np.asarray(cm["Q"], float); fp = float(cm["fpeak"]); r = np.asarray(json.loads(CURVE.read_text(encoding="utf-8"))["r"])
    for did in sys.argv[2:]:
        t0 = time.monotonic(); arr, meta = load(did); rec = []; sd = 20261114000 + DEV.index(did)
        for var, sig in (("O", arr["signal"]), ("F", tvfilt(arr["signal"], r))):
            for st, kind, x in (("res1", "res", 1.0), ("res4", "res", 4.0), ("gau4", "gau", 4.0)):
                D = sig + noise_of(arr, kind, x, sd); inv = snr_rows(D, TW, F10); wv = lib.wavelets(D, TW, cm); wr = libc.wavelets(D, cm, True)
                for m in MS:
                    o = methods(D, cm, Q, m * meta["c_eff_cum"], TW, 1.4, wv["w_out"], wr["w_in"], wr["w_out"], inv)
                    rec += [dict(state=var + st, method=k, m=m, S=score(Y, arr["ideal"])["S"], T2=t2(Y, arr["R"], fp)) for k, Y in o.items()]
                rec.append(dict(state=var + st, method="_bands", m=0, fs=wv["fmax_shallow"], fd=wv["fmax_target"], fs_r=wr["fmax_shallow"], fd_r=wr["fmax_target"]))
        (OUT / f"semi_{did}.json").write_text(json.dumps(rec), encoding="utf-8"); print(did, "%.0fs" % (time.monotonic() - t0), flush=True)
else:
    for h in fd.MODELS:
        D0, T, cm, tw, fp = fd.load_fd(h); rec = []; t0 = time.monotonic(); rk = fd.ricker(fp)
        for j, (lv, db) in enumerate((("30dB", 30.0), ("14dB", 14.0))):
            D = D0 + fd.noise(D0, db, 20261113900 + 10 * fd.MODELS.index(h) + j); inv = snr_rows(D, tw, F10)
            for m in MS:
                rec += [dict(state=lv, method=k, m=m, S=fd.score_fd(Y, T, tw)[0]) for k, Y in methods(D, cm, cm["Q"], m, tw, 0.0, rk, rk, rk, inv).items() if k in ("UQW", "H", "HWin", "VDB")]
        (OUT / f"fd_{h}.json").write_text(json.dumps(rec), encoding="utf-8"); print(h, "%.0fs" % (time.monotonic() - t0), flush=True)
