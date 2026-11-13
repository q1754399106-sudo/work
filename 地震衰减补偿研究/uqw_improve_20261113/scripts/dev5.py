"""Step 3, fifth probe (development data only): shaping toward a signal-to-noise weighted target.

dev4.py: half-strength shaping toward |W_out| undoes the noise suppression of UQ-W when noise dominates. Here the target in
each window is |W_out(f)| w(f), w = mean over the window rows of snr / (1 + snr) (the same coherence snr that UQ-W uses), so
that a correctly compensated UQ-W output is left alone. Methods: UQW, VDB, H (uqwh.half), HW50 / HW100 (weighted target,
strength 0.5 / 1). States: semi-synthetic development sets, calibrated and 4x noise, variants O and F (S, T2); FD H1-H6 at
30 and 14 dB (S). usage: dev5.py semi ID [ID ...] | dev5.py fd
"""
import sys, json, time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from probe2 import np, semi, score, sc11, TW, DT, snr_rows, CW, load, ROOT, t2, lib12, SEL, tvfilt, CURVE, MS, lib  # noqa: E402
import uqwh as H  # noqa: E402
fd = CW.fd; OUT = ROOT / "results" / "dev5"; OUT.mkdir(parents=True, exist_ok=True); F = np.fft.rfftfreq(1024, DT)


def half_w(X, w_out, inv, lam):
    """uqwh.half with the target |W_out| w; inv = 1/snr on (rows, rfftfreq(1024))."""
    n = X.shape[0]; L = int(round(H.WIN / DT)); step = int(round(H.HOP / DT)); h = np.hanning(L); f = np.fft.rfftfreq(H.NF, DT)
    live = np.abs(X).sum(axis=0) > 0; Wo = np.abs(np.fft.rfft(w_out, H.NF)); Wo = Wo / Wo.max(); band = (f >= H.LO) & (f <= H.HI)
    out = np.zeros_like(X); ws = np.zeros(n); W = 1.0 / (1.0 + inv)
    for s0 in range(-L // 2, n - L // 2, step):
        i0, i1 = max(s0, 0), min(s0 + L, n); hw = h[i0 - s0:i1 - s0]; Fx = np.fft.rfft(X[i0:i1] * hw[:, None], H.NF, axis=0)
        A = np.sqrt(lib.smooth(np.mean(np.abs(Fx[:, live]) ** 2, axis=1), H.NSM)); A = A / A[band].max()
        T = Wo * np.interp(f, F, (W[i0:i1] * hw[:, None]).sum(0) / hw.sum()); T = T / T[band].max()
        G0 = T * A / (A ** 2 + SEL["B3"]["eps"] * np.max(A ** 2)); G = np.where(band, G0 ** lam, G0)
        out[i0:i1] += np.fft.irfft(Fx * G[:, None], H.NF, axis=0)[:i1 - i0]; ws[i0:i1] += hw
    return out / np.maximum(ws, 1e-3 * ws.max())[:, None]


def five(D, cm, Q, c, tw, anchor, w_out):
    inv = snr_rows(D, tw, F); U = CW.uqw(D, cm, Q, c, H.SIG, tw, anchor, H.GAMMA, inv_snr=inv)
    return dict(UQW=U, VDB=lib12.tvw(fd.iq_t(D, cm, Q, None, c, 0, tw), w_out, SEL["B3"]["eps"], win=SEL["B3"]["win"]), H=H.half(U, w_out),
                HW50=half_w(U, w_out, inv, 0.5), HW100=half_w(U, w_out, inv, 1.0))


if sys.argv[1] == "semi":
    cm = semi.curtain_cm(); Q = np.asarray(cm["Q"], float); fp = float(cm["fpeak"]); r = np.asarray(json.loads(CURVE.read_text(encoding="utf-8"))["r"])
    for did in sys.argv[2:]:
        t0 = time.monotonic(); arr, meta = load(did); rec = []
        for var, sig in (("O", arr["signal"]), ("F", tvfilt(arr["signal"], r))):
            for st, x in (("cal", 1.0), ("nd", 4.0)):
                D = sig + x * arr["noise"]; w_out = sc11.prepare_sc(D, cm)["wv"]["w_out"]
                for m in MS:
                    rec += [dict(var=var, state=st, method=k, m=m, S=score(X, arr["ideal"])["S"], T2=t2(X, arr["R"], fp)) for k, X in five(D, cm, Q, m * meta["c_eff_cum"], TW, 1.4, w_out).items()]
        (OUT / f"semi_{did}.json").write_text(json.dumps(rec), encoding="utf-8"); print(did, "%.0fs" % (time.monotonic() - t0), flush=True)
else:
    for h in fd.MODELS:
        D0, T, cm, tw, fp = fd.load_fd(h); rec = []; t0 = time.monotonic()
        for j, (lv, db) in enumerate((("30dB", 30.0), ("14dB", 14.0))):
            D = D0 + fd.noise(D0, db, 20261113900 + 10 * fd.MODELS.index(h) + j)
            for m in MS:
                rec += [dict(state=lv, method=k, m=m, S=fd.score_fd(X, T, tw)[0]) for k, X in five(D, cm, cm["Q"], m, tw, 0.0, fd.ricker(fp)).items()]
        (OUT / f"fd_{h}.json").write_text(json.dumps(rec), encoding="utf-8"); print(h, "%.0fs" % (time.monotonic() - t0), flush=True)
