"""Step 3, sixth probe (development data only): shaping toward the anchor-window spectrum of the UQ-W output itself.

diag_nd.py: with 4x noise the wavelet estimate w_out collapses above about 45 Hz, so shaping toward |W_out| low-passes the
UQ-W output (this, not a re-amplification of noise, is what dev4.py / dev5.py showed). Here the target is the section-mean
amplitude spectrum of the UQ-W output in the anchor window 1.2-1.6 s (11-bin smoothing), where little compensation has been
applied; no separate wavelet estimate is needed.
  HA   G = (T / A)^0.5 inside 8-72 Hz, cosine taper to zero outside (dev3.taper)
  HAW  the same with T multiplied by w(t, f) / w(anchor, f), w = snr / (1 + snr), ratio clipped to [0.25, 4]
Methods also: UQW, VDB, H (uqwh.half). Semi-synthetic development sets, calibrated and 4x noise, variants O and F.
usage: dev6.py ID [ID ...]
"""
import sys, json, time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from probe2 import np, semi, score, sc11, TW, DT, snr_rows, CW, load, ROOT, t2, lib12, SEL, tvfilt, CURVE, MS, lib  # noqa: E402
import uqwh as H  # noqa: E402
from dev3 import taper  # noqa: E402
fd = CW.fd; OUT = ROOT / "results" / "dev6"; OUT.mkdir(parents=True, exist_ok=True); F = np.fft.rfftfreq(1024, DT); ANC = (1.2, 1.6)


def anchor_shape(X, tw, inv=None, lam=0.5):
    n = X.shape[0]; L = int(round(H.WIN / DT)); step = int(round(H.HOP / DT)); h = np.hanning(L); f = np.fft.rfftfreq(H.NF, DT)
    live = np.abs(X).sum(axis=0) > 0; band = (f >= H.LO) & (f <= H.HI); tp = taper(f); out = np.zeros_like(X); ws = np.zeros(n)
    amp = lambda seg: np.sqrt(lib.smooth(np.mean(np.abs(np.fft.rfft(seg, H.NF, axis=0)[:, live]) ** 2, axis=1), H.NSM))
    ka = (tw >= ANC[0]) & (tw <= ANC[1]); T0 = amp(X[ka] * np.hanning(ka.sum())[:, None]); T0 = T0 / T0[band].max()
    W = None if inv is None else 1.0 / (1.0 + inv); wa = None if inv is None else np.interp(f, F, W[ka].mean(0))
    for s0 in range(-L // 2, n - L // 2, step):
        i0, i1 = max(s0, 0), min(s0 + L, n); hw = h[i0 - s0:i1 - s0]; seg = X[i0:i1] * hw[:, None]; A = amp(seg); A = A / A[band].max(); T = T0
        if inv is not None:
            T = T0 * np.clip(np.interp(f, F, (W[i0:i1] * hw[:, None]).sum(0) / hw.sum()) / wa, 0.25, 4.0); T = T / T[band].max()
        G = (T / np.maximum(A, 1e-3)) ** lam * tp
        out[i0:i1] += np.fft.irfft(np.fft.rfft(seg, H.NF, axis=0) * G[:, None], H.NF, axis=0)[:i1 - i0]; ws[i0:i1] += hw
    return out / np.maximum(ws, 1e-3 * ws.max())[:, None]


if __name__ == "__main__":
    cm = semi.curtain_cm(); Q = np.asarray(cm["Q"], float); fp = float(cm["fpeak"]); r = np.asarray(json.loads(CURVE.read_text(encoding="utf-8"))["r"])
    for did in sys.argv[1:]:
        t0 = time.monotonic(); arr, meta = load(did); rec = []
        for var, sig in (("O", arr["signal"]), ("F", tvfilt(arr["signal"], r))):
            for st, x in (("cal", 1.0), ("nd", 4.0)):
                D = sig + x * arr["noise"]; w_out = sc11.prepare_sc(D, cm)["wv"]["w_out"]; inv = snr_rows(D, TW, F)
                for m in MS:
                    c = m * meta["c_eff_cum"]; U = CW.uqw(D, cm, Q, c, H.SIG, TW, 1.4, H.GAMMA, inv_snr=inv)
                    o = dict(UQW=U, VDB=lib12.tvw(fd.iq_t(D, cm, Q, None, c, 0, TW), w_out, SEL["B3"]["eps"], win=SEL["B3"]["win"]), H=H.half(U, w_out),
                             HA=anchor_shape(U, TW), HAW=anchor_shape(U, TW, inv))
                    rec += [dict(var=var, state=st, method=k, m=m, S=score(X, arr["ideal"])["S"], T2=t2(X, arr["R"], fp)) for k, X in o.items()]
        (OUT / f"{did}.json").write_text(json.dumps(rec), encoding="utf-8"); print(did, "%.0fs" % (time.monotonic() - t0), flush=True)
