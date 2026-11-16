"""Diagnostic (development set ZA_c0.75_r0, m = 1): why does shaping undo UQ-W when noise dominates?
Spectra in 2.3-2.8 s of the target wavelet, the snr weight, the ideal section, the UQ-W output and its signal / noise parts."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from probe2 import np, semi, score, sc11, TW, DT, snr_rows, CW, load, t2, lib  # noqa: E402
import uqwh as H  # noqa: E402
cm = semi.curtain_cm(); Q = np.asarray(cm["Q"], float); fp = float(cm["fpeak"]); F = np.fft.rfftfreq(1024, DT); f = np.fft.rfftfreq(512, DT)
arr, meta = load("ZA_c0.75_r0"); c = meta["c_eff_cum"]; k = (TW >= 2.3) & (TW <= 2.8)
spec = lambda X: np.sqrt(lib.smooth(np.mean(np.abs(np.fft.rfft(X[k] * np.hanning(k.sum())[:, None], 512, axis=0)) ** 2, axis=1), 11))
fr = [int(np.argmin(np.abs(f - x))) for x in (10, 20, 30, 40, 50, 60, 70, 90, 110)]; row = lambda a: " ".join("%6.2f" % v for v in a[fr]); nrm = lambda a: a / a[fr[:7]].max()
for st, x in (("cal", 1.0), ("nd", 4.0)):
    D = arr["signal"] + x * arr["noise"]; w_out = sc11.prepare_sc(D, cm)["wv"]["w_out"]; inv = snr_rows(D, TW, F)
    U = CW.uqw(D, cm, Q, c, H.SIG, TW, 1.4, H.GAMMA, inv_snr=inv); Us = CW.uqw(arr["signal"], cm, Q, c, H.SIG, TW, 1.4, H.GAMMA, inv_snr=inv); Hh = H.half(U, w_out)
    Wo = np.abs(np.fft.rfft(w_out, 512)); i = int(np.argmin(np.abs(TW - 2.55)))
    print(st, "f (Hz)                 ", " ".join("%6d" % f[j] for j in fr))
    print(st, "target |W_out|         ", row(nrm(Wo))); print(st, "weight snr/(1+snr)     ", row(np.interp(f, F, 1 / (1 + inv[i]))))
    print(st, "true snr amp (input)   ", row(spec(arr["signal"]) / spec(x * arr["noise"])))
    print(st, "ideal                  ", row(nrm(spec(arr["ideal"])))); print(st, "UQ-W output            ", row(nrm(spec(U))))
    print(st, "UQ-W signal/noise parts", row(spec(Us) / spec(U - Us))); print(st, "H output               ", row(nrm(spec(Hh))))
    for name, X in (("UQW", U), ("H", Hh), ("UQW signal part", Us)):
        print(st, "%-16s S %.3f T2 %.3f" % (name, score(X, arr["ideal"])["S"], t2(X, arr["R"], fp)))
