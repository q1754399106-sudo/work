"""Diagnostic (post hoc, no well tie, no score): where does the energy near Nyquist in UQ-W / IQ70 come from?
Prints, for one inline, the input power, 1/snr of the coherence estimate and the UQ-W gain |H| by frequency at three times.

usage: diag_nyq.py [INLINE]
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import run as R  # noqa: E402
from run import np, L, lib, semi, TW, CW, F_OP  # noqa: E402
sys.path.insert(0, str(R.BASE / "qunknown_comp_20261109" / "scripts"))
import uqw as U  # noqa: E402

ln = int(sys.argv[1]) if len(sys.argv) > 1 else 10105; cm = semi.curtain_cm(); Q = np.asarray(cm["Q"], float); D = L.read_inline(ln); c = L.c_sr(ln)
inv = CW.snr_rows(D, TW, F_OP); x = lib.absorption(dict(cm, Q=Q * c), TW); dx = x - lib.absorption(dict(cm, Q=Q * c), [1.4])[0]; s = semi.gain_sigma2(U.G_IQ)
lf = np.log(np.maximum(F_OP, np.finfo(float).tiny) / float(cm["fref"])); H = np.abs(U.rows_w(x, dx, 0.447, F_OP, lf, s, inv, 0.0))
print("inline %d  DT %.4f  n_t %d  c_sr %.3f" % (ln, R.DT, len(TW), c)); fs = [10, 30, 50, 70, 90, 100, 110, 118, 122, 124, 125]
for t in (1.6, 2.2, 2.8):
    i = int(np.argmin(np.abs(TW - t))); k = (TW >= t - 0.2) & (TW <= t + 0.2); P = np.mean(np.abs(np.fft.rfft(D[k] * np.hanning(k.sum())[:, None], 1024, axis=0)) ** 2, axis=1); P = 10 * np.log10(P / P.max() + 1e-30)
    print("t %.1f s  t* %.4f" % (t, x[i]))
    for name, v in (("input power dB", P), ("snr (1/inv)   ", 1 / inv[i]), ("gain |H| dB   ", 20 * np.log10(H[i] + 1e-30))):
        print("   %s " % name + " ".join("%8.1f" % v[int(np.argmin(np.abs(F_OP - f)))] for f in fs))
print("   f (Hz)         " + " ".join("%8d" % f for f in fs))
