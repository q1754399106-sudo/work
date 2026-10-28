"""Shared code of experiment 12: the improved pipeline of experiment 11 (SR + SC-a with c_SR x 1.3 + phase fix
to c_SR) against baselines that keep one Q for amplitude and phase and get their caution elsewhere.

Imported unchanged from cache/inputs: experiment 11 sc11.py (SC-a, phase_fix, M1) and everything it imports
(experiments 7-10). Their paths are redirected here. New code: iq (Wang-type stabilised inverse Q with separate
amplitude and phase factors), tvw (zero-phase time-varying Wiener shaping, as van der Baan 2012).
"""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
IN = ROOT / "cache" / "inputs"
sys.path.insert(0, str(ROOT.parent / "novelty_evidence_20261025" / "cache" / "deps"))   # scipy
for p in ("novelty_evidence_20261025/scripts", "field_application_20261028/scripts", "semisynthetic_volve_20261028/scripts",
          "specslope_q_volve_20261028/scripts", "multiwell_tie_volve_20261028/scripts", "srsc_volve_20261028/scripts"):
    sys.path.insert(0, str(IN / p))
import sc11  # noqa: E402
from sc11 import np, lib, semi, sq, gen, run, TW, DT  # noqa: E402

for mod in (sc11, sq, semi):
    mod.ROOT, mod.IN = ROOT, IN
sc11.wq.ROOT, sc11.wq.IN, sc11.wq.E9 = ROOT, IN, IN / "multiwell_tie_volve_20261028"
E10 = sc11.E10 = IN / "specslope_q_volve_20261028"
E11 = IN / "srsc_volve_20261028"
DEV_E10 = ROOT.parent / "specslope_q_volve_20261028" / "cache" / "out"
DEV_E11 = ROOT.parent / "srsc_volve_20261028" / "cache" / "dev"


def iq(D, cm, Q, c_amp, c_ph, g_db, nfft=1024):
    """semi.inverse_q with the amplitude term on t*(Q c_amp) and the dispersion term on t*(Q c_ph);
    c_amp=None gives the phase-only filter. c_amp = c_ph reproduces semi.inverse_q."""
    Q = np.asarray(Q, float); n = len(TW); f = np.fft.rfftfreq(nfft, DT); lag = np.arange(n)
    tp = lib.absorption(dict(cm, Q=Q * c_ph), TW)
    ta = None if c_amp is None else lib.absorption(dict(cm, Q=Q * c_amp), TW)
    s = None if c_amp is None else semi.gain_sigma2(g_db); M = np.empty((n, n))
    lf = np.log(np.maximum(f, np.finfo(float).tiny) / float(cm["fref"]))
    for i in range(n):
        H = np.exp(-2j * f * tp[i] * lf)
        if ta is not None:
            beta = np.exp(-np.pi * f * ta[i]); H = H * (beta + s) / (beta ** 2 + s)
        H[0] = H[0].real; H[-1] = H[-1].real
        M[i] = np.fft.irfft(H, nfft)[(i - lag) % nfft]
    return M @ D


def tvw(X, w_out, eps, win=0.4, hop=0.1, nfft=512, nsm=11):
    """Zero-phase time-varying Wiener shaping: in each Hann window the section-mean amplitude spectrum A(f)
    (live traces, smoothed over nsm bins) is replaced by the output wavelet spectrum |W_out(f)|, filter
    |W_out| A / (A^2 + eps max A^2); filtered windows are summed and divided by the summed window."""
    n, nx = X.shape; L = int(round(win / DT)); step = int(round(hop / DT)); h = np.hanning(L)
    live = np.abs(X).sum(axis=0) > 0; Wo = np.abs(np.fft.rfft(w_out, nfft)); Wo = Wo / Wo.max()
    starts = list(range(-L // 2, n - L // 2, step)); out = np.zeros_like(X); wsum = np.zeros(n)
    for s0 in starts:
        i0, i1 = max(s0, 0), min(s0 + L, n); hw = h[i0 - s0:i1 - s0]
        seg = X[i0:i1] * hw[:, None]; F = np.fft.rfft(seg, nfft, axis=0)
        A = np.sqrt(lib.smooth(np.mean(np.abs(F[:, live]) ** 2, axis=1), nsm))
        G = Wo * A / (A ** 2 + eps * np.max(A ** 2))
        out[i0:i1] += np.fft.irfft(F * G[:, None], nfft, axis=0)[:i1 - i0]; wsum[i0:i1] += hw
    return out / np.maximum(wsum, 1e-3 * wsum.max())[:, None]


def score(o, ideal):
    s = sc11.m1(o, ideal); s["S"] = 0.5 * (s["M1"] + s["M1_lag"]); return s
