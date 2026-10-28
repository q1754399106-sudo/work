"""Shared code of experiment 11: SR + SC-a (Q factor from the spectral decay, compensation by SC-a) against
QScan + SC-a, the current recommended pipeline.

Imported unchanged from cache/inputs: experiment 7 lib.py (SC solver, QScan), experiment 8 semi.py (generator,
inverse Q), experiment 10 sq.py / gen.py / run.py (SR, generator, M1 scores), experiment 9 wq.py / synth.py (Gardner
reflectivity of the long logs, semi-synthetic wells). Their paths are redirected here.
"""
import sys, json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
IN = ROOT / "cache" / "inputs"
sys.path.insert(0, str(ROOT.parent / "novelty_evidence_20261025" / "cache" / "deps"))   # scipy
for p in ("novelty_evidence_20261025/scripts", "field_application_20261028/scripts", "semisynthetic_volve_20261028/scripts",
          "specslope_q_volve_20261028/scripts", "multiwell_tie_volve_20261028/scripts"):
    sys.path.insert(0, str(IN / p))
import sq, gen, run, wq, synth  # noqa: E402
from sq import np, lib, semi, bm, TW, DT  # noqa: E402

sq.ROOT, sq.IN = ROOT, IN
semi.ROOT, semi.IN = ROOT, IN
wq.ROOT, wq.IN, wq.E9 = ROOT, IN, IN / "multiwell_tie_volve_20261028"
E10 = IN / "specslope_q_volve_20261028"
DEV_OUT = ROOT.parent / "specslope_q_volve_20261028" / "cache" / "out"


def haar_sigma(x):
    return run.haar_sigma(x)


def prepare_sc(D, cm):
    """Everything SC-a and QScan need from one section."""
    live = np.abs(D).sum(axis=0) > 0; K, _, _ = bm.predictor(D, live)
    sigma = np.array([haar_sigma(D[:, j]) for j in range(D.shape[1])]); wv = lib.wavelets(D, TW, cm)
    return dict(live=live, K=K, sigma=sigma, wv=wv)


def sca(D, cm, Q, c, p):
    """SC-a with Q = c x model Q."""
    return lib.solve(cm, Q * c, p["wv"], D, p["sigma"], p["K"], p["live"], *lib.SEL_SC)[0]


def qscan(D, cm, Q, p):
    fq, info = lib.qscan(cm, Q, p["wv"], D, p["sigma"], p["K"], p["live"]); return fq, info


def phase_fix(X, cm, Q, c_from, c_to, nfft=1024):
    """Time-variant phase-only operator that turns the dispersion correction of Q x c_from into that of Q x c_to
    (the phase term of semi.inverse_q with t* = t*(Q c_to) - t*(Q c_from); amplitude unchanged)."""
    a = dict(cm, Q=np.asarray(Q, float) * c_to); b = dict(cm, Q=np.asarray(Q, float) * c_from)
    dts = lib.absorption(a, TW) - lib.absorption(b, TW); n = len(TW); f = np.fft.rfftfreq(nfft, DT)
    M = np.empty((n, n)); lag = np.arange(n); lf = np.log(np.maximum(f, np.finfo(float).tiny) / float(cm["fref"]))
    for i in range(n):
        H = np.exp(-2j * f * dts[i] * lf); H[0] = 1.0; H[-1] = H[-1].real
        M[i] = np.fft.irfft(H, nfft)[(i - lag) % nfft]
    return M @ X


def m1(out, ideal):
    return run.m1(out, ideal)
