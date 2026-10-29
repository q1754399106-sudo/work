"""Shared code of experiment 13: the experiment-12 methods, with their frozen parameters, on real Volve inlines.

Imported unchanged from cache/inputs: experiment 12 lib12.py (iq, tvw) and everything it imports (sc11 of
experiment 11, lib/semi/sq of experiments 7-10), experiment 9 tie.py/wq.py (well logs, well paths, SR factor per
inline, tie score). Their paths are redirected here. New code: reading a whole inline at its own number (no
replacement), the QScan factor per inline, the method sets.
"""
import sys, json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
IN = ROOT / "cache" / "inputs"
sys.path.insert(0, str(ROOT.parent / "novelty_evidence_20261025" / "cache" / "deps"))   # scipy
for p in ("novelty_evidence_20261025/scripts", "field_application_20261028/scripts", "semisynthetic_volve_20261028/scripts",
          "specslope_q_volve_20261028/scripts", "multiwell_tie_volve_20261028/scripts", "srsc_volve_20261028/scripts",
          "baselines_volve_20261029/scripts"):
    sys.path.insert(0, str(IN / p))
import lib12, sc11, tie, wq  # noqa: E402
from sc11 import np, lib, semi, sq, TW, DT  # noqa: E402

for mod in (lib12, sc11, sq, semi, wq):
    mod.ROOT, mod.IN = ROOT, IN
wq.E9 = IN / "multiwell_tie_volve_20261028"
sc11.E10 = lib12.E10 = IN / "specslope_q_volve_20261028"

SEL = json.loads((IN / "baselines_volve_20261029" / "analysis" / "dev_selection.json").read_text(encoding="utf-8"))
A_TUNED = SEL["B2"]["a"]                       # 2^-13
XL0, XL1, XSTEP = 2100, 2500, 2                # read_lines window: 201 xlines
QSCAN_LINES = [10043, 10087, 10113, 10129, 10141, 10165, 10183, 10193, 10207, 10209]   # well mid inlines + 10141
WIN_T = semi.win(np.arange(1126) * DT, semi.T0, semi.T1)
_CUBE = None


def cube():
    global _CUBE
    if _CUBE is None:
        _CUBE = tie.cube()
    return _CUBE


def read_inline(L):
    """Inline L, xlines XL0..XL1 every XSTEP, 1.2-3.0 s; missing traces are zero (dead)."""
    h, x, y, il, xl = cube(); D = np.zeros((int(WIN_T.sum()), (XL1 - XL0) // XSTEP + 1))
    for j in np.where((il == L) & (xl >= XL0) & (xl <= XL1) & ((xl - XL0) % XSTEP == 0))[0]:
        D[:, (xl[j] - XL0) // XSTEP] = semi.ibm(h[j, 240:].copy().view(">u4"))[WIN_T]
    return D


def c_sr(L):
    """Experiment-10 SR factor interpolated in inline number (as experiment 9)."""
    return tie.sr_factor(L)


def c_qscan(L, qs):
    """QScan factor of inline L: log-linear interpolation between the QScan lines (held constant beyond them)."""
    k = sorted((int(a), v) for a, v in qs.items())
    return float(np.exp(np.interp(L, [a for a, _ in k], [np.log(v) for _, v in k])))


def methods_a(D, cm, Q, c, p):
    """SR-based methods (PLAN_ZH.md section 3)."""
    out = dict(Input=D, IQ30=lib12.iq(D, cm, Q, c, c, 30), B1_IQ70=lib12.iq(D, cm, Q, c, c, SEL["B1"]["g_db"]))
    out["B2_SCa"] = lib.solve(cm, Q * c, p["wv"], D, p["sigma"], p["K"], p["live"], A_TUNED, lib.SEL_SC[1])[0]
    out["B3_VDB"] = lib12.tvw(lib12.iq(D, cm, Q, None, c, 0), p["wv"]["w_out"], SEL["B3"]["eps"], win=SEL["B3"]["win"])
    out["PROP"] = sc11.phase_fix(sc11.sca(D, cm, Q, c * 1.3, p), cm, Q, c * 1.3, c)
    return out


def methods_b(D, cm, Q, fq, p):
    """The user's pipeline with the QScan factor fq: as is (a = 0.125) and with the tuned weight."""
    return dict(QSCa=sc11.sca(D, cm, Q, fq, p),
                QSCa_tuned=lib.solve(cm, Q * fq, p["wv"], D, p["sigma"], p["K"], p["live"], A_TUNED, lib.SEL_SC[1])[0])


def bands(o):
    """Mean coherent SNR (2.0-3.0 s) in 15 Hz bands, and the 5-95% energy band."""
    f, snr, _ = lib.coherence_snr(o, TW, 2.0, 3.0)
    b = {f"{a}-{a + 15}": float(snr[(f >= a) & (f < a + 15)].mean()) for a in (10, 25, 40, 55, 70)}
    b["band05_95"] = [float(v) for v in lib.spectral(o, o, TW)[1:]]
    return b
