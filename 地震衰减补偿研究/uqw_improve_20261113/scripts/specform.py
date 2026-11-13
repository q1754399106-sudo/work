"""Step 3 diagnostic (no wells, no compensation): is the spectral change from the shallow window (1.2-1.6 s) to the target
window (2.3-2.8 s) of constant-Q form, ln[A_deep / A_shallow] = a - pi f dt*, on the field lines and on the development sets?

Prints the mean log amplitude ratio by frequency, its straight-line fit over 10-55 Hz (dt*), the rms residual and the
quadratic coefficient, per source (field: every 7th analysis line; semi: the 24 development sets; also their noise-free signal).

usage: specform.py
"""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]; BASE = ROOT.parent
sys.path.insert(0, str(BASE / "field_qperturb_20261113" / "scripts"))
import run as R  # noqa: E402
from run import np, L, lib, semi, TW, DT  # noqa: E402
sys.path.insert(0, str(BASE / "qunknown_comp_20261109" / "scripts"))
from stage0 import DEV, load  # noqa: E402

NF = 512; F = np.fft.rfftfreq(NF, DT); FS = np.arange(8, 73, 4); FIT = (F >= 10) & (F <= 55)


def amp(D, a, b):
    k = (TW >= a) & (TW <= b); live = np.abs(D).sum(0) > 0
    return np.sqrt(lib.smooth(np.mean(np.abs(np.fft.rfft(D[k][:, live] * np.hanning(k.sum())[:, None], NF, axis=0)) ** 2, axis=1), 5))


def report(name, sections, cm):
    r = np.array([np.log(amp(D, 2.3, 2.8) / amp(D, 1.2, 1.6)) for D in sections]); m = r.mean(0); p1 = np.polyfit(F[FIT], m[FIT], 1); p2 = np.polyfit(F[FIT], m[FIT], 2)
    res = m[FIT] - np.polyval(p1, F[FIT]); slopes = np.array([np.polyfit(F[FIT], q[FIT], 1)[0] for q in r]); x = lib.absorption(cm, [1.4, 2.55])
    print("%s (n %d): dt* fit %.4f s (per-section sd %.4f; curtain model %.4f), rms residual %.3f, quadratic %+.2e /Hz^2 (curvature over 10-55 Hz: %+.2f)" % (
        name, len(r), -p1[0] / np.pi, slopes.std() / np.pi, x[1] - x[0], np.sqrt(np.mean(res ** 2)), p2[0], p2[0] * (45 / 2) ** 2))
    print("   ln ratio minus fit at " + " ".join("%d" % f for f in FS) + " Hz:"); print("   " + " ".join("%+.2f" % (m[int(np.argmin(np.abs(F - f)))] - np.polyval(p1, f)) for f in FS))
    print("   ln ratio (re 10 Hz):   " + " ".join("%+.2f" % (m[int(np.argmin(np.abs(F - f)))] - m[int(np.argmin(np.abs(F - 10)))]) for f in FS), flush=True)


if __name__ == "__main__":
    cm = semi.curtain_cm(); lines = R.lines()[::7]
    report("field", [L.read_inline(ln) for ln in lines], cm)
    A = [load(d)[0] for d in DEV]
    report("semi D", [a["D"] for a in A], cm); report("semi signal", [a["signal"] for a in A], cm); report("semi ideal", [a["ideal"] for a in A], cm)
