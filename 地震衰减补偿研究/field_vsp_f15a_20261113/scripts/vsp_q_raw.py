"""Check of vsp_q.py on the raw stack (VSPNI_COMPUTED_1.SEGY, Z, 200 levels): first arrival taken at its own peak,
two window lengths; t* slope against the peak time in three bands; spectra (dB re 30 Hz) of shallow / middle / deep groups.

usage: vsp_q_raw.py RAW_DIR
"""
import sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
from vsp_q import read, F, NF, DT

if __name__ == "__main__":
    D, md, tvd = read(Path(sys.argv[1]) / "VSPNI_COMPUTED_1.SEGY"); o = np.argsort(md); D, md, tvd = D[o], md[o], tvd[o]; pk = np.argmax(np.abs(D), axis=1); t = pk * DT
    print("levels %d; peak time %.3f-%.3f s; non-monotonic steps %d" % (len(md), t.min(), t.max(), int((np.diff(t) < -0.005).sum())))
    for pre, post in ((0.02, 0.06), (0.03, 0.10)):
        n0, n1 = int(pre / DT), int(post / DT); w = np.hanning(n0 + n1); ok = (pk > n0) & (tvd > 1000)
        A = np.array([np.convolve(np.abs(np.fft.rfft(D[i, pk[i] - n0:pk[i] + n1] * w, NF)), np.ones(3) / 3, mode="same") for i in np.where(ok)[0]]); tt, zz = t[ok], tvd[ok]
        ref = np.median(A[:10] / A[:10].max(axis=1, keepdims=True), axis=0)
        for lo, hi in ((10, 60), (10, 80), (15, 100)):
            b = (F >= lo) & (F <= hi); ts = np.array([-np.polyfit(F[b], np.log(a[b] / ref[b]), 1)[0] / np.pi for a in A]); p = np.polyfit(tt, ts, 1)
            print("window -%.0f/+%.0f ms, band %d-%d Hz: 1000/Q %+6.2f, residual sd %.1f ms (t* range expected for Q 100: %.1f ms)" % (1e3 * pre, 1e3 * post, lo, hi, 1e3 * p[0], 1e3 * (ts - np.polyval(p, tt)).std(), 10 * (tt.max() - tt.min())))
    g = np.array_split(np.arange(len(A)), 3); fs = [10, 20, 30, 40, 50, 60, 80, 100]; print("spectra, dB re 30 Hz, at", fs, "Hz")
    for name, k in zip(("shallow", "middle", "deep"), g):
        s = np.median(A[k] / A[k][:, [np.argmin(abs(F - 30))]], axis=0); print("  %-7s TVD %4.0f-%4.0f m: %s" % (name, zz[k].min(), zz[k].max(), " ".join("%+5.1f" % (20 * np.log10(s[np.argmin(abs(F - f))])) for f in fs)))
