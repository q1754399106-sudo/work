"""Q from the downgoing wave of the 15/9-F-15A vertical-incidence VSP (spectral ratio; exploratory first look).

Input (raw/, not in the repository): VSPNI_COMPUTED_4.SEGY, "downwave input to deconvolution", Z, 191 levels, 1 ms,
aligned at the first break; receiver MD / TVD (RKB) in bytes 37 / 41 (scaled 1e-4); first-break times from
VSPNI_COMPUTED_STKLIST.ASC (station lines, microseconds, one-way from the source above the receiver).
Per level: amplitude spectrum of the first arrival in a tapered window; slope of ln(A / A_ref) over BAND = -pi dt*.
t*(level) relative to the reference group (median spectrum of the REF shallowest usable levels).

usage: vsp_q.py RAW_DIR
"""
import sys, re, json
from pathlib import Path
import numpy as np

NS, DT, NF = 5000, 0.001, 1024; F = np.fft.rfftfreq(NF, DT); BANDS = [(10.0, 60.0), (10.0, 80.0), (15.0, 100.0)]; PRE, POST = 0.03, 0.10


def ibm(b):
    u = b.view(">u4").astype(np.int64)
    return np.where(u >> 31, -1.0, 1.0) * (u & 0xffffff) / 16777216.0 * 16.0 ** (((u >> 24) & 0x7f) - 64)


def read(path):
    b = np.fromfile(path, np.uint8)[3600:].reshape(-1, 240 + NS * 4); H = np.ascontiguousarray(b[:, :240]); i4 = lambda o: H[:, o - 1:o + 3].copy().view(">i4").ravel()
    return ibm(np.ascontiguousarray(b[:, 240:])).reshape(len(b), NS), i4(37) * 1e-4, i4(41) * 1e-4


def stations(path):
    """MD -> first break (s) from the station lines of the stack list."""
    out = {}
    for line in Path(path).read_text(errors="replace").splitlines():
        m = re.match(r"^\s*(\d+)\s+(\d+)\s+(\d+)\s+(\d+\.\d)\s+(\d+\.\d)\s+(-?\d+)\s+(-?\d+)\s+(-?\d+)\s+(-?\d+)\s*$", line)
        if m:
            out[round(float(m.group(4)), 1)] = int(m.group(3)) * 1e-6
    return out


if __name__ == "__main__":
    raw = Path(sys.argv[1]); D, md, tvd = read(raw / "VSPNI_COMPUTED_4.SEGY"); st = stations(raw / "VSPNI_COMPUTED_STKLIST.ASC")
    t1 = np.array([st.get(round(float(m), 1), np.nan) for m in md]); o = np.argsort(md); D, md, tvd, t1 = D[o], md[o], tvd[o], t1[o]
    print("levels %d, MD %.0f-%.0f m, TVD %.0f-%.0f m, first break %.3f-%.3f s (missing %d)" % (len(md), md[0], md[-1], tvd[0], tvd[-1], np.nanmin(t1), np.nanmax(t1), np.isnan(t1).sum()))
    pk = np.argmax(np.abs(D), axis=1); i0 = int(np.median(pk)); k = slice(i0 - int(PRE / DT), i0 + int(POST / DT)); w = np.hanning(k.stop - k.start)
    A = np.abs(np.fft.rfft(D[:, k] * w, NF, axis=1)); A = np.array([np.convolve(a, np.ones(3) / 3, mode="same") for a in A])
    cen = (A * F).sum(1) / A.sum(1); print("alignment sample %d; centroid frequency by depth (every 15th level): %s" % (i0, " ".join("%.0f" % c for c in cen[::15])))
    ref = np.median(A[:8] / A[:8].max(axis=1, keepdims=True), axis=0); res = {}
    for lo, hi in BANDS:
        b = (F >= lo) & (F <= hi); ts = np.array([-np.polyfit(F[b], np.log(a[b] / ref[b]), 1)[0] / np.pi for a in A]); res["%g-%g" % (lo, hi)] = ts.tolist()
        ok = ~np.isnan(t1); p = np.polyfit(t1[ok], ts[ok], 1); print("band %g-%g Hz: overall slope dt*/dt = %.5f -> Q %.0f" % (lo, hi, p[0], 1 / p[0]))
        for a, c in [(0.1, 0.6), (0.6, 0.9), (0.9, 1.2), (1.2, 1.5)]:
            s = ok & (t1 >= a) & (t1 < c)
            if s.sum() >= 6:
                q = np.polyfit(t1[s], ts[s], 1)[0]; r = ts[s] - np.polyval(np.polyfit(t1[s], ts[s], 1), t1[s])
                print("   one-way %.1f-%.1f s (TVD %4.0f-%4.0f m, n %3d): 1000/Q %+6.2f (Q %6.0f), residual sd %.4f s" % (a, c, tvd[s].min(), tvd[s].max(), s.sum(), 1000 * q, 1 / q if q else np.inf, r.std()))
    b = (F >= 10) & (F <= 80); print("t* (ms, 10-80 Hz) at every 10th level:"); print("  " + "  ".join("%.0fm/%.2fs:%+.1f" % (tvd[i], t1[i], 1000 * res["10-80"][i]) for i in range(0, len(md), 10)))
    out = Path(__file__).resolve().parents[1] / "results"; out.mkdir(exist_ok=True)
    (out / "vsp_q.json").write_text(json.dumps(dict(md=md.tolist(), tvd=tvd.tolist(), t_oneway=t1.tolist(), centroid=cen.tolist(), tstar=res)), encoding="utf-8")
