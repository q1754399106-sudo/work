"""VSP Q, second pass: raw stack (VSPNI_COMPUTED_1.SEGY), first arrival = largest peak within +-40 ms of the listed
first break, window -20/+60 ms; t* per level against a shallow reference (10 shallowest levels below 1000 m TVD);
slope against one-way time with a bootstrap over levels (90% interval), per band, whole column and two halves.

usage: vsp_q2.py RAW_DIR
"""
import sys, json
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
from vsp_q import read, stations, F, NF, DT

if __name__ == "__main__":
    raw = Path(sys.argv[1]); D, md, tvd = read(raw / "VSPNI_COMPUTED_1.SEGY"); st = stations(raw / "VSPNI_COMPUTED_STKLIST.ASC")
    t1 = np.array([st.get(round(float(m), 1), np.nan) for m in md]); ok = ~np.isnan(t1) & (tvd > 1000); D, md, tvd, t1 = D[ok], md[ok], tvd[ok], t1[ok]; o = np.argsort(md); D, md, tvd, t1 = D[o], md[o], tvd[o], t1[o]
    A, tp = [], []
    for d, t in zip(D, t1):
        i = int(round(t / DT)); p = i - 40 + int(np.argmax(np.abs(d[i - 40:i + 41]))); tp.append(p * DT)
        A.append(np.convolve(np.abs(np.fft.rfft(d[p - 20:p + 60] * np.hanning(80), NF)), np.ones(3) / 3, mode="same"))
    A = np.array(A); tp = np.array(tp); ref = np.median(A[:10] / A[:10].max(axis=1, keepdims=True), axis=0); rng = np.random.default_rng(1); out = dict(tvd=tvd.tolist(), t=t1.tolist())
    print("levels %d, TVD %.0f-%.0f m, one-way %.3f-%.3f s; peak minus listed first break: median %+.0f ms, range %+.0f..%+.0f" % (len(tvd), tvd[0], tvd[-1], t1[0], t1[-1], *(1e3 * np.array([np.median(tp - t1), (tp - t1).min(), (tp - t1).max()]))))
    for lo, hi in ((10, 60), (10, 80), (15, 100), (20, 100), (30, 100)):
        b = (F >= lo) & (F <= hi); ts = np.array([-np.polyfit(F[b], np.log(a[b] / ref[b]), 1)[0] / np.pi for a in A]); out["%d-%d" % (lo, hi)] = ts.tolist(); row = []
        for name, s in (("all", np.ones(len(ts), bool)), ("above 2600 m", tvd < 2600), ("below 2600 m", tvd >= 2600)):
            k = np.where(s)[0]; q = np.polyfit(t1[k], ts[k], 1)[0]; bs = [np.polyfit(t1[j], ts[j], 1)[0] for j in (rng.choice(k, len(k)) for _ in range(500))]
            row.append("%s %+6.1f [%+6.1f, %+6.1f]" % (name, 1e3 * q, *(1e3 * np.percentile(bs, [5, 95]))))
        print("band %3d-%3d Hz, 1000/Q: %s" % (lo, hi, " | ".join(row)))
    (Path(__file__).resolve().parents[1] / "results" / "vsp_q2.json").write_text(json.dumps(out), encoding="utf-8")
