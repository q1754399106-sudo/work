"""Reflectivity of the Volve wells in two-way time (1 ms) for the colour stress group and the colour statistics.

For every well of experiment 9's inventory with AC and DEN: two-way time along hole from the sonic, stretched
piecewise-linearly so that each pick (Well_picks_Volve_v1.dat, MD/TWT) falls at its TWT; only the MD range
covered by picks is kept; gaps (missing AC or DEN) end a segment. Impedance -> RC = diff(ln Z)/2 at 1 ms.
"""
import sys, re
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import sq  # noqa: E402
from sq import np  # noqa: E402

E9 = sq.IN / "multiwell_tie_volve_20261028"
WELLS = ["F-1", "F-1_A", "F-1_B", "F-10", "F-11_A", "F-11_T2", "F-12", "F-14", "F-15", "F-15_A", "F-15_B", "F-4", "F-5"]
MIN_SEG = 0.1          # s, shortest kept segment


def picks():
    out, spans = {}, None
    for line in (E9 / "raw" / "picks" / "Well_picks_Volve_v1.dat").read_text(encoding="latin-1").splitlines():
        if line.strip().startswith("-----"):
            spans = [(m.start(), m.end()) for m in re.finditer(r"-+", line)]; continue
        if spans is None or not line.strip() or line.startswith("Well NO") or line.lstrip().startswith("Well name"):
            continue
        cols = [line[a:b + 1].strip() if i < len(spans) - 1 else line[a:].strip() for i, (a, b) in enumerate(spans)]
        try:
            md, twt = float(cols[4]), float(cols[7])
        except (ValueError, IndexError):
            continue
        out.setdefault(cols[0].replace("NO ", ""), []).append((md, twt / 1000))
    return out


def well_rc(name, P):
    """List of 1 ms RC segments (numpy arrays) of one well."""
    A = np.loadtxt(E9 / "inventory" / "logs" / f"15_9-{name}.csv", delimiter=",", skiprows=1)
    md, ac, den = A[:, 0], A[:, 1], A[:, 2]
    pk = sorted(P.get("15/9-" + name.replace("_", " "), []))
    if len(pk) < 2:
        return []
    pm, pt = np.array([p[0] for p in pk]), np.array([p[1] for p in pk])
    ok = np.isfinite(ac) & np.isfinite(den) & (md >= pm.min()) & (md <= pm.max())
    segs, cur = [], []
    idx = np.where(ok)[0]
    for k in idx:                                       # split on gaps
        if cur and k != cur[-1] + 1:
            segs.append(cur); cur = []
        cur.append(k)
    if cur:
        segs.append(cur)
    out = []
    for s in segs:
        s = np.array(s); m, v, d = md[s], 0.3048 / (ac[s] * 1e-6), den[s]
        ts = np.concatenate([[0], np.cumsum(2 * np.diff(m) / (0.5 * (v[1:] + v[:-1])))])
        tp = np.interp(m, pm, pt)                        # pick-based TWT (coarse)
        inside = (pm > m[0]) & (pm < m[-1]); km = np.r_[m[0], pm[inside], m[-1]]
        kt_s = np.interp(km, m, ts); kt_p = np.interp(km, m, tp)
        t = np.interp(ts, kt_s, kt_p)                    # sonic detail, pick anchors
        if t[-1] - t[0] < MIN_SEG:
            continue
        tt = np.arange(t[0], t[-1], 0.001); z = np.log(np.interp(tt, t, v * d))
        out.append(0.5 * np.diff(z))
    return out


def all_segments():
    P = picks(); return {w: well_rc(w, P) for w in WELLS}


def colour(rc, fs=1000.0, lo=10, hi=60):
    """Slope (per Hz) of 0.5 ln power of a Hann-tapered RC window over lo-hi Hz."""
    X = np.fft.rfft(rc * np.hanning(len(rc)), 4096); f = np.fft.rfftfreq(4096, 1 / fs)
    p = np.convolve(np.abs(X) ** 2, np.ones(41) / 41, mode="same"); s = (f >= lo) & (f <= hi)
    return float(np.polyfit(f[s], 0.5 * np.log(p[s]), 1)[0])


if __name__ == "__main__":
    import json
    S = all_segments(); rows = {}
    for w, segs in S.items():
        L = [len(s) / 1000 for s in segs]
        ks = [colour(s[i:i + 200]) for s in segs for i in range(0, len(s) - 199, 100)]
        rows[w] = dict(n_segments=len(segs), seconds=L, colour_0p2s=ks)
        print(w, len(segs), np.round(L, 3), "kappa", np.round(ks, 4), flush=True)
    allk = np.array([k for r in rows.values() for k in r["colour_0p2s"]])
    t_log, rc = sq.semi.log_reflectivity()
    ref = [colour(rc[i:i + 200]) for i in range(0, len(rc) - 199, 100)]
    summ = dict(n_windows=int(len(allk)), mean=float(allk.mean()), sd=float(allk.std(ddof=1)),
                exp8_log_F11T2_posthoc=ref)
    print(summ)
    (sq.ROOT / "analysis" / "wells_colour.json").write_text(json.dumps(dict(wells=rows, summary=summ), indent=1), encoding="utf-8")
