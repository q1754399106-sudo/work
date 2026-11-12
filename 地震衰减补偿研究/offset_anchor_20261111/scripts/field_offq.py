"""Offset spectral ratio on Volve ST10010 CIP gathers (exploratory first look, NOTE_ZH.md).

Data (raw/, not in the repository): moveout-corrected PSDM CIP gathers in time, one inline per file, offsets -4950..4950
step 100 m, 4 ms, 0-3.4 s; phase-only inverse Q (120) was applied in processing, amplitudes are not Q-compensated.
RMS stacking velocities of the 4D processing on the same grid (436 inlines x 800 crosslines).

Per block of crosslines: the moveout is restored with the RMS velocity (t_x = sqrt(t0^2 + x^2 / v^2), 8 x Fourier
oversampling then linear), spectra are taken in windows that follow the moveout curve, and for every offset class
  y = -slope over frequency of 0.5 ln(P_x / P_ref),  x = pi (t_x - t_ref)  at the window centre
so that y = x / Q_eff for a constant-Q overburden (extra t* = extra traveltime / Q). Q_eff = sum(x^2) / sum(x y).
"flat": the same ratio on the delivered (stretched) gathers, for comparison.

usage: field_offq.py GATHER_SEGY VELOCITY_SEGY [XL0 XL1 BLOCK]
"""
import sys, os, json
from pathlib import Path
import numpy as np

NS, DT, UP, NF = 851, 0.004, 8, 512
T = np.arange(NS) * DT; TU = np.arange(NS * UP) * DT / UP; F = np.fft.rfftfreq(NF, DT)
WIN = [(a, a + 0.5) for a in np.arange(0.9, 2.71, 0.3)]; BAND = (8.0, 60.0); REF = (50, 350); OFFS = np.arange(450, 3251, 200); MUTE = 1.25; FLOOR = 0.02


def ibm(b):
    u = b.view(">u4").astype(np.int64)
    return np.where(u >> 31, -1.0, 1.0) * (u & 0xffffff) / 16777216.0 * 16.0 ** (((u >> 24) & 0x7f) - 64)


def read(path, ibm_format=True):
    tl = 240 + NS * 4; n = (os.path.getsize(path) - 3600) // tl; m = np.memmap(path, dtype=np.uint8, mode="r", offset=3600, shape=(n, tl))
    return m, lambda o: np.ascontiguousarray(m[:, o - 1:o + 3]).view(">i4").ravel()


def traces(m, k):
    return ibm(np.ascontiguousarray(m[k, 240:])).reshape(len(k), NS)


def resample(g, tq):
    """g (n_trace, n_t) read at times tq (n_t), zero outside the record."""
    G = np.fft.rfft(g, axis=1); Z = np.zeros((g.shape[0], NS * UP // 2 + 1), complex); Z[:, :G.shape[1]] = G; u = np.fft.irfft(Z, NS * UP, axis=1) * UP
    return np.stack([np.interp(tq, TU, r, left=0.0, right=0.0) for r in u])


def power(y, a, b):
    k = (T >= a) & (T < b); seg = y[:, k]; seg = seg[(np.abs(seg) > 0).mean(axis=1) > 0.95]
    if len(seg) < 20:
        return None
    return np.convolve(np.mean(np.abs(np.fft.rfft(seg * np.hanning(k.sum()), NF, axis=1)) ** 2, axis=0), np.ones(5) / 5, mode="same")


def band(P0, Px):
    ok = (F >= BAND[0]) & (F <= BAND[1]) & (P0 > FLOOR * P0.max()) & (Px > FLOOR * Px.max()); i = np.where(F >= BAND[0])[0][0]; j = i
    while j < len(F) and ok[j]:
        j += 1
    return slice(i, j)


def block(m, off, xl, v, x0, x1):
    """rows: one per (window, offset class) with x, y (moveout restored) and y_flat."""
    sel = (xl >= x0) & (xl < x1); ao = np.abs(off); rows = []
    k = np.where(sel & (ao >= REF[0]) & (ao <= REF[1]))[0]; xr = float(np.sqrt(np.mean(ao[k].astype(float) ** 2))); ref_flat = traces(m, k)
    ref = resample(ref_flat, np.interp(T, np.sqrt(T ** 2 + (xr / v) ** 2), T, left=0.0))
    for o in OFFS:
        k = np.where(sel & (ao >= o - 100) & (ao < o + 100))[0]
        if len(k) < 20:
            continue
        flat = traces(m, k); xo = float(np.sqrt(np.mean(ao[k].astype(float) ** 2))); tx = np.sqrt(T ** 2 + (xo / v) ** 2); rec = resample(flat, np.interp(T, tx, T, left=0.0))
        for a, b in WIN:
            tc = 0.5 * (a + b); mv = lambda x, t: float(np.sqrt(t ** 2 + (x / np.interp(t, T, v)) ** 2))
            if mv(xo, tc) / tc > MUTE:
                continue
            P0, Px = power(ref, mv(xr, a), mv(xr, b)), power(rec, mv(xo, a), mv(xo, b)); Q0, Qx = power(ref_flat, a, b), power(flat, a, b)
            if P0 is None or Px is None or Q0 is None or Qx is None:
                continue
            s = band(P0, Px); sf = band(Q0, Qx)
            if s.stop - s.start < 12 or sf.stop - sf.start < 12:
                continue
            rows.append(dict(win=a, off=xo, x=np.pi * (mv(xo, tc) - mv(xr, tc)), fmax=float(F[s.stop - 1]), y=float(-np.polyfit(F[s], 0.5 * np.log(Px[s] / P0[s]), 1)[0]),
                             y_flat=float(-np.polyfit(F[sf], 0.5 * np.log(Qx[sf] / Q0[sf]), 1)[0])))
    return rows


if __name__ == "__main__":
    gpath, vpath = sys.argv[1], sys.argv[2]; X0, X1, BL = (int(q) for q in sys.argv[3:6]) if len(sys.argv) > 5 else (2160, 2480, 40)
    m, h = read(gpath); off, xl, il = h(37), h(193), int(h(189)[0]); mv_, hv = read(vpath); vil, vxl = hv(189), hv(193); out = []
    for x0 in range(X0, X1, BL):
        kv = np.where((vil == il) & (vxl >= x0) & (vxl < x0 + BL))[0]; v = ibm(np.ascontiguousarray(mv_[kv, 240:])).reshape(len(kv), NS).mean(axis=0)
        rows = block(m, off, xl, v, x0, x0 + BL); out += [dict(r, inline=il, xl0=x0) for r in rows]
        print("inline %d crosslines %d-%d  vrms(1.0, 2.0, 3.0 s) %.0f %.0f %.0f" % (il, x0, x0 + BL - 1, *np.interp([1.0, 2.0, 3.0], T, v)))
        for a, _ in WIN:
            g = [r for r in rows if r["win"] == a]
            if len(g) < 3:
                continue
            x, y, yf = (np.array([r[q] for r in g]) for q in ("x", "y", "y_flat")); q = np.sum(x * x) / np.sum(x * y); cc = np.corrcoef(x, y)[0, 1]
            print("  %.1f-%.1f s  n %2d  fmax %2.0f  1000/Q %+7.2f (Q %7.1f)  corr(x, y) %+.2f | flat: 1000/Q %+7.2f | y*1000 by offset: %s" % (
                a, a + 0.5, len(g), np.median([r["fmax"] for r in g]), 1000 / q, q, cc, 1000 * np.sum(x * yf) / np.sum(x * x), " ".join("%+.1f" % (1000 * w) for w in y)), flush=True)
    res = Path(__file__).resolve().parents[1] / "results"; res.mkdir(exist_ok=True); (res / f"field_offq_{il}.json").write_text(json.dumps(out), encoding="utf-8")
