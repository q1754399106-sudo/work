"""PLAN_FIELD_I_ZH.md: offset spectral ratio on the Volve seq-79 / cable-5081 PZ shot gathers (ST0202R08, not NMO-corrected,
top-muted) around the Volve F platform, against the F-15A zero-offset VSP interval Q of TVD 1005-2600 m (two-way 1.02-2.40 s).

Geometry along the line (azimuth 284 deg, platform at 0): CMP = (s + g) / 2, CMPs within +-1000 m, 8 blocks of 250 m.
Moveout and regressors from 1-D ray tracing in the VTI-2008 vertical Vp (mean over the CMP range): source at 6 m, receiver
at the seabed (median water depth), same ray parameter on both legs. Windows: recorded zero-offset time centres 0.6-2.4 s
(step 0.2 s), 0.4 s long; far windows have the same length, centred at the ray-traced time of the window centre (frozen
A1b rule EQLEN, no spectral smoothing). Reference class |x| <= 350 m, offset classes 500-3900 m (+-100 m). A far window
is used only if it starts >= 0.1 s after the mute line t = -0.13 + |x| / 1216 and the moveout stretch is <= 1.3.
Band 8-60 Hz, both powers > 2 % of their in-band peak, contiguous from 8 Hz, >= 10 Hz wide (no noise estimate exists).
y = -d/df 0.5 ln(P_x / P_ref) = pi (dtau_A / Q_A + dtau_B / Q_B): dtau_L = extra time spent in layer L by the far ray
(A: seabed to two-way 1.02 s, B: 1.02-2.40 s; the water legs do not attenuate). Least squares without intercept.
Primary: 1000 / Q_B from all blocks pooled, 90 % interval from a bootstrap over the 8 blocks (2000 resamples, seed 79).

usage: field79.py synth   self-check on synthetic gathers with the field geometry and mute (Q_A 60, Q_B 90)
       field79.py run [vscale]   field data (vscale 1 = primary; 0.95 / 1.05 descriptive)  -> results/field79[_v..].json"""
import sys, json, math, time
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
MAIN = Path("D:/ChatGPT/研究/地震衰减补偿研究/uq_rtm_20261117/cache")
SEGY = MAIN / "volve_seq079" / "seq079_dline5081.segy"; VP = MAIN / "fld" / "vp.bin"; META = MAIN / "fld" / "meta.json"
DT = 0.004; NT = 800; T = np.arange(NT) * DT; NF = 512; F = np.fft.rfftfreq(NF, DT)
CENT = np.round(np.arange(0.6, 2.21, 0.2), 3); WL = 0.4; REF = 350.0; OFFS = np.arange(500, 3901, 200); HW = 100.0   # amendment: centres <= 2.2 s (2.4 s windows reached below layer B)
CMP_LIM = 1000.0; NBLK = 8; BAND = (8.0, 60.0); FLOOR = 0.02; MINBW = 10.0; STRETCH = 1.3
MUTE = lambda x: -0.13 + np.abs(x) / 1216.0; MUTE_PAD = 0.1; TB = (1.02, 2.40); ZS = 6.0


# ---------------- data ----------------
def read_headers():
    raw = np.memmap(SEGY, dtype=np.uint8, mode="r"); ns = int(np.frombuffer(raw[3220:3222].tobytes(), ">i2")[0])
    tsz = 240 + 4 * ns; n = (raw.size - 3600) // tsz; H = np.ndarray((n, tsz), np.uint8, raw[3600:3600 + n * tsz])
    h = lambda o, t=">i4": np.frombuffer(np.ascontiguousarray(H[:, o:o + (4 if t == ">i4" else 2)]).tobytes(), t).astype(np.float64)
    sc = h(70, ">i2"); f = np.where(sc < 0, -1.0 / sc, np.where(sc > 0, sc, 1.0)); A = math.radians(14); X0, Y0 = 435050.02, 6478563.53
    il = lambda x, y: -(x - X0) * math.cos(A) + (y - Y0) * math.sin(A)
    s = il(h(72) * f, h(76) * f); g = il(h(80) * f, h(84) * f)
    return H, dict(cmp=(s + g) / 2, off=np.abs(h(36)), wd=h(64))


def ibm2f(b):
    b = b.astype(np.uint32); s = np.where(b >> 31, -1.0, 1.0); e = ((b >> 24) & 0x7F).astype(np.int32) - 64
    return (s * (b & 0x00FFFFFF) / float(1 << 24) * np.power(16.0, e)).astype(np.float32)


def traces(H, k):
    return ibm2f(np.frombuffer(np.ascontiguousarray(H[k, 240:240 + 4 * NT]).tobytes(), ">u4").reshape(len(k), NT))


# ---------------- 1-D rays ----------------
class Rays:
    def __init__(self, vscale, zw):
        m = json.loads(META.read_text()); nz, nx = m["nz"], m["nx"]; V = np.fromfile(VP, np.float32).reshape(nz, nx)
        x = m["xmin"] + (np.arange(nx) - m["npml"]) * m["dx"]; z = (np.arange(nz) - m["npml"]) * m["dx"]
        k = z >= 0; v = V[np.ix_(k, (x >= -CMP_LIM) & (x <= CMP_LIM))].mean(axis=1).astype(float) * vscale
        self.dz = 1.0; self.z = np.arange(0, z[k][-1], self.dz) + 0.5 * self.dz; self.v = np.interp(self.z, z[k], v)
        self.T0 = 2 * np.cumsum(self.dz / self.v)                                 # vertical two-way time from the sea surface
        self.zw = zw; self.zA = self.z[np.searchsorted(self.T0, TB[0])]; self.zB = self.z[np.searchsorted(self.T0, TB[1])]
        self.sed = self.z > zw; self.lay = [(self.z > zw) & (self.z <= self.zA), (self.z > self.zA) & (self.z <= self.zB)]
        self.legs = ((self.z > ZS).astype(float) + (self.z > zw).astype(float))   # number of legs crossing each 1 m cell
        self.vmax = None

    def shoot(self, zr, p):
        """recorded time, offset and time in each layer (A, B, all sediment) for reflector depth zr, ray parameters p."""
        k = self.z < zr; v = self.v[k]; w = self.legs[k] * self.dz; q = np.sqrt(np.maximum(1 - (p[:, None] * v[None]) ** 2, 1e-9))
        dt = w / (v[None] * q); x = (w * v)[None] * p[:, None] / q
        lay = [L[k] for L in self.lay] + [self.sed[k]]
        return dt.sum(1), x.sum(1), [dt[:, L].sum(1) for L in lay]

    def table(self, tc):
        """for a recorded zero-offset time tc: functions of offset -> (recorded time, layer times A, B, sediment)."""
        t0 = lambda zr: self.shoot(zr, np.zeros(1))[0][0]
        lo, hi = self.zw + 1, self.z[-1] - 1
        for _ in range(60):
            mid = 0.5 * (lo + hi); lo, hi = (mid, hi) if t0(mid) < tc else (lo, mid)
        zr = 0.5 * (lo + hi); vm = self.v[self.z < zr].max(); p = np.linspace(0, 0.999 / vm, 4000)
        t, x, L = self.shoot(zr, p); ok = np.r_[True, np.diff(x) > 0]; t, x = t[ok], x[ok]; L = [l[ok] for l in L]
        return lambda xq: (np.interp(xq, x, t, right=np.nan), *[np.interp(xq, x, l, right=np.nan) for l in L])


# ---------------- estimator ----------------
def power(y, a, b):
    k = (T >= a) & (T < b)
    if k.sum() < int(0.9 * WL / DT):
        return None
    seg = y[:, k]; seg = seg[(np.abs(seg) > 0).mean(axis=1) > 0.95]
    if len(seg) < 10:
        return None
    return np.mean(np.abs(np.fft.rfft(seg * np.hanning(k.sum()), NF, axis=1)) ** 2, axis=0)


def band(P0, Px):
    ib = (F >= BAND[0]) & (F <= BAND[1])
    ok = ib & (P0 > FLOOR * P0[ib].max()) & (Px > FLOOR * Px[ib].max()); i = np.where(F >= BAND[0])[0][0]; j = i
    while j < len(F) and ok[j]:
        j += 1
    return slice(i, j) if F[j - 1] - F[i] >= MINBW else None


def block_rows(G, off, tabs):
    """G: traces of one CMP block (n, NT), off: |offset|. One row per (window, offset class)."""
    rows = []; kr = off <= REF
    if kr.sum() < 10:
        return rows
    xr = float(np.sqrt(np.mean(off[kr] ** 2)))
    for tc, tab in zip(CENT, tabs):
        r0 = tab(np.array([xr])); tr = float(r0[0][0]); P0 = power(G[kr], tr - WL / 2, tr + WL / 2)
        if P0 is None:
            continue
        for o in OFFS:
            k = (off >= o - HW) & (off < o + HW)
            if k.sum() < 10:
                continue
            xo = float(np.sqrt(np.mean(off[k] ** 2))); tx, la, lb, ls = (float(v[0]) for v in tab(np.array([xo])))
            if not np.isfinite(tx) or tx / tc > STRETCH or tx - WL / 2 < MUTE(o + HW) + MUTE_PAD or tx + WL / 2 > T[-1]:
                continue
            Px = power(G[k], tx - WL / 2, tx + WL / 2)
            if Px is None:
                continue
            s = band(P0, Px)
            if s is None:
                continue
            y = float(-np.polyfit(F[s], 0.5 * np.log(Px[s] / P0[s]), 1)[0])
            rows.append(dict(tc=float(tc), off=xo, y=y, fmax=float(F[s.stop - 1]), dA=la - float(r0[1][0]), dB=lb - float(r0[2][0]),
                             dS=ls - float(r0[3][0])))
    return rows


def fit(rows):
    if len(rows) < 3:
        return [np.nan, np.nan]
    X = np.pi * np.array([[r["dA"], r["dB"]] for r in rows]); y = np.array([r["y"] for r in rows])
    return list(1000 * np.linalg.lstsq(X, y, rcond=None)[0])


def fit1(rows):
    """amendment primary: one effective 1000/Q of the whole sediment path (regressor dS)."""
    if len(rows) < 3:
        return np.nan
    x = np.pi * np.array([r["dS"] for r in rows]); y = np.array([r["y"] for r in rows]); return float(1000 * (x @ y) / (x @ x))


def summarise(blocks):
    allr = [r for b in blocks for r in b]; pt = fit(allr); rng = np.random.default_rng(79); bs = []; b1 = []
    for _ in range(2000):
        pick = rng.integers(0, NBLK, NBLK); sub = [r for i in pick for r in blocks[i]]; bs.append(fit(sub)); b1.append(fit1(sub))
    bs = np.array(bs); lo, hi = np.nanpercentile(bs, [5, 95], axis=0)
    eff = {}
    for tc in CENT:
        g = [r for r in allr if r["tc"] == tc]
        if len(g) >= 3:
            x = np.pi * np.array([r["dS"] for r in g]); y = np.array([r["y"] for r in g])
            eff["%.1f" % tc] = dict(n=len(g), inv_q_eff_x1000=float(1000 * (x @ y) / (x @ x)), corr=float(np.corrcoef(x, y)[0, 1]) if len(g) > 2 else None,
                                    fmax_median=float(np.median([r["fmax"] for r in g])))
    l1, h1 = np.nanpercentile(b1, [5, 95])
    return dict(invQ_eff=fit1(allr), ci90_eff=[float(l1), float(h1)], per_block_invQ_eff=[fit1(b) for b in blocks], n_rows=len(allr), rows_per_block=[len(b) for b in blocks], rows_layerB_windows_per_block=[sum(r["tc"] > TB[0] for r in b) for b in blocks],
                invQ_A=pt[0], invQ_B=pt[1], ci90_A=[float(lo[0]), float(hi[0])], ci90_B=[float(lo[1]), float(hi[1])],
                per_block_invQ_B=[fit(b)[1] for b in blocks], q_eff_by_window=eff)


def blocks_of(cmp):
    e = np.linspace(-CMP_LIM, CMP_LIM, NBLK + 1); return [np.where((cmp >= e[i]) & (cmp < e[i + 1]))[0] for i in range(NBLK)]


def run(vscale=1.0):
    t0 = time.monotonic(); H, h = read_headers(); sel = np.abs(h["cmp"]) < CMP_LIM; zw = float(np.median(h["wd"][sel]))
    R = Rays(vscale, zw); tabs = [R.table(tc) for tc in CENT]; blocks = []
    for k in blocks_of(h["cmp"]):
        G = traces(H, k); blocks.append(block_rows(G, h["off"][k], tabs)); print("block rows", len(blocks[-1]), "%.0fs" % (time.monotonic() - t0), flush=True)
    out = dict(vscale=vscale, water_depth=zw, zA=float(R.zA), zB=float(R.zB), **summarise(blocks), rows=[dict(r, block=i) for i, b in enumerate(blocks) for r in b])
    tag = "" if vscale == 1.0 else "_v%03d" % round(100 * vscale)
    (ROOT / "results" / ("field79%s.json" % tag)).write_text(json.dumps(out, indent=1))
    print(json.dumps({k: v for k, v in out.items() if k != "rows"}, indent=1))


QS = (60.0, 90.0, 100.0)   # synthetic Q_A, Q_B, Q below B
SYN_SEED = 600_000          # amendment self-check; 500 000 (first self-check and diagnosis) used


def synth_block(args):
    ib, snr = args
    """one CMP block of synthetic traces with the field geometry and mute: sparse reflectivity per 12.5 m CMP bin (seeded by
    the bin), Ricker 25 Hz, amplitude decay exp(-pi f t*), white noise at 24 dB."""
    H, h = read_headers(); sel = np.abs(h["cmp"]) < CMP_LIM; zw = float(np.median(h["wd"][sel])); R = Rays(1.0, zw)
    tk = np.arange(0.25, 2.95, 0.02); xs = np.arange(0, 4101, 25.0)
    grid = np.array([np.stack(R.table(t)(xs), 1) for t in tk])                  # (tk, xs, 4)
    tf = np.arange(0.25, 2.95, DT); fr = np.fft.rfftfreq(2 * NT, DT); W = (fr / 25.0) ** 2 * np.exp(-(fr / 25.0) ** 2)
    k = blocks_of(h["cmp"])[ib]; G = np.zeros((len(k), NT), np.float32); rng = np.random.default_rng(SYN_SEED + 1000 + ib)
    for j, i in enumerate(k):
        b = int(np.round(h["cmp"][i] / 12.5)); rb = np.random.default_rng(SYN_SEED + b); r = (rb.random(tf.size) < 0.3) * rb.standard_normal(tf.size)
        x = h["off"][i]; ix = min(int(round(x / 25.0)), xs.size - 1); v = np.array([np.interp(tf, tk, grid[:, ix, q]) for q in range(4)])
        ok = np.isfinite(v[0]) & (r != 0); tt, la, lb, ls = v[:, ok]; tst = la / QS[0] + lb / QS[1] + (ls - la - lb) / QS[2]
        S = (r[ok][None] * np.exp(-np.pi * fr[:, None] * tst[None] - 2j * np.pi * fr[:, None] * tt[None])).sum(1) * W
        tr = np.fft.irfft(S, 2 * NT)[:NT]; tr += (np.sqrt(np.mean(tr ** 2)) / 10 ** (snr / 20) if np.isfinite(snr) else 0.0) * rng.standard_normal(NT)
        tr[T < MUTE(x)] = 0.0; G[j] = tr
    tabs = [R.table(tc) for tc in CENT]; rows = block_rows(G, h["off"][k], tabs); print("synth block", ib, "rows", len(rows), flush=True)
    for r in rows:   # true y of the row (diagnostic; the amendment's truth = the same single fit on y_true)
        r["y_true"] = float(np.pi * (r["dA"] / QS[0] + r["dB"] / QS[1] + (r["dS"] - r["dA"] - r["dB"]) / QS[2]))
    return rows


def synth(snr=24.0):
    """Self-check: must recover 1000/Q_B = 11.1 (PLAN_FIELD_I_ZH.md section 3)."""
    from multiprocessing import Pool
    with Pool(4) as pool:
        blocks = pool.map(synth_block, [(i, snr) for i in range(NBLK)])
    out = summarise(blocks); out["truth"] = dict(invQ_A=1000 / QS[0], invQ_B=1000 / QS[1]); out["snr"] = snr; out["seed"] = SYN_SEED
    allr = [r for b in blocks for r in b]; x = np.pi * np.array([r["dS"] for r in allr]); yt = np.array([r["y_true"] for r in allr])
    out["truth"]["invQ_eff"] = float(1000 * (x @ yt) / (x @ x)); e = out["invQ_eff"]; tt = out["truth"]["invQ_eff"]
    out["selfcheck_pass"] = bool(abs(e / tt - 1) <= 0.15 and out["ci90_eff"][0] <= tt <= out["ci90_eff"][1])
    out["rows"] = [dict(r, block=i) for i, b in enumerate(blocks) for r in b]
    (ROOT / "results" / ("field79_synth_s%d%s.json" % (SYN_SEED, "" if snr == 24.0 else "_snr%s" % snr))).write_text(json.dumps(out, indent=1)); print(json.dumps({k: v for k, v in out.items() if k != "rows"}, indent=1))


if __name__ == "__main__":
    if sys.argv[1] == "synth":
        synth(float(sys.argv[2]) if len(sys.argv) > 2 else 24.0)
    elif sys.argv[1] == "run":
        run(float(sys.argv[2]) if len(sys.argv) > 2 else 1.0)
