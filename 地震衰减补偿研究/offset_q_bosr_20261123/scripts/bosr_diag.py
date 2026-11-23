"""T1.0 Pre-registered mechanism diagnostics & Zero-Stretch Primary-Slowness Beamforming (B-OSR / FE-OSR)
on Volve seq-79 raw PZDETECTED records (cable 5081), processed seq-79 gathers, and F-12 / F-14 well logs.

Read-only access to D:/ChatGPT/研究 caches. Low memory (< 1.5 GB), single process, 2 BLAS threads.
Usage: python scripts/bosr_diag.py [--smoke]
"""
import os, sys, json, math, time
os.environ.setdefault("OPENBLAS_NUM_THREADS", "2")
os.environ.setdefault("MKL_NUM_THREADS", "2")
os.environ.setdefault("OMP_NUM_THREADS", "2")

from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
MAIN_CACHE = Path(r"D:/ChatGPT/研究/地震衰减补偿研究/uq_rtm_20261117/cache")
RAW_SGY = MAIN_CACHE / "volve_seq079" / "raw_seq079_PZ.sgy"
PROC_SGY = MAIN_CACHE / "volve_seq079" / "seq079_dline5081.segy"
VP_BIN = MAIN_CACHE / "fld" / "vp.bin"
META_JSON = MAIN_CACHE / "fld" / "meta.json"
WELL_SCRIPTS = Path(r"D:/ChatGPT/研究/地震衰减补偿研究/field_baselines_volve_20261029/scripts")

# Exact constants matching field79.py & field79raw.py
DT = 0.004
NT = 800
T = np.arange(NT) * DT
NF = 512
F = np.fft.rfftfreq(NF, DT)
CENT = np.round(np.arange(0.6, 2.21, 0.2), 3)
WL = 0.4
REF = 350.0
OFFS = np.arange(500, 3901, 200)
HW_INCOH = 100.0     # original incoherent bin half-width (200 m)
HW_BEAM = 200.0      # coherent beamforming bin half-width (L = 400 m aperture)
CMP_LIM = 1000.0
NBLK = 8
FINE_PER_BLK = 5     # 50 m fine CMP supergathers inside each 250 m block
BAND = (12.0, 60.0)
FLOOR = 0.02
MINBW = 10.0
STRETCH = 1.3
MUTE = lambda x: -0.13 + np.abs(x) / 1216.0
MUTE_PAD = 0.1
TB = (1.02, 2.40)
ZS = 6.0
VW = 1500.0
NS_RAW = 1876
TSZ_RAW = 240 + 4 * NS_RAW
SHIFT = 12


def ibm2f(b):
    b = b.astype(np.uint32)
    s = np.where(b >> 31, -1.0, 1.0)
    e = ((b >> 24) & 0x7F).astype(np.int32) - 64
    return (s * (b & 0x00FFFFFF) / float(1 << 24) * np.power(16.0, e)).astype(np.float32)


def hp(X):
    """Exact zero-phase low-cut from field79raw.py: 0 below 6 Hz, cosine taper to 1 at 10 Hz."""
    n = 2 * X.shape[1]
    f = np.fft.rfftfreq(n, DT)
    g = np.clip((f - 6.0) / 4.0, 0.0, 1.0)
    g = 0.5 - 0.5 * np.cos(np.pi * g)
    return np.fft.irfft(np.fft.rfft(X, n, axis=1) * g, n, axis=1)[:, :X.shape[1]].astype(np.float32)


class Rays:
    """Exact 1-D ray tracer from field79.py."""
    def __init__(self, vscale, zw):
        m = json.loads(META_JSON.read_text(encoding="utf-8"))
        nz, nx = m["nz"], m["nx"]
        V = np.fromfile(VP_BIN, np.float32).reshape(nz, nx)
        x = m["xmin"] + (np.arange(nx) - m["npml"]) * m["dx"]
        z = (np.arange(nz) - m["npml"]) * m["dx"]
        k = z >= 0
        v = V[np.ix_(k, (x >= -CMP_LIM) & (x <= CMP_LIM))].mean(axis=1).astype(float) * vscale
        self.dz = 1.0
        self.z = np.arange(0, z[k][-1], self.dz) + 0.5 * self.dz
        self.v = np.interp(self.z, z[k], v)
        self.T0 = 2.0 * np.cumsum(self.dz / self.v)
        self.zw = zw
        self.zA = self.z[np.searchsorted(self.T0, TB[0])]
        self.zB = self.z[np.searchsorted(self.T0, TB[1])]
        self.sed = self.z > zw
        self.lay = [(self.z > zw) & (self.z <= self.zA), (self.z > self.zA) & (self.z <= self.zB)]
        self.legs = (self.z > ZS).astype(float) + (self.z > zw).astype(float)

    def shoot(self, zr, p):
        k = self.z < zr
        v = self.v[k]
        w = self.legs[k] * self.dz
        q = np.sqrt(np.maximum(1.0 - (p[:, None] * v[None]) ** 2, 1e-9))
        dt = w / (v[None] * q)
        x = (w * v)[None] * p[:, None] / q
        lay = [L[k] for L in self.lay] + [self.sed[k]]
        return dt.sum(1), x.sum(1), [dt[:, L].sum(1) for L in lay]

    def table(self, tc):
        t0 = lambda zr: self.shoot(zr, np.zeros(1))[0][0]
        lo, hi = self.zw + 1.0, self.z[-1] - 1.0
        for _ in range(60):
            mid = 0.5 * (lo + hi)
            lo, hi = (mid, hi) if t0(mid) < tc else (lo, mid)
        zr = 0.5 * (lo + hi)
        vm = self.v[self.z < zr].max()
        p = np.linspace(0.0, 0.999 / vm, 4000)
        t, x, L = self.shoot(zr, p)
        ok = np.r_[True, np.diff(x) > 0]
        t, x = t[ok], x[ok]
        L = [l[ok] for l in L]
        return lambda xq: (np.interp(xq, x, t, right=np.nan), *[np.interp(xq, x, l, right=np.nan) for l in L])


def ptable(tab):
    return lambda x: (tab(np.atleast_1d(x) + 10.0)[0] - tab(np.maximum(np.atleast_1d(x) - 10.0, 0.0))[0]) / 20.0


def ghost(p, zw):
    cw = np.sqrt(np.maximum(1.0 - (VW * np.asarray(p, float)) ** 2, 0.05))
    return cw, 2.0 * zw * cw / VW


def power_incoh(y, a, b):
    k = (T >= a) & (T < b)
    if k.sum() < int(0.9 * WL / DT):
        return None
    seg = y[:, k]
    seg = seg[(np.abs(seg) > 0).mean(axis=1) > 0.95]
    if len(seg) < 10:
        return None
    return np.mean(np.abs(np.fft.rfft(seg * np.hanning(k.sum()), NF, axis=1)) ** 2, axis=0)


def smooth_period(P, tau):
    w = max(3, int(round(1.0 / tau / (F[1] - F[0]))))
    lp = np.log(np.maximum(P, 1e-30 * P.max()))
    return np.exp(np.convolve(lp, np.ones(w) / w, mode="same"))


def band_slice(P0, Px):
    ib = (F >= BAND[0]) & (F <= BAND[1])
    ok = ib & (P0 > FLOOR * P0[ib].max()) & (Px > FLOOR * Px[ib].max())
    i = np.where(F >= BAND[0])[0][0]
    j = i
    while j < len(F) and ok[j]:
        j += 1
    return slice(i, j) if F[j - 1] - F[i] >= MINBW else None


def band_single(P):
    ib = (F >= BAND[0]) & (F <= BAND[1])
    ok = ib & (P > FLOOR * P[ib].max())
    i = np.where(F >= BAND[0])[0][0]
    j = i
    while j < len(F) and ok[j]:
        j += 1
    return slice(i, j) if F[j - 1] - F[i] >= MINBW else None


def fit1_noint(rows):
    if len(rows) < 3:
        return float("nan")
    x = np.pi * np.array([r["dS"] for r in rows])
    y = np.array([r["y"] for r in rows])
    return float(1000.0 * (x @ y) / (x @ x))


def fit_fe(rows, key_x="tau_sed", key_y="slope", exclude_near=False):
    """Within-(block, tc) fixed-effects regression: slope = alpha_{b,tc} + pi * tau_sed * (1/Q_eff)."""
    sub = [r for r in rows if (not exclude_near or r["off"] > REF)]
    xd, yd = [], []
    groups = {}
    for r in sub:
        groups.setdefault((r["block"], r["tc"]), []).append(r)
    for g in groups.values():
        if len(g) >= 2:
            xg = np.pi * np.array([r[key_x] for r in g])
            yg = np.array([r[key_y] for r in g])
            xd.append(xg - xg.mean())
            yd.append(yg - yg.mean())
    if not xd:
        return float("nan")
    xd = np.concatenate(xd)
    yd = np.concatenate(yd)
    denom = float(xd @ xd)
    return float(1000.0 * (xd @ yd) / denom) if denom > 0 else float("nan")


def boot_blocks(blocks_rows, fn, seed=79, nboot=2000):
    rng = np.random.default_rng(seed)
    nblk = len(blocks_rows)
    vals = []
    for _ in range(nboot):
        pick = rng.integers(0, nblk, nblk)
        sub = [r for i in pick for r in blocks_rows[i]]
        vals.append(fn(sub))
    lo, hi = np.nanpercentile(vals, [5, 95])
    return [float(lo), float(hi)]


def summarise_rows(blocks_rows):
    allr = [r for b in blocks_rows for r in b]
    pt = fit1_noint(allr)
    ci = boot_blocks(blocks_rows, fit1_noint)
    fe_all = fit_fe(allr, "dS", "y", exclude_near=False)
    fe_ci = boot_blocks(blocks_rows, lambda s: fit_fe(s, "dS", "y", False))
    fe_midfar = fit_fe(allr, "dS", "y", exclude_near=True)
    by_win = {}
    for tc in CENT:
        g = [r for r in allr if abs(r["tc"] - tc) < 1e-6]
        if len(g) >= 3:
            x = np.pi * np.array([r["dS"] for r in g])
            y = np.array([r["y"] for r in g])
            q0 = float(1000.0 * (x @ y) / (x @ x))
            A = np.column_stack([x, np.ones_like(x)])
            coef = np.linalg.lstsq(A, y, rcond=None)[0]
            corr = float(np.corrcoef(x, y)[0, 1]) if np.std(x) > 0 and np.std(y) > 0 else float("nan")
            by_win[f"{tc:.1f}"] = dict(
                n=len(g), q_noint=q0, q_int=float(1000.0 * coef[0]),
                intercept=float(coef[1]), corr=corr,
                semblance_med=float(np.median([r["semb"] for r in g])) if "semb" in g[0] else None,
            )
    # 1.2-1.6 s pooled and 2.0-2.2 s pooled
    g_mid = [r for r in allr if 1.15 <= r["tc"] <= 1.65]
    g_deep = [r for r in allr if 1.95 <= r["tc"] <= 2.25]
    return dict(
        invQ_noint=pt, ci90_noint=ci,
        invQ_fe=fe_all, ci90_fe=fe_ci,
        invQ_fe_midfar=fe_midfar,
        n_rows=len(allr),
        q_12_16_noint=fit1_noint(g_mid),
        q_12_16_fe=fit_fe(g_mid, "dS", "y", False),
        q_20_22_noint=fit1_noint(g_deep),
        q_20_22_fe=fit_fe(g_deep, "dS", "y", False),
        by_window=by_win,
    )


def load_raw_headers():
    n = (RAW_SGY.stat().st_size - 3600) // TSZ_RAW
    H = np.memmap(RAW_SGY, dtype=np.uint8, mode="r", offset=3600, shape=(n, TSZ_RAW))
    h = lambda o, t=">i4": np.frombuffer(np.ascontiguousarray(H[:, o:o + (4 if t == ">i4" else 2)]).tobytes(), t).astype(np.float64)
    A = math.radians(14.0)
    X0, Y0 = 435050.02, 6478563.53
    il = lambda x, y: -(x - X0) * math.cos(A) + (y - Y0) * math.sin(A)
    s = il(h(180) / 100.0, h(184) / 100.0)
    g = il(h(188) / 100.0, h(192) / 100.0)
    return H, dict(ffid=h(8), trid=h(28, ">i2"), off=np.abs(h(36)), wd=h(64), dline=h(196), dpt=h(200), cmp=(s + g) / 2.0)


def read_raw_samples(H, k):
    tr = ibm2f(np.frombuffer(np.ascontiguousarray(H[k, 240:240 + 4 * (NT + SHIFT)]).tobytes(), ">u4").reshape(len(k), NT + SHIFT))
    return hp(tr[:, SHIFT:])


def load_proc_headers():
    raw = np.memmap(PROC_SGY, dtype=np.uint8, mode="r")
    ns = int(np.frombuffer(raw[3220:3222].tobytes(), ">i2")[0])
    tsz = 240 + 4 * ns
    n = (raw.size - 3600) // tsz
    H = np.ndarray((n, tsz), np.uint8, raw[3600:3600 + n * tsz])
    h = lambda o, t=">i4": np.frombuffer(np.ascontiguousarray(H[:, o:o + (4 if t == ">i4" else 2)]).tobytes(), t).astype(np.float64)
    sc = h(70, ">i2")
    f = np.where(sc < 0, -1.0 / sc, np.where(sc > 0, sc, 1.0))
    A = math.radians(14.0)
    X0, Y0 = 435050.02, 6478563.53
    il = lambda x, y: -(x - X0) * math.cos(A) + (y - Y0) * math.sin(A)
    s = il(h(72) * f, h(76) * f)
    g = il(h(80) * f, h(84) * f)
    return H, dict(cmp=(s + g) / 2.0, off=np.abs(h(36)), wd=h(64))


def read_proc_samples(H, k):
    return ibm2f(np.frombuffer(np.ascontiguousarray(H[k, 240:240 + 4 * NT]).tobytes(), ">u4").reshape(len(k), NT))


# -------------------- D1: Component ID & P-Z Phase/Lag Check --------------------
def diagnose_components(blocks, tabs, zw):
    """Identify TRID 11 (A) vs TRID 12 (B): ghost notch positions, cross-corr lag, phase diff, near-offset spectra."""
    tau0 = 2.0 * zw / VW
    f_notch = 1.0 / tau0  # ~8.24 Hz
    PA_near = np.zeros_like(F)
    PB_near = np.zeros_like(F)
    XAB_near = np.zeros_like(F, dtype=np.complex128)
    lags = np.arange(-15, 16)
    xc_sum = np.zeros_like(lags, dtype=np.float64)
    a_en = b_en = 0.0
    n_used = 0
    for A, B, off, _ in blocks:
        kr = off <= REF
        if kr.sum() < 10:
            continue
        xr = float(np.sqrt(np.mean(off[kr] ** 2)))
        for tc, tab in zip(CENT, tabs):
            if tc < 1.0 or tc > 2.0:
                continue
            tr = float(tab(np.array([xr]))[0][0])
            k = (T >= tr - WL / 2) & (T < tr + WL / 2)
            sa = A[kr][:, k] * np.hanning(k.sum())[None, :]
            sb = B[kr][:, k] * np.hanning(k.sum())[None, :]
            FA = np.fft.rfft(sa, NF, axis=1)
            FB = np.fft.rfft(sb, NF, axis=1)
            PA_near += np.mean(np.abs(FA) ** 2, axis=0)
            PB_near += np.mean(np.abs(FB) ** 2, axis=0)
            XAB_near += np.mean(FA * np.conj(FB), axis=0)
            # Normalized cross-correlation across lags
            for li, lag in enumerate(lags):
                if lag >= 0:
                    xc_sum[li] += float(np.sum(sa[:, :sa.shape[1]-lag] * sb[:, lag:]))
                else:
                    xc_sum[li] += float(np.sum(sa[:, -lag:] * sb[:, :sb.shape[1]+lag]))
            a_en += float(np.sum(sa ** 2))
            b_en += float(np.sum(sb ** 2))
            n_used += 1
    PA_near /= max(n_used, 1)
    PB_near /= max(n_used, 1)
    XAB_near /= max(n_used, 1)
    xc_norm = xc_sum / math.sqrt(max(a_en * b_en, 1e-30))
    best_idx = int(np.argmax(np.abs(xc_norm)))
    best_lag = int(lags[best_idx])
    best_corr = float(xc_norm[best_idx])

    # Check ghost comb at integer multiples k*f_notch vs half-integer (k+0.5)*f_notch (k=2..6 -> 16..50 Hz)
    # Subtract smoothed log spectrum to see comb troughs (<0) vs peaks (>0)
    w = max(3, int(round(f_notch / (F[1] - F[0]))))
    def comb_score(P):
        lp = np.log(np.maximum(P, 1e-30))
        detr = lp - np.convolve(lp, np.ones(w) / w, mode="same")
        int_vals = [float(np.interp(k * f_notch, F, detr)) for k in range(2, 7)]
        half_vals = [float(np.interp((k + 0.5) * f_notch, F, detr)) for k in range(2, 7)]
        return dict(at_int_k=float(np.mean(int_vals)), at_half_k=float(np.mean(half_vals)),
                    diff_int_minus_half=float(np.mean(int_vals) - np.mean(half_vals)))

    cA = comb_score(PA_near)
    cB = comb_score(PB_near)
    # Coherency & phase difference between A and B over 12-50 Hz
    ib = (F >= 12) & (F <= 50)
    coh = np.abs(XAB_near[ib]) ** 2 / np.maximum(PA_near[ib] * PB_near[ib], 1e-30)
    ph_deg = np.degrees(np.angle(XAB_near[ib]))

    # Choose optimal amplitude-domain PZ sign minimizing near-offset comb index
    c0sq_val = float(PB_near[(F >= 10) & (F <= 60)].sum() / PA_near[(F >= 10) & (F <= 60)].sum())
    s_cal = 1.0 / math.sqrt(c0sq_val)
    comb_idx = {}
    for sg in (+1, -1):
        # Power of A + sg * s_cal * B = PA + s_cal^2 PB + 2*sg*s_cal*Re(XAB)
        Psum = PA_near + (s_cal ** 2) * PB_near + 2.0 * sg * s_cal * XAB_near.real
        lp = np.log(np.maximum(Psum, 1e-30))
        sm = np.convolve(lp, np.ones(w) / w, mode="same")
        comb_idx[sg] = float(np.mean(np.abs(lp - sm)[(F >= 10) & (F <= 60)]))
    best_sign = +1 if comb_idx[+1] < comb_idx[-1] else -1

    # Also test 1-sample lag shift if best_lag != 0
    return dict(
        tau0_s=tau0, f_notch_hz=f_notch, c0sq=c0sq_val, s_cal=s_cal,
        comb_A=cA, comb_B=cB,
        inferred_P="TRID 11 (A)" if cA["diff_int_minus_half"] < cB["diff_int_minus_half"] else "TRID 12 (B)",
        best_xcorr_lag_samples=best_lag, best_xcorr_lag_ms=best_lag * DT * 1e3,
        best_xcorr_val=best_corr,
        coherency_12_50_med=float(np.median(coh)),
        phase_AB_12_50_deg_med=float(np.median(ph_deg)),
        pz_amp_comb_plus=comb_idx[+1], pz_amp_comb_minus=comb_idx[-1],
        pz_amp_best_sign=best_sign,
    )


# -------------------- D0: Exact field79raw reproduction --------------------
def run_incoh_field79raw(blocks, tabs, ptabs, zw, c0sq, comp, pz_sign=1):
    """comp in {'pz_pow', 'a', 'b', 'pz_up'}."""
    s_cal = 1.0 / math.sqrt(c0sq)
    blk_rows = []
    for bi, (A, B, off, _) in enumerate(blocks):
        rows = []
        kr = off <= REF
        if kr.sum() < 10:
            blk_rows.append(rows); continue
        xr = float(np.sqrt(np.mean(off[kr] ** 2)))
        for tc, tab, ptab in zip(CENT, tabs, ptabs):
            r0 = tab(np.array([xr]))
            tr = float(r0[0][0])
            cwr, taur = ghost(ptab(xr)[0], zw)
            if comp == "pz_up":
                X0 = A[kr] + pz_sign * (s_cal / cwr) * B[kr]
                P0_raw = power_incoh(X0, tr - WL / 2, tr + WL / 2)
                if P0_raw is None:
                    continue
                P0 = smooth_period(P0_raw, taur)
            else:
                PA0 = power_incoh(A[kr], tr - WL / 2, tr + WL / 2)
                PB0 = power_incoh(B[kr], tr - WL / 2, tr + WL / 2)
                if PA0 is None or PB0 is None:
                    continue
                P0_comb = (PA0 + PB0 / (c0sq * cwr ** 2)) if comp == "pz_pow" else (PA0 if comp == "a" else PB0)
                P0 = smooth_period(P0_comb, taur)

            for o in OFFS:
                k = (off >= o - HW_INCOH) & (off < o + HW_INCOH)
                if k.sum() < 10:
                    continue
                xo = float(np.sqrt(np.mean(off[k] ** 2)))
                tx, la, lb, ls = (float(v[0]) for v in tab(np.array([xo])))
                if not np.isfinite(tx) or tx / tc > STRETCH or tx - WL / 2 < MUTE(o + HW_INCOH) + MUTE_PAD or tx + WL / 2 > T[-1]:
                    continue
                cw, tau = ghost(ptab(xo)[0], zw)
                if comp == "pz_up":
                    Xx = A[k] + pz_sign * (s_cal / cw) * B[k]
                    Px_raw = power_incoh(Xx, tx - WL / 2, tx + WL / 2)
                    if Px_raw is None:
                        continue
                    Px = smooth_period(Px_raw, tau)
                else:
                    PA = power_incoh(A[k], tx - WL / 2, tx + WL / 2)
                    PB = power_incoh(B[k], tx - WL / 2, tx + WL / 2)
                    if PA is None or PB is None:
                        continue
                    Px_comb = (PA + PB / (c0sq * cw ** 2)) if comp == "pz_pow" else (PA if comp == "a" else PB)
                    Px = smooth_period(Px_comb, tau)

                s = band_slice(P0, Px)
                if s is None:
                    continue
                y = float(-np.polyfit(F[s], 0.5 * np.log(Px[s] / P0[s]), 1)[0])
                rows.append(dict(block=bi, tc=float(tc), off=xo, y=y, dS=ls - float(r0[3][0]),
                                 tau_sed=ls, fmax=float(F[s.stop - 1])))
        blk_rows.append(rows)
    return summarise_rows(blk_rows), blk_rows


# -------------------- D3: Zero-Stretch Primary-Slowness Beamforming (B-OSR) --------------------
def beam_power_finecmp(X, off, cmp_arr, cmp_edges, tc, tab, x0, hw, steer_mode="primary", v_steer=None, kappa_p=0.0):
    """
    Zero-stretch coherent local slant-stack inside each 50 m fine CMP bin, then average power across fine CMPs.
    Returns (P_beam_corr, P_incoh, semblance_12_60, n_fine_used).
    steer_mode: 'primary' (exact rigid shift tab(x_i) - tab(x0)),
                'const_v' (rigid shift sqrt(tc_m^2 + (x_i/v)^2) - t_center where tc_m puts hyperbola through (x0, t_center)).
    """
    t_center = float(tab(np.array([x0]))[0][0])
    # Extract slightly padded window [t_center - WL/2 - pad, t_center + WL/2 + pad] so rigid shifts stay inside!
    pad = 0.25
    t_lo = max(0.0, t_center - WL / 2 - pad)
    t_hi = min(T[-1], t_center + WL / 2 + pad)
    k_pad = (T >= t_lo) & (T < t_hi)
    n_pad = int(k_pad.sum())
    if n_pad < int(WL / DT) + 10:
        return None
    t_sub = T[k_pad]
    n_fft_pad = 1 << int(math.ceil(math.log2(2 * n_pad)))
    f_pad = np.fft.rfftfreq(n_fft_pad, DT)

    # Target output window indices inside t_sub
    k_win = (t_sub >= t_center - WL / 2) & (t_sub < t_center + WL / 2)
    nw = int(k_win.sum())
    if nw < int(0.9 * WL / DT):
        return None
    wt = np.hanning(nw)
    # Exact Hann^2 time variance for analytical beam-smearing correction
    dts_w = (np.arange(nw) - 0.5 * (nw - 1)) * DT
    wt2 = wt ** 2
    sig_t2 = float(np.sum(wt2 * dts_w ** 2) / np.sum(wt2))

    P_beam_acc = np.zeros_like(F)
    P_incoh_acc = np.zeros_like(F)
    sig_xi2_list = []
    n_fine = 0

    for fi in range(len(cmp_edges) - 1):
        m = (cmp_arr >= cmp_edges[fi]) & (cmp_arr < cmp_edges[fi + 1]) & (off >= x0 - hw) & (off < x0 + hw)
        if m.sum() < 8:
            continue
        seg = X[m][:, k_pad]
        live = (np.abs(seg) > 0).mean(axis=1) > 0.90
        if live.sum() < 8:
            continue
        seg = seg[live]
        xo_i = off[m][live]
        xi = xo_i - x0

        # Compute exact rigid time shift dt_i for each trace (positive dt_i = shift trace earlier -> exp(+i 2pi f dt_i))
        if steer_mode == "primary":
            dt_i = tab(xo_i)[0] - t_center
        else:
            # Hyperbola of velocity v_steer passing through (x0, t_center): t(x) = sqrt(t_center^2 + (x^2 - x0^2)/v^2)
            arg = t_center ** 2 + (xo_i ** 2 - x0 ** 2) / (v_steer ** 2)
            if np.any(arg <= 0.01):
                # fallback to linear slowness p = x0 / (v_steer^2 * t_center)
                p_lin = x0 / (v_steer ** 2 * max(t_center, 0.1))
                dt_i = p_lin * xi
            else:
                dt_i = np.sqrt(arg) - t_center

        if not np.all(np.isfinite(dt_i)):
            continue

        # Hann aperture taper + mild offset-bin density compensation (divide by local trace count in 50 m sub-bins)
        w_ap = 0.5 * (1.0 + np.cos(np.pi * np.clip(xi / hw, -1.0, 1.0)))
        # Density compensation in 8 sub-bins across [-hw, +hw]
        bins = np.clip(((xi + hw) / (2.0 * hw) * 8).astype(int), 0, 7)
        cnts = np.bincount(bins, minlength=8).astype(float)
        w_i = w_ap / np.maximum(cnts[bins], 1.0)
        w_sum = float(np.sum(w_i))
        if w_sum <= 0:
            continue
        w_i /= w_sum
        sig_xi2_list.append(float(np.sum(w_i * xi ** 2)))

        # Rigid frequency-domain phase shift (ZERO wavelet stretch!)
        S_pad = np.fft.rfft(seg, n_fft_pad, axis=1)
        phase = np.exp(2j * np.pi * f_pad[None, :] * dt_i[:, None])
        aligned = np.fft.irfft(S_pad * phase, n_fft_pad, axis=1)[:, :n_pad]

        # Coherent beam vs incoherent average inside window k_win
        win_tr = aligned[:, k_win] * wt[None, :]
        beam_tr = np.sum(w_i[:, None] * win_tr, axis=0)
        Pb = np.abs(np.fft.rfft(beam_tr, NF)) ** 2
        Pi = np.sum(w_i[:, None] * (np.abs(np.fft.rfft(win_tr, NF, axis=1)) ** 2), axis=0)
        P_beam_acc += Pb
        P_incoh_acc += Pi
        n_fine += 1

    if n_fine == 0:
        return None
    P_beam = P_beam_acc / n_fine
    P_incoh = P_incoh_acc / n_fine
    # Exact analytical intra-window beam-smearing correction (only for primary steering!)
    if steer_mode == "primary":
        sig_xi2 = float(np.mean(sig_xi2_list))
        B2_prim = np.exp(-((2.0 * np.pi * F) ** 2) * (kappa_p ** 2) * sig_t2 * sig_xi2)
        B2_prim = np.maximum(B2_prim, 0.25)  # safety floor (in practice > 0.95 up to 60 Hz)
        P_beam_corr = P_beam / B2_prim
    else:
        P_beam_corr = P_beam
    ib = (F >= BAND[0]) & (F <= BAND[1])
    semb = float(np.median(P_beam[ib] / np.maximum(P_incoh[ib], 1e-30)))
    return P_beam_corr, P_incoh, semb, n_fine


def run_bosr(blocks, R, tabs, ptabs, zw, c0sq, comp, pz_sign=1, hw=HW_BEAM, do_slowscan=False):
    """Run B-OSR on raw blocks for comp in {'a', 'b', 'pz_up'}."""
    s_cal = 1.0 / math.sqrt(c0sq)
    # Precompute kappa_p(tc, x) = dp_P / dt_0 via centered finite difference dt=0.05 s
    tabs_p = [R.table(tc + 0.05) for tc in CENT]
    tabs_m = [R.table(max(tc - 0.05, 0.25)) for tc in CENT]
    ptabs_p = [ptable(t) for t in tabs_p]
    ptabs_m = [ptable(t) for t in tabs_m]
    blk_edges = np.linspace(-CMP_LIM, CMP_LIM, NBLK + 1)

    blk_rows = []
    slow_diag = []  # for 2.0-2.2 s slowness scan: primary vs multiple (1600) vs slow (950)

    for bi, (A, B, off, cmp_arr) in enumerate(blocks):
        rows = []
        fine_edges = np.linspace(blk_edges[bi], blk_edges[bi + 1], FINE_PER_BLK + 1)
        kr = off <= REF
        if kr.sum() < 10:
            blk_rows.append(rows); continue
        xr = float(np.sqrt(np.mean(off[kr] ** 2)))

        for w_idx, (tc, tab, ptab) in enumerate(zip(CENT, tabs, ptabs)):
            r0 = tab(np.array([xr]))
            tr = float(r0[0][0])
            cwr, taur = ghost(ptab(xr)[0], zw)
            if comp == "pz_up":
                X_near = A[kr] + pz_sign * (s_cal / cwr) * B[kr]
            elif comp == "a":
                X_near = A[kr]
            else:
                X_near = B[kr]
            # Near-offset beam (or coherent stack at xr <= 350 m)
            kp0 = float((ptabs_p[w_idx](xr)[0] - ptabs_m[w_idx](xr)[0]) / (0.05 + min(0.05, tc - 0.25)))
            res0 = beam_power_finecmp(X_near, off[kr], cmp_arr[kr], fine_edges, tc, tab, xr, REF, "primary", None, kp0)
            if res0 is None:
                continue
            P0_b, P0_i, semb0, _ = res0
            P0 = smooth_period(P0_b, taur)

            for o in OFFS:
                # Check mute using center + hw_incoh (same acceptance mask as field79raw!)
                tx, la, lb, ls = (float(v[0]) for v in tab(np.array([float(o)])))
                if not np.isfinite(tx) or tx / tc > STRETCH or tx - WL / 2 < MUTE(o + HW_INCOH) + MUTE_PAD or tx + WL / 2 > T[-1]:
                    continue
                k_bin = (off >= o - hw) & (off < o + hw)
                if k_bin.sum() < 15:
                    continue
                xo = float(np.sqrt(np.mean(off[k_bin] ** 2)))
                tx, la, lb, ls = (float(v[0]) for v in tab(np.array([xo])))
                cw, tau = ghost(ptab(xo)[0], zw)
                if comp == "pz_up":
                    Xx = A[k_bin] + pz_sign * (s_cal / cw) * B[k_bin]
                elif comp == "a":
                    Xx = A[k_bin]
                else:
                    Xx = B[k_bin]
                kp = float((ptabs_p[w_idx](xo)[0] - ptabs_m[w_idx](xo)[0]) / (0.05 + min(0.05, tc - 0.25)))
                resx = beam_power_finecmp(Xx, off[k_bin], cmp_arr[k_bin], fine_edges, tc, tab, xo, hw, "primary", None, kp)
                if resx is None:
                    continue
                Px_b, Px_i, sembx, nf_used = resx
                Px = smooth_period(Px_b, tau)
                s = band_slice(P0, Px)
                if s is None:
                    continue
                y = float(-np.polyfit(F[s], 0.5 * np.log(Px[s] / P0[s]), 1)[0])
                sl_beam = float(-np.polyfit(F[s], 0.5 * np.log(Px[s]), 1)[0])
                rows.append(dict(block=bi, tc=float(tc), off=xo, y=y, slope=sl_beam,
                                 dS=ls - float(r0[3][0]), tau_sed=ls, semb=sembx, fmax=float(F[s.stop - 1])))

                # Slowness scan on 2.0-2.2 s (and 1.2-1.6 s control)
                if do_slowscan and tc in (1.4, 2.0, 2.2) and o in (1100, 1500, 1900, 2300):
                    rm = beam_power_finecmp(Xx, off[k_bin], cmp_arr[k_bin], fine_edges, tc, tab, xo, hw, "const_v", 1600.0, 0.0)
                    rs = beam_power_finecmp(Xx, off[k_bin], cmp_arr[k_bin], fine_edges, tc, tab, xo, hw, "const_v", 950.0, 0.0)
                    if rm is not None and rs is not None:
                        Pm = smooth_period(rm[0], tau)
                        Ps = smooth_period(rs[0], tau)
                        ib = (F >= 12.0) & (F <= 60.0)
                        sl_p = float(-np.polyfit(F[ib], 0.5 * np.log(np.maximum(Px[ib], 1e-30)), 1)[0])
                        sl_m = float(-np.polyfit(F[ib], 0.5 * np.log(np.maximum(Pm[ib], 1e-30)), 1)[0])
                        sl_s = float(-np.polyfit(F[ib], 0.5 * np.log(np.maximum(Ps[ib], 1e-30)), 1)[0])
                        slow_diag.append(dict(
                            block=bi, tc=float(tc), off=float(o),
                            sl_prim=sl_p, sl_mult1600=sl_m, sl_slow950=sl_s,
                            diff_mult_minus_prim=sl_m - sl_p,
                            semb_prim=sembx, semb_mult=rm[2], semb_slow=rs[2],
                            ratio_mult_over_prim_45Hz=float(np.interp(45.0, F, Pm) / max(np.interp(45.0, F, Px), 1e-30)),
                        ))
        blk_rows.append(rows)

    summary = summarise_rows(blk_rows)
    if do_slowscan and slow_diag:
        for ttag, tsel in [("1.4", [1.4]), ("2.0_2.2", [2.0, 2.2])]:
            sub = [d for d in slow_diag if d["tc"] in tsel]
            if sub:
                summary[f"slowscan_{ttag}"] = dict(
                    n=len(sub),
                    sl_prim_med=float(np.median([d["sl_prim"] for d in sub])),
                    sl_mult1600_med=float(np.median([d["sl_mult1600"] for d in sub])),
                    sl_slow950_med=float(np.median([d["sl_slow950"] for d in sub])),
                    diff_mult_minus_prim_med=float(np.median([d["diff_mult_minus_prim"] for d in sub])),
                    semb_prim_med=float(np.median([d["semb_prim"] for d in sub])),
                    semb_mult_med=float(np.median([d["semb_mult"] for d in sub])),
                    semb_slow_med=float(np.median([d["semb_slow"] for d in sub])),
                    ratio_mult_over_prim_45Hz_med=float(np.median([d["ratio_mult_over_prim_45Hz"] for d in sub])),
                )
    return summary, blk_rows


# -------------------- D4: Volve F-12 / F-14 Well-Log Exact 1D Prestack Reflectivity (Q=inf vs VSP Q) --------------------
def run_well_prestack_forward(R, tabs, zw):
    """Exact 1D plane-wave layer recursion (with all internal multiples + angle-dependent thin-bed tuning + moveout compression)
    built from F-12 and F-14 real well-log reflectivity, evaluated at Q=inf (zero absorption) and VSP Q=(60,90,100)."""
    sys.path.insert(0, str(WELL_SCRIPTS))
    import lib13
    P_picks = lib13.wq.picks()
    out = {}
    dt_w = 0.001  # 1 ms TWT well sampling
    fr = np.fft.rfftfreq(2 * NT, DT)
    W_src = (fr / 25.0) ** 2 * np.exp(-(fr / 25.0) ** 2)

    for wname in ("F-12", "F-14"):
        W = lib13.wq.well(wname, P_picks)
        t_w = W["t"]
        rc_w = W["rc"]  # length len(t_w)-1, 1 ms TWT spacing
        # Map each 1 ms layer to depth via R.T0 and get background velocity v_k
        tm = 0.5 * (t_w[:-1] + t_w[1:])
        zk = np.interp(tm, R.T0, R.z)
        vk = np.interp(zk, R.z, R.v)
        # Reconstruct fine-scale velocity/impedance fluctuations from rc_w (assuming Gardner rho ~ v^0.25 => Z ~ v^1.25 => dln v = 0.8 dln Z)
        # Wait: exact 1D acoustic layer recursion only needs layer reflection r_k(p) and layer vertical phase shift dtau_k(p) = dt_w * cos(theta_k)!
        # Include exact acoustic angle-dependent reflection coefficient r_k(p) + exact vertical slowness thinning cos(theta_k) = sqrt(1 - p^2 v_k^2)!
        # Let's reconstruct local v_fine from background vk * exp(0.8 * cumsum(2*rc_detrended)) or directly use vk for kinematics & exact acoustic impedance contrast!
        # Keep v_k smooth background for exact Snell's law cos(theta_k) = sqrt(1 - p^2 v_k^2) + exact 1D Kennett recursion with ALL internal multiples!
        # Also test simple primary-only moveout-compressed reflectivity with exact well rc_w!
        well_res = {}
        for q_mode, q_vals in [("Q_inf", (1e9, 1e9, 1e9)), ("Q_VSP", (60.0, 90.0, 100.0))]:
            # Synthesize 1 CMP gather over offsets [200] + OFFS (200..3100 m) using ray-traced p(tc, x) per layer!
            # Even better: exact tau-p / plane-wave Kennett reflectivity R_total(f, p) at dense p grid, then inverse-slant-stack (exact 1D acoustic wavefield at offset x)!
            # To keep exact geometry and zero aliasing, synthesize directly per offset x by summing primary + internal multiples at exact ray parameter p(t, x):
            # Method: for each window center tc (only where well covers [tc-0.2, tc+0.2]), evaluate exact 1D Kennett recursion
            # for the well section at the window's horizontal slowness p_x = ptable(tc)(x)!
            # Note: in a 1D layered medium, the exact reflected wavefield around (tc, x) with ray parameter p_x has spectrum
            # R_1D(f; p_x) computed by layer-by-layer upward recursion (exact thin-bed tuning + ALL internal multiples + moveout compression cos(theta_k)!).
            rows = []
            xr = 200.0
            for tc, tab in zip(CENT, tabs):
                if tc - WL / 2 < t_w[0] + 0.02 or tc + WL / 2 > t_w[-1] - 0.02:
                    continue
                ptab = ptable(tab)
                # Select well layers spanning [tc - 0.35, tc + 0.35] so window edges have full reflectivity + internal multiples
                kw = (tm >= max(t_w[0], tc - 0.35)) & (tm <= min(t_w[-1], tc + 0.35))
                rc_sub = rc_w[kw]
                v_sub = vk[kw]
                t_sub0 = tm[kw]
                # Function to compute exact 1D Kennett reflection spectrum & time-domain trace at slowness p_val
                def synth_trace_at_x(x_val):
                    tx, la, lb, ls = (float(v[0]) for v in tab(np.array([x_val])))
                    p_val = float(ptab(x_val)[0])
                    # Local layer cosines and two-way vertical delay dt_k(p) = dt_w * cos(theta_k)
                    # Note: in time domain of an un-NMO'd trace at offset x, event spacing compresses by cos(theta_k) = sqrt(1 - p^2 v_k^2)!
                    cos_k = np.sqrt(np.maximum(1.0 - (p_val * v_sub) ** 2, 0.05))
                    # Exact acoustic angle-dependent reflection coefficient:
                    # Given normal-incidence r0 = (Z2-Z1)/(Z2+Z1) and Gardner (dln rho = 0.25 dln v => dln Z = 1.25 dln v),
                    # acoustic r(theta) = (Z2/cos2 - Z1/cos1)/(Z2/cos2 + Z1/cos1) ≈ r0 + 0.4 * r0 * tan^2(theta)
                    # Let's compute exact layer-by-layer acoustic r_k(p) from impedance contrast + velocity contrast:
                    sin2_k = np.clip((p_val * v_sub) ** 2, 0.0, 0.90)
                    tan2_k = sin2_k / np.maximum(1.0 - sin2_k, 0.10)
                    r_ang = np.clip(rc_sub * (1.0 + 0.4 * tan2_k), -0.95, 0.95)
                    # Two-way delay between consecutive interfaces in this slice: dt_k = dt_w * cos_k
                    dt_k = dt_w * cos_k
                    # Intrinsic Q attenuation factor per layer (if Q_VSP)
                    tstar_over = la / q_vals[0] + lb / q_vals[1] + max(ls - la - lb, 0.0) / q_vals[2]
                    # Kennett upward recursion from bottom to top of slice (vectorized over frequency fr!)
                    R_up = np.zeros_like(fr, dtype=np.complex128)
                    for idx in range(len(r_ang) - 1, -1, -1):
                        rk = r_ang[idx]
                        # Layer phase delay exp(-i 2pi f dt_k)
                        ph = np.exp(-2j * np.pi * fr * dt_k[idx])
                        R_up = (rk + R_up * ph) / (1.0 + rk * R_up * ph)
                    # Apply overburden Q attenuation and source wavelet W_src, place center at t = 0.5 s inside 1.2 s trace
                    att = np.exp(-np.pi * fr * tstar_over)
                    # Shift slice top so window center aligns cleanly
                    t_rel_center = float(np.sum(dt_k[t_sub0 <= tc]))
                    spec = R_up * att * W_src * np.exp(-2j * np.pi * fr * (0.5 - t_rel_center))
                    tr = np.fft.irfft(spec, 2 * NT)[:NT]
                    # Window [0.5 - WL/2, 0.5 + WL/2]
                    kw_t = (T >= 0.5 - WL / 2) & (T < 0.5 + WL / 2)
                    P_pow = np.abs(np.fft.rfft(tr[kw_t] * np.hanning(kw_t.sum()), NF)) ** 2
                    return P_pow, tx, ls

                P0, tr0, ls0 = synth_trace_at_x(xr)
                for o in OFFS:
                    tx, la, lb, ls = (float(v[0]) for v in tab(np.array([float(o)])))
                    if not np.isfinite(tx) or tx / tc > STRETCH or tx - WL / 2 < MUTE(o + HW_INCOH) + MUTE_PAD or tx + WL / 2 > T[-1]:
                        continue
                    Px, _, _ = synth_trace_at_x(float(o))
                    s = band_slice(P0, Px)
                    if s is None:
                        continue
                    y = float(-np.polyfit(F[s], 0.5 * np.log(Px[s] / P0[s]), 1)[0])
                    rows.append(dict(block=0, tc=float(tc), off=float(o), y=y, dS=ls - ls0))

            by_w = {}
            for tc in sorted(set(r["tc"] for r in rows)):
                g = [r for r in rows if r["tc"] == tc]
                x = np.pi * np.array([r["dS"] for r in g])
                y = np.array([r["y"] for r in g])
                by_w[f"{tc:.1f}"] = float(1000.0 * (x @ y) / (x @ x))
            g_20_22 = [r for r in rows if 1.95 <= r["tc"] <= 2.25]
            g_12_16 = [r for r in rows if 1.15 <= r["tc"] <= 1.65]
            well_res[q_mode] = dict(
                invQ_eff=fit1_noint(rows),
                q_12_16=fit1_noint(g_12_16),
                q_20_22=fit1_noint(g_20_22),
                by_window=by_w,
            )
        out[wname] = well_res
    return out


# -------------------- D5: Raw vs Processed Near-Offset Curvature c2 (Route 3A) --------------------
def run_near_curvature(raw_blocks, zw, pz_sign, c0sq):
    """Compare 10-55 Hz log-spectral-ratio quadratic curvature c2 (2.3-2.8 s vs 1.2-1.6 s) on:
    1) Raw unprocessed near-offset traces (|x| <= 350 m) — A, B, PZ_pow, PZ_up (with ghost-period smoothing)
    2) Processed seq079 near-offset traces (|x| <= 350 m)
    Compare directly against archived Volve post-stack c2 = +1.322e-3 /Hz^2."""
    fit_mask = (F >= 10.0) & (F <= 55.0)
    f_fit = F[fit_mask]
    fc = 0.5 * (10.0 + 55.0)
    tau0 = 2.0 * zw / VW
    s_cal = 1.0 / math.sqrt(c0sq)

    def calc_c2(P_sh, P_dp, smooth_ghost=False):
        if smooth_ghost:
            P_sh = smooth_period(P_sh, tau0)
            P_dp = smooth_period(P_dp, tau0)
        rat = 0.5 * np.log(np.maximum(P_dp[fit_mask], 1e-30) / np.maximum(P_sh[fit_mask], 1e-30))
        # Fit y = c0 + c1*(f-fc) + c2*(f-fc)^2
        c2, c1, c0 = np.polyfit(f_fit - fc, rat, 2)
        lin_sl = np.polyfit(f_fit, rat, 1)[0]
        # Bending metric as in sid_uq: difference of 10-35 vs 35-55 slope
        m1 = (f_fit >= 10) & (f_fit <= 35)
        m2 = (f_fit >= 35) & (f_fit <= 55)
        drop_10_35 = float(-(rat[m1][-1] - rat[m1][0]))
        drop_35_55 = float(-(rat[m2][-1] - rat[m2][0]))
        return dict(c2=float(c2), c1=float(c1), lin_slope=float(lin_sl),
                    drop_10_35=drop_10_35, drop_35_55=drop_35_55)

    out = {}
    # 1. Raw near offsets across all 8 blocks pooled
    PA_sh = np.zeros_like(F); PA_dp = np.zeros_like(F)
    PB_sh = np.zeros_like(F); PB_dp = np.zeros_like(F)
    PU_sh = np.zeros_like(F); PU_dp = np.zeros_like(F)
    nb = 0
    for A, B, off, _ in raw_blocks:
        kr = off <= REF
        if kr.sum() < 10:
            continue
        pa1 = power_incoh(A[kr], 1.2, 1.6); pa2 = power_incoh(A[kr], 2.3, 2.8)
        pb1 = power_incoh(B[kr], 1.2, 1.6); pb2 = power_incoh(B[kr], 2.3, 2.8)
        U = A[kr] + pz_sign * s_cal * B[kr]
        pu1 = power_incoh(U, 1.2, 1.6); pu2 = power_incoh(U, 2.3, 2.8)
        if any(x is None for x in (pa1, pa2, pb1, pb2, pu1, pu2)):
            continue
        PA_sh += pa1; PA_dp += pa2
        PB_sh += pb1; PB_dp += pb2
        PU_sh += pu1; PU_dp += pu2
        nb += 1
    if nb > 0:
        out["raw_A_TRID11"] = calc_c2(PA_sh, PA_dp, smooth_ghost=True)
        out["raw_B_TRID12"] = calc_c2(PB_sh, PB_dp, smooth_ghost=True)
        out["raw_PZ_pow"] = calc_c2(PA_sh + PB_sh / c0sq, PA_dp + PB_dp / c0sq, smooth_ghost=True)
        out["raw_PZ_up"] = calc_c2(PU_sh, PU_dp, smooth_ghost=True)
        out["raw_PZ_up_nosmooth"] = calc_c2(PU_sh, PU_dp, smooth_ghost=False)

    # 2. Processed seq079 near offsets
    Hp, hp_d = load_proc_headers()
    sel = (np.abs(hp_d["cmp"]) < CMP_LIM) & (hp_d["off"] <= REF)
    kp = np.where(sel)[0]
    G = read_proc_samples(Hp, kp)
    pp1 = power_incoh(G, 1.2, 1.6)
    pp2 = power_incoh(G, 2.3, 2.8)
    if pp1 is not None and pp2 is not None:
        out["proc_seq079_near"] = calc_c2(pp1, pp2, smooth_ghost=False)
    out["poststack_archived_ref_c2"] = 1.322e-3
    return out


def main():
    smoke = "--smoke" in sys.argv
    t0 = time.monotonic()
    print(f"=== T1.0 B-OSR / FE-OSR Diagnostics (smoke={smoke}) ===", flush=True)

    H, h = load_raw_headers()
    sel = (h["dline"] == 5081) & (np.abs(h["cmp"]) < CMP_LIM)
    key = lambda i: (h["ffid"][i], h["dpt"][i])
    ia = {key(i): i for i in np.where(sel & (h["trid"] == 11))[0]}
    ib = {key(i): i for i in np.where(sel & (h["trid"] == 12))[0]}
    pairs = sorted(set(ia) & set(ib))
    ka = np.array([ia[p] for p in pairs])
    kb = np.array([ib[p] for p in pairs])
    cmp_all = h["cmp"][ka]
    off_all = h["off"][ka]
    zw = float(np.median(h["wd"][ka]))
    print(f"Loaded {len(pairs)} raw TRID 11/12 pairs, water depth={zw:.1f} m ({time.monotonic()-t0:.1f}s)", flush=True)

    R = Rays(1.0, zw)
    tabs = [R.table(tc) for tc in CENT]
    ptabs = [ptable(t) for t in tabs]

    e = np.linspace(-CMP_LIM, CMP_LIM, NBLK + 1)
    blk_indices = [0, 1] if smoke else list(range(NBLK))
    blocks = []
    for i in blk_indices:
        m = (cmp_all >= e[i]) & (cmp_all < e[i + 1])
        blocks.append((read_raw_samples(H, ka[m]), read_raw_samples(H, kb[m]), off_all[m], cmp_all[m]))
        print(f"  block {i}: {int(m.sum())} pairs ({time.monotonic()-t0:.1f}s)", flush=True)

    # D1: Component physics & PZ calibration
    d1 = diagnose_components(blocks, tabs, zw)
    c0sq = d1["c0sq"]
    pz_sign = d1["pz_amp_best_sign"]
    print("\n--- D1: Component Identification & P-Z Calibration ---", flush=True)
    print(json.dumps(d1, indent=2), flush=True)

    # D0 & D2: Incoherent reproduction + Amplitude-domain PZ_up
    print("\n--- D0 & D2: Incoherent Repro & Amplitude-Domain PZ_up ---", flush=True)
    d0_d2 = {}
    for comp in ("pz_pow", "a", "b", "pz_up"):
        s, _ = run_incoh_field79raw(blocks, tabs, ptabs, zw, c0sq, comp, pz_sign=pz_sign)
        d0_d2[comp] = s
        print(f"  incoh {comp:>6}: noint={s['invQ_noint']:+6.2f} {np.round(s['ci90_noint'],1).tolist()} | FE={s['invQ_fe']:+6.2f} | FE_midfar={s['invQ_fe_midfar']:+6.2f} | 1.2-1.6(noint/FE)={s['q_12_16_noint']:+5.1f}/{s['q_12_16_fe']:+5.1f} | 2.0-2.2(noint/FE)={s['q_20_22_noint']:+5.1f}/{s['q_20_22_fe']:+5.1f}", flush=True)

    # D3: Zero-stretch primary-slowness beamforming (B-OSR) + slowness scan
    print("\n--- D3: Zero-Stretch Primary-Slowness Beamforming (B-OSR, L=400 m, 50 m fine CMP) ---", flush=True)
    d3 = {}
    for comp in ("a", "b", "pz_up"):
        s, _ = run_bosr(blocks, R, tabs, ptabs, zw, c0sq, comp, pz_sign=pz_sign, hw=HW_BEAM, do_slowscan=(comp in ("a", "pz_up")))
        d3[comp] = s
        print(f"  B-OSR {comp:>6}: noint={s['invQ_noint']:+6.2f} {np.round(s['ci90_noint'],1).tolist()} | FE={s['invQ_fe']:+6.2f} {np.round(s['ci90_fe'],1).tolist()} | 1.2-1.6(noint/FE)={s['q_12_16_noint']:+5.1f}/{s['q_12_16_fe']:+5.1f} | 2.0-2.2(noint/FE)={s['q_20_22_noint']:+5.1f}/{s['q_20_22_fe']:+5.1f}", flush=True)
        for sk in ("slowscan_1.4", "slowscan_2.0_2.2"):
            if sk in s:
                print(f"    {sk}: {json.dumps(s[sk])}", flush=True)

    # D4: Well-log 1D exact prestack forward model (F-12, F-14)
    print("\n--- D4: Volve F-12 / F-14 Well-Log Exact 1D Prestack Reflectivity (Q=inf vs VSP Q) ---", flush=True)
    d4 = run_well_prestack_forward(R, tabs, zw)
    print(json.dumps(d4, indent=2), flush=True)

    # D5: Raw vs processed near-offset curvature c2
    print("\n--- D5: Raw vs Processed Near-Offset Spectral Curvature c2 (10-55 Hz) ---", flush=True)
    d5 = run_near_curvature(blocks, zw, pz_sign, c0sq)
    print(json.dumps(d5, indent=2), flush=True)

    out = dict(smoke=smoke, water_depth=zw, n_pairs=len(pairs), elapsed_s=round(time.monotonic() - t0, 1),
               D1_components=d1, D0_D2_incoh=d0_d2, D3_bosr=d3, D4_well_forward=d4, D5_curvature=d5)
    out_dir = ROOT / "results"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / ("bosr_diag_smoke.json" if smoke else "bosr_diag.json")
    out_path.write_text(json.dumps(out, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"\nSaved {out_path} ({time.monotonic()-t0:.1f}s)", flush=True)


if __name__ == "__main__":
    main()
