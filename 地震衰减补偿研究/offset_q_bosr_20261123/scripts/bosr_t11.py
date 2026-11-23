"""T1.1 Pre-registered Aperture-Invariance, Primary-Corridor Autofocus (B-OSR-Auto),
Processed Volve Control, and Brand-New Synthetic Self-Check (seed=900000).

Read-only access to D:/ChatGPT/研究 caches. Low RAM (< 1.5 GB), single process, 2 BLAS threads.
Usage: python scripts/bosr_t11.py [--smoke]
"""
import os, sys, json, math, time
os.environ.setdefault("OPENBLAS_NUM_THREADS", "2")
os.environ.setdefault("MKL_NUM_THREADS", "2")
os.environ.setdefault("OMP_NUM_THREADS", "2")

from pathlib import Path
import numpy as np

# Import shared constants & utilities from T1.0 bosr_diag.py
sys.path.insert(0, str(Path(__file__).resolve().parent))
import bosr_diag as bd

ROOT = bd.ROOT
DT, NT, T, NF, F = bd.DT, bd.NT, bd.T, bd.NF, bd.F
CENT, WL, REF, OFFS = bd.CENT, bd.WL, bd.REF, bd.OFFS
HW_INCOH, CMP_LIM, NBLK, FINE_PER_BLK = bd.HW_INCOH, bd.CMP_LIM, bd.NBLK, bd.FINE_PER_BLK
BAND, FLOOR, MINBW, STRETCH, MUTE, MUTE_PAD = bd.BAND, bd.FLOOR, bd.MINBW, bd.STRETCH, bd.MUTE, bd.MUTE_PAD
ZS, VW = bd.ZS, bd.VW

# Primary-corridor slowness scaling grid: alpha = p / p_prim in [0.84, 1.18] (+/- 8% velocity error)
# Strictly excludes water/shallow multiples (v ~ 1600 m/s -> alpha in [1.5, 2.7]!)
ALPHA_GRID = np.round(np.linspace(0.84, 1.18, 18), 4)
AL_BAND = (F >= 12.0) & (F <= 35.0)  # Low-to-mid primary energy band for autofocus (zero HF overfitting)
BAND_LO = (12.0, 30.0)
BAND_HI = (30.0, 55.0)


class RaysLocal(bd.Rays):
    """1-D ray tracer using local 2-D Volve vp.bin averaged over [x_lo, x_hi]."""
    def __init__(self, x_lo, x_hi, zw, vscale=1.0):
        m = json.loads(bd.META_JSON.read_text(encoding="utf-8"))
        nz, nx = m["nz"], m["nx"]
        V = np.fromfile(bd.VP_BIN, np.float32).reshape(nz, nx)
        x = m["xmin"] + (np.arange(nx) - m["npml"]) * m["dx"]
        z = (np.arange(nz) - m["npml"]) * m["dx"]
        k = z >= 0
        xm = (x >= x_lo) & (x <= x_hi)
        v = V[np.ix_(k, xm)].mean(axis=1).astype(float) * vscale
        self.dz = 1.0
        self.z = np.arange(0, z[k][-1], self.dz) + 0.5 * self.dz
        self.v = np.interp(self.z, z[k], v)
        self.T0 = 2.0 * np.cumsum(self.dz / self.v)
        self.zw = zw
        self.zA = self.z[np.searchsorted(self.T0, bd.TB[0])]
        self.zB = self.z[np.searchsorted(self.T0, bd.TB[1])]
        self.sed = self.z > zw
        self.lay = [(self.z > zw) & (self.z <= self.zA), (self.z > self.zA) & (self.z <= self.zB)]
        self.legs = (self.z > ZS).astype(float) + (self.z > zw).astype(float)


def compute_c0sq_allwin(blocks, tabs):
    """Exact field79raw.py c0sq across ALL windows tc in CENT (0.6-2.2 s), 10-60 Hz."""
    ib = (F >= 10.0) & (F <= 60.0)
    num = den = 0.0
    for A, B, off, _ in blocks:
        kr = off <= REF
        if kr.sum() < 10:
            continue
        xr = float(np.sqrt(np.mean(off[kr] ** 2)))
        for tc, tab in zip(CENT, tabs):
            tr = float(tab(np.array([xr]))[0][0])
            PA = bd.power_incoh(A[kr], tr - WL / 2, tr + WL / 2)
            PB = bd.power_incoh(B[kr], tr - WL / 2, tr + WL / 2)
            if PA is None or PB is None:
                continue
            num += float(PB[ib].sum())
            den += float(PA[ib].sum())
    return num / max(den, 1e-30)


def fit_subband_and_curv(P0, Px, s_full):
    """Compute full-band slope, low-band (12-30) slope, high-band (30-55) slope, and quadratic curvature c2."""
    rat_full = 0.5 * np.log(np.maximum(Px, 1e-30) / np.maximum(P0, 1e-30))
    f_full = F[s_full]
    y_full = float(-np.polyfit(f_full, rat_full[s_full], 1)[0])
    # Quadratic curvature c2 over the valid band s_full (only if bandwidth >= 20 Hz)
    if f_full[-1] - f_full[0] >= 20.0:
        fc = 0.5 * (f_full[0] + f_full[-1])
        c2 = float(np.polyfit(f_full - fc, rat_full[s_full], 2)[0])
    else:
        c2 = float("nan")
    # Sub-bands restricted inside s_full valid range
    fmax_v = float(f_full[-1])
    mlo = (F >= BAND_LO[0]) & (F <= min(BAND_LO[1], fmax_v))
    mhi = (F >= BAND_HI[0]) & (F <= min(BAND_HI[1], fmax_v))
    y_lo = float(-np.polyfit(F[mlo], rat_full[mlo], 1)[0]) if mlo.sum() >= 5 and (F[mlo][-1] - F[mlo][0] >= 10.0) else float("nan")
    y_hi = float(-np.polyfit(F[mhi], rat_full[mhi], 1)[0]) if mhi.sum() >= 5 and (F[mhi][-1] - F[mhi][0] >= 10.0) else float("nan")
    return y_full, y_lo, y_hi, c2


def fast_boot_ci(blocks_rows, mode="noint", seed=79, nboot=2000):
    """Exact algebraic sufficient-statistic block bootstrap (1000x faster, identical to bd.boot_blocks to 1e-15)."""
    nblk = len(blocks_rows)
    if nblk == 0:
        return [float("nan"), float("nan")]
    num = np.zeros(nblk, dtype=np.float64)
    den = np.zeros(nblk, dtype=np.float64)
    for bi, rows in enumerate(blocks_rows):
        if not rows:
            continue
        if mode == "noint":
            x = np.pi * np.array([r["dS"] for r in rows], dtype=np.float64)
            y = np.array([r["y"] for r in rows], dtype=np.float64)
            num[bi] = float(x @ y)
            den[bi] = float(x @ x)
        else:
            groups = {}
            for r in rows:
                groups.setdefault(r["tc"], []).append(r)
            for g in groups.values():
                if len(g) >= 2:
                    xg = np.pi * np.array([r["dS"] for r in g], dtype=np.float64)
                    yg = np.array([r["y"] for r in g], dtype=np.float64)
                    xd = xg - xg.mean()
                    yd = yg - yg.mean()
                    num[bi] += float(xd @ yd)
                    den[bi] += float(xd @ xd)
    if np.sum(den) <= 0:
        return [float("nan"), float("nan")]
    rng = np.random.default_rng(seed)
    picks = rng.integers(0, nblk, size=(nboot, nblk))
    s_num = num[picks].sum(axis=1)
    s_den = den[picks].sum(axis=1)
    vals = np.where(s_den > 0, 1000.0 * s_num / np.maximum(s_den, 1e-30), np.nan)
    if not np.any(np.isfinite(vals)):
        return [float("nan"), float("nan")]
    lo, hi = np.nanpercentile(vals, [5, 95])
    return [float(lo), float(hi)]


def summarise_t11(blocks_rows):
    """Fast extended summary including sub-bands (12-30 vs 30-55 Hz), curvature c2/dS, and alpha diagnostics."""
    allr = [r for b in blocks_rows for r in b]
    pt = bd.fit1_noint(allr)
    ci = fast_boot_ci(blocks_rows, "noint")
    fe_all = bd.fit_fe(allr, "dS", "y", exclude_near=False)
    fe_ci = fast_boot_ci(blocks_rows, "fe")
    g_mid = [r for r in allr if 1.15 <= r["tc"] <= 1.65]
    g_deep = [r for r in allr if 1.95 <= r["tc"] <= 2.25]
    base = dict(
        invQ_noint=pt, ci90_noint=ci,
        invQ_fe=fe_all, ci90_fe=fe_ci,
        invQ_fe_midfar=fe_all,
        n_rows=len(allr),
        q_12_16_noint=bd.fit1_noint(g_mid),
        q_12_16_fe=bd.fit_fe(g_mid, "dS", "y", False),
        q_20_22_noint=bd.fit1_noint(g_deep),
        q_20_22_fe=bd.fit_fe(g_deep, "dS", "y", False),
    )

    def fit_key(rows, yk):
        sub = [r for r in rows if np.isfinite(r.get(yk, float("nan")))]
        if len(sub) < 3:
            return float("nan"), float("nan")
        x = np.pi * np.array([r["dS"] for r in sub])
        y = np.array([r[yk] for r in sub])
        q_noint = float(1000.0 * (x @ y) / (x @ x))
        q_fe = bd.fit_fe(sub, "dS", yk, exclude_near=False)
        return q_noint, q_fe

    qlo_ni, qlo_fe = fit_key(allr, "y_lo")
    qhi_ni, qhi_fe = fit_key(allr, "y_hi")
    c2_norm = [r["c2"] / r["dS"] for r in allr if np.isfinite(r.get("c2", float("nan"))) and r["dS"] > 0.05]
    alphas = [r["alpha"] for r in allr if "alpha" in r and np.isfinite(r["alpha"])]
    base.update(dict(
        q_lo_12_30_noint=qlo_ni, q_lo_12_30_fe=qlo_fe,
        q_hi_30_55_noint=qhi_ni, q_hi_30_55_fe=qhi_fe,
        c2_over_dS_med=float(np.median(c2_norm)) if c2_norm else float("nan"),
        alpha_med=float(np.median(alphas)) if alphas else 1.0,
        alpha_p10=float(np.percentile(alphas, 10)) if alphas else 1.0,
        alpha_p90=float(np.percentile(alphas, 90)) if alphas else 1.0,
        v_eff_ratio_med=float(1.0 / math.sqrt(np.median(alphas))) if alphas else 1.0,
    ))
    return base


def extract_fine_cmp_pack(X, off, cmp_arr, fine_edges, tab, x0, hw):
    """Pre-extract padded trace packs for each fine CMP inside [x0-hw, x0+hw] so Fixed/Auto share work."""
    t_center = float(tab(np.array([x0]))[0][0])
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
    k_win = (t_sub >= t_center - WL / 2) & (t_sub < t_center + WL / 2)
    nw = int(k_win.sum())
    if nw < int(0.9 * WL / DT):
        return None
    wt = np.hanning(nw)
    dts_w = (np.arange(nw) - 0.5 * (nw - 1)) * DT
    wt2 = wt ** 2
    sig_t2 = float(np.sum(wt2 * dts_w ** 2) / np.sum(wt2))

    packs = []
    for fi in range(len(fine_edges) - 1):
        m = (cmp_arr >= fine_edges[fi]) & (cmp_arr < fine_edges[fi + 1]) & (off >= x0 - hw) & (off < x0 + hw)
        if m.sum() < 6:
            continue
        seg = X[m][:, k_pad]
        live = (np.abs(seg) > 0).mean(axis=1) > 0.90
        if live.sum() < 6:
            continue
        seg = seg[live]
        xo_i = off[m][live]
        xi = xo_i - x0
        dt0 = tab(xo_i)[0] - t_center
        if not np.all(np.isfinite(dt0)):
            continue
        w_ap = 0.5 * (1.0 + np.cos(np.pi * np.clip(xi / hw, -1.0, 1.0)))
        bins = np.clip(((xi + hw) / (2.0 * hw) * 8).astype(int), 0, 7)
        cnts = np.bincount(bins, minlength=8).astype(float)
        w_i = w_ap / np.maximum(cnts[bins], 1.0)
        ws = float(np.sum(w_i))
        if ws <= 0:
            continue
        w_i /= ws
        sig_xi2 = float(np.sum(w_i * xi ** 2))
        S_pad = np.fft.rfft(seg, n_fft_pad, axis=1)
        # Base alignment at alpha = 1.0 for ultra-fast residual-shift autofocus scoring
        ph0 = np.exp(2j * np.pi * f_pad[None, :] * dt0[:, None])
        al0 = np.fft.irfft(S_pad * ph0, n_fft_pad, axis=1)[:, :n_pad]
        win0 = al0[:, k_win] * wt[None, :]
        W0 = np.fft.rfft(win0, NF, axis=1)  # (ntr, nF)
        packs.append(dict(S_pad=S_pad, dt0=dt0, w_i=w_i, sig_xi2=sig_xi2, W0=W0, win0=win0))

    if not packs:
        return None
    return dict(packs=packs, n_pad=n_pad, n_fft_pad=n_fft_pad, f_pad=f_pad, k_win=k_win, wt=wt, sig_t2=sig_t2)


def score_alphas_on_pack(pack_obj):
    """Vectorized 12-35 Hz coherent beam energy across ALPHA_GRID using frequency-domain residual shift."""
    f_al = F[AL_BAND]  # (nf_al,)
    d_alphas = ALPHA_GRID - 1.0  # (na,)
    e_curve = np.zeros(len(ALPHA_GRID), dtype=np.float64)
    for p in pack_obj["packs"]:
        wW = p["w_i"][:, None] * p["W0"][:, AL_BAND]  # (ntr, nf_al)
        dt0 = p["dt0"]  # (ntr,)
        # Broadcast (na, ntr, nf_al): fast and tiny (~18 x 40 x 15 complex128 = 170 KB)
        ang = (2.0 * np.pi) * (d_alphas[:, None, None] * dt0[None, :, None] * f_al[None, None, :])
        B = np.sum(wW[None, :, :] * (np.cos(ang) + 1j * np.sin(ang)), axis=1)  # (na, nf_al)
        e_curve += np.sum(B.real ** 2 + B.imag ** 2, axis=1)
    return e_curve / len(pack_obj["packs"])


def eval_beam_at_alpha(pack_obj, alpha, kappa_p):
    """Exact padded rigid-shift coherent beam power at chosen alpha, with analytical intra-window smear correction."""
    P_beam_acc = np.zeros_like(F)
    P_incoh_acc = np.zeros_like(F)
    sig_xi2_list = []
    f_pad = pack_obj["f_pad"]
    n_fft_pad = pack_obj["n_fft_pad"]
    n_pad = pack_obj["n_pad"]
    k_win = pack_obj["k_win"]
    wt = pack_obj["wt"]
    for p in pack_obj["packs"]:
        if abs(alpha - 1.0) < 1e-12:
            win_tr = p["win0"]
            W_tr = p["W0"]
        else:
            dt_i = alpha * p["dt0"]
            ph = np.exp(2j * np.pi * f_pad[None, :] * dt_i[:, None])
            al = np.fft.irfft(p["S_pad"] * ph, n_fft_pad, axis=1)[:, :n_pad]
            win_tr = al[:, k_win] * wt[None, :]
            W_tr = np.fft.rfft(win_tr, NF, axis=1)
        beam_W = np.sum(p["w_i"][:, None] * W_tr, axis=0)
        P_beam_acc += beam_W.real ** 2 + beam_W.imag ** 2
        P_incoh_acc += np.sum(p["w_i"][:, None] * (W_tr.real ** 2 + W_tr.imag ** 2), axis=0)
        sig_xi2_list.append(p["sig_xi2"])
    nf = len(pack_obj["packs"])
    P_beam = P_beam_acc / nf
    P_incoh = P_incoh_acc / nf
    sig_xi2 = float(np.mean(sig_xi2_list))
    kp_eff = alpha * kappa_p
    B2 = np.exp(-((2.0 * np.pi * F) ** 2) * (kp_eff ** 2) * pack_obj["sig_t2"] * sig_xi2)
    B2 = np.maximum(B2, 0.25)
    P_corr = P_beam / B2
    ib = (F >= BAND[0]) & (F <= BAND[1])
    semb = float(np.median(P_beam[ib] / np.maximum(P_incoh[ib], 1e-30)))
    return P_corr, P_incoh, semb


def run_suite_on_dataset(blocks, R_global, R_blocks, zw, c0sq, comp, pz_sign, hws,
                         do_ghost_smooth=True, run_incoh=False):
    """Run Fixed1D, Local2D, AutoBlock, AutoGlobal across apertures hws on a dataset."""
    out = {}
    s_cal = 1.0 / math.sqrt(max(c0sq, 1e-30))
    blk_edges = np.linspace(-CMP_LIM, CMP_LIM, NBLK + 1)

    # Precompute ray tables & kappa_p for global 1D and each block local 2D
    def make_ray_bundle(R_obj):
        tabs = [R_obj.table(tc) for tc in CENT]
        ptabs = [bd.ptable(t) for t in tabs]
        tp = [bd.ptable(R_obj.table(tc + 0.05)) for tc in CENT]
        tm = [bd.ptable(R_obj.table(max(tc - 0.05, 0.25))) for tc in CENT]
        return tabs, ptabs, tp, tm

    glob_bundle = make_ray_bundle(R_global)
    loc_bundles = [make_ray_bundle(Rb) for Rb in R_blocks]

    def get_sig(A, B, mask, cw):
        if comp == "pz_up":
            return A[mask] + pz_sign * (s_cal / cw) * B[mask]
        elif comp == "a":
            return A[mask]
        else:
            return B[mask]

    # Optional incoherent baseline (for processed seq079 and synthetics)
    if run_incoh:
        tabs, ptabs, _, _ = glob_bundle
        blk_rows = []
        for bi, (A, B, off, _) in enumerate(blocks):
            rows = []
            kr = off <= REF
            if kr.sum() < 8:
                blk_rows.append(rows); continue
            xr = float(np.sqrt(np.mean(off[kr] ** 2)))
            for tc, tab, ptab in zip(CENT, tabs, ptabs):
                r0 = tab(np.array([xr]))
                tr = float(r0[0][0])
                cwr, taur = bd.ghost(ptab(xr)[0], zw)
                P0_r = bd.power_incoh(get_sig(A, B, kr, cwr), tr - WL / 2, tr + WL / 2)
                if P0_r is None:
                    continue
                P0 = bd.smooth_period(P0_r, taur) if do_ghost_smooth else P0_r
                for o in OFFS:
                    k = (off >= o - HW_INCOH) & (off < o + HW_INCOH)
                    if k.sum() < 8:
                        continue
                    xo = float(np.sqrt(np.mean(off[k] ** 2)))
                    tx, la, lb, ls = (float(v[0]) for v in tab(np.array([xo])))
                    if not np.isfinite(tx) or tx / tc > STRETCH or tx - WL / 2 < MUTE(o + HW_INCOH) + MUTE_PAD or tx + WL / 2 > T[-1]:
                        continue
                    cw, tau = bd.ghost(ptab(xo)[0], zw)
                    Px_r = bd.power_incoh(get_sig(A, B, k, cw), tx - WL / 2, tx + WL / 2)
                    if Px_r is None:
                        continue
                    Px = bd.smooth_period(Px_r, tau) if do_ghost_smooth else Px_r
                    s = bd.band_slice(P0, Px)
                    if s is None:
                        continue
                    yf, ylo, yhi, c2 = fit_subband_and_curv(P0, Px, s)
                    rows.append(dict(block=bi, tc=float(tc), off=xo, y=yf, y_lo=ylo, y_hi=yhi, c2=c2,
                                     dS=ls - float(r0[3][0]), tau_sed=ls, alpha=1.0))
            blk_rows.append(rows)
        out["Incoh"] = summarise_t11(blk_rows)

    # Beamforming suite across apertures
    for hw in hws:
        # 1) Fixed1D (using global 1D R_global, alpha=1.0)
        tabs_g, ptabs_g, tp_g, tm_g = glob_bundle
        rows_f1d = []
        for bi, (A, B, off, cmp_arr) in enumerate(blocks):
            rows = []
            fine_edges = np.linspace(blk_edges[bi], blk_edges[bi + 1], FINE_PER_BLK + 1)
            kr = off <= REF
            if kr.sum() < 8:
                rows_f1d.append(rows); continue
            xr = float(np.sqrt(np.mean(off[kr] ** 2)))
            for wi, (tc, tab, ptab) in enumerate(zip(CENT, tabs_g, ptabs_g)):
                r0 = tab(np.array([xr]))
                cwr, taur = bd.ghost(ptab(xr)[0], zw)
                pk0 = extract_fine_cmp_pack(get_sig(A, B, kr, cwr), off[kr], cmp_arr[kr], fine_edges, tab, xr, REF)
                if pk0 is None:
                    continue
                kp0 = float((tp_g[wi](xr)[0] - tm_g[wi](xr)[0]) / (0.05 + min(0.05, tc - 0.25)))
                P0_b, _, _ = eval_beam_at_alpha(pk0, 1.0, kp0)
                P0 = bd.smooth_period(P0_b, taur) if do_ghost_smooth else P0_b
                for o in OFFS:
                    tx0, _, _, _ = (float(v[0]) for v in tab(np.array([float(o)])))
                    if not np.isfinite(tx0) or tx0 / tc > STRETCH or tx0 - WL / 2 < MUTE(o + HW_INCOH) + MUTE_PAD or tx0 + WL / 2 > T[-1]:
                        continue
                    kb = (off >= o - hw) & (off < o + hw)
                    if kb.sum() < 12:
                        continue
                    xo = float(np.sqrt(np.mean(off[kb] ** 2)))
                    tx, la, lb, ls = (float(v[0]) for v in tab(np.array([xo])))
                    cw, tau = bd.ghost(ptab(xo)[0], zw)
                    pkx = extract_fine_cmp_pack(get_sig(A, B, kb, cw), off[kb], cmp_arr[kb], fine_edges, tab, xo, hw)
                    if pkx is None:
                        continue
                    kp = float((tp_g[wi](xo)[0] - tm_g[wi](xo)[0]) / (0.05 + min(0.05, tc - 0.25)))
                    Px_b, _, sembx = eval_beam_at_alpha(pkx, 1.0, kp)
                    Px = bd.smooth_period(Px_b, tau) if do_ghost_smooth else Px_b
                    s = bd.band_slice(P0, Px)
                    if s is None:
                        continue
                    yf, ylo, yhi, c2 = fit_subband_and_curv(P0, Px, s)
                    rows.append(dict(block=bi, tc=float(tc), off=xo, y=yf, y_lo=ylo, y_hi=yhi, c2=c2,
                                     dS=ls - float(r0[3][0]), tau_sed=ls, semb=sembx, alpha=1.0))
            rows_f1d.append(rows)

        # 2) Local2D + AutoBlock + AutoGlobal (share extracted fine-CMP packs!)
        # Store lightweight packs per (bi, wi, o) to compute AutoGlobal across all blocks cleanly
        rows_loc = [[] for _ in range(len(blocks))]
        rows_ab = [[] for _ in range(len(blocks))]
        rows_ag = [[] for _ in range(len(blocks))]
        # First pass: collect per-(wi, o) pooled energy curves for AutoGlobal, and compute Local2D & AutoBlock on the fly
        glob_e = {}  # (wi, o) -> energy curve sum across blocks
        cache_bins = []  # store only needed scalars/packs if doing AutoGlobal
        # To keep RAM minimal (< 300 MB), store P0_loc, P0_ab and per-(bi, wi, o) score curves + re-eval or hold compact packs?
        # Even simpler: run per-window wi (only 9 windows!), so RAM holds 1 window at a time (< 40 MB)!
        for wi, tc in enumerate(CENT):
            # Extract near-offset packs for all blocks for window wi
            near_info = {}
            far_info = {}  # (bi, o) -> dict(pk, xo, ls, ls0, kp, tau, e_curve)
            e_pool = {}    # o -> sum of e_curve over blocks
            for bi, (A, B, off, cmp_arr) in enumerate(blocks):
                tabs_b, ptabs_b, tp_b, tm_b = loc_bundles[bi]
                tab, ptab = tabs_b[wi], ptabs_b[wi]
                fine_edges = np.linspace(blk_edges[bi], blk_edges[bi + 1], FINE_PER_BLK + 1)
                kr = off <= REF
                if kr.sum() < 8:
                    continue
                xr = float(np.sqrt(np.mean(off[kr] ** 2)))
                r0 = tab(np.array([xr]))
                cwr, taur = bd.ghost(ptab(xr)[0], zw)
                pk0 = extract_fine_cmp_pack(get_sig(A, B, kr, cwr), off[kr], cmp_arr[kr], fine_edges, tab, xr, REF)
                if pk0 is None:
                    continue
                kp0 = float((tp_b[wi](xr)[0] - tm_b[wi](xr)[0]) / (0.05 + min(0.05, tc - 0.25)))
                ec0 = score_alphas_on_pack(pk0)
                a0_ab = float(ALPHA_GRID[int(np.argmax(ec0))])
                P0_loc_b, _, _ = eval_beam_at_alpha(pk0, 1.0, kp0)
                P0_ab_b, _, _ = eval_beam_at_alpha(pk0, a0_ab, kp0)
                P0_loc = bd.smooth_period(P0_loc_b, taur) if do_ghost_smooth else P0_loc_b
                P0_ab = bd.smooth_period(P0_ab_b, taur) if do_ghost_smooth else P0_ab_b
                near_info[bi] = dict(pk0=pk0, kp0=kp0, taur=taur, ls0=float(r0[3][0]),
                                     P0_loc=P0_loc, P0_ab=P0_ab, ec0=ec0)
                e_pool["NEAR"] = e_pool.get("NEAR", np.zeros_like(ec0)) + ec0

                for o in OFFS:
                    tx0, _, _, _ = (float(v[0]) for v in tab(np.array([float(o)])))
                    if not np.isfinite(tx0) or tx0 / tc > STRETCH or tx0 - WL / 2 < MUTE(o + HW_INCOH) + MUTE_PAD or tx0 + WL / 2 > T[-1]:
                        continue
                    kb = (off >= o - hw) & (off < o + hw)
                    if kb.sum() < 12:
                        continue
                    xo = float(np.sqrt(np.mean(off[kb] ** 2)))
                    tx, la, lb, ls = (float(v[0]) for v in tab(np.array([xo])))
                    cw, tau = bd.ghost(ptab(xo)[0], zw)
                    pkx = extract_fine_cmp_pack(get_sig(A, B, kb, cw), off[kb], cmp_arr[kb], fine_edges, tab, xo, hw)
                    if pkx is None:
                        continue
                    kp = float((tp_b[wi](xo)[0] - tm_b[wi](xo)[0]) / (0.05 + min(0.05, tc - 0.25)))
                    ecx = score_alphas_on_pack(pkx)
                    ax_ab = float(ALPHA_GRID[int(np.argmax(ecx))])
                    e_pool[int(o)] = e_pool.get(int(o), np.zeros_like(ecx)) + ecx

                    # Local2D (alpha=1.0)
                    Px_loc_b, _, semb_loc = eval_beam_at_alpha(pkx, 1.0, kp)
                    Px_loc = bd.smooth_period(Px_loc_b, tau) if do_ghost_smooth else Px_loc_b
                    s = bd.band_slice(P0_loc, Px_loc)
                    if s is not None:
                        yf, ylo, yhi, c2 = fit_subband_and_curv(P0_loc, Px_loc, s)
                        rows_loc[bi].append(dict(block=bi, tc=float(tc), off=xo, y=yf, y_lo=ylo, y_hi=yhi, c2=c2,
                                                 dS=ls - float(r0[3][0]), tau_sed=ls, semb=semb_loc, alpha=1.0))

                    # AutoBlock (block-optimal alpha in primary corridor)
                    Px_ab_b, _, semb_ab = eval_beam_at_alpha(pkx, ax_ab, kp)
                    Px_ab = bd.smooth_period(Px_ab_b, tau) if do_ghost_smooth else Px_ab_b
                    s = bd.band_slice(P0_ab, Px_ab)
                    if s is not None:
                        yf, ylo, yhi, c2 = fit_subband_and_curv(P0_ab, Px_ab, s)
                        rows_ab[bi].append(dict(block=bi, tc=float(tc), off=xo, y=yf, y_lo=ylo, y_hi=yhi, c2=c2,
                                                dS=ls - float(r0[3][0]), tau_sed=ls, semb=semb_ab, alpha=ax_ab))

                    far_info[(bi, int(o))] = dict(pkx=pkx, xo=xo, ls=ls, kp=kp, tau=tau)

            # AutoGlobal: 1 shared alpha per (tc, offset) pooled across ALL blocks!
            if "NEAR" in e_pool:
                a0_ag = float(ALPHA_GRID[int(np.argmax(e_pool["NEAR"]))])
                P0_ag_map = {}
                for bi, ni in near_info.items():
                    P0_b, _, _ = eval_beam_at_alpha(ni["pk0"], a0_ag, ni["kp0"])
                    P0_ag_map[bi] = bd.smooth_period(P0_b, ni["taur"]) if do_ghost_smooth else P0_b
                for (bi, oi), fi_d in far_info.items():
                    if bi not in P0_ag_map:
                        continue
                    ax_ag = float(ALPHA_GRID[int(np.argmax(e_pool[oi]))])
                    Px_b, _, semb_ag = eval_beam_at_alpha(fi_d["pkx"], ax_ag, fi_d["kp"])
                    Px_ag = bd.smooth_period(Px_b, fi_d["tau"]) if do_ghost_smooth else Px_b
                    s = bd.band_slice(P0_ag_map[bi], Px_ag)
                    if s is not None:
                        yf, ylo, yhi, c2 = fit_subband_and_curv(P0_ag_map[bi], Px_ag, s)
                        rows_ag[bi].append(dict(block=bi, tc=float(tc), off=fi_d["xo"], y=yf, y_lo=ylo, y_hi=yhi, c2=c2,
                                                dS=fi_d["ls"] - near_info[bi]["ls0"], tau_sed=fi_d["ls"],
                                                semb=semb_ag, alpha=ax_ag))

        out[f"H{int(hw)}"] = dict(
            Fixed1D=summarise_t11(rows_f1d),
            Local2D=summarise_t11(rows_loc),
            AutoBlock=summarise_t11(rows_ab),
            AutoGlobal=summarise_t11(rows_ag),
        )
    return out


# -------------------- Part D: Brand-New Synthetic Self-Check (seed=900000) --------------------
def synthesize_volve_like(R_model, zw, q_vals, seed=900000, n_blocks=NBLK):
    """
    Generate realistic 8-block raw synthetic CMP gathers with:
    - True velocity error relative to R_model: +2.5% systematic + U(-2%, +2%) per block
    - Causal constant-Q attenuation (amplitude + dispersion) along exact sediment ray path
    - 6 m source ghost + 91 m seabed hydrophone ghost (P)
    - Crossing blue-spectrum water/shallow peg-leg multiples (v in [1520, 1680] m/s, Q_m = 500)
    - Near-offset low-frequency Scholte/direct contamination (x <= 350 m, 8-18 Hz)
    - 18 dB broadband white noise
    Returns (blocks, true_invQ_eff, v_err_blocks).
    """
    rng = np.random.default_rng(seed)
    fr = np.fft.rfftfreq(2 * NT, DT)
    f_ref = 25.0
    W_src = (fr / 25.0) ** 2 * np.exp(-(fr / 25.0) ** 2)
    # Near-offset Scholte low-freq filter (8-18 Hz)
    W_sch = np.exp(-0.5 * ((fr - 12.0) / 3.5) ** 2)

    blk_edges = np.linspace(-CMP_LIM, CMP_LIM, NBLK + 1)
    # Dense zero-offset traveltime grid for fast interpolation of exact rays per block
    tc_grid = np.linspace(0.35, 2.95, 27)
    x_grid = np.linspace(0.0, 3800.0, 96)

    # Precompute analytical true effective 1000/Q under R_model acceptance mask
    true_num = true_den = 0.0
    tabs_m = [R_model.table(tc) for tc in CENT]
    for tc, tab in zip(CENT, tabs_m):
        r0 = tab(np.array([200.0]))
        t0_a, t0_b, t0_s = float(r0[1][0]), float(r0[2][0]), float(r0[3][0])
        for o in OFFS:
            tx, la, lb, ls = (float(v[0]) for v in tab(np.array([float(o)])))
            if not np.isfinite(tx) or tx / tc > STRETCH or tx - WL / 2 < MUTE(o + HW_INCOH) + MUTE_PAD or tx + WL / 2 > T[-1]:
                continue
            dS = ls - t0_s
            dts = (la - t0_a) / q_vals[0] + (lb - t0_b) / q_vals[1] + max((ls - la - lb) - (t0_s - t0_a - t0_b), 0.0) / q_vals[2]
            true_num += dS * dts
            true_den += dS * dS
    true_invQ = 1000.0 * true_num / max(true_den, 1e-30)

    blocks = []
    v_errs = []
    # Realistic 20 m offset spacing per fine CMP (180 traces spanning 60..3640 m -> 10-20 traces per H100/H200 bin!)
    base_offs = np.linspace(60.0, 3640.0, 180)

    for bi in range(n_blocks):
        v_err = float(1.025 + rng.uniform(-0.020, 0.020))
        v_errs.append(v_err)
        R_true = bd.Rays(v_err, zw)
        # Precompute 2D lookup tables (tc_grid x x_grid) for t, la, lb, ls, p
        T_lut = np.zeros((len(tc_grid), len(x_grid)))
        LA_lut = np.zeros_like(T_lut)
        LB_lut = np.zeros_like(T_lut)
        LS_lut = np.zeros_like(T_lut)
        P_lut = np.zeros_like(T_lut)
        for ti, tcg in enumerate(tc_grid):
            tb = R_true.table(tcg)
            ptb = bd.ptable(tb)
            t_v, la_v, lb_v, ls_v = tb(x_grid)
            T_lut[ti] = t_v; LA_lut[ti] = la_v; LB_lut[ti] = lb_v; LS_lut[ti] = ls_v
            P_lut[ti] = ptb(x_grid)

        fine_edges = np.linspace(blk_edges[bi], blk_edges[bi + 1], FINE_PER_BLK + 1)
        A_list, off_list, cmp_list = [], [], []
        for fi in range(FINE_PER_BLK):
            offs = np.clip(base_offs + rng.uniform(-6.0, 6.0, size=len(base_offs)), 40.0, 3750.0)
            cmps = rng.uniform(fine_edges[fi] + 1.0, fine_edges[fi + 1] - 1.0, size=len(offs))
            ntr = len(offs)
            spec = np.zeros((ntr, len(fr)), dtype=np.complex128)

            # 1) 42 primary reflectors in [0.45, 2.75] s
            t0_refs = np.linspace(0.48, 2.72, 42) + rng.uniform(-0.02, 0.02, size=42)
            amps_ref = rng.normal(0.0, 1.0, size=42) * np.exp(-0.20 * (t0_refs - 0.5))
            for t0_r, amp_r in zip(t0_refs, amps_ref):
                idx = int(np.clip(np.searchsorted(tc_grid, t0_r) - 1, 0, len(tc_grid) - 2))
                w_t = (t0_r - tc_grid[idx]) / (tc_grid[idx + 1] - tc_grid[idx])
                t_row = (1 - w_t) * T_lut[idx] + w_t * T_lut[idx + 1]
                la_row = (1 - w_t) * LA_lut[idx] + w_t * LA_lut[idx + 1]
                lb_row = (1 - w_t) * LB_lut[idx] + w_t * LB_lut[idx + 1]
                ls_row = (1 - w_t) * LS_lut[idx] + w_t * LS_lut[idx + 1]
                p_row = (1 - w_t) * P_lut[idx] + w_t * P_lut[idx + 1]

                t_x = np.interp(offs, x_grid, t_row)
                la_x = np.interp(offs, x_grid, la_row)
                lb_x = np.interp(offs, x_grid, lb_row)
                ls_x = np.interp(offs, x_grid, ls_row)
                p_x = np.interp(offs, x_grid, p_row)
                ok = np.isfinite(t_x) & (t_x < T[-1] + 0.2)
                if not np.any(ok):
                    continue
                tstar = la_x[ok] / q_vals[0] + lb_x[ok] / q_vals[1] + np.maximum(ls_x[ok] - la_x[ok] - lb_x[ok], 0.0) / q_vals[2]
                cw_x, tau_r = bd.ghost(p_x[ok], zw)
                tau_s = 2.0 * ZS * cw_x / VW
                lnf = np.log(np.maximum(fr, 1.0) / f_ref)
                att = np.exp((-np.pi * fr[None, :] + 2j * fr[None, :] * lnf[None, :]) * tstar[:, None])
                g_s = 1.0 - np.exp(-2j * np.pi * fr[None, :] * tau_s[:, None])
                g_r = 1.0 - np.exp(-2j * np.pi * fr[None, :] * tau_r[:, None])
                ph = np.exp(-2j * np.pi * fr[None, :] * t_x[ok, None])
                spec[ok] += amp_r * W_src[None, :] * g_s * g_r * att * ph

            # 2) 22 crossing blue-spectrum water/shallow peg-leg multiples (v in 1520..1680 m/s, Q_m = 500)
            t0_m = np.linspace(0.75, 2.55, 22) + rng.uniform(-0.03, 0.03, size=22)
            v_m = rng.uniform(1520.0, 1680.0, size=22)
            amp_m = rng.normal(0.0, 0.55, size=22) * (0.5 + 0.6 * (t0_m / 2.2))
            for tm0, vm, am in zip(t0_m, v_m, amp_m):
                t_mx = np.sqrt(tm0 ** 2 + (offs / vm) ** 2)
                p_mx = offs / (vm ** 2 * t_mx)
                cw_m, tau_rm = bd.ghost(p_mx, zw)
                tau_sm = 2.0 * ZS * cw_m / VW
                tstar_m = t_mx / 500.0  # mostly water / very shallow -> high Q (bluer than deep primaries!)
                att_m = np.exp(-np.pi * fr[None, :] * tstar_m[:, None])
                g_sm = 1.0 - np.exp(-2j * np.pi * fr[None, :] * tau_sm[:, None])
                g_rm = 1.0 - np.exp(-2j * np.pi * fr[None, :] * tau_rm[:, None])
                ph_m = np.exp(-2j * np.pi * fr[None, :] * t_mx[:, None])
                spec += am * W_src[None, :] * g_sm * g_rm * att_m * ph_m

            tr = np.fft.irfft(spec, 2 * NT, axis=1)[:, :NT].astype(np.float32)

            # 3) Near-offset Scholte/direct low-frequency contamination (x <= 350 m)
            kn = offs <= REF
            if np.any(kn):
                rms_n = float(np.sqrt(np.mean(tr[kn] ** 2))) + 1e-12
                raw_sch = rng.normal(0.0, 1.0, size=(int(kn.sum()), 2 * NT))
                sch = np.fft.irfft(np.fft.rfft(raw_sch, axis=1) * W_sch[None, :], 2 * NT, axis=1)[:, :NT]
                sch *= (0.45 * rms_n / (float(np.sqrt(np.mean(sch ** 2))) + 1e-12))
                tr[kn] += sch.astype(np.float32)

            # 4) Broadband white noise (~18 dB SNR) + exact hp filter
            rms_all = float(np.sqrt(np.mean(tr ** 2))) + 1e-12
            tr += rng.normal(0.0, 0.10 * rms_all, size=tr.shape).astype(np.float32)
            tr = bd.hp(tr)
            A_list.append(tr); off_list.append(offs); cmp_list.append(cmps)

        A_blk = np.vstack(A_list)
        blocks.append((A_blk, np.zeros_like(A_blk), np.concatenate(off_list), np.concatenate(cmp_list)))
    return blocks, float(true_invQ), v_errs


def main():
    smoke = "--smoke" in sys.argv
    t0 = time.monotonic()
    print(f"=== T1.1 Aperture Invariance, B-OSR-Auto & Synthetic Self-Check (smoke={smoke}) ===", flush=True)

    # Load raw Volve headers & blocks
    H, h = bd.load_raw_headers()
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
    e = np.linspace(-CMP_LIM, CMP_LIM, NBLK + 1)
    blk_indices = [0, 1] if smoke else list(range(NBLK))
    raw_blocks = []
    for i in blk_indices:
        m = (cmp_all >= e[i]) & (cmp_all < e[i + 1])
        raw_blocks.append((bd.read_raw_samples(H, ka[m]), bd.read_raw_samples(H, kb[m]), off_all[m], cmp_all[m]))
    print(f"Loaded raw Volve ({len(raw_blocks)} blocks, zw={zw:.1f} m) in {time.monotonic()-t0:.1f}s", flush=True)

    R_glob = bd.Rays(1.0, zw)
    R_loc_blocks = [RaysLocal(e[i], e[i + 1], zw, 1.0) for i in blk_indices]
    tabs_g = [R_glob.table(tc) for tc in CENT]
    ptabs_g = [bd.ptable(t) for t in tabs_g]

    # Part A: Exact all-window c0sq & M0 reproduction
    c0sq_all = compute_c0sq_allwin(raw_blocks, tabs_g)
    d1_tmp = bd.diagnose_components(raw_blocks, tabs_g, zw)
    pz_sign = d1_tmp["pz_amp_best_sign"]
    print(f"\n--- Part A: Exact all-window c0sq={c0sq_all:.4f} (vs field79raw 2.114), pz_sign={pz_sign} ---", flush=True)
    part_a = dict(c0sq_allwin=c0sq_all, pz_sign=pz_sign)
    for comp in ("pz_pow", "a", "b"):
        s, _ = bd.run_incoh_field79raw(raw_blocks, tabs_g, ptabs_g, zw, c0sq_all, comp, pz_sign=pz_sign)
        part_a[f"incoh_{comp}"] = s
        print(f"  incoh {comp:>6}: noint={s['invQ_noint']:+6.2f} {np.round(s['ci90_noint'],1).tolist()} | FE={s['invQ_fe']:+6.2f}", flush=True)

    # Part B: Processed seq079_dline5081.segy control
    print(f"\n--- Part B: Processed Volve seq079_dline5081 Control ({time.monotonic()-t0:.1f}s) ---", flush=True)
    Hp, hp_d = bd.load_proc_headers()
    proc_blocks = []
    for i in blk_indices:
        m = (hp_d["cmp"] >= e[i]) & (hp_d["cmp"] < e[i + 1])
        idx = np.where(m)[0]
        G = bd.read_proc_samples(Hp, idx)
        proc_blocks.append((G, np.zeros_like(G), hp_d["off"][idx], hp_d["cmp"][idx]))
    hws_proc = [100.0, 150.0, 200.0] if not smoke else [100.0, 200.0]
    part_b = run_suite_on_dataset(proc_blocks, R_glob, R_loc_blocks, zw, 1.0, "a", 1, hws_proc,
                                  do_ghost_smooth=False, run_incoh=True)
    print(f"  Proc Incoh: noint={part_b['Incoh']['invQ_noint']:+6.2f} | FE={part_b['Incoh']['invQ_fe']:+6.2f}", flush=True)
    for hw in hws_proc:
        hk = f"H{int(hw)}"
        for mname in ("Fixed1D", "Local2D", "AutoBlock", "AutoGlobal"):
            r = part_b[hk][mname]
            print(f"  Proc {hk} {mname:>10}: noint={r['invQ_noint']:+6.2f} | FE={r['invQ_fe']:+6.2f} {np.round(r['ci90_fe'],1).tolist()} | lo/hi(FE)={r['q_lo_12_30_fe']:+5.1f}/{r['q_hi_30_55_fe']:+5.1f} | c2/dS={r['c2_over_dS_med']:+.2e} | a_med={r['alpha_med']:.3f}", flush=True)
    del proc_blocks

    # Part C: Raw Volve aperture invariance & autofocus (a, b, pz_up)
    print(f"\n--- Part C: Raw Volve Aperture Invariance & B-OSR-Auto ({time.monotonic()-t0:.1f}s) ---", flush=True)
    hws_raw = [75.0, 100.0, 150.0, 200.0] if not smoke else [100.0, 200.0]
    part_c = {}
    for comp in ("a", "b", "pz_up"):
        print(f"  [Raw comp={comp}]", flush=True)
        res_c = run_suite_on_dataset(raw_blocks, R_glob, R_loc_blocks, zw, c0sq_all, comp, pz_sign, hws_raw,
                                     do_ghost_smooth=True, run_incoh=False)
        part_c[comp] = res_c
        for hw in hws_raw:
            hk = f"H{int(hw)}"
            for mname in ("Fixed1D", "Local2D", "AutoBlock", "AutoGlobal"):
                r = res_c[hk][mname]
                print(f"    {hk} {mname:>10}: noint={r['invQ_noint']:+6.2f} | FE={r['invQ_fe']:+6.2f} {np.round(r['ci90_fe'],1).tolist()} | midfar={r['invQ_fe_midfar']:+5.1f} | lo/hi(FE)={r['q_lo_12_30_fe']:+5.1f}/{r['q_hi_30_55_fe']:+5.1f} | c2/dS={r['c2_over_dS_med']:+.2e} | a_med={r['alpha_med']:.3f}", flush=True)
    del raw_blocks

    # Part D: Brand-New Synthetic Self-Check (seed=900000)
    print(f"\n--- Part D: Synthetic Self-Check (seed=900000, +2.5% +/- 2% v_err + multiples + ghosts + near noise) ({time.monotonic()-t0:.1f}s) ---", flush=True)
    hws_syn = [100.0, 150.0, 200.0] if not smoke else [100.0, 200.0]
    n_syn_blk = 2 if smoke else NBLK
    part_d = {}
    for sname, qv, sd in [("S0_Qinf", (1e9, 1e9, 1e9), 900000), ("S1_QVSP", (60.0, 90.0, 100.0), 900001)]:
        syn_blocks, true_iq, v_errs = synthesize_volve_like(R_glob, zw, qv, seed=sd, n_blocks=n_syn_blk)
        R_syn_loc = [R_glob for _ in range(n_syn_blk)]  # Processing only knows wrong model R_glob!
        res_s = run_suite_on_dataset(syn_blocks, R_glob, R_syn_loc, zw, 1.0, "a", 1, hws_syn,
                                     do_ghost_smooth=True, run_incoh=True)
        res_s["true_invQ_eff"] = true_iq
        res_s["v_err_blocks"] = v_errs
        part_d[sname] = res_s
        print(f"  [{sname}] true 1000/Q={true_iq:+.2f}, v_errs={np.round(v_errs, 3).tolist()}", flush=True)
        print(f"    Incoh: noint={res_s['Incoh']['invQ_noint']:+6.2f} | FE={res_s['Incoh']['invQ_fe']:+6.2f}", flush=True)
        for hw in hws_syn:
            hk = f"H{int(hw)}"
            for mname in ("Fixed1D", "AutoBlock", "AutoGlobal"):
                r = res_s[hk][mname]
                print(f"    {hk} {mname:>10}: noint={r['invQ_noint']:+6.2f} | FE={r['invQ_fe']:+6.2f} {np.round(r['ci90_fe'],1).tolist()} | lo/hi(FE)={r['q_lo_12_30_fe']:+5.1f}/{r['q_hi_30_55_fe']:+5.1f} | c2/dS={r['c2_over_dS_med']:+.2e} | a_med={r['alpha_med']:.3f}", flush=True)

    out = dict(
        smoke=smoke, water_depth=zw, elapsed_s=round(time.monotonic() - t0, 1),
        PartA_M0_repro=part_a, PartB_Proc_Ctrl=part_b,
        PartC_Raw_Volve=part_c, PartD_Synthetic=part_d,
    )
    out_dir = ROOT / "results"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / ("bosr_t11_smoke.json" if smoke else "bosr_t11.json")
    out_path.write_text(json.dumps(out, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"\nSaved {out_path} ({time.monotonic()-t0:.1f}s)", flush=True)


if __name__ == "__main__":
    main()
