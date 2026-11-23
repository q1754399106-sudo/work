"""T1.1 Pre-registered Full 8-Block Aperture-Invariance, Primary-Corridor Autofocus (B-OSR-Auto),
Constant-Wavenumber Zero-Stretch Beamforming (CW-BOSR: f*H(f)=C0), Processed Volve Control,
and Brand-New 8-Block Synthetic Self-Check (seed=900000/900001).

Read-only access to D:/ChatGPT/研究 caches. Low RAM (< 1.2 GB), single process, 2 BLAS threads.
Usage: python scripts/bosr_t11.py [--smoke]
"""
import os, sys, json, math, time
os.environ.setdefault("OPENBLAS_NUM_THREADS", "2")
os.environ.setdefault("MKL_NUM_THREADS", "2")
os.environ.setdefault("OMP_NUM_THREADS", "2")

from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import bosr_diag as bd

ROOT = bd.ROOT
DT, NT, T, NF, F = bd.DT, bd.NT, bd.T, bd.NF, bd.F
CENT, WL, REF, OFFS = bd.CENT, bd.WL, bd.REF, bd.OFFS
HW_INCOH, CMP_LIM, NBLK, FINE_PER_BLK = bd.HW_INCOH, bd.CMP_LIM, bd.NBLK, bd.FINE_PER_BLK
BAND, FLOOR, MINBW, STRETCH, MUTE, MUTE_PAD = bd.BAND, bd.FLOOR, bd.MINBW, bd.STRETCH, bd.MUTE, bd.MUTE_PAD
ZS, VW = bd.ZS, bd.VW

ALPHA_GRID = np.round(np.linspace(0.84, 1.18, 18), 4)
AL_BAND = (F >= 12.0) & (F <= 35.0)
BAND_LO = (12.0, 30.0)
BAND_HI = (30.0, 55.0)


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
    rat_full = 0.5 * np.log(np.maximum(Px, 1e-30) / np.maximum(P0, 1e-30))
    f_full = F[s_full]
    y_full = float(-np.polyfit(f_full, rat_full[s_full], 1)[0])
    if f_full[-1] - f_full[0] >= 20.0:
        fc = 0.5 * (f_full[0] + f_full[-1])
        c2 = float(np.polyfit(f_full - fc, rat_full[s_full], 2)[0])
    else:
        c2 = float("nan")
    fmax_v = float(f_full[-1])
    mlo = (F >= BAND_LO[0]) & (F <= min(BAND_LO[1], fmax_v))
    mhi = (F >= BAND_HI[0]) & (F <= min(BAND_HI[1], fmax_v))
    y_lo = float(-np.polyfit(F[mlo], rat_full[mlo], 1)[0]) if mlo.sum() >= 5 and (F[mlo][-1] - F[mlo][0] >= 10.0) else float("nan")
    y_hi = float(-np.polyfit(F[mhi], rat_full[mhi], 1)[0]) if mhi.sum() >= 5 and (F[mhi][-1] - F[mhi][0] >= 10.0) else float("nan")
    return y_full, y_lo, y_hi, c2


def fast_boot_ci(blocks_rows, mode="noint", seed=79, nboot=2000):
    """Exact algebraic sufficient-statistic block bootstrap (identical to bd.boot_blocks to 1e-15)."""
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
    allr = [r for b in blocks_rows for r in b]
    pt = bd.fit1_noint(allr)
    ci = fast_boot_ci(blocks_rows, "noint")
    fe_all = bd.fit_fe(allr, "dS", "y", exclude_near=False)
    fe_ci = fast_boot_ci(blocks_rows, "fe")
    g_mid = [r for r in allr if 1.15 <= r["tc"] <= 1.65]
    g_deep = [r for r in allr if 1.95 <= r["tc"] <= 2.25]

    def fit_key(rows, yk):
        sub = [r for r in rows if np.isfinite(r.get(yk, float("nan")))]
        if len(sub) < 3:
            return float("nan"), float("nan")
        x = np.pi * np.array([r["dS"] for r in sub])
        y = np.array([r[yk] for r in sub])
        return float(1000.0 * (x @ y) / (x @ x)), bd.fit_fe(sub, "dS", yk, exclude_near=False)

    qlo_ni, qlo_fe = fit_key(allr, "y_lo")
    qhi_ni, qhi_fe = fit_key(allr, "y_hi")
    c2_norm = [r["c2"] / r["dS"] for r in allr if np.isfinite(r.get("c2", float("nan"))) and r["dS"] > 0.05]
    alphas = [r["alpha"] for r in allr if "alpha" in r and np.isfinite(r["alpha"])]
    return dict(
        invQ_noint=pt, ci90_noint=ci,
        invQ_fe=fe_all, ci90_fe=fe_ci,
        n_rows=len(allr),
        q_12_16_noint=bd.fit1_noint(g_mid), q_12_16_fe=bd.fit_fe(g_mid, "dS", "y", False),
        q_20_22_noint=bd.fit1_noint(g_deep), q_20_22_fe=bd.fit_fe(g_deep, "dS", "y", False),
        q_lo_12_30_noint=qlo_ni, q_lo_12_30_fe=qlo_fe,
        q_hi_30_55_noint=qhi_ni, q_hi_30_55_fe=qhi_fe,
        c2_over_dS_med=float(np.median(c2_norm)) if c2_norm else float("nan"),
        alpha_med=float(np.median(alphas)) if alphas else 1.0,
        v_eff_ratio_med=float(1.0 / math.sqrt(np.median(alphas))) if alphas else 1.0,
    )


# -------------------- Unified Single-Pass Pack Extractor for Fixed + Auto + CW-BOSR --------------------
def align_window_finecmps(X, off, cmp_arr, fine_edges, tab, x0, hw_max):
    """Extract and rigidly zero-stretch align ALL traces in [x0 - hw_max, x0 + hw_max] ONCE per fine CMP!
    Returns lightweight frequency-domain spectra W_tr (ntr, nF) + xi + dt0 so Fixed(H=100,150,200),
    AutoBlock(H=100,150,200), AND CW-BOSR(C0=3600,4800) ALL share the exact same single FFT pass!"""
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
        m = (cmp_arr >= fine_edges[fi]) & (cmp_arr < fine_edges[fi + 1]) & (off >= x0 - hw_max) & (off < x0 + hw_max)
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
        # Sort by xi once for exact trapezoidal density quadrature
        ord_i = np.argsort(xi)
        seg = seg[ord_i]; xi = xi[ord_i]; dt0 = dt0[ord_i]
        S_pad = np.fft.rfft(seg, n_fft_pad, axis=1)
        al0 = np.fft.irfft(S_pad * np.exp(2j * np.pi * f_pad[None, :] * dt0[:, None]), n_fft_pad, axis=1)[:, :n_pad]
        W0 = np.fft.rfft(al0[:, k_win] * wt[None, :], NF, axis=1)  # (ntr, nF)
        packs.append(dict(xi=xi, dt0=dt0, W0=W0))

    if not packs:
        return None
    return dict(packs=packs, sig_t2=sig_t2)


def trapz_dxi(xi, hw_cap):
    n = len(xi)
    if n == 1:
        return np.ones(1, dtype=np.float64)
    dxi = np.empty(n, dtype=np.float64)
    dxi[1:-1] = 0.5 * (xi[2:] - xi[:-2])
    dxi[0] = xi[1] - xi[0]
    dxi[-1] = xi[-1] - xi[-2]
    return np.clip(dxi, 1e-3, 0.5 * hw_cap)


def eval_fixed_and_auto(pack_obj, hw, kappa_p):
    """Evaluate BOTH Fixed1D (alpha=1) AND AutoBlock (12-35 Hz primary-corridor optimal alpha)
    at aperture half-width `hw` purely in frequency domain (zero extra FFTs!)."""
    f_al = F[AL_BAND]
    d_alphas = ALPHA_GRID - 1.0
    sub_packs = []
    e_curve = np.zeros(len(ALPHA_GRID), dtype=np.float64)

    for p in pack_obj["packs"]:
        m = np.abs(p["xi"]) < hw
        if m.sum() < 6:
            continue
        xi = p["xi"][m]; dt0 = p["dt0"][m]; W0 = p["W0"][m]
        w_ap = 0.5 * (1.0 + np.cos(np.pi * np.clip(xi / hw, -1.0, 1.0)))
        w_raw = w_ap * trapz_dxi(xi, hw)
        ws = float(np.sum(w_raw))
        if ws <= 0:
            continue
        w_i = w_raw / ws
        sig_xi2 = float(np.sum(w_i * xi ** 2))
        sub_packs.append((w_i, dt0, W0, sig_xi2))
        wW_al = w_i[:, None] * W0[:, AL_BAND]
        ang = (2.0 * np.pi) * (d_alphas[:, None, None] * dt0[None, :, None] * f_al[None, None, :])
        B_al = np.sum(wW_al[None, :, :] * (np.cos(ang) + 1j * np.sin(ang)), axis=1)
        e_curve += np.sum(B_al.real ** 2 + B_al.imag ** 2, axis=1)

    if not sub_packs:
        return None
    alpha_opt = float(ALPHA_GRID[int(np.argmax(e_curve))])
    da_opt = alpha_opt - 1.0

    Pb_fix = np.zeros_like(F)
    Pb_auto = np.zeros_like(F)
    sig_xi2_m = 0.0
    for w_i, dt0, W0, sxi2 in sub_packs:
        sig_xi2_m += sxi2
        Bf = np.sum(w_i[:, None] * W0, axis=0)
        Pb_fix += Bf.real ** 2 + Bf.imag ** 2
        if abs(da_opt) < 1e-12:
            Pb_auto += Bf.real ** 2 + Bf.imag ** 2
        else:
            ph = np.exp((2j * np.pi * da_opt) * (dt0[:, None] * F[None, :]))
            Ba = np.sum(w_i[:, None] * W0 * ph, axis=0)
            Pb_auto += Ba.real ** 2 + Ba.imag ** 2

    nf = len(sub_packs)
    Pb_fix /= nf; Pb_auto /= nf; sig_xi2_m /= nf
    B2_f = np.maximum(np.exp(-((2.0 * np.pi * F) ** 2) * (kappa_p ** 2) * pack_obj["sig_t2"] * sig_xi2_m), 0.25)
    B2_a = np.maximum(np.exp(-((2.0 * np.pi * F) ** 2) * ((alpha_opt * kappa_p) ** 2) * pack_obj["sig_t2"] * sig_xi2_m), 0.25)
    return Pb_fix / B2_f, Pb_auto / B2_a, alpha_opt


def eval_cw_bosr(pack_obj, c0_m_hz, f_min=12.0):
    """Constant-Wavenumber Zero-Stretch Beamformer: H(f) = c0_m_hz / f across [12, 60] Hz
    (padded cleanly on [6, 75] Hz so ghost boxcar smoothing has zero boundary drop)."""
    hw_max = c0_m_hz / f_min
    ib = np.where((F >= 6.0) & (F <= 75.0))[0]
    f_use = F[ib]
    H_f = np.minimum(c0_m_hz / np.maximum(f_use, f_min), hw_max)

    Pb_acc = np.zeros_like(F)
    n_fine = 0
    for p in pack_obj["packs"]:
        m = np.abs(p["xi"]) < hw_max
        if m.sum() < 6:
            continue
        xi = p["xi"][m]
        W_sub = p["W0"][m][:, ib]
        dxi = trapz_dxi(xi, hw_max)
        u = xi[:, None] / H_f[None, :]
        w_ap = np.where(np.abs(u) < 1.0, 0.5 * (1.0 + np.cos(np.pi * np.clip(u, -1.0, 1.0))), 0.0)
        w_raw = w_ap * dxi[:, None]
        w_sum = np.sum(w_raw, axis=0)
        ok_f = w_sum > 0
        if not np.all(ok_f[(f_use >= 12.0) & (f_use <= 55.0)]):
            continue
        w_norm = np.where(ok_f[None, :], w_raw / np.maximum(w_sum[None, :], 1e-30), 0.0)
        B = np.sum(w_norm * W_sub, axis=0)
        Pb_acc[ib] += B.real ** 2 + B.imag ** 2
        n_fine += 1

    if n_fine == 0:
        return None
    return Pb_acc / n_fine


def run_unified_suite(blocks, R_glob, zw, c0sq, comp, pz_sign, hws, c0_list, do_ghost_smooth=True, run_incoh=False):
    """Single-pass extraction computing Incoh (opt), Fixed1D(hws), AutoBlock(hws), AND CW-BOSR(c0_list)!"""
    out = {}
    s_cal = 1.0 / math.sqrt(max(c0sq, 1e-30))
    blk_edges = np.linspace(-CMP_LIM, CMP_LIM, NBLK + 1)
    tabs = [R_glob.table(tc) for tc in CENT]
    ptabs = [bd.ptable(t) for t in tabs]
    tp = [bd.ptable(R_glob.table(tc + 0.05)) for tc in CENT]
    tm = [bd.ptable(R_glob.table(max(tc - 0.05, 0.25))) for tc in CENT]

    def get_sig(A, B, mask, cw):
        if comp == "pz_up":
            return A[mask] + pz_sign * (s_cal / cw) * B[mask]
        elif comp == "a":
            return A[mask]
        return B[mask]

    # Optional incoherent baseline
    if run_incoh:
        blk_rows = []
        for bi, (A, B, off, _) in enumerate(blocks):
            rows = []
            kr = off <= REF
            if kr.sum() < 8:
                blk_rows.append(rows); continue
            xr = float(np.sqrt(np.mean(off[kr] ** 2)))
            for tc, tab, ptab in zip(CENT, tabs, ptabs):
                r0 = tab(np.array([xr])); tr = float(r0[0][0])
                cwr, taur = bd.ghost(ptab(xr)[0], zw)
                P0_r = bd.power_incoh(get_sig(A, B, kr, cwr), tr - WL / 2, tr + WL / 2)
                if P0_r is None: continue
                P0 = bd.smooth_period(P0_r, taur) if do_ghost_smooth else P0_r
                for o in OFFS:
                    k = (off >= o - HW_INCOH) & (off < o + HW_INCOH)
                    if k.sum() < 8: continue
                    xo = float(np.sqrt(np.mean(off[k] ** 2)))
                    tx, _, _, ls = (float(v[0]) for v in tab(np.array([xo])))
                    if not np.isfinite(tx) or tx / tc > STRETCH or tx - WL / 2 < MUTE(o + HW_INCOH) + MUTE_PAD or tx + WL / 2 > T[-1]:
                        continue
                    cw, tau = bd.ghost(ptab(xo)[0], zw)
                    Px_r = bd.power_incoh(get_sig(A, B, k, cw), tx - WL / 2, tx + WL / 2)
                    if Px_r is None: continue
                    Px = bd.smooth_period(Px_r, tau) if do_ghost_smooth else Px_r
                    s = bd.band_slice(P0, Px)
                    if s is None: continue
                    yf, ylo, yhi, c2 = fit_subband_and_curv(P0, Px, s)
                    rows.append(dict(block=bi, tc=float(tc), off=xo, y=yf, y_lo=ylo, y_hi=yhi, c2=c2,
                                     dS=ls - float(r0[3][0]), tau_sed=ls, alpha=1.0))
            blk_rows.append(rows)
        out["Incoh"] = summarise_t11(blk_rows)

    # Single-pass unified beamforming extraction (max aperture = max(max(hws), max(c0_list)/12.0))
    hw_max_all = max(max(hws), max(c0_list) / 12.0)
    rows_fix = {hw: [[] for _ in range(len(blocks))] for hw in hws}
    rows_auto = {hw: [[] for _ in range(len(blocks))] for hw in hws}
    rows_cw = {c0: [[] for _ in range(len(blocks))] for c0 in c0_list}

    for bi, (A, B, off, cmp_arr) in enumerate(blocks):
        fine_edges = np.linspace(blk_edges[bi], blk_edges[bi + 1], FINE_PER_BLK + 1)
        kr = off <= REF
        xr = float(np.sqrt(np.mean(off[kr] ** 2))) if kr.sum() >= 8 else float("nan")

        for wi, (tc, tab, ptab) in enumerate(zip(CENT, tabs, ptabs)):
            # Near-offset pack (used by Fixed & AutoBlock)
            P0_fix = {}; P0_auto = {}; ls0_near = float("nan")
            if np.isfinite(xr):
                r0 = tab(np.array([xr])); ls0_near = float(r0[3][0])
                cwr, taur = bd.ghost(ptab(xr)[0], zw)
                pk0 = align_window_finecmps(get_sig(A, B, kr, cwr), off[kr], cmp_arr[kr], fine_edges, tab, xr, REF)
                if pk0 is not None:
                    kp0 = float((tp[wi](xr)[0] - tm[wi](xr)[0]) / (0.05 + min(0.05, tc - 0.25)))
                    r_n = eval_fixed_and_auto(pk0, REF, kp0)
                    if r_n is not None:
                        pf0, pa0, _ = r_n
                        P0_fix_val = bd.smooth_period(pf0, taur) if do_ghost_smooth else pf0
                        P0_auto_val = bd.smooth_period(pa0, taur) if do_ghost_smooth else pa0
                        for hw in hws:
                            P0_fix[hw] = P0_fix_val; P0_auto[hw] = P0_auto_val

            # Offset loop: ONE pack extraction per (bi, wi, o) feeds ALL hws AND ALL c0_list!
            cw_bins = {c0: [] for c0 in c0_list}
            for o in OFFS:
                tx0, _, _, ls0_o = (float(v[0]) for v in tab(np.array([float(o)])))
                if not np.isfinite(tx0) or tx0 / tc > STRETCH or tx0 - WL / 2 < MUTE(o + HW_INCOH) + MUTE_PAD or tx0 + WL / 2 > T[-1]:
                    continue
                kb = (off >= o - hw_max_all) & (off < o + hw_max_all)
                if kb.sum() < 12:
                    continue
                cw_o, tau_o = bd.ghost(ptab(float(o))[0], zw)
                pkx = align_window_finecmps(get_sig(A, B, kb, cw_o), off[kb], cmp_arr[kb], fine_edges, tab, float(o), hw_max_all)
                if pkx is None:
                    continue
                kp = float((tp[wi](float(o))[0] - tm[wi](float(o))[0]) / (0.05 + min(0.05, tc - 0.25)))

                # Fixed1D & AutoBlock at each hw
                if P0_fix:
                    for hw in hws:
                        res_fa = eval_fixed_and_auto(pkx, hw, kp)
                        if res_fa is None:
                            continue
                        pfx, pax, aopt = res_fa
                        Pxf = bd.smooth_period(pfx, tau_o) if do_ghost_smooth else pfx
                        Pxa = bd.smooth_period(pax, tau_o) if do_ghost_smooth else pax
                        sf = bd.band_slice(P0_fix[hw], Pxf)
                        if sf is not None:
                            yf, ylo, yhi, c2 = fit_subband_and_curv(P0_fix[hw], Pxf, sf)
                            rows_fix[hw][bi].append(dict(block=bi, tc=float(tc), off=float(o), y=yf, y_lo=ylo, y_hi=yhi,
                                                         c2=c2, dS=ls0_o - ls0_near, tau_sed=ls0_o, alpha=1.0))
                        sa = bd.band_slice(P0_auto[hw], Pxa)
                        if sa is not None:
                            yf, ylo, yhi, c2 = fit_subband_and_curv(P0_auto[hw], Pxa, sa)
                            rows_auto[hw][bi].append(dict(block=bi, tc=float(tc), off=float(o), y=yf, y_lo=ylo, y_hi=yhi,
                                                          c2=c2, dS=ls0_o - ls0_near, tau_sed=ls0_o, alpha=aopt))

                # CW-BOSR at each c0 (symmetric reference at smallest valid offset o >= 500 m!)
                for c0 in c0_list:
                    pcw = eval_cw_bosr(pkx, c0)
                    if pcw is None or pcw.max() <= 0:
                        continue
                    P_sm = bd.smooth_period(np.maximum(pcw, 1e-30 * pcw.max()), tau_o) if do_ghost_smooth else pcw
                    cw_bins[c0].append((float(o), ls0_o, P_sm))

            for c0 in c0_list:
                bl = cw_bins[c0]
                if len(bl) < 2:
                    continue
                o0, ls0_c, P0_c = bl[0]
                for ox, lsx, Px_c in bl[1:]:
                    sc = bd.band_slice(P0_c, Px_c)
                    if sc is None:
                        continue
                    yf, ylo, yhi, c2 = fit_subband_and_curv(P0_c, Px_c, sc)
                    rows_cw[c0][bi].append(dict(block=bi, tc=float(tc), off=ox, y=yf, y_lo=ylo, y_hi=yhi, c2=c2,
                                                dS=lsx - ls0_c, tau_sed=lsx, alpha=1.0))

    for hw in hws:
        out[f"H{int(hw)}"] = dict(Fixed1D=summarise_t11(rows_fix[hw]), AutoBlock=summarise_t11(rows_auto[hw]))
    for c0 in c0_list:
        out[f"CW{int(c0)}"] = summarise_t11(rows_cw[c0])
    return out


# -------------------- Part D: Realistic 8-Block Synthetic Self-Check (seed=900000 / 900001) --------------------
def synthesize_volve_v2(R_model, zw, q_vals, q_mult, amp_m_scale, seed, n_blocks=NBLK):
    rng = np.random.default_rng(seed)
    fr = np.fft.rfftfreq(2 * NT, DT)
    f_ref = 25.0
    W_src = (fr / 25.0) ** 2 * np.exp(-(fr / 25.0) ** 2)
    W_sch = np.exp(-0.5 * ((fr - 12.0) / 3.5) ** 2)
    blk_edges = np.linspace(-CMP_LIM, CMP_LIM, NBLK + 1)
    tc_grid = np.linspace(0.35, 2.95, 27)
    x_grid = np.linspace(0.0, 3800.0, 96)

    # Analytical true effective 1000/Q under acceptance mask
    tnum = tden = 0.0
    tabs_m = [R_model.table(tc) for tc in CENT]
    for tc, tab in zip(CENT, tabs_m):
        r0 = tab(np.array([500.0]))
        for o in OFFS:
            tx, la, lb, ls = (float(v[0]) for v in tab(np.array([float(o)])))
            if not np.isfinite(tx) or tx / tc > STRETCH or tx - WL / 2 < MUTE(o + HW_INCOH) + MUTE_PAD or tx + WL / 2 > T[-1]:
                continue
            dS = ls - float(r0[3][0])
            dts = (la - float(r0[1][0])) / q_vals[0] + (lb - float(r0[2][0])) / q_vals[1] + max((ls - la - lb) - (float(r0[3][0]) - float(r0[1][0]) - float(r0[2][0])), 0.0) / q_vals[2]
            tnum += dS * dts; tden += dS * dS
    true_iq = 1000.0 * tnum / max(tden, 1e-30)

    base_offs = np.linspace(60.0, 3700.0, 180)  # 20 m realistic spacing
    blocks = []
    v_errs = []
    for bi in range(n_blocks):
        v_err = float(1.025 + rng.uniform(-0.020, 0.020))
        v_errs.append(v_err)
        R_true = bd.Rays(v_err, zw)
        T_lut = np.zeros((len(tc_grid), len(x_grid)))
        LA_lut = np.zeros_like(T_lut); LB_lut = np.zeros_like(T_lut); LS_lut = np.zeros_like(T_lut); P_lut = np.zeros_like(T_lut)
        for ti, tcg in enumerate(tc_grid):
            tb = R_true.table(tcg); ptb = bd.ptable(tb)
            T_lut[ti], LA_lut[ti], LB_lut[ti], LS_lut[ti] = tb(x_grid)
            P_lut[ti] = ptb(x_grid)

        fine_edges = np.linspace(blk_edges[bi], blk_edges[bi + 1], FINE_PER_BLK + 1)
        A_list, off_list, cmp_list = [], [], []
        for fi in range(FINE_PER_BLK):
            offs = np.clip(base_offs + rng.uniform(-5.0, 5.0, size=len(base_offs)), 40.0, 3750.0)
            cmps = rng.uniform(fine_edges[fi] + 1.0, fine_edges[fi + 1] - 1.0, size=len(offs))
            spec = np.zeros((len(offs), len(fr)), dtype=np.complex128)
            t0_refs = np.linspace(0.48, 2.72, 42) + rng.uniform(-0.02, 0.02, size=42)
            amps_ref = rng.normal(0.0, 1.0, size=42)
            for t0_r, amp_r in zip(t0_refs, amps_ref):
                idx = int(np.clip(np.searchsorted(tc_grid, t0_r) - 1, 0, len(tc_grid) - 2))
                wt_i = (t0_r - tc_grid[idx]) / (tc_grid[idx + 1] - tc_grid[idx])
                t_x = np.interp(offs, x_grid, (1 - wt_i) * T_lut[idx] + wt_i * T_lut[idx + 1])
                la_x = np.interp(offs, x_grid, (1 - wt_i) * LA_lut[idx] + wt_i * LA_lut[idx + 1])
                lb_x = np.interp(offs, x_grid, (1 - wt_i) * LB_lut[idx] + wt_i * LB_lut[idx + 1])
                ls_x = np.interp(offs, x_grid, (1 - wt_i) * LS_lut[idx] + wt_i * LS_lut[idx + 1])
                p_x = np.interp(offs, x_grid, (1 - wt_i) * P_lut[idx] + wt_i * P_lut[idx + 1])
                ok = np.isfinite(t_x) & np.isfinite(p_x) & np.isfinite(la_x) & np.isfinite(ls_x) & (t_x < T[-1] + 0.2)
                if not np.any(ok): continue
                tstar = la_x[ok] / q_vals[0] + lb_x[ok] / q_vals[1] + np.maximum(ls_x[ok] - la_x[ok] - lb_x[ok], 0.0) / q_vals[2]
                cw_x, tau_r = bd.ghost(p_x[ok], zw); tau_s = 2.0 * ZS * cw_x / VW
                lnf = np.log(np.maximum(fr, 1.0) / f_ref)
                att = np.exp((-np.pi * fr[None, :] + 2j * fr[None, :] * lnf[None, :]) * tstar[:, None])
                g_r = 1.0 - 0.85 * np.exp(-2j * np.pi * fr[None, :] * tau_r[:, None])
                g_s = 1.0 - 0.90 * np.exp(-2j * np.pi * fr[None, :] * tau_s[:, None])
                spec[ok] += amp_r * W_src[None, :] * g_s * g_r * att * np.exp(-2j * np.pi * fr[None, :] * t_x[ok, None])

            t0_m = np.linspace(0.75, 2.55, 22) + rng.uniform(-0.03, 0.03, size=22)
            v_m = rng.uniform(1520.0, 1680.0, size=22)
            amp_m = rng.normal(0.0, amp_m_scale, size=22) * (0.4 + 0.8 * (t0_m / 2.2))
            for tm0, vm, am in zip(t0_m, v_m, amp_m):
                t_mx = np.sqrt(tm0 ** 2 + (offs / vm) ** 2)
                p_mx = offs / (vm ** 2 * t_mx)
                cw_m, tau_rm = bd.ghost(p_mx, zw); tau_sm = 2.0 * ZS * cw_m / VW
                att_m = np.exp(-np.pi * fr[None, :] * (t_mx / q_mult)[:, None])
                g_rm = 1.0 - 0.85 * np.exp(-2j * np.pi * fr[None, :] * tau_rm[:, None])
                g_sm = 1.0 - 0.90 * np.exp(-2j * np.pi * fr[None, :] * tau_sm[:, None])
                spec += am * W_src[None, :] * g_sm * g_rm * att_m * np.exp(-2j * np.pi * fr[None, :] * t_mx[:, None])

            tr = np.fft.irfft(spec, 2 * NT, axis=1)[:, :NT].astype(np.float32)
            kn = offs <= REF
            if np.any(kn):
                rms_n = float(np.sqrt(np.mean(tr[kn] ** 2))) + 1e-12
                sch = np.fft.irfft(np.fft.rfft(rng.normal(0, 1, (int(kn.sum()), 2 * NT)), axis=1) * W_sch[None, :], 2 * NT, axis=1)[:, :NT]
                tr[kn] += (0.40 * rms_n / (float(np.sqrt(np.mean(sch ** 2))) + 1e-12) * sch).astype(np.float32)
            rms_dp = float(np.sqrt(np.mean(tr[:, int(1.6 / DT):int(2.4 / DT)] ** 2))) + 1e-12
            tr += rng.normal(0.0, 0.08 * rms_dp, size=tr.shape).astype(np.float32)
            tr = bd.hp(tr)
            A_list.append(tr); off_list.append(offs); cmp_list.append(cmps)
        A_blk = np.vstack(A_list)
        blocks.append((A_blk, A_blk, np.concatenate(off_list), np.concatenate(cmp_list)))
    return blocks, float(true_iq), v_errs


def fmt_row(tag, r):
    return (f"  {tag:14s}: noint={r['invQ_noint']:+6.2f} {np.round(r['ci90_noint'],1).tolist()} | "
            f"FE={r['invQ_fe']:+6.2f} {np.round(r['ci90_fe'],1).tolist()} | "
            f"12-16/20-22(FE)={r['q_12_16_fe']:+5.1f}/{r['q_20_22_fe']:+5.1f} | "
            f"lo/hi(FE)={r['q_lo_12_30_fe']:+5.1f}/{r['q_hi_30_55_fe']:+5.1f} | "
            f"c2/dS={r['c2_over_dS_med']:+.2e} | a={r['alpha_med']:.3f}")


def main():
    smoke = "--smoke" in sys.argv
    t0 = time.monotonic()
    print(f"=== T1.1 Full Aperture + Auto + CW-BOSR & Synthetic Self-Check (smoke={smoke}) ===", flush=True)

    H, h = bd.load_raw_headers()
    sel = (h["dline"] == 5081) & (np.abs(h["cmp"]) < CMP_LIM)
    key = lambda i: (h["ffid"][i], h["dpt"][i])
    ia = {key(i): i for i in np.where(sel & (h["trid"] == 11))[0]}
    ib = {key(i): i for i in np.where(sel & (h["trid"] == 12))[0]}
    pairs = sorted(set(ia) & set(ib))
    ka = np.array([ia[p] for p in pairs]); kb = np.array([ib[p] for p in pairs])
    cmp_all = h["cmp"][ka]; off_all = h["off"][ka]
    zw = float(np.median(h["wd"][ka]))
    e = np.linspace(-CMP_LIM, CMP_LIM, NBLK + 1)
    blk_indices = [0, 1] if smoke else list(range(NBLK))
    raw_blocks = []
    for i in blk_indices:
        m = (cmp_all >= e[i]) & (cmp_all < e[i + 1])
        raw_blocks.append((bd.read_raw_samples(H, ka[m]), bd.read_raw_samples(H, kb[m]), off_all[m], cmp_all[m]))
    print(f"Loaded raw Volve ({len(raw_blocks)} blocks, zw={zw:.1f} m) in {time.monotonic()-t0:.1f}s", flush=True)

    R_glob = bd.Rays(1.0, zw)
    tabs_g = [R_glob.table(tc) for tc in CENT]; ptabs_g = [bd.ptable(t) for t in tabs_g]

    # Part A: Exact all-window c0sq & M0 reproduction
    c0sq_all = compute_c0sq_allwin(raw_blocks, tabs_g)
    d1_tmp = bd.diagnose_components(raw_blocks, tabs_g, zw)
    pz_sign = d1_tmp["pz_amp_best_sign"]
    print(f"\n--- Part A: All-window c0sq={c0sq_all:.4f} (field79raw=2.114), pz_sign={pz_sign} ---", flush=True)
    part_a = dict(c0sq_allwin=c0sq_all, pz_sign=pz_sign)
    for comp in ("pz_pow", "a", "b"):
        s, _ = bd.run_incoh_field79raw(raw_blocks, tabs_g, ptabs_g, zw, c0sq_all, comp, pz_sign=pz_sign)
        part_a[f"incoh_{comp}"] = s
        print(f"  incoh {comp:>6}: noint={s['invQ_noint']:+6.2f} {np.round(s['ci90_noint'],1).tolist()} | FE={s['invQ_fe']:+6.2f}", flush=True)

    hws = [100.0, 150.0, 200.0]
    c0s = [3600.0, 4800.0]

    # Part B: Processed Volve seq079_dline5081
    print(f"\n--- Part B: Processed Volve seq079_dline5081 ({time.monotonic()-t0:.1f}s) ---", flush=True)
    Hp, hp_d = bd.load_proc_headers()
    proc_blocks = []
    for i in blk_indices:
        idx = np.where((hp_d["cmp"] >= e[i]) & (hp_d["cmp"] < e[i + 1]))[0]
        G = bd.read_proc_samples(Hp, idx)
        proc_blocks.append((G, G, hp_d["off"][idx], hp_d["cmp"][idx]))
    part_b = run_unified_suite(proc_blocks, R_glob, zw, 1.0, "a", 1, hws, c0s, do_ghost_smooth=False, run_incoh=True)
    print(fmt_row("Proc Incoh", part_b["Incoh"]), flush=True)
    for hw in hws:
        print(fmt_row(f"Proc H{int(hw)} Fix", part_b[f"H{int(hw)}"]["Fixed1D"]), flush=True)
        print(fmt_row(f"Proc H{int(hw)} Auto", part_b[f"H{int(hw)}"]["AutoBlock"]), flush=True)
    for c0 in c0s:
        print(fmt_row(f"Proc CW{int(c0)}", part_b[f"CW{int(c0)}"]), flush=True)
    del proc_blocks

    # Part C: Raw Volve (a, b, pz_up)
    print(f"\n--- Part C: Raw Volve Fixed / Auto / CW-BOSR ({time.monotonic()-t0:.1f}s) ---", flush=True)
    part_c = {}
    for comp in ("a", "b", "pz_up"):
        print(f"  [Raw comp={comp}]", flush=True)
        rc = run_unified_suite(raw_blocks, R_glob, zw, c0sq_all, comp, pz_sign, hws, c0s, do_ghost_smooth=True, run_incoh=False)
        part_c[comp] = rc
        for hw in hws:
            print(fmt_row(f"H{int(hw)} Fix", rc[f"H{int(hw)}"]["Fixed1D"]), flush=True)
            print(fmt_row(f"H{int(hw)} Auto", rc[f"H{int(hw)}"]["AutoBlock"]), flush=True)
        for c0 in c0s:
            print(fmt_row(f"CW{int(c0)}", rc[f"CW{int(c0)}"]), flush=True)
    del raw_blocks

    # Part D: Brand-New 8-Block Synthetic Self-Check (seed=900000 / 900001)
    print(f"\n--- Part D: Synthetic Self-Check (seed=900000/900001, +2.5% +/- 2% v_err + mult + ghosts + near) ({time.monotonic()-t0:.1f}s) ---", flush=True)
    n_syn = 2 if smoke else NBLK
    part_d = {}
    for sname, qv, qm, am, sd in [
        ("S0_Qinf", (1e9, 1e9, 1e9), 1e9, 0.12, 900000),
        ("S1_QVSP", (60.0, 90.0, 100.0), 400.0, 0.12, 900001),
    ]:
        sblks, truth, verrs = synthesize_volve_v2(R_glob, zw, qv, qm, am, sd, n_blocks=n_syn)
        rs = run_unified_suite(sblks, R_glob, zw, 1.0, "a", 1, hws, c0s, do_ghost_smooth=True, run_incoh=True)
        rs["true_invQ_eff"] = truth; rs["v_err_blocks"] = verrs
        part_d[sname] = rs
        print(f"  [{sname}] truth={truth:+.2f}, v_errs={np.round(verrs,3).tolist()}", flush=True)
        print(fmt_row("Incoh", rs["Incoh"]), flush=True)
        for hw in hws:
            print(fmt_row(f"H{int(hw)} Fix", rs[f"H{int(hw)}"]["Fixed1D"]), flush=True)
            print(fmt_row(f"H{int(hw)} Auto", rs[f"H{int(hw)}"]["AutoBlock"]), flush=True)
        for c0 in c0s:
            print(fmt_row(f"CW{int(c0)}", rs[f"CW{int(c0)}"]), flush=True)

    out = dict(smoke=smoke, water_depth=zw, elapsed_s=round(time.monotonic() - t0, 1),
               PartA_M0_repro=part_a, PartB_Proc_Ctrl=part_b, PartC_Raw_Volve=part_c, PartD_Synthetic=part_d)
    out_dir = ROOT / "results"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / ("bosr_t11_smoke.json" if smoke else "bosr_t11.json")
    out_path.write_text(json.dumps(out, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"\nSaved {out_path} ({time.monotonic()-t0:.1f}s)", flush=True)


if __name__ == "__main__":
    main()
