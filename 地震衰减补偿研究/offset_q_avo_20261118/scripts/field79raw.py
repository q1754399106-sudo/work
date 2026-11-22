"""PLAN_FIELD_I_RAW_ZH.md: offset spectral ratio on the RAW Volve seq-79 PZDETECTED records (cable 5081), same estimator,
geometry, windows and mute guard as field79.py (amendment 1: single effective Q), with:
- recording-filter delay 48 ms removed (12 samples), first 3.2 s;
- P / Z summation per receiver: X = A + sigma * s * B, A = TRID 11, B = TRID 12 of the same shot and receiver point,
  s = rms(A) / rms(B) over 0.5-3.0 s (frequency independent, so it does not change any spectral slope); one sign sigma for
  the whole line, the one with the smaller ghost-comb index D of the near-offset spectra (windows 1.0-2.2 s);
  D = mean over 10-60 Hz of |ln P - ln P smoothed over 8 Hz|;
- band 10-60 Hz (swell noise below 10 Hz in raw records).
usage: field79raw.py synth           self-check: field geometry, receiver ghost on both components, swell noise
       field79raw.py run [comp]      comp = pz (primary), a (TRID 11 only), b (TRID 12 only)   -> results/field79raw[_a|_b].json"""
import sys, json, math, time
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
import field79 as F  # noqa: E402

F.BAND = (12.0, 60.0)   # raw: swell noise, low-cut below (hp)
RAW = F.MAIN / "volve_seq079" / "raw_seq079_PZ.sgy"; NS = 1876; TSZ = 240 + 4 * NS; SHIFT = 12; NT = F.NT; DT = F.DT; T = F.T
VW = 1500.0; SYN_SEED = 800_000; R_SB = 0.3   # seabed reflection coefficient of the synthetic water-layer reverberation
GHOST = True; SWELL = True; NOISE = True   # synthetic switches (diagnosis)


def headers():
    n = (RAW.stat().st_size - 3600) // TSZ; H = np.memmap(RAW, dtype=np.uint8, mode="r", offset=3600, shape=(n, TSZ))
    h = lambda o, t=">i4": np.frombuffer(np.ascontiguousarray(H[:, o:o + (4 if t == ">i4" else 2)]).tobytes(), t).astype(np.float64)
    A = math.radians(14); X0, Y0 = 435050.02, 6478563.53; il = lambda x, y: -(x - X0) * math.cos(A) + (y - Y0) * math.sin(A)
    s = il(h(180) / 100, h(184) / 100); g = il(h(188) / 100, h(192) / 100)
    return H, dict(ffid=h(8), trid=h(28, ">i2"), off=np.abs(h(36)), wd=h(64), dline=h(196), dpt=h(200), cmp=(s + g) / 2)


def samples(H, k):
    tr = F.ibm2f(np.frombuffer(np.ascontiguousarray(H[k, 240:240 + 4 * (NT + SHIFT)]).tobytes(), ">u4").reshape(len(k), NT + SHIFT))
    return hp(tr[:, SHIFT:])


def hp(X):
    """zero-phase low-cut applied to whole traces before any window: 0 below 6 Hz, cosine taper to 1 at 10 Hz."""
    n = 2 * X.shape[1]; f = np.fft.rfftfreq(n, DT); g = np.clip((f - 6.0) / 4.0, 0, 1); g = 0.5 - 0.5 * np.cos(np.pi * g)
    return np.fft.irfft(np.fft.rfft(X, n, axis=1) * g, n, axis=1)[:, :X.shape[1]].astype(np.float32)


def comb_index(P):
    ib = (F.F >= 10) & (F.F <= 60); w = max(3, int(round(8.0 / (F.F[1] - F.F[0])))); lp = np.log(np.maximum(P, 1e-30))
    sm = np.convolve(lp, np.ones(w) / w, mode="same"); return float(np.mean(np.abs(lp - sm)[ib]))


def near_D(Xs, offs, tabs):
    """ghost-comb index of the near-offset (|x| <= 350 m) mean spectrum, averaged over the windows 1.0-2.2 s."""
    k = offs <= F.REF; xr = float(np.sqrt(np.mean(offs[k] ** 2))); out = []
    for tc, tab in zip(F.CENT, tabs):
        if tc < 1.0:
            continue
        tr = float(tab(np.array([xr]))[0][0]); P = F.power(Xs[k], tr - F.WL / 2, tr + F.WL / 2)
        if P is not None:
            out.append(comb_index(P))
    return float(np.mean(out))


def pz(A, B, sig):
    k = (T >= 0.5) & (T < 3.0); s = np.sqrt(np.mean(A[:, k] ** 2, 1) / np.maximum(np.mean(B[:, k] ** 2, 1), 1e-30))
    return A + sig * s[:, None] * B


def estimate(blocksAB, tabs, comp):
    """blocksAB: per block (A, B, offsets). Returns summary with the chosen sign."""
    offs = np.concatenate([o for _, _, o in blocksAB]); info = {}
    if comp == "pz":
        Aall = np.concatenate([a for a, _, _ in blocksAB]); Ball = np.concatenate([b for _, b, _ in blocksAB])
        D = {sg: near_D(pz(Aall, Ball, sg), offs, tabs) for sg in (1, -1)}; sig = min(D, key=D.get)
        info = dict(sign=sig, D_plus=D[1], D_minus=D[-1], D_A=near_D(Aall, offs, tabs), D_B=near_D(Ball, offs, tabs))
    rows = []
    for a, b, o in blocksAB:
        X = pz(a, b, info["sign"]) if comp == "pz" else (a if comp == "a" else b); rows.append(F.block_rows(X, o, tabs))
    return dict(**info, **F.summarise(rows)), rows


def ptable(R, tc):
    """horizontal slowness p(x) at the receiver for the reflection of recorded zero-offset time tc."""
    tab = R.table(tc); return lambda x: (tab(np.atleast_1d(x) + 10.0)[0] - tab(np.maximum(np.atleast_1d(x) - 10.0, 0.0))[0]) / 20.0


def ghost(p, zw):
    cw = np.sqrt(np.maximum(1 - (VW * np.asarray(p, float)) ** 2, 0.05)); return cw, 2 * zw * cw / VW


def smooth_period(P, tau):
    """boxcar over one receiver-ghost period 1/tau of ln P: removes any response periodic in f with period 1/tau."""
    w = max(3, int(round(1.0 / tau / (F.F[1] - F.F[0])))); lp = np.log(np.maximum(P, 1e-30 * P.max()))
    return np.exp(np.convolve(lp, np.ones(w) / w, mode="same"))


def rows_pz(A, B, off, tabs, ptabs, zw, c0sq, comp):
    """as field79.block_rows, with the power-domain P/Z combination P_A + P_B / (c0^2 cos^2) per row (angle from the ray
    table at the receiver) and ghost-period smoothing of the near and far spectra; comp 'pz', 'a' or 'b'."""
    rows = []; kr = off <= F.REF
    if kr.sum() < 10:
        return rows
    xr = float(np.sqrt(np.mean(off[kr] ** 2)))
    comb = lambda PA, PB, cw: PA + PB / (c0sq * cw ** 2) if comp == "pz" else (PA if comp == "a" else PB)
    for tc, tab, ptab in zip(F.CENT, tabs, ptabs):
        r0 = tab(np.array([xr])); tr = float(r0[0][0]); cwr, taur = ghost(ptab(xr)[0], zw)
        PA0, PB0 = F.power(A[kr], tr - F.WL / 2, tr + F.WL / 2), F.power(B[kr], tr - F.WL / 2, tr + F.WL / 2)
        if PA0 is None or PB0 is None:
            continue
        P0 = smooth_period(comb(PA0, PB0, cwr), taur)
        for o in F.OFFS:
            k = (off >= o - F.HW) & (off < o + F.HW)
            if k.sum() < 10:
                continue
            xo = float(np.sqrt(np.mean(off[k] ** 2))); tx, la, lb, ls = (float(v[0]) for v in tab(np.array([xo])))
            if not np.isfinite(tx) or tx / tc > F.STRETCH or tx - F.WL / 2 < F.MUTE(o + F.HW) + F.MUTE_PAD or tx + F.WL / 2 > T[-1]:
                continue
            PA, PB = F.power(A[k], tx - F.WL / 2, tx + F.WL / 2), F.power(B[k], tx - F.WL / 2, tx + F.WL / 2)
            if PA is None or PB is None:
                continue
            cw, tau = ghost(ptab(xo)[0], zw); Px = smooth_period(comb(PA, PB, cw), tau); s = F.band(P0, Px)
            if s is None:
                continue
            y = float(-np.polyfit(F.F[s], 0.5 * np.log(Px[s] / P0[s]), 1)[0])
            rows.append(dict(tc=float(tc), off=xo, y=y, fmax=float(F.F[s.stop - 1]), cw=float(cw), dA=la - float(r0[1][0]), dB=lb - float(r0[2][0]),
                             dS=ls - float(r0[3][0])))
    return rows


def c0sq_near(blocks, tabs):
    """relative B / A sensitivity from the near-offset reference windows (cos ~ 1): pooled power ratio over 10-60 Hz."""
    ib = (F.F >= 10) & (F.F <= 60); sa = sb = 0.0
    for A, B, off in blocks:
        kr = off <= F.REF
        if kr.sum() < 10:
            continue
        xr = float(np.sqrt(np.mean(off[kr] ** 2)))
        for tc, tab in zip(F.CENT, tabs):
            tr = float(tab(np.array([xr]))[0][0]); PA, PB = F.power(A[kr], tr - F.WL / 2, tr + F.WL / 2), F.power(B[kr], tr - F.WL / 2, tr + F.WL / 2)
            if PA is not None and PB is not None:
                sa += PA[ib].sum(); sb += PB[ib].sum()
    return sb / sa


def estimate2(blocks, R, zw, comp):
    tabs = [R.table(tc) for tc in F.CENT]; ptabs = [ptable(R, tc) for tc in F.CENT]; c0sq = c0sq_near(blocks, tabs)
    rows = [rows_pz(A, B, off, tabs, ptabs, zw, c0sq, comp) for A, B, off in blocks]
    return dict(c0sq=c0sq, **F.summarise(rows)), rows


def run(comp="pz"):
    t0 = time.monotonic(); H, h = headers(); sel = (h["dline"] == 5081) & (np.abs(h["cmp"]) < F.CMP_LIM)
    key = lambda i: (h["ffid"][i], h["dpt"][i]); ia = {key(i): i for i in np.where(sel & (h["trid"] == 11))[0]}
    ib = {key(i): i for i in np.where(sel & (h["trid"] == 12))[0]}; pairs = sorted(set(ia) & set(ib))
    ka = np.array([ia[p] for p in pairs]); kb = np.array([ib[p] for p in pairs]); cmp = h["cmp"][ka]; off = h["off"][ka]
    zw = float(np.median(h["wd"][ka])); R = F.Rays(1.0, zw); tabs = [R.table(tc) for tc in F.CENT]; e = np.linspace(-F.CMP_LIM, F.CMP_LIM, F.NBLK + 1)
    print("pairs", len(pairs), "unpaired A", len(ia) - len(pairs), "B", len(ib) - len(pairs), "water depth", zw, flush=True)
    blocks = []
    for i in range(F.NBLK):
        m = (cmp >= e[i]) & (cmp < e[i + 1]); blocks.append((samples(H, ka[m]), samples(H, kb[m]), off[m]))
        print("block", i, m.sum(), "%.0fs" % (time.monotonic() - t0), flush=True)
    out, rows = estimate2(blocks, R, zw, comp); out.update(comp=comp, n_pairs=len(pairs), water_depth=zw)
    out["rows"] = [dict(r, block=i) for i, b in enumerate(rows) for r in b]
    (F.ROOT / "results" / ("field79raw%s.json" % ("" if comp == "pz" else "_" + comp))).write_text(json.dumps(out, indent=1))
    print(json.dumps({k: v for k, v in out.items() if k not in ("rows", "q_eff_by_window")}, indent=1))


def synth_block(ib):
    """as field79.synth_block, plus a first-order receiver ghost (sea surface, R = -1 for pressure) on component A,
    (1 + ghost) x cos(theta_w) on component B, and swell noise (2-8 Hz, 3 x the signal rms) on both."""
    H, h = F.read_headers(); sel = np.abs(h["cmp"]) < F.CMP_LIM; zw = float(np.median(h["wd"][sel])); R = F.Rays(1.0, zw)
    tk = np.arange(0.25, 2.95, 0.02); xs = np.arange(0, 4101, 25.0)
    grid = np.array([np.stack(R.table(t)(xs), 1) for t in tk]); pg = np.gradient(grid[:, :, 0], xs, axis=1)    # horizontal slowness
    tf = np.arange(0.25, 2.95, DT); fr = np.fft.rfftfreq(2 * NT, DT); W = (fr / 25.0) ** 2 * np.exp(-(fr / 25.0) ** 2)
    sw = ((fr >= 2) & (fr <= 8)).astype(float)
    k = F.blocks_of(h["cmp"])[ib]; A = np.zeros((len(k), NT), np.float32); B = np.zeros_like(A); rng = np.random.default_rng(SYN_SEED + 1000 + ib)
    for j, i in enumerate(k):
        b = int(np.round(h["cmp"][i] / 12.5)); rb = np.random.default_rng(SYN_SEED + b); r = (rb.random(tf.size) < 0.3) * rb.standard_normal(tf.size)
        x = h["off"][i]; ix = min(int(round(x / 25.0)), xs.size - 1); v = np.array([np.interp(tf, tk, grid[:, ix, q]) for q in range(4)])
        p = np.interp(tf, tk, pg[:, ix]); ok = np.isfinite(v[0]) & (r != 0); tt, la, lb, ls = v[:, ok]; tst = la / F.QS[0] + lb / F.QS[1] + (ls - la - lb) / F.QS[2]
        cw = np.sqrt(np.maximum(1 - (VW * p[ok]) ** 2, 0.05)); tau = 2 * zw * cw / VW
        U = r[ok][None] * np.exp(-np.pi * fr[:, None] * tst[None] - 2j * np.pi * fr[:, None] * tt[None]); gh = np.exp(-2j * np.pi * fr[:, None] * tau[None]) * GHOST
        for C, fac in ((A, (1 - gh) / (1 + R_SB * gh)), (B, (1 + gh) / (1 + R_SB * gh) * cw[None])):
            tr = np.fft.irfft((U * fac).sum(1) * W, 2 * NT)[:NT]; rms = np.sqrt(np.mean(tr ** 2))
            tr += NOISE * rms / 10 ** (24 / 20) * rng.standard_normal(NT)
            nsw = np.fft.irfft(np.fft.rfft(rng.standard_normal(2 * NT)) * sw, 2 * NT)[:NT]; tr += SWELL * 3 * rms * nsw / np.sqrt(np.mean(nsw ** 2))
            C[j] = tr
    allr = []
    print("synth block", ib, flush=True)
    return hp(A), hp(B), h["off"][k]


def synth():
    from multiprocessing import Pool
    with Pool(4) as pool:
        blocks = pool.map(synth_block, range(F.NBLK))
    H, h = F.read_headers(); sel = np.abs(h["cmp"]) < F.CMP_LIM; zw = float(np.median(h["wd"][sel])); R = F.Rays(1.0, zw)
    res = {}
    for comp in ("pz", "a", "b"):
        out, rows = estimate2(blocks, R, zw, comp)
        allr = [r for b in rows for r in b]; x = np.pi * np.array([r["dS"] for r in allr])
        yt = np.pi * np.array([r["dA"] / F.QS[0] + r["dB"] / F.QS[1] + (r["dS"] - r["dA"] - r["dB"]) / F.QS[2] for r in allr])
        tt = float(1000 * (x @ yt) / (x @ x)); out["truth_invQ_eff"] = tt
        out["selfcheck_pass"] = bool(abs(out["invQ_eff"] / tt - 1) <= 0.15 and out["ci90_eff"][0] <= tt <= out["ci90_eff"][1])
        res[comp] = {k: v for k, v in out.items() if k != "q_eff_by_window"}
        print(comp, json.dumps({k: res[comp][k] for k in ("invQ_eff", "ci90_eff", "truth_invQ_eff", "selfcheck_pass") + ("c0sq",)}), flush=True)
    (F.ROOT / "results" / "field79raw_synth.json").write_text(json.dumps(dict(seed=SYN_SEED, **res), indent=1))


if __name__ == "__main__":
    if sys.argv[1] == "synth":
        synth()
    elif sys.argv[1] == "run":
        run(sys.argv[2] if len(sys.argv) > 2 else "pz")
