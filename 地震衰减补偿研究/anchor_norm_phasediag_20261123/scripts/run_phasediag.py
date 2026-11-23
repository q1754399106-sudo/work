"""T3 Pre-registered Anchor Theory Normalization Phase Diagram (PLAN_ZH.md):
3 L1 normalization rules (R_med, R_fix, R_unit) x 3 wavelet depths (tw = 1.4, 2.0, 2.5 s)
x 3 Q-profiles (S, D, N3) on brand-new seeds 9931 & 9932.

Usage: python scripts/run_phasediag.py [--smoke]
Resumable via results/<job>.json + .claim. Low RAM (< 250 MB per worker, 1 BLAS thread).
"""
import os, sys, json, time
for _k in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_k, "1")

from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parents[1]
DATA_ROOT = Path(r"D:\ChatGPT\研究\地震衰减补偿研究")
ANCH_CODE = DATA_ROOT / "qscan_anchor_theory_20261030" / "scripts"
sys.path.insert(0, str(ANCH_CODE))
import qprofile as qp  # noqa: E402
from env15 import lib, semi, sc11, bm, gen  # noqa: E402

LAM = 0.5
EPS = 3e-3
TOL_CF = 1e-3
bm.MAXITER = 20000
ANCHORS = (1.0, 0.7, 1.4)
DVA = 1.0 / 0.7 - 1.0 / 1.4
VG = 1.0 / np.asarray(lib.QGRID, float)

PROFILES = {
    "S":  [30.0, 200.0, 200.0, 200.0],
    "D":  [250.0, 250.0, 40.0, 40.0],
    "N3": [100.0, 300.0, 40.0, 300.0],
}
qp.PROFILES.clear()
qp.PROFILES.update(PROFILES)

RULES = ("R_med", "R_fix", "R_unit")
TWS = (1.4, 2.0, 2.5)
OUT = HERE / "results"


def job_list():
    jobs = []
    # Core 3D phase diagram: seed 9931, 3 profiles x 3 rules x 3 tw = 27 cells
    for tw in TWS:
        for prof in ("S", "D", "N3"):
            for rule in RULES:
                jobs.append((prof, rule, tw, 9931))
    # Independent seed replication at tw=1.4: seed 9932, 3 profiles x 3 rules = 9 cells
    for prof in ("S", "D", "N3"):
        for rule in RULES:
            jobs.append((prof, rule, 1.4, 9932))
    return jobs


def freq(cm):
    N = int(cm["nfft"])
    f = np.arange(N // 2 + 1) / (N * float(cm["dt"]))
    b = -np.pi * f + 2j * f * np.log(np.maximum(f, np.finfo(float).tiny) / float(cm["fref"]))
    b[0] = 0
    return N, f, b


def gmat0(cm, aeff):
    t = bm.vec(cm["t"])
    n = len(t)
    N, _, b = freq(cm)
    E = np.exp(np.outer(b, aeff)) * np.exp(-2j * np.pi * np.outer(np.arange(N // 2 + 1), np.arange(n)) / N)
    E[-1] = E[-1].real
    return np.fft.irfft(E, n=N, axis=0)[:n]


def wavelets_tw(D, cm, va, t_w, raw=False, nfft=512):
    """Anchor wavelet extracted from [t_w - 0.20, t_w + 0.20] s using absorption a(t_w) * va."""
    Q = np.asarray(cm["Q"], float)
    c_eff = dict(cm, Q=Q / va)
    f, snr_s, P_s = lib.coherence_snr(D, semi.TW, t_w - 0.20, t_w + 0.20, nfft)
    _, snr_d, _ = lib.coherence_snr(D, semi.TW, 2.00, 3.00, nfft)
    amp = np.sqrt(P_s) / np.exp(-np.pi * f * lib.absorption(c_eff, [t_w])[0])
    fpk = float(f[np.argmax(np.sqrt(P_s) * (f >= 5))])
    fs = lib.snr_band(f, snr_s, fpk)
    fd = lib.snr_band(f, snr_d, fpk)

    def cut(fc):
        g = np.ones_like(f)
        r = (f > fc - 5) & (f < fc)
        g[r] = 0.5 * (1 + np.cos(np.pi * (f[r] - (fc - 5)) / 5))
        g[f >= fc] = 0
        g[f < 3] = 0
        return g

    a_in = amp * cut(fs)
    if raw:
        return np.roll(np.fft.irfft(a_in, nfft), 100)[:201]

    a_out = a_in * cut(min(fd, fs))
    def zp(a):
        w = np.fft.irfft(a, nfft)
        w = np.roll(w, 100)[:201]
        return w / np.max(np.abs(w))

    w_in, w_out = zp(a_in), zp(a_out)
    w_out *= np.max(np.abs(np.fft.rfft(w_in, nfft))) / np.max(np.abs(np.fft.rfft(w_out, nfft)))
    return dict(w_in=w_in, w_out=w_out, fmax_win=float(fs), fmax_target=float(fd))


def unnorm_dict(cm, v, w_in):
    a = lib.absorption(cm, bm.vec(cm["t"]))
    G = gmat0(cm, a * v)
    Win = lib.toeplitz_wavelet(w_in, len(a))
    return G @ Win


def pen_vec_y(Bt, rule, s_ref):
    """Per-column L1 penalty vector p_k(v, va) on unnormalized y = X / scale (excluding sigma factor)."""
    n = Bt.shape[0]
    s = np.linalg.norm(Bt, axis=0)
    base = LAM * np.sqrt(2.0 * np.log(n))
    if rule == "R_med":
        return np.full(n, base * float(np.median(s)))
    elif rule == "R_fix":
        return np.full(n, base * float(s_ref))
    elif rule == "R_unit":
        return base * s
    raise ValueError(rule)


def solve_rule(D, Bt, pp, rule, s_ref, start_y=None, tol=None):
    """Solve L1 problem under `rule` using bm.solve via exact algebraic mapping; return (y, J, info)."""
    n = D.shape[0]
    scale = np.linalg.norm(Bt, axis=0)
    B = Bt / scale
    if rule == "R_med":
        a_eff = LAM
        sc_arg = scale
    elif rule == "R_fix":
        a_eff = LAM * (float(s_ref) / float(np.median(scale)))
        sc_arg = scale
    elif rule == "R_unit":
        a_eff = LAM
        sc_arg = np.ones(n, dtype=float)
    else:
        raise ValueError(rule)

    old_tol = bm.TOL
    if tol is not None:
        bm.TOL = tol
    try:
        X, info = bm.solve(
            B, np.eye(n), D, pp["sigma"], sc_arg, pp["K"], a_eff, 0.0, pp["live"],
            start=None if start_y is None else start_y * scale[:, None],
        )
    finally:
        bm.TOL = old_tol

    y = X / scale[:, None]
    sigma, live = pp["sigma"], pp["live"]
    s2 = np.where(live, sigma ** 2, np.inf)
    p_y = pen_vec_y(Bt, rule, s_ref)
    res = D - B @ X
    J = float(np.sum((0.5 * np.sum(res ** 2, axis=0) + np.sum((p_y[:, None] * sigma[None, :]) * np.abs(y), axis=0))[live] / s2[live]))
    return y, J, dict(converged=bool(info["converged"]), iterations=int(info["iterations"]), kkt=float(info["kkt_ratio"]))


def polish_vec(D, Bt, y, pp, p_y):
    """Exact fixed-active-set polish for arbitrary per-column penalty vector p_y."""
    live = np.where(pp["live"])[0]
    sigma = pp["sigma"]
    A, S, Y = {}, {}, np.zeros_like(y)
    bad_sign = 0
    nact = 0
    kkt = []
    for j in live:
        Aj = np.flatnonzero(y[:, j] != 0)
        s = np.sign(y[Aj, j])
        BA = Bt[:, Aj]
        if len(Aj):
            yA = np.linalg.solve(BA.T @ BA, BA.T @ D[:, j] - sigma[j] * (p_y[Aj] * s))
            bad_sign += int(np.sum(np.sign(yA) != s))
            nact += len(Aj)
            Y[Aj, j] = yA
            r = D[:, j] - BA @ yA
        else:
            r = D[:, j]
        A[j], S[j] = Aj, s
        g = np.abs(Bt.T @ r) / np.maximum(sigma[j] * p_y, 1e-30)
        g[Aj] = 0.0
        kkt.append(float(g.max()))
    return Y, A, S, dict(n_active=nact, bad_sign=bad_sign, kkt_inactive_max=float(np.max(kkt)))


def jfix_vec(D, Bt, A, S, pp, p_y):
    """Universal exact closed-form fixed-active-set objective J for R_med, R_fix, and R_unit."""
    J = 0.0
    for j, Aj in A.items():
        d = D[:, j]
        sg = pp["sigma"][j]
        if len(Aj):
            BA = Bt[:, Aj]
            c = BA.T @ d - sg * (p_y[Aj] * S[j])
            F = 0.5 * float(d @ d) - 0.5 * float(c @ np.linalg.solve(BA.T @ BA, c))
        else:
            F = 0.5 * float(d @ d)
        J += F / (sg ** 2)
    return float(J)


def scan_anchor_rule(D, cm, pp, ca, rule, t_w, s_ref):
    """Fixed-anchor QScan under normalization rule `rule` and wavelet window `t_w`."""
    wv = wavelets_tw(D, cm, 1.0 / ca, t_w, raw=False)
    J = np.full(len(lib.QGRID), np.nan)
    C = np.zeros(J.shape, bool)
    prev_y = None
    for j in range(len(lib.QGRID) - 1, -1, -1):
        c = lib.QGRID[j]
        if np.isfinite(J).any() and j <= int(np.nanargmin(J)) - 3:
            break
        Bt = unnorm_dict(cm, 1.0 / c, wv["w_in"])
        y, Jv, info = solve_rule(D, Bt, pp, rule, s_ref, start_y=prev_y)
        prev_y = y
        J[j] = Jv
        C[j] = info["converged"]
    jm = int(np.nanargmin(J))
    edge = jm in (0, len(J) - 1)
    v = float(VG[jm])
    if not edge:
        idx = [jm - 1, jm, jm + 1]
        if all(np.isfinite(J[idx])):
            p2 = np.polyfit(VG[idx], J[idx], 2)
            if p2[0] > 0:
                v = float(np.clip(-p2[1] / (2 * p2[0]), VG[jm + 1], VG[jm - 1]))
    trust = bool(all(C[max(0, jm - 1):jm + 2]))
    return dict(v_grid=float(VG[jm]), v_sel=v, edge=bool(edge), neighbourhood_converged=trust,
                J_min=float(J[jm]))


def closed_form_rule(D, cm, pp, v0, rule, t_w, s_ref):
    """Single-solve universal fixed-active-set closed-form alpha_F for R_med, R_fix, and R_unit."""
    t = bm.vec(cm["t"]).astype(float)
    a = lib.absorption(cm, t)
    a_tw = float(lib.absorption(cm, [t_w])[0])
    # Use raw (unnormalized) wavelet so WR[va] is smooth in va; note peak-normalized w_in = WR[va] / max|WR[va]|.
    # For R_med and R_unit the problem is strictly invariant to scalar rescaling of w_in; for R_fix (where scan_anchor_rule
    # fixes s_ref relative to peak-normalized w_in), rescaling into raw coordinates multiplies s_ref by max|WR[va]|.
    WR = {va: wavelets_tw(D, cm, va, t_w, raw=True) for va in (1.0 - EPS, 1.0, 1.0 + EPS)}
    WMAX = {va: float(np.max(np.abs(WR[va]))) for va in WR}

    Bt0 = unnorm_dict(cm, v0, WR[1.0])
    s_ref_0 = s_ref * WMAX[1.0]
    y0, J0, info0 = solve_rule(D, Bt0, pp, rule, s_ref_0, tol=TOL_CF)
    p0 = pen_vec_y(Bt0, rule, s_ref_0)
    Y, A, S, pinfo = polish_vec(D, Bt0, y0, pp, p0)

    def JF(dv=0.0, dva=0.0):
        va = 1.0 + dva
        Bt = unnorm_dict(cm, v0 + dv, WR[va])
        p_y = pen_vec_y(Bt, rule, s_ref * WMAX[va])
        return jfix_vec(D, Bt, A, S, pp, p_y)

    Jc = JF(0.0, 0.0)
    e = EPS
    Jvv = (JF(e, 0.0) - 2.0 * Jc + JF(-e, 0.0)) / (e ** 2)
    Jvx = (JF(e, e) - JF(e, -e) - JF(-e, e) + JF(-e, -e)) / (4.0 * e * e)
    alpha_F = float(-Jvx / Jvv) if abs(Jvv) > 1e-30 else float("nan")
    rel_closure = abs(Jc - J0) / max(abs(J0), 1e-30)

    kmed = int(np.argsort(np.linalg.norm(Bt0, axis=0))[len(t) // 2])
    U = float(a_tw * np.mean(a) / np.mean(a ** 2))
    med_ratio = float(a_tw / a[kmed])
    return dict(
        v0=v0, a_tw=a_tw, J0=J0, Jc=Jc, rel_closure=float(rel_closure),
        Jvv=float(Jvv), Jvx=float(Jvx), alpha_F=alpha_F,
        t_med=float(t[kmed]), U=U, med_ratio=med_ratio,
        polish=pinfo, solve0=info0,
    )


def run_one(prof, rule, t_w, seed, base):
    t0 = time.monotonic()
    ctx = qp.context(base, prof)
    g = gen.signal(ctx, 1.0, seed)
    D = g["x"] + qp.MULT * g["level"] * ctx["noise"][0]
    cm = ctx["cm"]
    pp = sc11.prepare_sc(D, cm)

    # Fixed reference column norm at baseline (v=1.0, ca=1.0) for R_fix
    wv_ref = wavelets_tw(D, cm, 1.0, t_w, raw=False)
    Bt_ref = unnorm_dict(cm, 1.0, wv_ref["w_in"])
    s_ref = float(np.median(np.linalg.norm(Bt_ref, axis=0)))

    scans = {}
    for ca in ANCHORS:
        scans[f"a{ca}"] = scan_anchor_rule(D, cm, pp, ca, rule, t_w, s_ref)
    t_scans = time.monotonic() - t0

    v0 = scans["a1.0"]["v_sel"]
    cf = closed_form_rule(D, cm, pp, v0, rule, t_w, s_ref)
    khat = float((scans["a0.7"]["v_sel"] - scans["a1.4"]["v_sel"]) / DVA)
    valid = bool(all(not x["edge"] and x["neighbourhood_converged"] for x in scans.values()) and cf["solve0"]["converged"])
    return dict(
        profile=prof, rule=rule, tw=t_w, seed=seed, c_true=1.0, lam=LAM,
        Q_profile=PROFILES[prof], s_ref=s_ref,
        scans=scans, v_sel_10=v0, kappa_hat=khat, valid=valid,
        seconds_scans=round(t_scans, 1), seconds=round(time.monotonic() - t0, 1),
        **cf,
    )


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    base = gen.prepare(qp.STRUCT, qp.NOISE)
    if "--smoke" in sys.argv:
        for rule in RULES:
            rec = run_one("D", rule, 1.4, 9919, base)
            print(f"SMOKE D {rule:6s} tw=1.4: v10={rec['v_sel_10']:.4f} khat={rec['kappa_hat']:.4f} aF={rec['alpha_F']:.4f} diff={rec['alpha_F']-rec['kappa_hat']:+.4f} clos={rec['rel_closure']:.2e} valid={rec['valid']} ({rec['seconds']}s)", flush=True)
        return

    for prof, rule, tw, seed in job_list():
        job = f"{prof}_{rule}_tw{tw:.1f}_s{seed}"
        path = OUT / f"{job}.json"
        if path.exists():
            continue
        claim = OUT / f"{job}.claim"
        try:
            os.close(os.open(claim, os.O_CREAT | os.O_EXCL | os.O_WRONLY))
        except FileExistsError:
            continue
        try:
            rec = run_one(prof, rule, tw, seed, base)
            path.write_text(json.dumps(rec, indent=1), encoding="utf-8")
            print(f"{job:26s} {'OK' if rec['valid'] else 'EDGE/NC':7s} v10={rec['v_sel_10']:.4f} khat={rec['kappa_hat']:.4f} aF={rec['alpha_F']:.4f} diff={rec['alpha_F']-rec['kappa_hat']:+.4f} clos={rec['rel_closure']:.1e} ({rec['seconds']}s)", flush=True)
        finally:
            if claim.exists():
                claim.unlink()


if __name__ == "__main__":
    main()
