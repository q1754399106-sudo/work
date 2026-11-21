"""Step 0 of direction A (PLAN_ZH.md, frozen in e3df339): closed-form optimum of the sat -> fs design slope g against
Monte Carlo with real (Volve log), thin-bed and sparse-random reflectivity spectra. Amplitude only, noise free.

  ln Gamma(x; gamma) = ln int p_x^gamma(f) rho_d(f) exp[-pi (u - 1) a(x) f] df,   p_x^gamma(f) = r^g p0(r^g f),  r = a(x)/a(0)
  dGA = least-squares slope of ln Gamma over s = sin^2(theta) (11 offsets)

usage: tradeoff.py run [nproc]  -> results/tradeoff.json ;  tradeoff.py wells -> results/wells.json (window counts only)"""
import os, sys, json, re, time
from pathlib import Path
os.environ.setdefault("OMP_NUM_THREADS", "2"); os.environ.setdefault("OPENBLAS_NUM_THREADS", "2")
import numpy as np

HERE = Path(__file__).resolve().parent; ROOT = HERE.parent; RES = ROOT / "results"
LEG = ROOT.parent / "offset_q_avo_20261118" / "scripts" / "legacy2b"
sys.path.insert(0, str(LEG))
import gate1 as g1  # noqa: E402  (geometry, a0_ax, SIN2; unchanged)

# Volve logs: the experiment-9 inventory (local cache, not in git); same files as multiwell_tie_volve_20261028
E9 = Path(os.environ.get("E9", r"D:\ChatGPT\研究\地震衰减补偿研究\baselines_volve_20261029\cache\inputs\multiwell_tie_volve_20261028"))
TIE_WELLS = ["F-1", "F-1_A", "F-11_A", "F-11_T2", "F-12", "F-14", "F-15", "F-15_A", "F-15_B", "F-4", "F-5"]
GAP_MAX, GAP_TIE, MIN_SPAN, FINE = 40.0, 5.0, 0.15, 0.001

MODELS = {"M1": [40, 50, 70, 100], "M12": [55, 40, 100, 70], "M24": [42, 70, 55, 95], "M25": [75, 45, 30, 65],
          "M26": [35, 40, 85, 50], "M27": [95, 60, 45, 35]}
FPS = (40.0, 60.0); SIGS = (0.447, 0.09); SOURCES = ("A", "B", "C", "L")   # L = log-linear algebra self-check
GAM = np.round(np.arange(0, 1.5001, 0.025), 4); NDRAW = 20000; NBOOT = 200; SEED = 62_000_000
DF = 0.25; DTR = 0.001; NFFT = int(round(1 / (DF * DTR)))      # 4000-point FFT on the 1 ms grid -> 0.25 Hz bins
EV = np.round(g1.C.event_t / g1.DT).astype(int)


# ------------------------------------------------------------------ wells (copied from multiwell_tie_volve_20261028 wq.py / tie.py)
def picks():
    out, spans = {}, None

    def fl(x):
        try:
            return float(x)
        except ValueError:
            return np.nan
    for line in (E9 / "raw" / "picks" / "Well_picks_Volve_v1.dat").read_text(encoding="latin-1").splitlines():
        if line.strip().startswith("-----"):
            spans = [(m.start(), m.end()) for m in re.finditer(r"-+", line)]; continue
        if spans is None or not line.strip() or line.startswith("Well NO") or line.lstrip().startswith("Well name"):
            continue
        cols = [line[a:b + 1].strip() if i < len(spans) - 1 else line[a:].strip() for i, (a, b) in enumerate(spans)]
        if len(cols) < 12:
            continue
        md, tvd, twt, e, n = fl(cols[4]), fl(cols[5]), fl(cols[7]) / 1000, fl(cols[10]), fl(cols[11])
        if np.isfinite(md) and np.isfinite(tvd):
            out.setdefault(cols[0].replace("NO ", ""), []).append((cols[1], md, tvd, twt, e, n))
    return out


def running_median(x, k=5):
    h = k // 2; xp = np.pad(x, h, mode="edge")
    return np.median(np.stack([xp[i:i + len(x)] for i in range(k)]), axis=0)


def interp_x(x, xp, fp):
    y = np.interp(x, xp, fp); lo, hi = x < xp[0], x > xp[-1]
    y[lo] = fp[0] + (x[lo] - xp[0]) * (fp[1] - fp[0]) / (xp[1] - xp[0]); y[hi] = fp[-1] + (x[hi] - xp[-1]) * (fp[-1] - fp[-2]) / (xp[-1] - xp[-2])
    return y


def well(name, P):
    A = np.loadtxt(E9 / "inventory" / "logs" / f"15_9-{name}.csv", delimiter=",", skiprows=1)
    md, ac = A[:, 0], A[:, 1]; ok = np.isfinite(ac); md, ac = md[ok], ac[ok]
    pieces = np.split(np.arange(len(md)), np.where(np.diff(md) > GAP_MAX)[0] + 1)
    s = max(pieces, key=lambda p: md[p[-1]] - md[p[0]]); md, ac = md[s], running_median(ac[s])
    pk = sorted(P["15/9-" + name.replace("_", " ")], key=lambda r: r[1])
    pm = np.array([r[1] for r in pk]); pv = np.array([r[2] for r in pk]); pt = np.array([r[3] for r in pk])
    _, u = np.unique(pm, return_index=True); pm, pv, pt = pm[u], pv[u], pt[u]
    tvd = interp_x(md, pm, pv); v = 0.3048 / (ac * 1e-6)
    ts = np.r_[0, np.cumsum(2 * np.diff(tvd) / (0.5 * (v[1:] + v[:-1])))]
    a = np.isfinite(pt) & (pm >= md[0]) & (pm <= md[-1]); km, kt = pm[a], pt[a]
    ks = np.interp(km, md, ts); t = ts + np.interp(ts, ks, kt - ks)
    return dict(md=md, t=t, anchors=np.c_[km, kt])


def tie_log(name, P):
    W = well(name, P); am, at = W["anchors"][:, 0], W["anchors"][:, 1]
    A = np.loadtxt(E9 / "inventory" / "logs" / f"15_9-{name}.csv", delimiter=",", skiprows=1)
    md, ac, den = A.T; ok = np.isfinite(ac) & np.isfinite(den); md, ac, den = md[ok], ac[ok], den[ok]
    pieces = np.split(np.arange(len(md)), np.where(np.diff(md) > GAP_TIE)[0] + 1)
    ov = [max(0.0, min(md[p[-1]], am.max()) - max(md[p[0]], am.min())) for p in pieces]
    s = pieces[int(np.argmax(ov))]; md, ac, den = md[s], running_median(ac[s]), running_median(den[s])
    t = np.interp(md, W["md"], W["t"]); lo, hi = max(t[0], at.min()), min(t[-1], at.max())
    tt = np.round(np.arange(np.ceil(lo / FINE), np.floor(hi / FINE) + 1) * FINE, 6)
    if len(tt) < 3:
        return None
    z = np.log(np.interp(tt, t, 0.3048 / (ac * 1e-6) * den))
    return dict(name=name, rc=0.5 * np.diff(z), t0=float(tt[0]), t1=float(tt[-1]))


def log_windows(fp):
    """Hann windows of 6/fp s, step half a window, over every anchored span >= MIN_SPAN; |FFT|^2 on the 0.25 Hz grid."""
    P = picks(); out, info = [], {}
    for n in TIE_WELLS:
        try:
            g = tie_log(n, P)
        except Exception as e:  # noqa: BLE001
            info[n] = "skip: %s" % e; continue
        if g is None or g["t1"] - g["t0"] < MIN_SPAN:
            info[n] = "span < %.2f s" % MIN_SPAN; continue
        L = int(round(6 / fp / DTR)); st = L // 2; r = g["rc"]; k = 0
        for i0 in range(0, len(r) - L + 1, st):
            out.append(np.abs(np.fft.rfft(r[i0:i0 + L] * np.hanning(L), NFFT)) ** 2); k += 1
        info[n] = dict(span_s=round(g["t1"] - g["t0"], 3), windows=k)
    return np.array(out), info


def synth_windows(kind, fp, nd, rng):
    L = int(round(6 / fp / DTR)); w = np.hanning(L); c = L // 2; X = np.zeros((nd, L))
    if kind == "B":
        q = rng.uniform(-1, 1, nd); tau = rng.uniform(0.002, 0.040, nd); i1 = np.round(c - tau / 2 / DTR).astype(int); i2 = np.round(c + tau / 2 / DTR).astype(int)
        i2 = np.maximum(i2, i1 + 1); X[np.arange(nd), i1] = 1.0; X[np.arange(nd), i2] = q
    else:  # C: 30 % of a 2 ms grid, N(0, 1)
        on = rng.random((nd, L // 2)) < 0.3; X[:, ::2][:, :L // 2] = on * rng.standard_normal((nd, L // 2))
    return np.abs(np.fft.rfft(X * w, NFFT, axis=1)) ** 2


# ------------------------------------------------------------------ design and readout
def p0(f, fp):
    ric = f ** 2 * np.exp(-(f / fp) ** 2)  # Ricker amplitude spectrum (shape)
    bp = np.clip((f - 10.0 + 5) / 5, 0, 1) * np.clip((fp + 8 - f) / 8, 0, 1)
    return (ric * bp) ** 2


def setup(model, fp, j):
    Q = np.array(MODELS[model], float); a0, ax = g1.a0_ax(Q); k0 = EV[j]
    a = a0[k0] + ax[k0]; s = g1.SIN2[k0]; r = a / a[0]
    F = np.arange(1, int(round((fp + 8) / DF)) + 1) * DF
    P = np.stack([np.stack([rr ** gm * p0(rr ** gm * F, fp) for rr in r]) for gm in GAM])  # (G, X, F)
    P /= P.sum(axis=2, keepdims=True)
    fbar = (P * F).sum(axis=2)                                                             # (G, X)
    sc = s - s.mean(); slope = lambda y: (y * sc).sum(-1) / (sc ** 2).sum()                 # noqa: E731
    h = slope(np.log(a)); gg = slope(fbar); f0 = fbar[0, 0]
    return dict(a=a, s=s, sc=sc, F=F, P=P, h=float(h), g=gg, f0=float(f0), a0=float(a[0]))


def dga(S, rho, eps):
    """rho (N, F), eps (N,) -> dGA (G, N)."""
    u1 = np.exp(eps) - 1; out = np.zeros((len(GAM), len(eps), len(S["a"])))
    for i, ai in enumerate(S["a"]):
        RE = rho * np.exp(-np.pi * np.outer(u1 * ai, S["F"])); out[:, :, i] = (RE @ S["P"][:, i, :].T).T
    L = np.log(out); return (L * S["sc"]).sum(-1) / (S["sc"] ** 2).sum()


def beta_of(rho, S):
    w = S["P"][0, 0]; F = S["F"]; mf = (w * F).sum(); vf = (w * (F - mf) ** 2).sum()
    Er = rho @ w; cov = rho @ (w * (F - mf)); b = cov / (vf * Er)
    lr = np.log(np.maximum(rho, 1e-12 * rho.max(axis=1, keepdims=True))); bl = (lr @ (w * (F - mf))) / vf   # weighted LS slope of ln rho
    return b, bl


def vertex(m):
    i = int(np.argmin(m)); lo, hi = max(0, i - 5), min(len(GAM), i + 6)
    if hi - lo < 3:
        return float(GAM[i])
    c = np.polyfit(GAM[lo:hi], m[lo:hi], 2)
    v = -c[1] / (2 * c[0]) if c[0] > 0 else GAM[i]
    return float(np.clip(v, GAM[lo], GAM[hi - 1]))


def cell(args):
    src, model, fp, j, sig, spec = args
    S = setup(model, fp, j); F = S["F"]; rng = np.random.default_rng(SEED + 1000 * list(MODELS).index(model) + 100 * int(fp) // 10 + 10 * j + SIGS.index(sig) + SOURCES.index(src) * 7)
    nb = F.size
    if src == "L":
        betas = spec  # sigma_beta from source A, log-linear spectra
        bd = rng.normal(0, betas, NDRAW); rho = np.exp(np.outer(bd, F - S["f0"]))
    else:
        if spec is None:  # B, C: same population per fp in every cell, generated in the worker
            spec = synth_windows(src, fp, NDRAW, np.random.default_rng(SEED + (500 if src == "C" else 0) + int(fp)))
        R = spec[:, 1:nb + 1]
        rho = R[rng.integers(0, len(R), NDRAW)] if src == "A" else R
        rho = rho / rho.mean(axis=1, keepdims=True)
    eps = rng.normal(0, sig, NDRAW)
    # k from the design alone (beta = 0)
    one = np.ones((2, nb)); dq = dga(S, one, np.array([0.01, -0.01])); SQ = (dq[:, 0] - dq[:, 1]) / 0.02
    i5 = int(np.where(np.isclose(GAM, 0.5))[0][0]); k = float(-(SQ[i5] - SQ[0]) / (np.pi * S["a0"] * S["g"][i5]))
    hf0_first = S["h"] * S["f0"]; hf0_sens = float(-SQ[0] / (np.pi * S["a0"]))
    b, bl = beta_of(rho, S); sb2 = float(np.mean(b ** 2)); lam = np.pi ** 2 * S["a0"] ** 2 * sig ** 2 / sb2
    gt_cf = k ** 2 * lam / (1 + k ** 2 * lam); ratio_cf = 1 / (1 + k ** 2 * lam)
    D = dga(S, rho, eps); M = np.mean(D ** 2, axis=1)
    gtil = -k * S["g"] / (S["h"] * S["f0"])                     # fraction of the way to fs, per family gamma
    gmc = vertex(M); gt_mc = float(np.interp(gmc, GAM, gtil)); Mi = lambda x: float(np.interp(x, GAM, M))  # noqa: E731
    ratio_mc = Mi(gmc) / M[0]
    g_cf = float(np.interp(gt_cf, gtil, GAM)) if np.all(np.diff(gtil) > 0) else float("nan")
    red_cf = 1 - np.sqrt(Mi(g_cf) / M[0])
    bs = []
    br = np.random.default_rng(SEED + 7)
    for _ in range(NBOOT):
        ii = br.integers(0, NDRAW, NDRAW); bs.append(np.interp(vertex(np.mean(D[:, ii] ** 2, axis=1)), GAM, gtil))
    return dict(src=src, model=model, fp=fp, event=j, sig=sig, a0=S["a0"], h=S["h"], f0=S["f0"], k=k, hf0_first=hf0_first, hf0_sens=hf0_sens,
                g_fs=float(np.interp(1.0, GAM, S["g"])), sigma_beta=float(np.sqrt(sb2)), beta_mean=float(np.mean(b)), beta_sd=float(np.std(b)),
                beta_loglin_rms=float(np.sqrt(np.mean(bl ** 2))), lam=float(lam), gt_cf=float(gt_cf), gt_mc=gt_mc, gt_mc_se=float(np.std(bs)),
                ratio_cf=float(ratio_cf), gt_cf_H=float(gt_cf * hf0_sens / hf0_first), ratio_mc=float(ratio_mc), red_cf_in_mc=float(red_cf), red_mc=float(1 - np.sqrt(ratio_mc)),
                gamma_cf=g_cf, gamma_mc=gmc, rmse_sat=float(np.sqrt(M[0])), rmse_fs=float(np.sqrt(Mi(1.0))), gtil_at_gamma1=float(np.interp(1.0, GAM, gtil)))


def main(nproc):
    from multiprocessing import Pool
    t0 = time.monotonic(); specs, info = {}, {}
    for fp in FPS:
        W, inf = log_windows(fp); specs[("A", fp)] = W; info[str(fp)] = inf
        specs[("A", fp)] = W[:, :int(round((fp + 8) / DF)) + 2]; specs[("B", fp)] = None; specs[("C", fp)] = None
    print("log windows %.0fs" % (time.monotonic() - t0), {str(k): len(v) for k, v in specs.items() if v is not None}, flush=True)
    jobs = [(src, m, fp, j, sig, specs[(src, fp)]) for src in "ABC" for m in MODELS for fp in FPS for j in range(3) for sig in SIGS]
    with Pool(nproc) as pool:
        rows = pool.map(cell, jobs, chunksize=2)
    print("A/B/C %.0fs" % (time.monotonic() - t0), flush=True)
    sbA = {(r["model"], r["fp"], r["event"], r["sig"]): r["sigma_beta"] for r in rows if r["src"] == "A"}
    jobsL = [("L", m, fp, j, sig, sbA[(m, fp, j, sig)]) for m in MODELS for fp in FPS for j in range(3) for sig in SIGS]
    with Pool(nproc) as pool:
        rows += pool.map(cell, jobsL, chunksize=2)
    RES.mkdir(exist_ok=True)
    (RES / "tradeoff.json").write_text(json.dumps(dict(rows=rows, wells=info, ndraw=NDRAW), indent=0), encoding="utf-8")
    print("done %.0fs" % (time.monotonic() - t0), flush=True)


if __name__ == "__main__":
    if sys.argv[1] == "run":
        main(int(sys.argv[2]) if len(sys.argv) > 2 else 6)
    elif sys.argv[1] == "wells":
        for fp in FPS:
            W, inf = log_windows(fp); print(fp, len(W), json.dumps(inf))
