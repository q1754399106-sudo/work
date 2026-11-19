"""UQ-W prior / rho ablation (PLAN_ZH.md): mechanism diagnosis on the used development data only.

Solver side (function `compensate`) receives the input data D, the model (cm with the model Q profile), the declared
Q multiplier c, the time axis, the anchor, the 1/snr field estimated from D, and, for VDB, an output wavelet taken from
D (semi-synthetic, sc11.prepare_sc) or the known source Ricker (FD, as stage 2b). It never sees the ideal section, the
reflectivity, the FD target or the dataset metadata; the harness forms c = m x c_eff_cum (semi) or c = m (FD), which
is the experimental condition (declared Q = m x effective true Q), exactly as in stages 0, 2b, W, MMR, CB and CBL.

Variants (all gamma = 0, 70 dB floor s, lognormal law on u / u_hat = exp(sig e), e ~ N(0, 1)):
  V0      current UQ-W: uqw.uqw (41-point Gauss-Hermite), rho = s + E|G|^2 / snr            sig 0.447
  V1      as V0 but rho = s + |G(u=1)|^2 / snr  (|G(u=1)|^2 = exp(-2 pi f x_hat))           sig 0.447
  V2      prior truncated to u / u_hat in [0.5, 2], renormalised, rho formula as V0           sig 0.447
  V3      prior truncated to u / u_hat in [0.1, 2], renormalised, rho formula as V0           sig 0.447
  V0s25   V0 at sig 0.25 ; V1s25  V1 at sig 0.25                                             (V4)
  V0d     quadrature control: V0 with the dense Gauss-Legendre rule of V2 / V3 over the whole Gauss-Hermite span
  W0      sig 0 (no prior spread; descriptive reference)
  VDB, IQ70 controls (stage 0 / stage 2b definitions)
V2, V3, V0d use Gauss-Legendre nodes in e on the stated interval with weights w_GL phi(e), renormalised to sum 1
(exact truncation bounds; the Gauss-Hermite nodes would put the [0.5, 2] bounds at 0.519 / 1.926 with mass 0.916
instead of 0.879). Everything else (x = x_hat + dx_hat (u - 1), log-sum-exp, row filters, irfft, M @ D) is uqw.uqw.

usage: ablate.py smoke       code-path and quadrature checks only (no scores)
       ablate.py run K N     worker K of N over 24 semi + 36 FD jobs -> results/raw/<job>.json
"""
import os, sys, json, time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]; BASE = ROOT.parent
sys.path.insert(0, str(BASE / "qunknown_comp_20261109" / "scripts"))
import uqw as U  # noqa: E402
from uqw import np, lib, semi, lib12, TW, G_IQ, GH, snr_rows, DT  # noqa: E402
from uqiq import load  # noqa: E402
from stage0 import DEV, score, sc11, SEL  # noqa: E402
import fdcheck as fd  # noqa: E402

NFFT = 1024; SIG = 0.447; SIG2 = 0.25; MS = (0.5, 0.7, 1.0, 1.4, 2.0); CELLS = 1 << 22    # row block size: CELLS / (nf x nodes)
N_TR, N_FULL = 256, 1024                                  # Gauss-Legendre nodes: truncated priors / V0d span
E_LO, E_HI = float(GH[0].min()), float(GH[0].max())       # Gauss-Hermite span (|e| <= 11.6)
TRUNC = {"V2": (0.5, 2.0), "V3": (0.1, 2.0)}
OUT = ROOT / "results" / "raw"


def gl_rule(lo, hi, n):
    """Nodes e in [lo, hi], log weights of the standard normal restricted there (renormalised)."""
    z, w = np.polynomial.legendre.leggauss(n); e = 0.5 * (hi - lo) * z + 0.5 * (hi + lo)
    lw = np.log(w) - 0.5 * e ** 2; lw -= lw.max(); lw -= np.log(np.exp(lw).sum()); return e, lw


def trunc_rule(sig, a, b, n=N_TR):
    return gl_rule(np.log(a) / sig, np.log(b) / sig, n)


def expect_nodes(x_hat, dx_hat, sig, f, lf, e, lw):
    """uqw.expectations (lognormal branch) with explicit nodes e and log weights lw. Returns A, C, P, lG2."""
    k = 2 * f * lf; x_hat = np.asarray(x_hat, float)[:, None]; dx_hat = np.asarray(dx_hat, float)[:, None]
    xs = (x_hat + dx_hat * (np.exp(sig * e) - 1)[None, :])[:, None, :]
    l2 = -2 * np.pi * f[None, :, None] * xs + lw; mx = l2.max(-1, keepdims=True)
    lG2 = mx[..., 0] + np.log(np.exp(l2 - mx).sum(-1))
    A = np.exp((-np.pi * f - 1j * k)[None, :, None] * xs + lw - lG2[..., None]).sum(-1)
    P = np.exp(-1j * k[None, :, None] * xs + lw).sum(-1)
    return A, np.exp(-lG2), P, lG2


def gh_rule():
    e, w = GH; return e, np.log(w / w.sum())


def rows(x_hat, dx_hat, sig, f, lf, s, inv, rule, rho):
    """gamma = 0 rows: H = E[G*] / (E|G|^2 + rho), divided through by E|G|^2.
    rho 'scaled': s + E|G|^2 inv (uqw.rows_w with gamma 0); 'fixed': s + |G(u=1)|^2 inv, |G(u=1)|^2 = exp(-2 pi f x_hat)."""
    e, lw = rule; out = np.empty((len(x_hat), len(f)), complex); blk = max(1, CELLS // (len(f) * len(e)))
    for a in range(0, len(x_hat), blk):
        b = slice(a, a + blk); A, C, P, lG2 = expect_nodes(x_hat[b], dx_hat[b], sig, f, lf, e, lw)
        if rho == "scaled":
            out[b] = A / (1 + s * C + inv[b])
        else:
            r = np.exp(-2 * np.pi * f[None, :] * np.asarray(x_hat[b], float)[:, None] - lG2)   # |G(u=1)|^2 / E|G|^2
            out[b] = A / (1 + s * C + r * inv[b])
    return out


def apply_rows(H, D, nfft=NFFT):
    n = len(H); lag = np.arange(n); M = np.empty((n, n))
    for i in range(n):
        h = H[i].copy(); h[0] = h[0].real; h[-1] = h[-1].real
        M[i] = np.fft.irfft(h, nfft)[(i - lag) % nfft]
    return M @ D


def uqw_var(D, cm, Q, c, sig, tw, t_anchor, inv, rule, rho):
    """As uqw.uqw (gamma 0) with an explicit quadrature rule and rho form."""
    Q = np.asarray(Q, float); f = np.fft.rfftfreq(NFFT, DT); s = semi.gain_sigma2(G_IQ)
    x_hat = lib.absorption(dict(cm, Q=Q * c), tw); dx_hat = x_hat - lib.absorption(dict(cm, Q=Q * c), [t_anchor])[0]
    lf = np.log(np.maximum(f, np.finfo(float).tiny) / float(cm["fref"]))
    return apply_rows(rows(x_hat, dx_hat, sig, f, lf, s, inv, rule, rho), D)


def compensate(D, cm, Q, c, tw, anchor, inv, w_out):
    """All variants and controls for one declared c. Inputs are data-side only (see module docstring)."""
    gh = gh_rule()
    out = dict(V0=U.uqw(D, cm, Q, c, SIG, tw, anchor, 0.0, inv_snr=inv),
               V1=uqw_var(D, cm, Q, c, SIG, tw, anchor, inv, gh, "fixed"),
               V2=uqw_var(D, cm, Q, c, SIG, tw, anchor, inv, trunc_rule(SIG, *TRUNC["V2"]), "scaled"),
               V3=uqw_var(D, cm, Q, c, SIG, tw, anchor, inv, trunc_rule(SIG, *TRUNC["V3"]), "scaled"),
               V0s25=U.uqw(D, cm, Q, c, SIG2, tw, anchor, 0.0, inv_snr=inv),
               V1s25=uqw_var(D, cm, Q, c, SIG2, tw, anchor, inv, gh, "fixed"),
               V0d=uqw_var(D, cm, Q, c, SIG, tw, anchor, inv, gl_rule(E_LO, E_HI, N_FULL), "scaled"),
               W0=U.uqw(D, cm, Q, c, 0.0, tw, anchor, 0.0, inv_snr=inv),
               VDB=lib12.tvw(fd.iq_t(D, cm, Q, None, c, 0, tw), w_out, SEL["B3"]["eps"], win=SEL["B3"]["win"]),
               IQ70=fd.iq_t(D, cm, Q, c, c, G_IQ, tw))
    return out


# ---------------- harness (scoring side: reads the truth) ----------------
F = np.fft.rfftfreq(NFFT, DT); WS = semi.win(TW, 2.3, 2.8)


def t2(X, R, fp):
    """uqw_improve_20261113/scripts/probe.py t2 (same statements): field-style tie, 2.3-2.8 s, every 5th trace."""
    return float(np.mean([lib.tie(X[:, j], semi.ricker_ref(R[:, j], fp), TW, WS)[0] for j in range(2, X.shape[1], 5)]))


def run_semi(did, cm, Q):
    arr, meta = load(did); D = arr["D"]; w_out = sc11.prepare_sc(D, cm)["wv"]["w_out"]; inv = snr_rows(D, TW, F)
    ideal, R, fp = arr["ideal"], arr["R"], float(cm["fpeak"]); rec = []
    for m in MS:
        for k, X in compensate(D, cm, Q, m * meta["c_eff_cum"], TW, 1.4, inv, w_out).items():
            rec.append(dict(method=k, m=m, S=score(X, ideal)["S"], T2=t2(X, R, fp), finite=bool(np.isfinite(X).all())))
    return dict(id=did, kind="semi", model=meta["structure"], c0=meta["c0"], r=meta["r"], level="semi", records=rec)


def run_fd(h, lv, r):
    D0, T, cm, tw, fp = fd.load_fd(h); i = fd.MODELS.index(h) + 1
    seed = 20261113000 + 1000 * i + 10 * list(fd.LEVELS).index(lv) + r                     # stage 2b / W seeds
    D = D0 + fd.noise(D0, fd.LEVELS[lv], seed); inv = snr_rows(D, tw, F); rec = []
    for m in MS:
        for k, X in compensate(D, cm, cm["Q"], m, tw, 0.0, inv, fd.ricker(fp)).items():
            rec.append(dict(method=k, m=m, S=fd.score_fd(X, T, tw)[0], T2=None, finite=bool(np.isfinite(X).all())))
    return dict(id=f"FD_{h}_{lv}_r{r}", kind="fd", model=h, c0=None, r=r, level=lv, seed=seed, records=rec)


def smoke():
    t0 = time.monotonic(); cm = semi.curtain_cm(); Q = np.asarray(cm["Q"], float); arr, meta = load("ZA_c0.75_r0"); D = arr["D"]
    inv = snr_rows(D, TW, F); ok = True
    for c in (0.5 * meta["c_eff_cum"], 2.0 * meta["c_eff_cum"]):
        ref = U.uqw(D, cm, Q, c, SIG, TW, 1.4, 0.0, inv_snr=inv); den = np.abs(ref).max()
        a = uqw_var(D, cm, Q, c, SIG, TW, 1.4, inv, gh_rule(), "scaled"); d = np.abs(a - ref).max() / den
        print("semi c %.3f  generic GH path vs uqw.uqw rel %.1e" % (c, d)); ok &= d <= 1e-10
        z = np.zeros_like(inv)                                     # with 1/snr = 0 the rho form cannot matter
        d = np.abs(uqw_var(D, cm, Q, c, SIG, TW, 1.4, z, gh_rule(), "fixed") - U.uqw(D, cm, Q, c, SIG, TW, 1.4, 0.0, inv_snr=z)).max() / den
        print("   V1 vs V0 at 1/snr = 0 rel %.1e" % d); ok &= d <= 1e-10
        d = np.abs(uqw_var(D, cm, Q, c, 0.0, TW, 1.4, inv, gh_rule(), "fixed") - U.uqw(D, cm, Q, c, 0.0, TW, 1.4, 0.0, inv_snr=inv)).max() / den
        print("   V1 vs V0 at sig = 0 rel %.1e" % d); ok &= d <= 1e-10
        for name, rule2 in (("V2", lambda n: trunc_rule(SIG, *TRUNC["V2"], n)), ("V3", lambda n: trunc_rule(SIG, *TRUNC["V3"], n)),
                            ("V0d", lambda n: gl_rule(E_LO, E_HI, n))):
            n1 = N_FULL if name == "V0d" else N_TR
            x1 = uqw_var(D, cm, Q, c, SIG, TW, 1.4, inv, rule2(n1), "scaled"); x2 = uqw_var(D, cm, Q, c, SIG, TW, 1.4, inv, rule2(2 * n1), "scaled")
            d = np.abs(x1 - x2).max() / np.abs(x2).max(); print("   %s nodes %d vs %d rel %.1e" % (name, n1, 2 * n1, d)); ok &= d <= 1e-6
        x0d = uqw_var(D, cm, Q, c, SIG, TW, 1.4, inv, gl_rule(E_LO, E_HI, N_FULL), "scaled")
        print("   V0d vs V0 (quadrature difference, output rel) %.1e" % (np.abs(x0d - ref).max() / den))
    e, lw = trunc_rule(SIG, *TRUNC["V2"]); u = np.exp(SIG * e); print("V2 nodes u %.3f..%.3f, weights sum %.12f" % (u.min(), u.max(), np.exp(lw).sum()))
    e, lw = trunc_rule(SIG, *TRUNC["V3"]); u = np.exp(SIG * e); print("V3 nodes u %.3f..%.3f" % (u.min(), u.max()))
    D0, T, cmf, tw, fp = fd.load_fd("H1"); Df = D0 + fd.noise(D0, 14.0, 1); invf = snr_rows(Df, tw, F); t1 = time.monotonic()
    outs = compensate(Df, cmf, cmf["Q"], 0.5, tw, 0.0, invf, fd.ricker(fp))
    print("FD H1 (smoke noise seed 1) m 0.5: all finite %s, %.0fs for one c" % (all(np.isfinite(x).all() for x in outs.values()), time.monotonic() - t1))
    d = np.abs(uqw_var(Df, cmf, cmf["Q"], 0.5, SIG, tw, 0.0, invf, gh_rule(), "scaled") - outs["V0"]).max() / np.abs(outs["V0"]).max()
    print("FD generic GH path vs uqw.uqw rel %.1e" % d); ok &= d <= 1e-10
    print("smoke %s  %.0fs" % ("OK" if ok else "FAIL", time.monotonic() - t0)); return ok


if __name__ == "__main__":
    if sys.argv[1] == "smoke":
        sys.exit(0 if smoke() else 1)
    k, n = int(sys.argv[2]), int(sys.argv[3]); OUT.mkdir(parents=True, exist_ok=True)
    cm = semi.curtain_cm(); Q = np.asarray(cm["Q"], float)
    jobs = [("semi", d) for d in DEV] + [("fd", (h, lv, r)) for h in fd.MODELS for lv in fd.LEVELS for r in fd.R]
    for j, (kind, key) in enumerate(jobs):
        name = key if kind == "semi" else "FD_%s_%s_r%d" % key; path = OUT / f"{name}.json"
        if j % n != k or path.exists():
            continue
        t0 = time.monotonic(); rec = run_semi(key, cm, Q) if kind == "semi" else run_fd(*key)
        rec["seconds"] = time.monotonic() - t0; path.write_text(json.dumps(rec, indent=1), encoding="utf-8")
        print(name, len(rec["records"]), f"{rec['seconds']:.0f}s", flush=True)
