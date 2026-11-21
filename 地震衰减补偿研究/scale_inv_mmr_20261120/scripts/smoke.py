"""Smoke checks of si.py (no scores are computed): equivalences, monotone descent, ill-posedness of the original
form, timing. usage: smoke.py"""
import sys, time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import si  # noqa: E402
from si import np, semi, DT, NFFT, setup, cells, snr_rows, GH, SIG, U1  # noqa: E402
from mmr import h_of  # noqa: E402
from uqiq import TW, load  # noqa: E402
import confirmW as CW  # noqa: E402
import fdcheck as fd  # noqa: E402

t0 = time.monotonic(); cm = semi.curtain_cm(); Q = np.asarray(cm["Q"], float)
arr, meta = load("ZA_c0.75_r0"); D = arr["D"]; f = np.fft.rfftfreq(NFFT, DT); inv = snr_rows(D, TW, f)
w_out = CW.sc11.prepare_sc(D, cm)["wv"]["w_out"]; w = si.wspec(w_out); c = meta["c_eff_cum"]
print("prepare %.0fs, w_out len %d, w peak at %.1f Hz" % (time.monotonic() - t0, len(w_out), f[np.argmax(w)]))
fl, lf, xh, dx = setup(cm, Q, c, TW, 1.4); s = semi.gain_sigma2(si.G_IQ); rows = slice(200, 264)

# (1) one scenario: the scale-invariant optimum is the Wiener filter up to a row gain
R1 = si.Rows(xh[rows, None], fl, lf, s, inv[rows], w); H1, al1, i1 = R1.solve(np.ones(1))
Hw = R1.Gc[..., 0] / R1.den[..., 0]
cos = np.abs((H1 * np.conj(Hw)).sum(1)) / np.linalg.norm(H1, axis=1) / np.linalg.norm(Hw, axis=1)
print("(1) J=1: 1 - |cos(H, Wiener)| max %.1e, alternation its %d, regret max %.1e" % (1 - cos.min(), i1["its"], R1.regret(H1)[0].max()))

# (2) alpha = 1: first h-step is the Bayes form (UQ-W cells)
xs, pi = si.scenarios(xh[rows], dx[rows], "SIU"); RU = si.Rows(xs, fl, lf, s, inv[rows], w)
hb = h_of(pi, RU.Gc, RU.den); h0 = RU.h(pi, np.ones((xs.shape[0], len(pi))))
print("(2) alpha=1 h-step vs Bayes form rel %.1e" % (np.abs(h0 - hb).max() / np.abs(hb).max()))

# (3) monotone descent of the fixed-prior objective
al = np.ones((xs.shape[0], len(pi))); objs = []
for k in range(30):
    H = RU.h(pi, al); q, P = RU.qp(H); al = np.maximum(q, 0) / P; objs.append(((RU.E0 - np.maximum(q, 0) ** 2 / P) @ pi).mean())
d = np.diff(objs); print("(3) SIU objective: start %.6f -> %.6f after 30, max increase %.1e, its to 1e-9: %d" % (
    objs[0], objs[-1], d.max(), RU.solve(pi)[2]["its"]))
Hs, als, _ = RU.solve(pi)
print("    alpha quantiles over rows (5/50/95%%) at u = %s:" % np.round(np.exp(SIG * GH[0][[10, 20, 30]]), 2),
      np.round(np.quantile(als[:, [10, 20, 30]], [0.05, 0.5, 0.95], axis=0), 3).tolist())

# (4) ill-posedness of the original form: rho term not scaled with alpha -> objective decreases along H -> eps H
q, P = RU.qp(Hs); Pd = np.einsum("f,nfj->nj", w, np.abs(Hs[..., None]) ** 2 * np.abs(RU.Gc) ** 2)
Pr = np.einsum("f,nfj->nj", w, np.abs(Hs[..., None]) ** 2 * (RU.den - np.abs(RU.Gc) ** 2))
for eps in (1.0, 0.1, 0.01, 1e-4):
    orig = ((RU.E0 - np.maximum(q, 0) ** 2 / Pd + eps ** 2 * Pr) @ pi).mean()
    print("(4) original form at eps %.0e: %.6f" % (eps, orig))

# (5) timing, finiteness
for rule in ("SIU", "SIR"):
    t1 = time.monotonic(); X, info = si.design(D, cm, Q, c, TW, 1.4, rule, inv, w); dt = time.monotonic() - t1
    sm = {k: (np.round(v, 3).tolist() if k == "pi" else v) for k, v in info.items() if k not in ("regret", "alpha_q")}
    print("(5) semi %s %.0fs finite %s" % (rule, dt, bool(np.isfinite(X).all())), sm)
for m in (0.5, 2.0):
    t1 = time.monotonic(); X, info = si.design(D, cm, Q, m * c, TW, 1.4, "SIR", inv, w)
    print("    semi SIR m %.1f %.0fs gap %.1e its %d inner %d pi %s" % (m, time.monotonic() - t1, info["gap"], info["its"], info["inner"],
                                                                       np.round(info["pi"], 3).tolist()))
D0, T, cmf, tw, fp = fd.load_fd("H1"); Df = D0 + fd.noise(D0, 14.0, 1); invf = snr_rows(Df, tw, f); wf = si.wspec(fd.ricker(fp))
print("FD rows", len(tw), "fp", fp)
for rule in ("SIU", "SIR"):
    t1 = time.monotonic(); X, info = si.design(Df, cmf, cmf["Q"], 1.0, tw, 0.0, rule, invf, wf)
    print("(5) FD H1 %s %.0fs finite %s" % (rule, time.monotonic() - t1, bool(np.isfinite(X).all())),
          {k: (np.round(v, 3).tolist() if k == "pi" else v) for k, v in info.items() if k not in ("regret", "alpha_q")})
print("total %.0fs" % (time.monotonic() - t0))
