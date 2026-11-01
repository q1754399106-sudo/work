"""Segment perturbation (PLAN_SEGPHI_ZH.md): phi2 per time segment from the curvature of the re-optimized total J.

usage: segphi.py [--check]    worker; jobs claimed with results22/<job>.claim, resumable.
Base: the selection v0 at v_a = 1, lambda 1/2 (U1320A: results18; Volve: results19). Atom k gets the attenuation
a(t_k) (v0 + h) if t_k lies in segment b, else a(t_k) v0 (U1320A: plus the fixed water t*); h = -0.5, -0.25, 0, 0.25, 0.5; 8 equal segments over the
analysis window; plus the uniform perturbation (all atoms). Every point is a full lasso solve (warm start from h = 0).
--check: dict_a with a = a(t) v reproduces lib.dictionary (implementation check only).
-> results22/<job>.json  (J per segment and h, J uniform per h, t, a, a0, v0)
"""
import os, sys, json, time
from env15 import np, lib, semi, sc11, bm, gen, ROOT
import u1320a_kappa as uk
import qprofile as qp

OUT = ROOT / "results22"; bm.MAXITER = 20000; LAM = 0.5; NSEG = 8; H = (-0.5, -0.25, 0.25, 0.5)
qp.PROFILES.update({"D2": [150.0, 150.0, 60.0, 60.0], "M": [200.0, 45.0, 200.0, 200.0]})


def dict_a(cm, a, w_in, w_out):
    """lib.dictionary with the per-atom cumulative t* given directly."""
    t = bm.vec(cm["t"]); n = len(t); N = int(cm["nfft"]); f = np.arange(N // 2 + 1) / (N * float(cm["dt"]))
    b = -np.pi * f + 2j * f * np.log(np.maximum(f, np.finfo(float).tiny) / float(cm["fref"])); b[0] = 0
    G = np.empty((n, n))
    for k in range(n):
        Hk = np.exp(b * a[k]); Hk[-1] = Hk[-1].real
        G[:, k] = np.fft.irfft(Hk * np.exp(-2j * np.pi * np.arange(N // 2 + 1) * k / N), n=N)[:n]
    B = G @ lib.toeplitz_wavelet(w_in, n); scale = np.linalg.norm(B, axis=0)
    return B / scale, scale


def objective(D, B, scale, pp, start):
    sigma, K, live = pp["sigma"], pp["K"], pp["live"]; n = D.shape[0]; s2 = np.where(live, sigma ** 2, np.inf)
    X, info = bm.solve(B, np.eye(n), D, sigma, scale, K, LAM, 0.0, live, start=start)
    pen1 = np.sqrt(2 * np.log(n)) * sigma[None, :] * np.median(scale) / scale[:, None]; res = D - B @ X
    J = float(np.sum((0.5 * np.sum(res ** 2, axis=0) + LAM * np.sum(pen1 * np.abs(X), axis=0))[live] / s2[live]))
    return J, X, bool(info["converged"])


def perturb(D, cm, a1, v0, wv, pp, af=0.0):
    t = np.asarray(cm["t"]); edges = np.linspace(t[0], t[-1] + 1e-9, NSEG + 1)
    seg = np.clip(np.searchsorted(edges, t, side="right") - 1, 0, NSEG - 1)
    B, sc = dict_a(cm, af + a1 * v0, wv["w_in"], wv["w_out"]); J0, X0, c0 = objective(D, B, sc, pp, None)
    Js, conv = [], [c0]
    for b in range(NSEG):
        row = []
        for h in H:
            B, sc = dict_a(cm, af + a1 * (v0 + h * (seg == b)), wv["w_in"], wv["w_out"]); J, _, c = objective(D, B, sc, pp, X0)
            row.append(J); conv.append(c)
        Js.append(row)
    Ju = []
    for h in H:
        B, sc = dict_a(cm, af + a1 * (v0 + h), wv["w_in"], wv["w_out"]); J, _, c = objective(D, B, sc, pp, X0); Ju.append(J); conv.append(c)
    return dict(J0=J0, J_seg=Js, J_uniform=Ju, h=list(H), seg_edges=edges.tolist(), all_converged=all(conv))


def run_u(line, wname):
    D, cm = uk.load(line); Q = np.asarray(cm["Q"], float); pp = uk.prepare(D); t1, t2, tc = uk.WINS[wname]
    wv = uk.wavelets(D, cm["t"], dict(cm, Q=uk.model_q(Q, 1.0)), t1, t2, tc)
    v0 = json.loads((ROOT / "results18" / f"{line}__{wname}_va1.0.json").read_text())["v_sel"]
    af = lib.absorption(dict(cm, Q=np.where(Q >= 1e5, Q, np.inf)), cm["t"])       # water: fixed, not scaled by v
    a1 = lib.absorption(cm, cm["t"]) - af; t0 = time.monotonic(); r = perturb(D, cm, a1, v0, wv, pp, af)
    return dict(kind="U", line=line, window=wname, v0=v0, t=np.asarray(cm["t"]).tolist(), a=a1.tolist(),
                a0=float(np.interp(tc, cm["t"], a1)), seconds=time.monotonic() - t0, **r)


def run_v(prof, c, s, base):
    ctx = base if prof == "B" else qp.context(base, prof); cm = dict(ctx["cm"]); Q = np.asarray(cm["Q"], float)
    g = gen.signal(ctx, c, s); D = g["x"] + qp.MULT * g["level"] * ctx["noise"][0]; pp = sc11.prepare_sc(D, cm)
    wv = lib.wavelets(D, semi.TW, dict(cm, Q=Q * 1.0))
    v0 = json.loads((ROOT / "results19" / f"V_{prof}_c{c}_s{s}__a1.0.json").read_text())["v_sel"]
    a1 = lib.absorption(cm, cm["t"]); t0 = time.monotonic(); r = perturb(D, cm, a1, v0, wv, pp)
    return dict(kind="V", profile=prof, id=f"c{c}_s{s}", u=1.0 / c, v0=v0, t=np.asarray(cm["t"]).tolist(), a=a1.tolist(),
                a0=float(lib.absorption(cm, [1.4])[0]), seconds=time.monotonic() - t0, **r)


def jobs():
    J = [(f"U_{l}__{w}", ("u", l, w)) for w in ("w3.00", "w2.24") for l in uk.LINES]
    for p, s in (("B", 9971), ("S", 9971), ("D", 9971), ("D2", 9981), ("M", 9981)):
        J += [(f"V_{p}_c{c}_s{s}", ("v", p, c, s)) for c in qp.CS]
    return J


if __name__ == "__main__":
    if "--check" in sys.argv:
        D, cm = uk.load("3677-32"); Q = np.asarray(cm["Q"], float); wv = uk.wavelets(D, cm["t"], cm, 2.04, 2.44, 2.24)
        for v in (0.5, 2.0):
            B1, _, s1 = lib.dictionary(dict(cm, Q=uk.model_q(Q, v)), wv["w_in"], wv["w_out"])
            af = lib.absorption(dict(cm, Q=np.where(Q >= 1e5, Q, np.inf)), cm["t"]); a1 = lib.absorption(cm, cm["t"]) - af
            B2, s2 = dict_a(cm, af + a1 * v, wv["w_in"], wv["w_out"])
            print("v", v, "max |dB|", float(np.max(np.abs(B1 - B2))), "max rel dscale", float(np.max(np.abs(s1 / s2 - 1))))
        sys.exit()
    OUT.mkdir(exist_ok=True); base = None
    for job, spec in jobs():
        path = OUT / f"{job}.json"
        if path.exists():
            continue
        try:
            os.close(os.open(OUT / f"{job}.claim", os.O_CREAT | os.O_EXCL | os.O_WRONLY))
        except FileExistsError:
            continue
        if spec[0] == "u":
            rec = run_u(*spec[1:])
        else:
            base = base or gen.prepare(qp.STRUCT, qp.NOISE); rec = run_v(*spec[1:], base)
        path.write_text(json.dumps(rec), encoding="utf-8"); (OUT / f"{job}.claim").unlink()
        print(job, f"{rec['seconds']:.0f}s", "conv" if rec["all_converged"] else "NOT all converged", flush=True)
