"""POST HOC addendum (written after the main run started; not in PLAN_ZH.md section 3): the same curvature split on the
U1320A field lines, the known counter-example of the uniform formula (4508-32 deep window).

Pipeline of qscan_anchor_theory_20261030/scripts/u1320a_kappa.py: analysis window 2.04-3.70 s, water t* fixed (Q >= 1e5
not scaled), scaled t* a1(t) v, wavelet window [t1, t2] with absorption at tc, lambda 1/2. Expansion point v0 = archived
selection at v_a = 1 (results18). One lasso solve per case, then the fixed-active-set closed form (no re-solved FD).
Reference: two-anchor kappa-hat of PAPER_DRAFT_ZH.md table 1. -> results_u1320a/<line>__<win>.json
"""
import json, sys, time
from core import np, lib, bm, sc11, ANCH, gmats, solve, polish, jfix, mu_of, HERE, LAM, RESEARCH
import u1320a_kappa as uk

uk.SRC = RESEARCH / "field_vsp_u1320a_20261024" / "inputs"
OUT = HERE / "results_u1320a"; EPS = 3e-3; bm.MAXITER = 20000
KHAT = {("3677-32", "w2.24"): 0.211, ("3677-32", "w3.00"): 0.816, ("4508-32", "w2.24"): 0.255, ("4508-32", "w3.00"): 0.607}


def wavelet_raw(D, cm, Q, va, t1, t2, tc, nfft=512):
    c = dict(cm, Q=uk.model_q(Q, va)); f, snr_s, P_s = lib.coherence_snr(D, cm["t"], t1, t2, nfft)
    amp = np.sqrt(P_s) / np.exp(-np.pi * f * lib.absorption(c, [tc])[0])
    fpk = float(f[np.argmax(np.sqrt(P_s) * (f >= 5))]); fs = lib.snr_band(f, snr_s, fpk)
    g = np.ones_like(f); r = (f > fs - 5) & (f < fs)
    g[r] = 0.5 * (1 + np.cos(np.pi * (f[r] - (fs - 5)) / 5)); g[f >= fs] = 0; g[f < 3] = 0
    return np.roll(np.fft.irfft(amp * g, nfft), 100)[:201]


def run(line, wname):
    t0 = time.monotonic(); D, cm = uk.load(line); Q = np.asarray(cm["Q"], float); pp = uk.prepare(D)
    t1, t2, tc = uk.WINS[wname]; n = D.shape[0]; t = np.asarray(cm["t"], float)
    af = lib.absorption(dict(cm, Q=np.where(Q >= 1e5, Q, np.inf)), t)          # water, fixed
    a1 = lib.absorption(dict(cm, Q=uk.model_q(Q, 1.0)), t) - af                   # scaled part at v = 1
    a0 = float(lib.absorption(dict(cm, Q=uk.model_q(Q, 1.0)), [tc])[0] - lib.absorption(dict(cm, Q=np.where(Q >= 1e5, Q, np.inf)), [tc])[0])
    v0 = json.loads((ANCH / "results18" / f"{line}__{wname}_va1.0.json").read_text(encoding="utf-8"))["v_sel"]
    WR = {va: wavelet_raw(D, cm, Q, va, t1, t2, tc) for va in (1 - EPS, 1.0, 1 + EPS)}
    Win = {va: lib.toeplitz_wavelet(w, n) for va, w in WR.items()}

    def Bt(dv=0.0, dva=0.0):
        return gmats(cm, af + a1 * (v0 + dv))[0] @ Win[1.0 + dva]
    B0 = Bt(); y0, J0, info = solve(D, B0, pp, tol=1e-3); mu0 = mu_of(B0, n); Y, A, S, pinfo = polish(D, B0, y0, pp, mu0)

    def JF(dv=0.0, dva=0.0, mode="pipe"):
        B = Bt(dv, dva); mu = mu_of(B, n) if mode == "pipe" else (mu0 if mode == "fix" else mu0 * mode[1])
        return jfix(D, B, A, S, pp, mu)
    e = EPS; Jc = JF(mode="fix")

    def second(f):
        return dict(Jvv=(f(e, 0) - 2 * Jc + f(-e, 0)) / e ** 2, Jvx=(f(e, e) - f(e, -e) - f(-e, e) + f(-e, -e)) / (4 * e * e))
    F = dict(pipe=second(lambda a, b: JF(a, b)), fix=second(lambda a, b: JF(a, b, "fix")))
    for k in F.values():
        k["alpha"] = -k["Jvx"] / k["Jvv"]
    Jm = {k: JF(mode=("m", 1 + k * e)) for k in (-1, 1)}
    cross = lambda key: (JF(**{key: e}, mode=("m", 1 + e)) - JF(**{key: e}, mode=("m", 1 - e))
                         - JF(**{key: -e}, mode=("m", 1 + e)) + JF(**{key: -e}, mode=("m", 1 - e))) / (4 * e * e * mu0)
    mp = dict(J_mu=(Jm[1] - Jm[-1]) / (2 * e * mu0), J_mumu=(Jm[1] - 2 * Jc + Jm[-1]) / (e * mu0) ** 2,
              J_vmu=cross("dv"), J_vamu=cross("dva"))
    M = lambda dv, dva: mu_of(Bt(dv, dva), n)
    md = dict(mu_v=(M(e, 0) - M(-e, 0)) / (2 * e), mu_va=(M(0, e) - M(0, -e)) / (2 * e), mu_vv=(M(e, 0) - 2 * mu0 + M(-e, 0)) / e ** 2,
              mu_vva=(M(e, e) - M(e, -e) - M(-e, e) + M(-e, -e)) / (4 * e * e))
    C = dict(tt=F["fix"]["Jvv"], tm=2 * mp["J_vmu"] * md["mu_v"], mm=mp["J_mumu"] * md["mu_v"] ** 2 + mp["J_mu"] * md["mu_vv"])
    N = dict(tt=-F["fix"]["Jvx"], tm=-(mp["J_vmu"] * md["mu_va"] + mp["J_vamu"] * md["mu_v"]),
             mm=-(mp["J_mumu"] * md["mu_v"] * md["mu_va"] + mp["J_mu"] * md["mu_vva"]))
    kmed = int(np.argsort(np.linalg.norm(B0, axis=0))[n // 2])
    return dict(line=line, window=wname, v0=v0, kappa_hat=KHAT[(line, wname)], a0=a0, t_med=float(t[kmed]),
                a1_med=float(a1[kmed]), med=a0 / float(a1[kmed]) if a1[kmed] > 0 else None,
                uniform=float(a0 * a1.mean() / (a1 ** 2).mean()), F=F, C=C, N=N,
                alpha_block={k: N[k] / C[k] for k in C}, share_C={k: C[k] / sum(C.values()) for k in C},
                polish=pinfo, solve0=info, seconds=time.monotonic() - t0)


if __name__ == "__main__":
    OUT.mkdir(exist_ok=True)
    for line in uk.LINES:
        for wname in uk.WINS:
            r = run(line, wname); (OUT / f"{line}__{wname}.json").write_text(json.dumps(r, indent=1), encoding="utf-8")
            print(line, wname, "khat %.3f F %.3f fix %.3f med %s U %.3f t_med %.2f" % (
                r["kappa_hat"], r["F"]["pipe"]["alpha"], r["F"]["fix"]["alpha"], r["med"], r["uniform"], r["t_med"]),
                {k: round(v, 2) for k, v in r["share_C"].items()}, {k: round(v, 3) for k, v in r["alpha_block"].items()},
                f"{r['seconds']:.0f}s", flush=True)
