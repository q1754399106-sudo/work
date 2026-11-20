"""Per-dataset curvature of the anchor-point lasso objective (PLAN_ZH.md section 2).

usage: run.py [--smoke]    worker; jobs claimed with results/<job>.claim, resumable.
Expansion point: v0 = archived selection at anchor c_a = 1 (v_a = 1), lambda 1/2.
 R  re-solved finite differences of the scan objective J(v, v_a) (full lasso, warm start, KKT tol 1e-3) with steps h
    0.03 and 0.10: J_vv, J_v va -> alpha_R = -J_v va / J_vv (local two-anchor slope; the "truth").
 F  fixed active set (exact closed form, step 3e-3): alpha_F with the L1 weight mu following v and v_a as in the scan,
    alpha_F0 with mu following v_a only.
 H  analytic theta Hessian on the fixed active set (mu fixed): alpha_H = a0 a'H1 / a'Ha (delta model), direct part only
    (alpha_Hd), diagonals (alpha_diagH, alpha_diagHd); depth profiles and 8-segment / 0.1 s block sums.
 W  diagonal weights of the earlier tests: uniform, energy envelope of x^2, (|y| - pen)_+^2 with y = b_k'd, non-zero
    density (0.1 s moving average, as anchor_active_curv_20261119), wavefield energy z^2.
--smoke: development seed 9940, profile D, c 1.0 (runtime check only, not part of any table).
-> results/<job>.json
"""
import os, sys, json, time
from core import (np, lib, gen, qp, bm, sc11, HERE, datasets, load, v_selected, wavelet, unnorm_dict, mu_of, solve,
                  polish, jfix, theta_hessian, wavelet_raw, alpha_full, alpha_diag, smooth, LAM)

OUT = HERE / "results"; HR = 0.05; QS = (25.0, 50.0, 75.0); EPS = 3e-3; TOL = 1e-3; bm.MAXITER = 20000


def blocks(H, a, t, edges):
    seg = np.clip(np.searchsorted(edges, t, side="right") - 1, 0, len(edges) - 2); nb = len(edges) - 1
    P = np.zeros((nb, len(t)))
    P[seg, np.arange(len(t))] = a
    return (P @ H @ P.T).tolist()


def run(prof, did, c, s, base):
    t0 = time.monotonic(); D, cm = load(prof, did, c, s, base); pp = sc11.prepare_sc(D, cm); n = D.shape[0]
    t = bm.vec(cm["t"]).astype(float); a = lib.absorption(cm, t); a0 = float(lib.absorption(cm, [1.4])[0])
    v0 = 1.05 if did.endswith("9940") else v_selected(prof, did)
    WR = {va: wavelet_raw(D, cm, va) for va in (1.0 - EPS, 1.0, 1.0 + EPS, 1.0 - HR, 1.0 + HR)}
    one = np.ones(n); Bt0, _ = unnorm_dict(cm, v0, WR[1.0])
    y0, J0, info0 = solve(D, Bt0, pp, tol=TOL); mu0 = mu_of(Bt0, n)
    Y, A, S, pinfo = polish(D, Bt0, y0, pp, mu0)
    sc0 = np.linalg.norm(Bt0, axis=0); order = np.argsort(sc0)
    pct0 = {q: mu_of(Bt0, n, q) for q in QS}

    # F: fixed active set, closed form. mu modes: "pipe" (median column norm of the current dictionary, as the scan),
    # "fix" (mu0), ("q", p) (mu0 x percentile-p norm / its value at the base point: intervention), ("m", k) (mu0 x k)
    def JF(dv=0.0, dva=0.0, xth=0.0, xamp=0.0, mode="pipe"):
        Bt, _ = unnorm_dict(cm, v0 + dv, WR[1.0 + dva], theta=(-a0 * xth * one) if xth else None,
                            theta_amp=(-a0 * xamp * one) if xamp else None)
        if mode == "pipe":
            mu = mu_of(Bt, n)
        elif mode == "fix":
            mu = mu0
        elif mode[0] == "q":
            mu = mu0 * mu_of(Bt, n, mode[1]) / pct0[mode[1]]
        else:
            mu = mu0 * mode[1]
        return jfix(D, Bt, A, S, pp, mu)

    e = EPS; Jc = JF(mode="fix")

    def second(f):
        """J_vv, J_vx, J_v from f(dv, dx) on a +-e stencil."""
        fp, fm = f(e, 0), f(-e, 0)
        return dict(Jv=(fp - fm) / (2 * e), Jvv=(fp - 2 * Jc + fm) / e ** 2,
                    Jvx=(f(e, e) - f(e, -e) - f(-e, e) + f(-e, -e)) / (4 * e * e))
    F = {}
    F["pipe"] = second(lambda dv, dx: JF(dv, dva=dx))
    F["fix"] = second(lambda dv, dx: JF(dv, dva=dx, mode="fix"))
    F["fix_theta"] = second(lambda dv, dx: JF(dv, xth=dx, mode="fix"))
    F["fix_amp"] = second(lambda dv, dx: JF(dv, xamp=dx, mode="fix"))
    for q in QS:
        if q != 50.0:
            F[f"q{q:g}"] = second(lambda dv, dx, q=q: JF(dv, dva=dx, mode=("q", q)))
    for k in F.values():
        k["alpha"] = -k["Jvx"] / k["Jvv"]
    # mu partials at fixed dictionary (relative step e in mu) and the derivatives of the scan's mu(v, v_a)
    Jm = {k: JF(mode=("m", 1 + k * e)) for k in (-1, 1)}
    mpart = dict(J_mu=(Jm[1] - Jm[-1]) / (2 * e * mu0), J_mumu=(Jm[1] - 2 * Jc + Jm[-1]) / (e * mu0) ** 2,
                 J_vmu=(JF(e, mode=("m", 1 + e)) - JF(e, mode=("m", 1 - e)) - JF(-e, mode=("m", 1 + e))
                        + JF(-e, mode=("m", 1 - e))) / (4 * e * e * mu0),
                 J_vamu=(JF(dva=e, mode=("m", 1 + e)) - JF(dva=e, mode=("m", 1 - e)) - JF(dva=-e, mode=("m", 1 + e))
                         + JF(dva=-e, mode=("m", 1 - e))) / (4 * e * e * mu0))
    M = lambda dv, dva: mu_of(unnorm_dict(cm, v0 + dv, WR[1.0 + dva])[0], n)
    mder = dict(mu0=mu0, mu_v=(M(e, 0) - M(-e, 0)) / (2 * e), mu_va=(M(0, e) - M(0, -e)) / (2 * e),
                mu_vv=(M(e, 0) - 2 * mu0 + M(-e, 0)) / e ** 2,
                mu_vva=(M(e, e) - M(e, -e) - M(-e, e) + M(-e, -e)) / (4 * e * e))
    kmed = int(order[n // 2]); tq = {f"{q:g}": float(t[int(order[int(round(q / 100 * (n - 1)))])]) for q in QS}

    # R: re-solved finite differences (pipeline objective, normalised wavelets are equivalent with mu recomputed)
    R = {}; conv = [info0["converged"]]
    for h in (HR,):
        def JR(dv, dva):
            Bt, _ = unnorm_dict(cm, v0 + dv, WR[1.0 + dva]); _, J, inf = solve(D, Bt, pp, start=y0, tol=TOL)
            conv.append(inf["converged"]); return J
        Jp, Jm_ = JR(h, 0), JR(-h, 0)
        Jpp, Jpm, Jmp, Jmm = JR(h, h), JR(h, -h), JR(-h, h), JR(-h, -h)
        Jvv = (Jp - 2 * J0 + Jm_) / h ** 2; Jvva = (Jpp - Jpm - Jmp + Jmm) / (4 * h * h)
        R[f"{h:g}"] = dict(Jv=(Jp - Jm_) / (2 * h), Jvv=Jvv, Jvva=Jvva, alpha=-Jvva / Jvv)

    # H: analytic theta Hessian, mu fixed
    H, Hd, z2 = theta_hessian(D, cm, v0, WR[1.0], Y, A, pp)
    live = pp["live"]; X = Y * np.linalg.norm(Bt0, axis=0)[:, None]; XL = X[:, live]
    scale = np.linalg.norm(Bt0, axis=0); Bn = Bt0 / scale
    pen = LAM * np.sqrt(2 * np.log(n)) * pp["sigma"][None, live] * np.median(scale) / scale[:, None]
    ythr = (np.maximum(np.abs(Bn.T @ D[:, live]) - pen, 0) ** 2).mean(axis=1)
    prof_w = dict(e=(XL ** 2).mean(axis=1), n=(np.abs(XL) > 0).mean(axis=1), ythr=ythr, z2=z2)
    edges8 = np.linspace(t[0], t[-1] + 1e-9, 9); edges10 = np.arange(1.2, 3.0 + 1e-9 + 0.1, 0.1)
    rec = dict(profile=prof, id=did, c_true=c, seed=s, v0=v0, a0=a0, lam=LAM, t=t.tolist(), a=a.tolist(),
               n_live=int(live.sum()), solve0=info0, J0=J0, polish=pinfo, all_converged=bool(all(conv)),
               F=F, R=R, mu_partials=mpart, mu_derivs=mder, t_median_atom=float(t[kmed]), t_percentile_atom=tq,
               a_at=dict(zip(["t_med"] + [f"q{k}" for k in tq], [float(a[kmed])] + [float(np.interp(v, t, a)) for v in tq.values()])),
               alpha=dict(U=alpha_diag(one, a, a0), H=alpha_full(H, a, a0), Hd=alpha_full(Hd, a, a0),
                          diagH=alpha_diag(np.diag(H), a, a0), diagHd=alpha_diag(np.diag(Hd), a, a0),
                          E=alpha_diag(smooth(prof_w["e"]), a, a0), N=alpha_diag(smooth(prof_w["n"]), a, a0),
                          Y=alpha_diag(smooth(prof_w["ythr"]), a, a0), Z=alpha_diag(smooth(prof_w["z2"]), a, a0)),
               curv=dict(aHa=float(a @ H @ a), aH1=float(a @ H @ one), aHda=float(a @ Hd @ a), aHd1=float(a @ Hd @ one)),
               depth=dict(diagH=np.diag(H).tolist(), diagHd=np.diag(Hd).tolist(), Ha=(H @ a).tolist(), H1=(H @ one).tolist(),
                          Hda=(Hd @ a).tolist(), Hd1=(Hd @ one).tolist(), **{k: v.tolist() for k, v in prof_w.items()}),
               blocks8=dict(edges=edges8.tolist(), H=blocks(H, a, t, edges8), Hd=blocks(Hd, a, t, edges8)),
               blocks01=dict(edges=edges10.tolist(), H=blocks(H, a, t, edges10), Hd=blocks(Hd, a, t, edges10)),
               seconds=time.monotonic() - t0)
    return rec


if __name__ == "__main__":
    OUT.mkdir(exist_ok=True); base = gen.prepare(qp.STRUCT, qp.NOISE)
    if "--smoke" in sys.argv:
        r = run(sys.argv[2] if len(sys.argv) > 2 else "D", "c1.0_s9940", 1.0, 9940, base)
        print(json.dumps(dict(alpha=r["alpha"], F={k: v["alpha"] for k, v in r["F"].items()}, R=r["R"], mp=r["mu_partials"],
                              md=r["mu_derivs"], tq=r["t_percentile_atom"], a_at=r["a_at"], a0=r["a0"], polish=r["polish"],
                              conv=r["all_converged"], seconds=r["seconds"]), indent=1)); sys.exit()
    for prof, did, c, s in datasets():
        job = f"{prof}_{did}"; path = OUT / f"{job}.json"
        if path.exists():
            continue
        try:
            os.close(os.open(OUT / f"{job}.claim", os.O_CREAT | os.O_EXCL | os.O_WRONLY))
        except FileExistsError:
            continue
        rec = run(prof, did, c, s, base)
        path.write_text(json.dumps(rec), encoding="utf-8"); (OUT / f"{job}.claim").unlink()
        print(job, "conv" if rec["all_converged"] else "NOT-CONVERGED", "R %.3f F %.3f H %.3f" % (
            rec["R"]["0.05"]["alpha"], rec["F"]["pipe"]["alpha"], rec["alpha"]["H"]), f"{rec['seconds']:.0f}s", flush=True)
