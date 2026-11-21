"""Pre-registered confirmation (PLAN_ZH.md): single-solve fixed-active-set closed-form kappa vs the two-anchor scan.

usage: run.py [--smoke]    worker; jobs claimed with results/<job>.claim, resumable.
Per dataset (new seeds 9921 / 9922, 17 profiles, c 0.7 / 1.0 / 1.4):
 scans  fixed-anchor QScan at c_a = 1.0, 0.7, 1.4 (qprofile.scan, lambda 1/2): v_sel per anchor;
        truth kappa_hat = (v_sel(0.7) - v_sel(1.4)) / (1/0.7 - 1/1.4).
 F      one lasso solve at v0 = v_sel(1.0) (unnormalised wavelet, mu recomputed: same objective as the scan), fixed
        active set and signs, closed form J, step 3e-3: alpha_F = -J_v va / J_vv (method under test).
        Block split theta-theta / theta-mu / mu-mu (as anchor_active_hessian_20261120) for description.
 formulas  uniform U = a0 <a>/<a^2> and the post-hoc a0 / a(t_med) (t_med = atom of the median column norm).
--smoke: development seed 9919, profile N3, c 1.0 (runtime check only, not part of any table).
-> results/<job>.json
"""
import os, sys, json, time
from core import (np, lib, gen, qp, bm, sc11, HERE, datasets, load, scan_anchor, unnorm_dict, mu_of, solve, polish,
                  jfix, wavelet_raw, alpha_diag, ANCHORS, DEV_SEED, FAM_A)

OUT = HERE / "results"; EPS = 3e-3; TOL = 1e-3; bm.MAXITER = 20000
DVA = 1.0 / 0.7 - 1.0 / 1.4


def closed_form(D, cm, pp, v0):
    n = D.shape[0]; t = bm.vec(cm["t"]).astype(float); a = lib.absorption(cm, t); a0 = float(lib.absorption(cm, [1.4])[0])
    WR = {va: wavelet_raw(D, cm, va) for va in (1.0 - EPS, 1.0, 1.0 + EPS)}
    Bt0, _ = unnorm_dict(cm, v0, WR[1.0]); y0, J0, info0 = solve(D, Bt0, pp, tol=TOL); mu0 = mu_of(Bt0, n)
    Y, A, S, pinfo = polish(D, Bt0, y0, pp, mu0)

    def JF(dv=0.0, dva=0.0, mode="pipe"):
        Bt, _ = unnorm_dict(cm, v0 + dv, WR[1.0 + dva])
        mu = mu_of(Bt, n) if mode == "pipe" else (mu0 if mode == "fix" else mu0 * mode[1])
        return jfix(D, Bt, A, S, pp, mu)
    e = EPS; Jc = JF(mode="fix")

    def second(f):
        return dict(Jvv=(f(e, 0) - 2 * Jc + f(-e, 0)) / e ** 2, Jvx=(f(e, e) - f(e, -e) - f(-e, e) + f(-e, -e)) / (4 * e * e))
    F = dict(pipe=second(lambda x, y: JF(x, y)), fix=second(lambda x, y: JF(x, y, "fix")))
    for k in F.values():
        k["alpha"] = -k["Jvx"] / k["Jvv"]
    Jm = {k: JF(mode=("m", 1 + k * e)) for k in (-1, 1)}
    cross = lambda key: (JF(**{key: e}, mode=("m", 1 + e)) - JF(**{key: e}, mode=("m", 1 - e))
                         - JF(**{key: -e}, mode=("m", 1 + e)) + JF(**{key: -e}, mode=("m", 1 - e))) / (4 * e * e * mu0)
    mp = dict(J_mu=(Jm[1] - Jm[-1]) / (2 * e * mu0), J_mumu=(Jm[1] - 2 * Jc + Jm[-1]) / (e * mu0) ** 2,
              J_vmu=cross("dv"), J_vamu=cross("dva"))
    M = lambda dv, dva: mu_of(unnorm_dict(cm, v0 + dv, WR[1.0 + dva])[0], n)
    md = dict(mu_v=(M(e, 0) - M(-e, 0)) / (2 * e), mu_va=(M(0, e) - M(0, -e)) / (2 * e),
              mu_vv=(M(e, 0) - 2 * mu0 + M(-e, 0)) / e ** 2, mu_vva=(M(e, e) - M(e, -e) - M(-e, e) + M(-e, -e)) / (4 * e * e))
    C = dict(tt=F["fix"]["Jvv"], tm=2 * mp["J_vmu"] * md["mu_v"], mm=mp["J_mumu"] * md["mu_v"] ** 2 + mp["J_mu"] * md["mu_vv"])
    N = dict(tt=-F["fix"]["Jvx"], tm=-(mp["J_vmu"] * md["mu_va"] + mp["J_vamu"] * md["mu_v"]),
             mm=-(mp["J_mumu"] * md["mu_v"] * md["mu_va"] + mp["J_mu"] * md["mu_vva"]))
    kmed = int(np.argsort(np.linalg.norm(Bt0, axis=0))[n // 2])
    return dict(v0=v0, a0=a0, F=F, C=C, N=N, alpha_block={k: N[k] / C[k] for k in C},
                share_C={k: C[k] / sum(C.values()) for k in C},
                closure=dict(Jvv=sum(C.values()) / F["pipe"]["Jvv"] - 1, Jvx=-sum(N.values()) / F["pipe"]["Jvx"] - 1),
                t_med=float(t[kmed]), U=alpha_diag(np.ones(n), a, a0), med=a0 / float(a[kmed]),
                polish=pinfo, solve0=info0, J0=J0)


def run(prof, did, c, s, base):
    t0 = time.monotonic(); D, cm = load(prof, did, c, s, base); pp = sc11.prepare_sc(D, cm)
    scans = {}
    for ca in ANCHORS:
        r = scan_anchor(D, cm, pp, ca); r.pop("J"); r.pop("converged"); scans[f"a{ca}"] = r
    ts = time.monotonic() - t0
    cf = closed_form(D, cm, pp, scans["a1.0"]["v_sel"])
    khat = (scans["a0.7"]["v_sel"] - scans["a1.4"]["v_sel"]) / DVA
    valid = bool(all(not x["edge"] and x["neighbourhood_converged"] for x in scans.values()) and cf["solve0"]["converged"])
    return dict(profile=prof, family="A" if prof in FAM_A else "B", Q_profile=qp.PROFILES[prof], id=did, c_true=c,
                seed=s, lam=qp.LAM, scans=scans, kappa_hat=khat, valid=valid, seconds_scans=ts,
                seconds=time.monotonic() - t0, **cf)


if __name__ == "__main__":
    OUT.mkdir(exist_ok=True); base = gen.prepare(qp.STRUCT, qp.NOISE)
    if "--smoke" in sys.argv:
        r = run("N3", f"c1.0_s{DEV_SEED}", 1.0, DEV_SEED, base)
        print(json.dumps({k: r[k] for k in ("kappa_hat", "valid", "t_med", "U", "med", "alpha_block", "share_C", "closure",
                                             "polish", "seconds_scans", "seconds")}, indent=1),
              {k: v["v_sel"] for k, v in r["scans"].items()}, r["F"]["pipe"]["alpha"]); sys.exit()
    for prof, did, c, s in datasets():
        job = f"{prof}_{did}"; path = OUT / f"{job}.json"
        if path.exists():
            continue
        try:
            os.close(os.open(OUT / f"{job}.claim", os.O_CREAT | os.O_EXCL | os.O_WRONLY))
        except FileExistsError:
            continue
        rec = run(prof, did, c, s, base)
        path.write_text(json.dumps(rec, indent=1), encoding="utf-8"); (OUT / f"{job}.claim").unlink()
        print(job, "valid" if rec["valid"] else "INVALID", "khat %.3f F %.3f U %.3f med %.3f" % (
            rec["kappa_hat"], rec["F"]["pipe"]["alpha"], rec["U"], rec["med"]), f"{rec['seconds']:.0f}s", flush=True)
