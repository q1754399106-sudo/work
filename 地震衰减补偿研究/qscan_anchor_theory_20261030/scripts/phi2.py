"""phi2 from the objective (PLAN_PHI2_ZH.md): scans that record the lasso objective split by time sample.

usage: phi2.py [stage]    worker (stage 1 = U1320A, 2 = Volve, default both); jobs claimed with results19/<job>.claim.
Per grid point: J_t = 0.5 sum_traces r(t)^2/s2 + lam sum_traces pen1(t)|x(t)|/s2 over the live traces (atom k sits at
time t[k]); sum_t J_t = total J. Scans, wavelets, grids and early stopping are those of u1320a_kappa.py (stage 1) and
qprofile.py (stage 2). Analysis (curvature per sample, kappa_curv) in analyze_phi2.py.
-> results19/<job>.json  (v grid, total J, J_t residual and penalty parts, t, a(t), a0, v_sel)
"""
import os, sys, json, time
from env15 import np, lib, semi, sc11, bm, gen, ROOT
import u1320a_kappa as uk
import qprofile as qp

OUT = ROOT / "results19"; bm.MAXITER = 20000; LAM = 0.5
qp.PROFILES.update({"D2": [150.0, 150.0, 60.0, 60.0], "M": [200.0, 45.0, 200.0, 200.0]})


def split_J(D, B, X, sigma, scale, live):
    n = D.shape[0]; s2 = sigma[live] ** 2
    pen1 = np.sqrt(2 * np.log(n)) * sigma[None, :] * np.median(scale) / scale[:, None]
    res = (D - B @ X)[:, live]
    jr = np.sum(0.5 * res ** 2 / s2[None, :], axis=1); jp = np.sum(LAM * pen1[:, live] * np.abs(X[:, live]) / s2[None, :], axis=1)
    return jr, jp


def scan_up(D, pp, vgrid, build, past):
    """Scan v upwards (vgrid ascending), stop `past` points past the running minimum (U1320A 3, qprofile 2)."""
    sigma, K, live = pp["sigma"], pp["K"], pp["live"]; n = D.shape[0]
    rec = []; prev = None
    for j, v in enumerate(vgrid):
        if rec and j >= int(np.argmin([r[1] for r in rec])) + past + 1:
            break
        B, scale = build(v)
        X, info = bm.solve(B, np.eye(n), D, sigma, scale, K, LAM, 0.0, live, start=prev); prev = X
        jr, jp = split_J(D, B, X, sigma, scale, live)
        rec.append((float(v), float(jr.sum() + jp.sum()), jr, jp, bool(info["converged"])))
    J = np.array([r[1] for r in rec]); jm = int(np.argmin(J)); v = rec[jm][0]; vv = np.array([r[0] for r in rec])
    if 0 < jm < len(rec) - 1:
        a = np.polyfit(vv[jm - 1:jm + 2], J[jm - 1:jm + 2], 2)
        if a[0] > 0:
            v = float(np.clip(-a[1] / (2 * a[0]), vv[jm - 1], vv[jm + 1]))
    return dict(v_sel=v, v=vv.tolist(), J=J.tolist(), Jres=[np.round(r[2], 6).tolist() for r in rec], Jpen=[np.round(r[3], 6).tolist() for r in rec],
                converged=[r[4] for r in rec])


def run_u(line, wname, va):
    D, cm = uk.load(line); Q = np.asarray(cm["Q"], float); pp = uk.prepare(D); t1, t2, tc = uk.WINS[wname]
    wv = uk.wavelets(D, cm["t"], dict(cm, Q=uk.model_q(Q, va)), t1, t2, tc)

    def build(v):
        B, _, scale = lib.dictionary(dict(cm, Q=uk.model_q(Q, v)), wv["w_in"], wv["w_out"]); return B, scale

    t0 = time.monotonic(); r = scan_up(D, pp, uk.VG, build, 3)
    return dict(stage=1, line=line, window=wname, v_anchor=va, t=np.asarray(cm["t"]).tolist(),
                a=lib.absorption(cm, cm["t"]).tolist(), a0=float(lib.absorption(cm, [tc])[0]), seconds=time.monotonic() - t0, **r)


def run_v(prof, c, s, ca, base):
    ctx = base if prof == "B" else qp.context(base, prof)
    cm = dict(ctx["cm"]); Q = np.asarray(cm["Q"], float)          # B: the model Q of base (78.509/135.178)
    g = gen.signal(ctx, c, s); D = g["x"] + qp.MULT * g["level"] * ctx["noise"][0]; pp = sc11.prepare_sc(D, cm)
    wv = lib.wavelets(D, semi.TW, dict(cm, Q=Q * ca)); order = np.argsort(qp.VG)      # ascending v = descending c

    def build(v):
        B, _, scale = lib.dictionary(dict(cm, Q=Q / v), wv["w_in"], wv["w_out"]); return B, scale

    t0 = time.monotonic(); r = scan_up(D, pp, qp.VG[order], build, 2)
    return dict(stage=2, profile=prof, Q_profile=Q.tolist(), id=f"c{c}_s{s}", c_true=c, u=1.0 / c, seed=s, c_anchor=ca, v_anchor=1.0 / ca,
                t=np.asarray(cm["t"]).tolist(), a=lib.absorption(cm, cm["t"]).tolist(), a0=float(lib.absorption(cm, [1.4])[0]),
                seconds=time.monotonic() - t0, **r)


def jobs(stage):
    J = []
    if stage in (None, 1):
        J += [(f"U_{l}__{w}_va{va}", ("u", l, w, va)) for l in uk.LINES for w in uk.WINS for va in uk.VAS]
    if stage in (None, 2):
        for s in (9971, 9972, 9973, 9974):
            J += [(f"V_B_c{c}_s{s}__a{ca}", ("v", "B", c, s, ca)) for c in qp.CS for ca in (0.7, 1.0, 1.4)]
            J += [(f"V_{p}_c{c}_s{s}__a1.0", ("v", p, c, s, 1.0)) for p in ("S", "D") for c in qp.CS]
        for s in (9981, 9982, 9983, 9984):
            J += [(f"V_{p}_c{c}_s{s}__a1.0", ("v", p, c, s, 1.0)) for p in ("D2", "M") for c in qp.CS]
    return J


if __name__ == "__main__":
    OUT.mkdir(exist_ok=True); stage = int(sys.argv[1]) if len(sys.argv) > 1 else None; base = None
    for job, spec in jobs(stage):
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
        print(job, f"v {rec['v_sel']:.3f}", f"{rec['seconds']:.0f}s", flush=True)
