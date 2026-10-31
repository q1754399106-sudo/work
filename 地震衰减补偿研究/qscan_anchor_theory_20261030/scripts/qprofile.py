"""Q-profile test (PLAN_PROFILE_ZH.md): does kappa follow the shape of the absorption profile?

usage: qprofile.py [--smoke]    worker; jobs claimed with results7/<job>.claim, resumable.
Data: experiment-14 recipe (structure line 10227, noise lines 10149/10273/10355, noise x 0.5) with the model Q
profile replaced, for both generation and scanning: "S" = [30, 200, 200, 200] (lossy overburden, kappa 0.880) and
"D" = [250, 250, 40, 40] (lossy below 2.04 s, kappa 0.294); true Q = c x profile, c 0.7 / 1.0 / 1.4, new seeds
9971-9974. The source is re-estimated from the real line with the new profile (lib.wavelets), so that every
scenario reproduces the real shallow spectrum as in gen.signal.
Scans (lambda 1/2, full QGRID, warm start along the grid): fixed anchor c_a 0.7 / 1.0 / 1.4 (default wavelet
window 1.2-1.6 s); QScan-W (wavelet re-estimated for every c). Selection: grid minimum and the vertex of the
parabola in v = 1/c through it and its two neighbours.
-> results7/<profile>_<dataset>__<scan>.json
"""
import os, sys, json, time
from env15 import np, lib, semi, sc11, bm, gen, ROOT

STRUCT, NOISE = 10227, [10149, 10273, 10355]
MULT, LAM, MAXITER = 0.5, 0.5, 20000
PROFILES = {"S": [30.0, 200.0, 200.0, 200.0], "D": [250.0, 250.0, 40.0, 40.0]}
CS = (0.7, 1.0, 1.4); SEEDS = (9971, 9972, 9973, 9974); SCANS = ("a0.7", "a1.0", "a1.4", "w")
OUT = ROOT / "results7"; bm.MAXITER = MAXITER
VG = 1.0 / np.asarray(lib.QGRID, float)


def jobs():
    return [(p, f"c{c}_s{s}", c, s, sc) for p in PROFILES for s in SEEDS for c in CS for sc in SCANS]


def context(base, prof):
    cm = dict(base["cm"]); cm["Q"] = np.asarray(PROFILES[prof], float)
    return dict(base, cm=cm, w_in=lib.wavelets(base["D"], semi.TW, cm)["w_in"])


def scan(D, cm, Q, pp, wave):
    sigma, K, live = pp["sigma"], pp["K"], pp["live"]; n = D.shape[0]; s2 = np.where(live, sigma ** 2, np.inf)
    J = np.zeros(len(lib.QGRID)); C = np.zeros(J.shape, bool); prev = None
    for j, c in enumerate(lib.QGRID):
        cc = dict(cm); cc["Q"] = Q * c; wv = wave(cc)
        B, _, scale = lib.dictionary(cc, wv["w_in"], wv["w_out"])
        pen1 = np.sqrt(2 * np.log(n)) * sigma[None, :] * np.median(scale) / scale[:, None]
        X, info = bm.solve(B, np.eye(n), D, sigma, scale, K, LAM, 0.0, live, start=prev); prev = X
        res = D - B @ X
        J[j] = float(np.sum((0.5 * np.sum(res ** 2, axis=0) + LAM * np.sum(pen1 * np.abs(X), axis=0))[live] / s2[live]))
        C[j] = bool(info["converged"])
    jm = int(np.argmin(J)); edge = jm in (0, len(J) - 1); v = float(VG[jm])
    if not edge:
        idx = [jm - 1, jm, jm + 1]; a = np.polyfit(VG[idx], J[idx], 2)
        if a[0] > 0:
            v = float(np.clip(-a[1] / (2 * a[0]), VG[jm + 1], VG[jm - 1]))
    trust = bool(all(C[max(0, jm - 1):jm + 2]))
    return dict(v_grid=float(VG[jm]), v_sel=v, edge=bool(edge), neighbourhood_converged=trust,
                J=J.tolist(), converged=C.tolist())


def run(prof, did, c, s, sc, base):
    ctx = context(base, prof); cm = ctx["cm"]; Q = np.asarray(cm["Q"], float)
    g = gen.signal(ctx, c, s); D = g["x"] + MULT * g["level"] * ctx["noise"][0]; pp = sc11.prepare_sc(D, cm)
    t0 = time.monotonic()
    if sc == "w":
        ca, res = None, scan(D, cm, Q, pp, lambda cc: lib.wavelets(D, semi.TW, cc))
    else:
        ca = float(sc[1:]); wv = lib.wavelets(D, semi.TW, dict(cm, Q=Q * ca)); res = scan(D, cm, Q, pp, lambda cc: wv)
    return dict(profile=prof, Q_profile=PROFILES[prof], id=did, scan=sc, c_true=c, u=1.0 / c, seed=s, c_anchor=ca,
                v_anchor=(1.0 / ca if ca else None), noise_mult=MULT, lam=LAM, level=g["level"],
                seconds=time.monotonic() - t0, **res)


if __name__ == "__main__":
    OUT.mkdir(exist_ok=True); base = gen.prepare(STRUCT, NOISE)
    if "--smoke" in sys.argv:          # development seed 9940, not a test dataset; timing only
        rec = run("D", "c1.0_s9940", 1.0, 9940, "a1.0", base); rec.pop("J"); rec.pop("converged")
        print(json.dumps(rec, indent=1)); sys.exit()
    for prof, did, c, s, sc in jobs():
        job = f"{prof}_{did}__{sc}"; path = OUT / f"{job}.json"
        if path.exists():
            continue
        try:
            os.close(os.open(OUT / f"{job}.claim", os.O_CREAT | os.O_EXCL | os.O_WRONLY))
        except FileExistsError:
            continue
        rec = run(prof, did, c, s, sc, base)
        path.write_text(json.dumps(rec, indent=1), encoding="utf-8"); (OUT / f"{job}.claim").unlink()
        print(job, f"v {rec['v_sel']:.3f}", "edge" if rec["edge"] else "", f"{rec['seconds']:.0f}s", flush=True)
