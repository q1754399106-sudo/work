"""Development (before PLAN_BOUND_ZH.md is frozen): is kappa on profile D an amplitude (nonlinearity) effect?
Selection v_sel against the anchor v_a on one dataset per profile (development seed 9940, c_true = 1), anchors
c_a = 0.7, 0.8, 0.9, 1.0, 1.1, 1.25, 1.4. The coarse scan (qprofile.scan) is refined on a local grid in v (step 0.02,
+-0.16 around the coarse vertex) and the selection is the vertex of the parabola through the fine minimum.

usage: dev_anchor_curve.py <profile>      -> analysis/dev_anchor_curve_<profile>.json
"""
import sys, json, time
from env15 import np, lib, semi, sc11, bm, gen
from qprofile import scan

STRUCT, NOISE, MULT, LAM, SEED = 10227, [10149, 10273, 10355], 0.5, 0.5, 9940
PROF = {"base": None, "S": [30.0, 200.0, 200.0, 200.0], "D": [250.0, 250.0, 40.0, 40.0],
        "D2": [150.0, 150.0, 60.0, 60.0], "M": [200.0, 45.0, 200.0, 200.0]}
ANCHORS = (0.7, 0.8, 0.9, 1.0, 1.1, 1.25, 1.4)
bm.MAXITER = 20000


def objective(D, cm, Q, pp, wv, c, prev):
    sigma, K, live = pp["sigma"], pp["K"], pp["live"]; n = D.shape[0]; s2 = np.where(live, sigma ** 2, np.inf)
    cc = dict(cm); cc["Q"] = Q * c; B, _, scale = lib.dictionary(cc, wv["w_in"], wv["w_out"])
    pen1 = np.sqrt(2 * np.log(n)) * sigma[None, :] * np.median(scale) / scale[:, None]
    X, info = bm.solve(B, np.eye(n), D, sigma, scale, K, LAM, 0.0, live, start=prev); res = D - B @ X
    J = float(np.sum((0.5 * np.sum(res ** 2, axis=0) + LAM * np.sum(pen1 * np.abs(X), axis=0))[live] / s2[live]))
    return J, X, bool(info["converged"])


def fine(D, cm, Q, pp, wv, v0):
    vs = np.round(np.arange(v0 - 0.16, v0 + 0.1601, 0.02), 4); vs = vs[vs > 0.3]; J = []; conv = []; prev = None
    for v in vs[::-1]:
        j, prev, ok = objective(D, cm, Q, pp, wv, 1.0 / v, prev); J.append(j); conv.append(ok)
    J = np.array(J[::-1]); conv = conv[::-1]; m = int(np.argmin(J)); vsel = float(vs[m]); edge = m in (0, len(vs) - 1)
    if not edge:
        a = np.polyfit(vs[m - 1:m + 2], J[m - 1:m + 2], 2)
        if a[0] > 0:
            vsel = float(-a[1] / (2 * a[0]))
    return dict(v_fine=vs.tolist(), J_fine=J.tolist(), v_sel=vsel, edge=edge, converged=all(conv))


if __name__ == "__main__":
    p = sys.argv[1]; base = gen.prepare(STRUCT, NOISE); cm = dict(base["cm"])
    if PROF[p] is not None:
        cm["Q"] = np.asarray(PROF[p], float)
    ctx = dict(base, cm=cm, w_in=lib.wavelets(base["D"], semi.TW, cm)["w_in"]); Q = np.asarray(cm["Q"], float)
    g = gen.signal(ctx, 1.0, SEED); D = g["x"] + MULT * g["level"] * ctx["noise"][0]; pp = sc11.prepare_sc(D, cm)
    out = []
    for ca in ANCHORS:
        t0 = time.monotonic(); wv = lib.wavelets(D, semi.TW, dict(cm, Q=Q * ca)); co = scan(D, cm, Q, pp, lambda cc: wv)
        fr = fine(D, cm, Q, pp, wv, co["v_sel"]); rec = dict(c_anchor=ca, v_anchor=1 / ca, v_coarse=co["v_sel"], **fr,
                                                              seconds=time.monotonic() - t0)
        out.append(rec); print(p, ca, f"coarse {co['v_sel']:.3f} fine {fr['v_sel']:.3f}", "edge" if fr["edge"] else "",
                               f"{rec['seconds']:.0f}s", flush=True)
    va = np.array([r["v_anchor"] for r in out]); vs = np.array([r["v_sel"] for r in out])
    print("local slopes dv_sel/dv_a:", np.round(np.diff(vs) / np.diff(va), 3).tolist())
    json.dump(dict(profile=p, seed=SEED, c_true=1.0, rows=out), open(f"../analysis/dev_anchor_curve_{p}.json", "w"), indent=1)
