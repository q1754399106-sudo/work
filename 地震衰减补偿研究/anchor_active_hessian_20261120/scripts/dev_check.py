"""Implementation checks on the development seed 9940 (profile D, c 1.0; not part of any table).

1 unnormalised dictionary = lib.dictionary B x scale; 2 solve tolerance vs J; 3 polish consistency; 4 analytic theta
Hessian against finite differences of the fixed-set closed form (random direction and the v direction); 5 timing.
"""
import sys, time
from core import np, lib, gen, qp, bm, sc11, load, wavelet, unnorm_dict, mu_of, solve, polish, jfix, theta_hessian

prof = sys.argv[1] if len(sys.argv) > 1 else "D"
base = gen.prepare(qp.STRUCT, qp.NOISE)
D, cm = load(prof, "c1.0_s9940", 1.0, 9940, base); pp = sc11.prepare_sc(D, cm); n = D.shape[0]
print("D", D.shape, "live", int(pp["live"].sum()))
w1 = wavelet(D, cm, 1.0); v0 = 1.05
Bt, Win = unnorm_dict(cm, v0, w1)
B, _, scale = lib.dictionary(dict(cm, Q=np.asarray(cm["Q"]) / v0), w1, w1)
print("1 dict max rel diff", float(np.max(np.abs(B * scale - Bt)) / np.max(np.abs(Bt))))
for tol in (0.02, 0.005, 0.001):
    t0 = time.monotonic(); y, J, info = solve(D, Bt, pp, tol=tol); print("2 tol", tol, "J %.4f" % J, info, "%.0fs" % (time.monotonic() - t0))
mu = mu_of(Bt, n); Y, A, S, pinfo = polish(D, Bt, y, pp, mu)
print("3 polish", pinfo, "Jfix %.4f" % jfix(D, Bt, A, S, pp, mu), "J %.4f" % J)
t0 = time.monotonic(); H, Hd, z2 = theta_hessian(D, cm, v0, w1, Y, A, pp); print("5 hessian %.0fs" % (time.monotonic() - t0))
a = lib.absorption(cm, bm.vec(cm["t"])); rng = np.random.default_rng(1)
for name, dirn in (("random", rng.normal(size=n) * a.mean()), ("v (a)", a), ("deep half", a * (cm["t"] > 2.1))):
    for e in (1e-2, 3e-3):
        Jp = jfix(D, unnorm_dict(cm, v0, w1, e * dirn)[0], A, S, pp, mu)
        Jm = jfix(D, unnorm_dict(cm, v0, w1, -e * dirn)[0], A, S, pp, mu)
        J0 = jfix(D, Bt, A, S, pp, mu)
        print("4", name, e, "FD %.6g" % ((Jp - 2 * J0 + Jm) / e ** 2), "analytic %.6g" % (dirn @ H @ dirn),
              "direct %.6g" % (dirn @ Hd @ dirn))
