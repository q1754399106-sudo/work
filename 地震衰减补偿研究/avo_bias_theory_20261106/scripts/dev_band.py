"""Development on used M1: band-limited baseline. Noise-free sanity (no cut) then scores for limit rules; the rule that does best for it on E_G is frozen, to give it its best chance."""
import numpy as np
import step2b as s
from step2b import NEW, columns, metrics, spec
import gate1 as g
from gate1 import n, NX, T, TW, reflectors, synth
import theory as th, bandlim as bl

NEW["M1"] = np.array([40.0, 50, 70, 100]); s.SEEDS = 3; FP = 60.0; g.RIC = spec(FP); th.S2 = g.RIC ** 2; Q = NEW["M1"]
idx, R = reflectors("events_aligned", np.random.default_rng(0)); dp, tru = synth(idx, R, Q, False); O = bl.apply(dp[:n], Q, np.full(NX, 1e-9), nocut=True); kw = (T > .3) & (T < 1.0); Tb = s.bp(tru, (10.0, 60.0)); Ob = s.bp(O, (10.0, 60.0))
print("noise-free, no cut: relative error 10-60 Hz", round(float(np.linalg.norm((Ob - Tb)[kw]) / np.linalg.norm(Tb[kw])), 4), "gain near", round(float(np.sum(Ob[kw, 0] * Tb[kw, 0]) / np.sum(Tb[kw, 0] ** 2)), 4), "gain far", round(float(np.sum(Ob[kw, -1] * Tb[kw, -1]) / np.sum(Tb[kw, -1] ** 2)), 4))
K = ("ga", "eg", "egl", "mg", "cg", "wave")
for rule, kk in (("peak", 0), ("snr", 1.5), ("snr", 2.0), ("snr", 3.0), ("snr", 5.0)):
    print("\nrule", rule, kk)
    for sc, qf in (("events_aligned", 1.0), ("thin_aligned", 1.0), ("dense_aligned", 1.0), ("dense_aligned", 0.8), ("dense_aligned", 1.2), ("dense_nmo", 1.0)):
        D, Cl, Tr, sig, meta, keep = columns("M1", sc, FP); O = np.stack([bl.apply(D[:, :, c], Q * qf, sig[:, c], rule, kk) for c in range(D.shape[2])], 2); mt = metrics(O, O, Tr, keep, FP)
        print(f"  {sc:15s} Q x{qf} " + " ".join(f"{k} {np.mean(mt[k]):.3f}" for k in K))
