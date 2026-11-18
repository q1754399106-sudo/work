"""Development: noise-free M1, AVO off, c = 1, 4 seeds: bias of the estimator variants vs wavelet peak frequency."""
import sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
import oqa as o
Q = np.array(o.MODELS[sys.argv[1]], float); a0m, axm, P = o.rays(Q)
for fp in (30.0, 40.0, 60.0):
    out = {}
    for seed in range(4):
        rng = np.random.default_rng(777 + seed); A, Bg = o.reflectivity(rng, False); G = o.gathers(A, Bg, Q, 1.0, fp, P)
        for name, eq, sm, xw in (("base", False, 9, False), ("eq_ns", True, 1, False), ("xw", False, 9, True), ("eq_xw", True, 9, True), ("eq_ns_xw", True, 1, True)):
            o.EQLEN, o.SMOOTH, o.XWIN = eq, sm, xw; ce, _ = o.c_off(G, np.zeros(o.NX), a0m, axm, o.vrms); out.setdefault(name, []).append(np.log(ce))
    print(sys.argv[1], "fp %3.0f  " % fp + "  ".join("%s %+.3f" % (k, np.mean(v)) for k, v in out.items()), flush=True)
