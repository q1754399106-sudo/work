"""Development (M16-M19 are used after gate A1): base vs EQLEN vs EQLEN + no smoothing, 24 dB, AVO, N = 50 and 200,
seeds 0-1 with c 0.7/1.0/1.4, 40/60 Hz. usage: dev_fix.py <model>"""
import sys, json
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
import oqa as o

M = sys.argv[1]; Q = np.array(o.MODELS[M], float); a0m, axm, P = o.rays(Q); rows = []; mi = list(o.MODELS).index(M)
for fp in o.FPS:
    for c in o.CS:
        for seed in range(2):
            rng = np.random.default_rng(31_000_000 + mi * 100_000 + int(fp) * 1000 + int(c * 10) * 10 + seed)
            A, Bg = o.reflectivity(rng, True); G0 = o.gathers(A, Bg, Q, c, fp, P); rms = np.array([np.sqrt(np.mean(g ** 2)) for g in G0])
            nrng = np.random.default_rng(rng.integers(1 << 62) + 24); sig = rms / 10 ** (24 / 20); G = [g + s * nrng.standard_normal(g.shape) for g, s in zip(G0, sig)]
            for name, eq, sm, xw in (("base", False, 9, False), ("eq", True, 9, False), ("eq_ns", True, 1, False), ("eq_xw", True, 9, True), ("eq_ns_xw", True, 1, True), ("xw", False, 9, True)):
                o.EQLEN, o.SMOOTH, o.XWIN = eq, sm, xw
                for n in (50, 200):
                    for k in range(o.CMPS // n):
                        ce, _ = o.c_off([g[:, k * n:(k + 1) * n] for g in G], sig, a0m, axm, o.vrms)
                        rows.append(dict(model=M, fp=fp, c=c, seed=seed, est=name, n=n, l=float(np.log(ce / c)) if ce == ce and ce > 0 else None))
(o.ROOT / "results" / ("devfix_%s.json" % M)).write_text(json.dumps(rows))
