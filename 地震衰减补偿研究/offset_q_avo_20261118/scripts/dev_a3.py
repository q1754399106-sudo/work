"""Development for gate A3 (used models M24-M27, the gate-A2 data and c_hat): with the SAME calibrated Q (model Q x c_hat),
compare the saturated damping (gate A2 'cal', read back) with one-step Tikhonov, the split method, band-limited inverse Q
(Ratnett et al. 2019 idea, bandlim.py, SNR rule k = 3) and vertical-only Tikhonov (offset path ignored, a_x = 0).
Seeds 0-2, events/thin aligned, 24 and 30 dB. usage: dev_a3.py <model>"""
import sys, json, time
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent; sys.path.insert(0, str(HERE)); sys.path.insert(0, str(HERE / "legacy2b"))
import e2e  # noqa: E402
from e2e import s, np  # noqa: E402
import bandlim  # noqa: E402

M = sys.argv[1]; Qm = np.array(e2e.MODELS[M], float); mnum = int(M[1:]); SEEDS = range(3); rows = []; t0 = time.monotonic()
A2 = [r for fp in (60, 40) for r in json.loads((e2e.ROOT / "results" / ("a2_%s_f%d.json" % (M, fp))).read_text())]
CH = {(r["fp"], r["c"], r["snr"], r["seed"]): r["c_hat"] for r in A2 if r["method"] == "cal"}


def solve_generic(D, Cl, sig, keep, Qop, method):
    """one-step 'direct', 'split' (step 2b: offset part with LX then a0 with lam), 'tik0' (a0 only for every offset)."""
    nc = D.shape[2]; a0, ax = s.a0_ax(np.asarray(Qop, float)); RH = np.concatenate([D, Cl], axis=2); n = s.n; I = s.I
    Gd = [s.atten(a0 + ax[:, i]) for i in range(s.NX)]; target = np.sum(sig ** 2, axis=0) * n; ch = [None] * nc; out = np.zeros_like(RH)
    if method == "tik0":
        G0 = s.atten(a0); GtG0 = G0.T @ G0
    if method == "split":
        G0 = s.atten(a0); Gx = [s.atten(ax[:, i]) for i in range(s.NX)]
        Y = np.stack([np.linalg.solve(G.T @ G + s.LX * I, G.T @ RH[:, i, :]) for i, G in enumerate(Gx)], axis=1)
    else:
        GtG = [G.T @ G for G in Gd]; GtR = [Gd[i].T @ RH[:, i, :] for i in range(s.NX)]
    for lam in s.GRID:
        if method == "direct":
            O = np.stack([np.linalg.solve(GtG[i] + lam * I, GtR[i]) for i in range(s.NX)], axis=1)
        elif method == "tik0":
            O = np.stack([np.linalg.solve(GtG0 + lam * I, G0.T @ RH[:, i, :]) for i in range(s.NX)], axis=1)
        else:
            F0 = np.linalg.solve(G0.T @ G0 + lam * I, G0.T); O = np.stack([F0 @ Y[:, i, :] for i in range(s.NX)], axis=1)
        res = sum(np.sum(((Gd[i] @ O[:, i, :nc] - D[:, i, :]) * keep[:, i, None]) ** 2, axis=0) for i in range(s.NX))
        for c in range(nc):
            if ch[c] is None and res[c] <= target[c]:
                ch[c] = float(lam); out[:, :, c] = O[:, :, c]; out[:, :, nc + c] = O[:, :, nc + c]
        if all(v is not None for v in ch):
            break
    for c in range(nc):
        if ch[c] is None:
            out[:, :, c] = O[:, :, c]; out[:, :, nc + c] = O[:, :, nc + c]
    return out[:, :, :nc], out[:, :, nc:]


for fp in (60.0, 40.0):
    s.g.RIC = s.spec(fp); s.th.S2 = s.g.RIC ** 2
    for c in e2e.CS:
        for scene in ("events_aligned", "thin_aligned"):
            D, Cl, Tr, sig, meta, keep = e2e.columns(c * Qm, mnum, scene, fp, SEEDS); nc = D.shape[2]
            outs = {m: (np.zeros_like(D), np.zeros_like(D)) for m in ("direct", "split", "tik0", "band3")}
            for k, (snr, seed) in enumerate(meta):
                ce = CH[(fp, c, snr, seed)]; ce = ce if (ce == ce and ce > 0) else 1.0; Qc = ce * Qm; sl = slice(k, k + 1)
                for m in ("direct", "split", "tik0"):
                    a, b = solve_generic(D[:, :, sl], Cl[:, :, sl], sig[:, sl], keep, Qc, m); outs[m][0][:, :, k] = a[:, :, 0]; outs[m][1][:, :, k] = b[:, :, 0]
                outs["band3"][0][:, :, k] = bandlim.apply(D[:, :, k], Qc, sig[:, k], "snr", 3.0); outs["band3"][1][:, :, k] = bandlim.apply(Cl[:, :, k], Qc, sig[:, k], "snr", 3.0)
            for m, (O, Oc) in outs.items():
                q = s.metrics(O, Oc, Tr, keep, fp)
                for k, (snr, seed) in enumerate(meta):
                    rows.append(dict(model=M, fp=fp, scene=scene, c=c, snr=snr, seed=seed, method=m + "_cal", **{a: float(v[k]) for a, v in q.items()}))
            for r in A2:
                if r["fp"] == fp and r["c"] == c and r["scene"] == scene and r["seed"] in SEEDS:
                    rows.append(dict(r, method=r["method"] if r["method"] != "cal" else "sat_cal"))
            print(M, fp, c, scene, "%.0fs" % (time.monotonic() - t0), flush=True)
            (e2e.ROOT / "results" / ("dev_a3_%s.json" % M)).write_text(json.dumps(rows))
