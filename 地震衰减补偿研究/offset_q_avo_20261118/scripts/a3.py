"""Gate A3 (PLAN_GATE_A3_ZH.md): with the SAME calibrated Q (model Q x c_hat), do path-aware compensation operators tie on AVO,
and does vertical-only compensation ruin it? New models M28-M31, new c_hat seeds. Operators: saturated damping (e2e.solve 'sat'),
one-step Tikhonov, split, vertical-only Tikhonov (dev_a3.solve_generic, unchanged), band-limited inverse Q (bandlim.py, SNR rule
k = 3, unchanged). lambda by the residual principle (sigma known); band3 has no lambda.
usage: a3.py run <model> <fp>"""
import sys, json, time
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent; sys.path.insert(0, str(HERE)); sys.path.insert(0, str(HERE / "legacy2b"))
import e2e  # noqa: E402
from e2e import s  # noqa: E402
import bandlim  # noqa: E402

e2e.EST_BASE = 61_000_000; SEEDS = range(6); SCENES = ("events_aligned", "thin_aligned")
MODELS = {"M28": [55, 85, 40, 70], "M29": [33, 60, 95, 45], "M30": [70, 38, 50, 105], "M31": [48, 110, 36, 58]}


def solve_generic(D, Cl, sig, keep, Qop, method):
    """Copied from dev_a3.py (that file runs its development loop on import)."""
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


def run(M, fp):
    Qm = np.array(MODELS[M], float); mnum = int(M[1:]); rows = []; t0 = time.monotonic()
    s.g.RIC = s.spec(fp); s.th.S2 = s.g.RIC ** 2
    for c in e2e.CS:
        ch = {(snr, seed): e2e.c_hat(Qm, mnum, c, fp, snr, seed) for seed in SEEDS for snr in e2e.SNRS}
        for scene in SCENES:
            D, Cl, Tr, sig, meta, keep = e2e.columns(c * Qm, mnum, scene, fp, SEEDS)
            outs = {m: (np.zeros_like(D), np.zeros_like(D)) for m in ("sat", "direct", "split", "tik0", "band3")}
            for k, (snr, seed) in enumerate(meta):
                ce = ch[(snr, seed)]; ce = ce if (ce == ce and ce > 0) else 1.0; Qc = ce * Qm; sl = slice(k, k + 1)
                a, b, _ = e2e.solve(D[:, :, sl], Cl[:, :, sl], sig[:, sl], keep, Qc, "sat"); outs["sat"][0][:, :, k] = a[:, :, 0]; outs["sat"][1][:, :, k] = b[:, :, 0]
                for m in ("direct", "split", "tik0"):
                    a, b = solve_generic(D[:, :, sl], Cl[:, :, sl], sig[:, sl], keep, Qc, m); outs[m][0][:, :, k] = a[:, :, 0]; outs[m][1][:, :, k] = b[:, :, 0]
                outs["band3"][0][:, :, k] = bandlim.apply(D[:, :, k], Qc, sig[:, k], "snr", 3.0); outs["band3"][1][:, :, k] = bandlim.apply(Cl[:, :, k], Qc, sig[:, k], "snr", 3.0)
            for m, (O, Oc) in outs.items():
                q = s.metrics(O, Oc, Tr, keep, fp)
                for k, (snr, seed) in enumerate(meta):
                    rows.append(dict(model=M, fp=fp, scene=scene, c=c, snr=snr, seed=seed, method=m, c_hat=ch[(snr, seed)], **{a: float(v[k]) for a, v in q.items()}))
            print(M, fp, c, scene, "%.0fs" % (time.monotonic() - t0), flush=True)
            (e2e.ROOT / "results" / ("a3_%s_f%d.json" % (M, int(fp)))).write_text(json.dumps(rows))


if __name__ == "__main__":
    run(sys.argv[2], float(sys.argv[3]))
