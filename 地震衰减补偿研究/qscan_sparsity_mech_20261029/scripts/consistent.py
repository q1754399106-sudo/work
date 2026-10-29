"""Prototype 2 (exploratory): QScan with the input wavelet re-estimated for every scanned factor.

lib.wavelets estimates w_in as the shallow (1.2-1.6 s) spectrum divided by the model absorption at 1.4 s, with the
model Q (factor 1). lib.qscan builds every dictionary B(c) from that one w_in, so for c != 1 the atoms at the
shallow window no longer reproduce the shallow spectrum: the scan is anchored at c = 1. Hypothesis: this anchoring,
not the reflectivity, compresses the selection. Here w_in(c) = lib.wavelets(D, TW, cm with Q x c) for each c;
objective, solver, weights, grid and sigma/K are those of lib.qscan. Weights 0.5 (lib.qscan) and 0.125.
Data: the 20 held-out sections of experiment 12 (as the two pilots).

usage: consistent.py K N [lam,lam,...]    worker K of N; writes consistent/<id>.json, or with weights given
       consistent_lams/<id>.json (second run: 2^-4, 2^-6)
"""
import sys, json, time
from env13 import np, lib, semi, sc11, bm, ROOT, datasets

LAMS = [0.5, 0.125]
OUT = ROOT / "consistent"


def scan(D, cm, Q, pp):
    sigma, K, live = pp["sigma"], pp["K"], pp["live"]; n = D.shape[0]
    J = np.zeros((len(LAMS), len(lib.QGRID))); C = np.zeros(J.shape, bool)
    s2 = np.where(live, sigma ** 2, np.inf)
    for j, c in enumerate(lib.QGRID):
        cc = dict(cm); cc["Q"] = np.asarray(Q, float) * c
        wv = lib.wavelets(D, semi.TW, cc)
        B, _, scale = lib.dictionary(cc, wv["w_in"], wv["w_out"])
        pen1 = np.sqrt(2 * np.log(n)) * sigma[None, :] * np.median(scale) / scale[:, None]
        for i, lam in enumerate(LAMS):
            X, info = bm.solve(B, np.eye(n), D, sigma, scale, K, lam, 0.0, live)
            res = D - B @ X
            J[i, j] = float(np.sum((0.5 * np.sum(res ** 2, axis=0) + lam * np.sum(pen1 * np.abs(X), axis=0))[live] / s2[live]))
            C[i, j] = bool(info["converged"])
    return J, C


def main(k, n):
    global LAMS, OUT
    if len(sys.argv) > 3:
        LAMS = [float(eval(x)) for x in sys.argv[3].split(",")]; OUT = ROOT / "consistent_lams"
    OUT.mkdir(parents=True, exist_ok=True); cm = semi.curtain_cm(); Q = np.asarray(cm["Q"], float)
    for i, p in enumerate(datasets()):
        did = p.stem; path = OUT / f"{did}.json"
        if i % n != k or path.exists():
            continue
        t0 = time.monotonic(); D = np.load(p)["Input"].astype(float)
        pp = sc11.prepare_sc(D, cm); J, C = scan(D, cm, Q, pp)
        sel = {str(lam): float(lib.QGRID[int(np.argmin(J[i]))]) for i, lam in enumerate(LAMS)}
        rec = dict(id=did, c_true=float(did.split("_")[1][1:]), grid=lib.QGRID.tolist(), lams=LAMS, selected=sel,
                   J=J.tolist(), converged=C.tolist(), seconds=time.monotonic() - t0)
        path.write_text(json.dumps(rec, indent=1), encoding="utf-8"); print(did, f"{rec['seconds']:.0f}s", flush=True)


if __name__ == "__main__":
    main(int(sys.argv[1]), int(sys.argv[2]))
