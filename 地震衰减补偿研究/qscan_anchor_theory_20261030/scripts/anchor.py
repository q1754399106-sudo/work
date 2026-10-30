"""Anchor experiment (PLAN_ZH.md): QScan with the input wavelet estimated at a fixed Q factor c_a for all scanned c.

usage: anchor.py    worker; takes (dataset, anchor) jobs in order, claims each with results/<job>.claim, resumable.
Objective, sigma, K, grid as lib.qscan; weights 0.5 and 1/16; FISTA limit 20000, warm start along the grid per weight.
-> results/<id>__a<c_a>.json
"""
import os, json, time
from env15 import np, lib, semi, sc11, bm, ROOT, E12OUT

LAMS = [0.5, 2 ** -4]; ANCHORS = [0.7, 1.0, 1.4]; MAXITER = 20000
IDS = [f"T_c{c}_s{s}" for c in (0.65, 0.85, 1.05, 1.25, 1.5) for s in (9701, 9702)]
OUT = ROOT / "results"; bm.MAXITER = MAXITER


def scan(D, cm, Q, ca, pp):
    sigma, K, live = pp["sigma"], pp["K"], pp["live"]; n = D.shape[0]
    ca_cm = dict(cm); ca_cm["Q"] = Q * ca; wv = lib.wavelets(D, semi.TW, ca_cm)
    J = np.zeros((len(LAMS), len(lib.QGRID))); C = np.zeros(J.shape, bool); IT = np.zeros(J.shape, int)
    s2 = np.where(live, sigma ** 2, np.inf); prev = [None] * len(LAMS)
    for j, c in enumerate(lib.QGRID):
        cc = dict(cm); cc["Q"] = Q * c
        B, _, scale = lib.dictionary(cc, wv["w_in"], wv["w_out"])
        pen1 = np.sqrt(2 * np.log(n)) * sigma[None, :] * np.median(scale) / scale[:, None]
        for i, lam in enumerate(LAMS):
            X, info = bm.solve(B, np.eye(n), D, sigma, scale, K, lam, 0.0, live, start=prev[i]); prev[i] = X
            res = D - B @ X
            J[i, j] = float(np.sum((0.5 * np.sum(res ** 2, axis=0) + lam * np.sum(pen1 * np.abs(X), axis=0))[live] / s2[live]))
            C[i, j], IT[i, j] = bool(info["converged"]), int(info["iterations"])
    sel = {str(lam): float(lib.QGRID[int(np.argmin(J[i]))]) for i, lam in enumerate(LAMS)}
    return sel, dict(J=J.tolist(), converged=C.tolist(), iterations=IT.tolist())


if __name__ == "__main__":
    OUT.mkdir(exist_ok=True); cm = semi.curtain_cm(); Q = np.asarray(cm["Q"], float)
    for did in IDS:
        for ca in ANCHORS:
            job = f"{did}__a{ca}"; path = OUT / f"{job}.json"
            if path.exists():
                continue
            try:
                os.close(os.open(OUT / f"{job}.claim", os.O_CREAT | os.O_EXCL | os.O_WRONLY))
            except FileExistsError:
                continue
            t0 = time.monotonic(); D = np.load(E12OUT / f"{did}.npz")["Input"].astype(float); pp = sc11.prepare_sc(D, cm)
            sel, diag = scan(D, cm, Q, ca, pp)
            rec = dict(id=did, c_true=float(did.split("_")[1][1:]), c_anchor=ca, lams=LAMS, qscan_a=sel, diag=diag,
                       seconds=time.monotonic() - t0)
            path.write_text(json.dumps(rec, indent=1), encoding="utf-8"); (OUT / f"{job}.claim").unlink()
            print(job, sel, f"{rec['seconds']:.0f}s", flush=True)
