"""Pilot (exploratory, descriptive): does the Q factor QScan selects depend on the lasso weight of its objective?

QScan (lib.qscan of experiment 7, used unchanged in experiments 7-13) minimises the lasso objective with the fixed
weight SEL_LASSO = 0.5 over QGRID. Here the same objective is evaluated for several weights LAMS on the 20
held-out semi-synthetic sections of experiment 12 (true c known, 0.65-1.5), reusing their float32 copies in
baselines_volve_20261029/cache/out. LAMS includes 0.5, so the selection of lib.qscan is reproduced as a check.

Imported unchanged from field_baselines_volve_20261029/cache/inputs (the experiment-13 snapshot): lib (dictionary,
QGRID, SEL_LASSO), benchmark (lasso solver), sc11.prepare_sc (wavelets, sigma, predictor), semi.curtain_cm.

usage: pilot.py K N    worker K of N; writes results/<id>.json
"""
import sys, json, time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
E13 = ROOT.parent / "field_baselines_volve_20261029"; IN = E13 / "cache" / "inputs"
E12OUT = ROOT.parent / "baselines_volve_20261029" / "cache" / "out"
sys.path.insert(0, str(ROOT.parent / "novelty_evidence_20261025" / "cache" / "deps"))   # scipy
for p in ("novelty_evidence_20261025/scripts", "field_application_20261028/scripts", "semisynthetic_volve_20261028/scripts",
          "specslope_q_volve_20261028/scripts", "multiwell_tie_volve_20261028/scripts", "srsc_volve_20261028/scripts"):
    sys.path.insert(0, str(IN / p))
import sc11  # noqa: E402
from sc11 import np, lib, semi, sq, bm  # noqa: E402

for mod in (sc11, sq, semi):
    mod.ROOT, mod.IN = E13, IN
sc11.wq.ROOT, sc11.wq.IN, sc11.wq.E9 = E13, IN, IN / "multiwell_tie_volve_20261028"

LAMS = [2.0, 1.0, 0.5, 0.25, 0.125, 2 ** -5, 2 ** -7]          # 2^-9 does not converge in 5000 FISTA steps
OUT = ROOT / "results"


def scan(D, cm, Q, pp):
    """Lasso objective of lib.qscan for every weight in LAMS and every factor in QGRID (same formula)."""
    wv, sigma, K, live = pp["wv"], pp["sigma"], pp["K"], pp["live"]; n = D.shape[0]
    J = np.zeros((len(LAMS), len(lib.QGRID))); R = np.zeros_like(J); L1 = np.zeros_like(J); C = np.zeros(J.shape, bool)
    s2 = np.where(live, sigma ** 2, np.inf)
    for j, c in enumerate(lib.QGRID):
        cc = dict(cm); cc["Q"] = np.asarray(Q, float) * c
        B, P, scale = lib.dictionary(cc, wv["w_in"], wv["w_out"])
        for i, lam in enumerate(LAMS):
            X, info = bm.solve(B, np.eye(n), D, sigma, scale, K, lam, 0.0, live)
            res = D - B @ X
            w = np.sqrt(2 * np.log(n)) * sigma[None, :] * np.median(scale) / scale[:, None]
            r = float(np.sum((0.5 * np.sum(res ** 2, axis=0))[live] / s2[live]))
            l1 = float(np.sum(np.sum(w * np.abs(X), axis=0)[live] / s2[live]))
            R[i, j], L1[i, j], J[i, j], C[i, j] = r, l1, r + lam * l1, bool(info["converged"])
    return J, R, L1, C


def main(k, n):
    OUT.mkdir(parents=True, exist_ok=True); cm = semi.curtain_cm(); Q = np.asarray(cm["Q"], float)
    for i, p in enumerate(sorted(q for q in E12OUT.glob("T_*.npz") if "__" not in q.name)):
        did = p.stem; path = OUT / f"{did}.json"
        if i % n != k or path.exists():
            continue
        t0 = time.monotonic(); D = np.load(p)["Input"].astype(float)
        pp = sc11.prepare_sc(D, cm); J, R, L1, C = scan(D, cm, Q, pp)
        sel = {str(lam): float(lib.QGRID[int(np.argmin(J[i]))]) for i, lam in enumerate(LAMS)}
        rec = dict(id=did, c_true=float(did.split("_")[1][1:]), grid=lib.QGRID.tolist(), lams=LAMS, selected=sel,
                   J=J.tolist(), residual=R.tolist(), l1=L1.tolist(), converged=C.tolist(), seconds=time.monotonic() - t0)
        path.write_text(json.dumps(rec, indent=1), encoding="utf-8"); print(did, f"{rec['seconds']:.0f}s", flush=True)


if __name__ == "__main__":
    main(int(sys.argv[1]), int(sys.argv[2]))
