"""Exploratory baseline for pilot_mw.py: the literature multitrace method SC (structure-constrained sparse inversion,
benchmark.solve, imported unchanged) on the same noisy FD gathers, scored with the same S.

  SC_orig   a 0.125, mu_relative 10   (selection of novelty_evidence_20261025, used in sc_neutral_20261026)
  SC_tuned  a 2^-13, mu_relative 0.25 (experiments 12 and 17, tuned on the Volve semi-synthetic data)
Noise level per trace: standard deviation of the actual noise (favourable to SC; the methods under test estimate it).
Assumed Q = m x true Q, m in 0.5, 0.7, 1, 1.4, 2.

usage: pilot_sc.py run K N LEVEL [R ...]     realisations R (default 1)
"""
import sys, json, time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import uqp  # noqa: E402,F401  (path set-up of the experiment chain)
from uqp import np, lib, ROOT  # noqa: E402
import fdcheck as fd  # noqa: E402
from scipy.io import loadmat  # noqa: E402

bm = lib.bm; OUT = ROOT / "results" / "pilot_sc"
LEVELS = {"14dB": (14.0, 20261113000), "30dB": (30.0, 20261121000)}; M_SC = (0.5, 0.7, 1.0, 1.4, 2.0)
CONFIGS = {"SC_orig": (0.125, 10.0), "SC_tuned": (2.0 ** -13, 0.25)}


def run(name, level, r):
    D0, T, cm_, tw, fp = fd.load_fd(name); db, base = LEVELS[level]; seed = base + 1000 * (fd.MODELS.index(name) + 1) + r
    E = fd.noise(D0, db, seed); D = D0 + E; sig = E.std(axis=0); live = D.any(axis=0); K, _, _ = bm.predictor(D, live); rec = []
    cfg = dict(loadmat(fd.FD / "models" / f"{name}.mat", simplify_cells=True)["config"]); Q = np.asarray(cfg["Q"], float)
    for m in M_SC:
        c = dict(cfg); c["Q"] = Q * m; B, W, P, scale = bm.dictionary(c); d = dict(m=m)
        for lab, (a, mu) in CONFIGS.items():
            X, info = bm.solve(B, P, D, sig, scale, K, a, mu, live)
            d[lab] = float(fd.score_fd(X, T, tw)[0]); d[lab + "_conv"] = bool(info["converged"]); d[lab + "_s"] = float(info["seconds"])
        rec.append(d); print(name, level, r, d, flush=True)
    return dict(model=name, level=level, r=r, seed=seed, records=rec)


if __name__ == "__main__":
    k, n, level = int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]; rs = [int(x) for x in sys.argv[5:]] or [1]
    OUT.mkdir(parents=True, exist_ok=True); jobs = [(h, r) for r in rs for h in fd.MODELS]
    for j, (h, r) in enumerate(jobs):
        path = OUT / f"{h}_{level}_r{r}.json"
        if j % n != k or path.exists():
            continue
        t0 = time.monotonic(); out = run(h, level, r); out["seconds"] = time.monotonic() - t0
        path.write_text(json.dumps(out, indent=1), encoding="utf-8")
