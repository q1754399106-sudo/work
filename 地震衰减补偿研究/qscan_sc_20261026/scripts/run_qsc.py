"""QScan + SC on H1-H6 (PLAN_ZH.md "数据与方法").

For every model x scene x replicate, the archived QScan choice of the fd2d sweep (factor =
info.factors[best], index 0 = nominal Q) scales the used Q; SC is then solved exactly as in
sc_neutral_20261026/scripts/run_sc.py (benchmark.solve/predictor and common.dictionary of ff95a69,
a=0.125, mu_relative=10) on the sc_neutral rebuilt gather. Every local input is checked against the
SHA-256 recorded in cache/SNAPSHOT.json before use. Truth is not read.
Outputs: cache/solved/<model>_<scene>_r<r>_QSC.npz and .json.

usage: run_qsc.py H1 H2 ...   (primary r=1..6 first, then r=0; resumable)
"""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RESEARCH = ROOT.parent
SNAP = ROOT / "cache" / "inputs" / "novelty_evidence_20261025"
sys.path.insert(0, str(RESEARCH / "novelty_evidence_20261025" / "cache" / "deps"))
sys.path.insert(0, str(SNAP / "scripts"))
import benchmark as bm  # noqa: E402  (common.py pins BLAS to one thread before NumPy loads)
import hashlib, json, time  # noqa: E402,E401
import numpy as np  # noqa: E402
from scipy.io import loadmat  # noqa: E402

SEL = json.loads((SNAP / "analysis" / "selection.json").read_text(encoding="utf-8"))
assert (SEL["a"], SEL["mu_relative"]) == (0.125, 10.0)
SCENES = {"Q100_14dB": (1, 1.0), "Q100_08dB": (2, 1.0), "Q070_14dB": (1, 0.7)}
LOCAL = json.loads((ROOT / "cache" / "SNAPSHOT.json").read_text(encoding="utf-8"))["local_sha256"]


def checked(rel):
    p = RESEARCH / rel
    assert hashlib.sha256(p.read_bytes()).hexdigest() == LOCAL[rel], f"input changed: {rel}"
    return p


def main(models):
    out = ROOT / "cache" / "solved"
    out.mkdir(parents=True, exist_ok=True)
    cache, started = {}, time.monotonic()
    for rs in (range(1, 7), [0]):
        for s, (level, qs) in SCENES.items():
            for m in models:
                for r in rs:
                    f = out / f"{m}_{s}_r{r}_QSC"
                    if f.with_suffix(".npz").exists() and f.with_suffix(".json").exists():
                        continue
                    L = loadmat(checked(f"fd2d_20261025/cache/sweep/{m}_{s}/r{r:03d}.mat"), simplify_cells=True)["info"]
                    best = int(L["best"])
                    factor = float(L["factors"][best])
                    key = (m, qs, best)
                    if key not in cache:
                        cm = dict(loadmat(ROOT / "cache/inputs/fd2d_20261025/models" / f"{m}.mat", simplify_cells=True)["config"])
                        cm["Q"] = np.asarray(cm["Q"], float) * qs * factor
                        cache[key] = bm.dictionary(cm)
                    B, W, P, scale = cache[key]
                    g = loadmat(checked(f"sc_neutral_20261026/cache/noisy/{m}_L{level}_r{r}.mat"))
                    D = np.asarray(g["D"], float)
                    sig = bm.vec(g["sig"])
                    live = D.any(axis=0)
                    K, _, adj = bm.predictor(D, live)
                    X, info = bm.solve(B, P, D, sig, scale, K, SEL["a"], SEL["mu_relative"], live)
                    info.update(model=m, scene=s, r=r, label="QSC", qscan_best=best, qscan_factor=factor,
                                predictor_adjoint_error=adj)
                    np.savez_compressed(f.with_suffix(".npz"), out=X)
                    bm.writejson(f.with_suffix(".json"), info)
                    print(f"{m} {s} r{r} QSC factor {factor:.4f}: iter {info['iterations']} KKT {info['kkt_ratio']:.4f} "
                          f"conv {info['converged']} {info['seconds']:.0f}s (total {time.monotonic() - started:.0f}s)", flush=True)
    print("DONE", " ".join(models), flush=True)


if __name__ == "__main__":
    main(sys.argv[1:])
