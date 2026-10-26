"""SC and QScan + SC on the Z models (PLAN_ZH.md "方法"); same solver calls as
qscan_sc_round3_20261026/scripts/run_sc.py. SC uses the used Q = Q0 * (1 + e/100); QScan + SC uses
scanned.Q archived by run_combo3. Models are checked against logs/MODEL_HASHES.json; truth is not read.
Outputs: cache/solved/<model>_<noise>_Q<e>_<SC|QSC>.npz and .json.

usage: run_sc.py Z1_000 Z1_001 ...   (resumable; needs the run_combo3 results of these models)
"""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RESEARCH = ROOT.parent
SNAP = ROOT / "cache" / "inputs"
sys.path.insert(0, str(RESEARCH / "novelty_evidence_20261025" / "cache" / "deps"))
sys.path.insert(0, str(SNAP / "novelty_evidence_20261025" / "scripts"))
import benchmark as bm  # noqa: E402  (common.py pins BLAS to one thread before NumPy loads)
import hashlib, json, time  # noqa: E402,E401
import numpy as np  # noqa: E402
from scipy.io import loadmat  # noqa: E402

HASHES = json.loads((ROOT / "logs" / "MODEL_HASHES.json").read_text(encoding="utf-8"))
SEL = json.loads((SNAP / "novelty_evidence_20261025/analysis/selection.json").read_text(encoding="utf-8"))
assert (SEL["a"], SEL["mu_relative"]) == (0.125, 10.0)
NOISES = ["Moderate", "Strong"]
ERRORS = [-30, 0, 30]


def main(models):
    out = ROOT / "cache" / "solved"
    out.mkdir(parents=True, exist_ok=True)
    started = time.monotonic()
    for name in models:
        p = ROOT / "models" / f"{name}.mat"
        assert hashlib.sha256(p.read_bytes()).hexdigest() == HASHES[f"models/{name}.mat"], name
        dev = loadmat(p, simplify_cells=True)
        cm0 = dict(dev["config"])
        Q0 = np.asarray(cm0["Q"], float).reshape(-1)
        obs = dev["obs"] if isinstance(dev["obs"], list) else [dev["obs"]]
        dicts = {}
        for noise in NOISES:
            o = next(x for x in obs if x["name"] == noise)
            D = np.asarray(o["D"], float)
            sig = bm.vec(o["sigma"])
            live = D.any(axis=0)
            K, _, adj = bm.predictor(D, live)
            for e in ERRORS:
                rf = ROOT / "cache/runner/results" / f"{name}_Q{e:+d}_{noise}.mat"
                res = loadmat(rf, simplify_cells=True)
                used = Q0 * (1 + e / 100)
                factor = float(res["scan"]["factor"])
                scanned = np.asarray(res["scanned"]["Q"], float).reshape(-1)
                assert np.allclose(np.asarray(res["used"]["Q"], float).reshape(-1), used, rtol=1e-12)
                assert np.allclose(scanned, used * factor, rtol=1e-12)
                for label, q in (("SC", used), ("QSC", scanned)):
                    f = out / f"{name}_{noise}_Q{e:+d}_{label}"
                    if f.with_suffix(".npz").exists() and f.with_suffix(".json").exists():
                        continue
                    key = tuple(np.round(q, 12))
                    if key not in dicts:
                        cm = dict(cm0)
                        cm["Q"] = q
                        dicts[key] = bm.dictionary(cm)
                    B, W, P, scale = dicts[key]
                    X, info = bm.solve(B, P, D, sig, scale, K, SEL["a"], SEL["mu_relative"], live)
                    info.update(model=name, noise=noise, q_error=e, label=label, qscan_factor=factor,
                                q_used=[float(x) for x in q], predictor_adjoint_error=adj)
                    np.savez_compressed(f.with_suffix(".npz"), out=X)
                    bm.writejson(f.with_suffix(".json"), info)
                    print(f"{name} {noise} Q{e:+d} {label} factor {factor:.4f}: iter {info['iterations']} "
                          f"KKT {info['kkt_ratio']:.4f} conv {info['converged']} {info['seconds']:.0f}s "
                          f"(total {time.monotonic() - started:.0f}s)", flush=True)
    print("DONE", " ".join(models), flush=True)


if __name__ == "__main__":
    main(sys.argv[1:])
