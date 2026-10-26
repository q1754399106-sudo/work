"""SC and QScan + SC on F3/F4/V4/V5 x Moderate/Strong x 7 Q errors (PLAN_ZH.md "数据与方法").

SC uses the used Q = Q0 * (1 + e/100); QScan + SC uses the archived round3 choice scanned.Q
(= used Q * scan.factor, the factor QScan + Standard used). Both are solved exactly as in
novelty_evidence_20261025 (benchmark.solve/predictor and common.dictionary of ff95a69, a=0.125,
mu_relative=10). Every snapshot input is checked against cache/SNAPSHOT.json. Truth is not read.
Q+0 SC is solved first for each gather and must match the local novelty output (relative
difference <= 1e-8), otherwise the run stops.
Outputs: cache/solved/<set>_<noise>_Q<e>_<SC|QSC>.npz and .json.

usage: run_sc.py F3 F4 ...   (resumable)
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

S = json.loads((ROOT / "cache" / "SNAPSHOT.json").read_text(encoding="utf-8"))
NOISES = ["Moderate", "Strong"]
ERRORS = [0, -50, -30, -10, 10, 30, 50]          # Q+0 first: consistency check before anything else


def snap(rel):
    p = SNAP / rel
    assert hashlib.sha256(p.read_bytes()).hexdigest() == S["files"][rel]["sha256"], f"input changed: {rel}"
    return loadmat(p, simplify_cells=True)


def local(rel):
    p = RESEARCH / rel
    assert hashlib.sha256(p.read_bytes()).hexdigest() == S["local_sha256"][rel], f"input changed: {rel}"
    return p


SEL = json.loads((SNAP / "novelty_evidence_20261025/analysis/selection.json").read_text(encoding="utf-8"))
assert (SEL["a"], SEL["mu_relative"]) == (0.125, 10.0)


def main(sets):
    out = ROOT / "cache" / "solved"
    out.mkdir(parents=True, exist_ok=True)
    started = time.monotonic()
    for name in sets:
        dev = snap(f"round3_20261024/dev/{name}.mat")
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
                res = snap(f"round3_20261024/results/{name}_Q{e:+d}_{noise}.mat")
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
                    info.update(set=name, noise=noise, q_error=e, label=label, qscan_factor=factor,
                                q_used=[float(x) for x in q], predictor_adjoint_error=adj)
                    if label == "SC" and e == 0:
                        ref = np.load(local(f"novelty_evidence_20261025/cache/solved/{name}_{noise}_0.125_10.npz"))["out"]
                        diff = float(np.linalg.norm(X - ref) / max(np.linalg.norm(ref), 1e-30))
                        info["novelty_relative_difference"] = diff
                        if diff > 1e-8:
                            bm.writejson(f.with_suffix(".failed.json"), info)
                            raise SystemExit(f"consistency check failed: {name} {noise} Q+0 SC differs by {diff:.3e}")
                    np.savez_compressed(f.with_suffix(".npz"), out=X)
                    bm.writejson(f.with_suffix(".json"), info)
                    print(f"{name} {noise} Q{e:+d} {label} factor {factor:.4f}: iter {info['iterations']} "
                          f"KKT {info['kkt_ratio']:.4f} conv {info['converged']} {info['seconds']:.0f}s "
                          f"(total {time.monotonic() - started:.0f}s)", flush=True)
    print("DONE", " ".join(sets), flush=True)


if __name__ == "__main__":
    main(sys.argv[1:])
