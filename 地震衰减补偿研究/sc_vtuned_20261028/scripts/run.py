"""SC with V-family tuning (PLAN_ZH.md). Same solver as novelty_evidence_20261025/scripts/benchmark.py (unchanged).

  run.py dev [--extended]   V1-V3 x {Moderate, Strong} x Q{-30,-10,0,+10,+30}: every grid candidate (reads truth;
                            development data). Candidates: GRID, plus analysis/extension.json with --extended.
  run.py eval               W1-W6, X1-X3 x noise x Q{-30,0,+30}: SC and QScan + SC with analysis/selection.json
                            (and the secondary selection if it differs). Does not read truth.
Resumable; one result per solve in cache/. Two worker processes (2 cores), BLAS one thread each.
"""
import os
for _k in ("OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "OMP_NUM_THREADS"):
    os.environ[_k] = "1"
import sys, json, csv, hashlib, time  # noqa: E401
from pathlib import Path
from multiprocessing import Pool

ROOT = Path(__file__).resolve().parents[1]
RESEARCH = ROOT.parent
sys.path.insert(0, str(RESEARCH / "novelty_evidence_20261025" / "scripts"))
import benchmark as bm  # noqa: E402
import numpy as np  # noqa: E402
from scipy.io import loadmat  # noqa: E402

CODE_SHA = {"benchmark.py": "a2f593c07ffd32547c15d9629e4232c7e6ab1de8fb388c13d31b04d825e3de75",
            "common.py": "17a6f89faf0eddfc7469328ab80a5189e4ed3ecd4b3ad1ec58c70cd81f1ef5ab"}
for _f, _h in CODE_SHA.items():
    assert hashlib.sha256((RESEARCH / "novelty_evidence_20261025/scripts" / _f).read_bytes()).hexdigest() == _h, _f
assert (bm.MAXITER, bm.TOL) == (5000, .02)

A_GRID = [1 / 32, 1 / 16, 1 / 8, 1 / 4]
M_GRID = [1., 10., 100., 1000.]
GRID = [(a, m) for a in A_GRID for m in M_GRID]
NOISES = ["Moderate", "Strong"]
DEV_SCENES, DEV_ERRORS = ["V1", "V2", "V3"], [-30, -10, 0, 10, 30]
EVAL_SCENES, EVAL_ERRORS = [f"W{k}" for k in range(1, 7)] + [f"X{k}" for k in range(1, 4)], [-30, 0, 30]
VF = RESEARCH / "vfamily_20261027"
CS = RESEARCH / "coherence_soft_20261023"


def tag(a, m):
    return f"a{a:.6g}_m{m:g}"


def load(scene):
    if scene.startswith("V"):
        return loadmat(CS / "inputs/validation" / f"{scene}.mat", simplify_cells=True)
    p = VF / "models" / f"{scene}.mat"
    hashes = json.loads((VF / "logs/MODEL_HASHES.json").read_text(encoding="utf-8"))
    assert hashlib.sha256(p.read_bytes()).hexdigest() == hashes[f"models/{scene}.mat"], scene
    return loadmat(p, simplify_cells=True)


def factors():
    """QScan factor per (scene, noise, e) archived by vfamily run_combo3 (rows.csv)."""
    out = {}
    with open(VF / "analysis/rows.csv", encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            out[(r["scene"], r["noise"], int(r["q_error"]))] = float(r["qscan_factor"])
    return out


def task(job):
    """job = (mode, scene, noise, e, [(label, qmult, a, m), ...]); one dictionary per distinct Q."""
    mode, scene, noise, e, items = job
    d = load(scene)
    cm0 = dict(d["config"])
    Q0 = np.asarray(cm0["Q"], float).reshape(-1)
    obs = d["obs"] if isinstance(d["obs"], list) else [d["obs"]]
    o = next(x for x in obs if x["name"] == noise)
    D = np.asarray(o["D"], float)
    sig = bm.vec(o["sigma"])
    live = D.any(axis=0)
    K, _, adj = bm.predictor(D, live)
    T = None
    if mode == "dev":
        tr = loadmat(CS / "truth" / f"{scene}_truth.mat", simplify_cells=True)["tr"]
        T = np.asarray(next(t for t in tr if t["name"] == noise)["target"], float)
    dicts, log = {}, []
    for label, qmult, a, m in items:
        f = ROOT / "cache" / mode / f"{scene}_{noise}_Q{e:+d}_{label}_{tag(a, m)}"
        if f.with_suffix(".json").exists():
            continue
        q = Q0 * (1 + e / 100) * qmult
        key = tuple(np.round(q, 12))
        if key not in dicts:
            cm = dict(cm0)
            cm["Q"] = q
            dicts = {key: bm.dictionary(cm)}  # keep one dictionary in memory (V3/X: 851 x 851)
        B, W, P, scale = dicts[key]
        X, info = bm.solve(B, P, D, sig, scale, K, a, m, live)
        info.update(scene=scene, noise=noise, q_error=e, label=label, qscan_factor=qmult,
                    q_used=[float(x) for x in q], predictor_adjoint_error=adj)
        if T is not None:
            info["rel_error"] = float(np.linalg.norm(X - T) / max(np.linalg.norm(T), 1e-30))
        f.parent.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(f.with_suffix(".npz"), out=X)
        bm.writejson(f.with_suffix(".json"), info)
        log.append(f"{scene} {noise} Q{e:+d} {label} {tag(a, m)}: it {info['iterations']} KKT {info['kkt_ratio']:.4f} "
                   f"conv {info['converged']} {info['seconds']:.0f}s"
                   + (f" err {info['rel_error']:.5f}" if T is not None else ""))
    return log


def jobs(mode, extended):
    out = []
    if mode == "dev":
        cands = list(GRID)
        if extended:
            cands += [tuple(x) for x in json.loads((ROOT / "analysis/extension.json").read_text())["added"]]
        for s in DEV_SCENES:
            for e in DEV_ERRORS:
                for n in NOISES:
                    out.append(("dev", s, n, e, [("SC", 1.0, a, m) for a, m in cands]))
    else:
        sel = json.loads((ROOT / "analysis/selection.json").read_text(encoding="utf-8"))
        pars = [(sel["primary"]["a"], sel["primary"]["mu_relative"])]
        sec = (sel["secondary"]["a"], sel["secondary"]["mu_relative"])
        if sec != pars[0]:
            pars.append(sec)
        fac = factors()
        for s in EVAL_SCENES:
            for e in EVAL_ERRORS:
                for n in NOISES:
                    items = [("SC", 1.0, a, m) for a, m in pars] + [("QSC", fac[(s, n, e)], a, m) for a, m in pars]
                    out.append(("eval", s, n, e, items))
    # heavy (851-sample) scenes first so the two workers finish together
    return sorted(out, key=lambda j: j[1] not in ("V3", "X1", "X2", "X3"))


if __name__ == "__main__":
    mode = sys.argv[1]
    assert mode in ("dev", "eval")
    if mode == "eval":
        sel = json.loads((ROOT / "analysis/selection.json").read_text(encoding="utf-8"))
        assert sel["frozen"] is True
    started = time.monotonic()
    with Pool(2) as pool:
        for log in pool.imap_unordered(task, jobs(mode, "--extended" in sys.argv)):
            for line in log:
                print(f"[{time.monotonic() - started:6.0f}s] {line}", flush=True)
    print("DONE", mode, flush=True)
