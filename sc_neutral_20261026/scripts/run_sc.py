"""Solve the SC baseline on the rebuilt H1-H6 gathers (PLAN_ZH.md, 计算).

The solver, spatial predictor and dictionary are imported unchanged from the snapshot of
novelty_evidence_20261025/scripts (benchmark.solve, benchmark.predictor, common.dictionary), with
the parameters frozen in its analysis/selection.json: SC a=0.125, mu_relative=10, and the
no-structure ablation SC0 a=0.5, mu=0. Same noisy gather, sigma and Q as the fd2d Standard run.
Truth is not read. Outputs: cache/solved/<model>_<scene>_r<r>_<label>.npz and .json (solver info).

usage: run_sc.py H1 H2 ...   (models; tasks run in the priority order of PLAN_ZH.md, resumable)
"""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SNAP = ROOT / "cache" / "inputs" / "novelty_evidence_20261025"
sys.path.insert(0, str(ROOT.parent / "novelty_evidence_20261025" / "cache" / "deps"))
sys.path.insert(0, str(SNAP / "scripts"))
import benchmark as bm  # noqa: E402  (common.py pins BLAS to one thread before NumPy loads)
import json, time  # noqa: E402,E401
import numpy as np  # noqa: E402
from scipy.io import loadmat  # noqa: E402

SEL = json.loads((SNAP / "analysis" / "selection.json").read_text(encoding="utf-8"))
assert (SEL["a"], SEL["mu_relative"], SEL["no_structure_a"]) == (0.125, 10.0, 0.5)
CONFIGS = {"SC": (SEL["a"], SEL["mu_relative"]), "SC0": (SEL["no_structure_a"], 0.0)}
SCENES = {"Q100_14dB": (1, 1.0), "Q100_08dB": (2, 1.0), "Q070_14dB": (1, 0.7)}


def tasks(models):
    # Primary (Q0, r=1..6, SC) first, then the descriptive runs.
    order = [(s, rs, lab) for s in ("Q100_14dB", "Q100_08dB") for rs in [range(1, 7)] for lab in ["SC"]]
    order += [(s, range(1, 7), "SC0") for s in ("Q100_14dB", "Q100_08dB")]
    order += [("Q070_14dB", range(1, 7), lab) for lab in ("SC", "SC0")]
    order += [(s, [0], lab) for s in SCENES for lab in ("SC", "SC0")]
    for s, rs, lab in order:
        for m in models:
            for r in rs:
                yield m, s, r, lab


def main(models):
    cache = {}
    out = ROOT / "cache" / "solved"
    out.mkdir(parents=True, exist_ok=True)
    started = time.monotonic()
    for m, s, r, lab in tasks(models):
        f = out / f"{m}_{s}_r{r}_{lab}"
        if f.with_suffix(".npz").exists() and f.with_suffix(".json").exists():
            continue
        level, qs = SCENES[s]
        if (m, qs) not in cache:
            cm = dict(loadmat(ROOT / "cache/inputs/fd2d_20261025/models" / f"{m}.mat", simplify_cells=True)["config"])
            cm["Q"] = np.asarray(cm["Q"], float) * qs
            cache[(m, qs)] = bm.dictionary(cm)
        B, W, P, scale = cache[(m, qs)]
        g = loadmat(ROOT / "cache/noisy" / f"{m}_L{level}_r{r}.mat")
        D = np.asarray(g["D"], float)
        sig = bm.vec(g["sig"])
        live = D.any(axis=0)
        K, _, adj = bm.predictor(D, live)
        a, mu = CONFIGS[lab]
        X, info = bm.solve(B, P, D, sig, scale, K, a, mu, live)
        info.update(model=m, scene=s, r=r, label=lab, predictor_adjoint_error=adj)
        np.savez_compressed(f.with_suffix(".npz"), out=X)
        bm.writejson(f.with_suffix(".json"), info)
        print(f"{m} {s} r{r} {lab}: iter {info['iterations']} KKT {info['kkt_ratio']:.4f} "
              f"conv {info['converged']} {info['seconds']:.0f}s (total {time.monotonic() - started:.0f}s)", flush=True)
    print("DONE", " ".join(models), flush=True)


if __name__ == "__main__":
    main(sys.argv[1:])
