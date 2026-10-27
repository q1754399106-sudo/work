"""SC and QScan + SC on the field curtains (PLAN_ZH.md "方法").

The solver, spatial predictor and dictionary are imported unchanged from the snapshot of
novelty_evidence_20261025/scripts (benchmark.solve, benchmark.predictor, common.dictionary of ff95a69),
with SC a=0.125, mu_relative=10 from its analysis/selection.json. For each run of a well the curtain D,
the per-trace sigma and the dictionary model cm come from that well's curtain.mat; SC uses the run's Q,
QSC uses the run's Q times the QScan factor that combo_outputs selected in the same run.
  F11:  inputs/curtain.mat, results/<run>.mat                    (runs V0 Vm50 Vm30 Vp30 Vp50)
  F15A: the archived experiment-5 curtain.mat and results/<run>.mat (runs V0 Vm50 Vm30 Vp30 Vp50 R0)
Outputs: results_sc/<well>_<run>_<SC|QSC>.mat (X) and .json (solver info).

usage: run_sc.py F11|F15A   (resumable)
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
from scipy.io import loadmat, savemat  # noqa: E402

SEL = json.loads((SNAP / "analysis" / "selection.json").read_text(encoding="utf-8"))
assert (SEL["a"], SEL["mu_relative"]) == (0.125, 10.0)
F15A = ROOT / "cache" / "inputs" / "field_vsp_volve_20261025"
WELLS = {"F11": (ROOT / "inputs" / "curtain.mat", ROOT / "results", ["V0", "Vm50", "Vm30", "Vp30", "Vp50"]),
         "F15A": (F15A / "inputs" / "curtain.mat", F15A / "results", ["V0", "Vm50", "Vm30", "Vp30", "Vp50", "R0"])}


def main(well):
    curtain, results, runs = WELLS[well]
    c = loadmat(curtain, simplify_cells=True)
    D = np.asarray(c["D"], float)
    sig = bm.vec(c["sigma"])
    live = D.any(axis=0)
    K, _, adj = bm.predictor(D, live)
    out = ROOT / "results_sc"
    out.mkdir(exist_ok=True)
    started = time.monotonic()
    for run in runs:
        r = loadmat(results / f"{run}.mat", simplify_cells=True)
        Q = np.asarray(r["Q"], float).ravel()
        factor = float(r["info"]["qscan_factor"])
        for label, q in (("SC", Q), ("QSC", Q * factor)):
            f = out / f"{well}_{run}_{label}"
            if f.with_suffix(".mat").exists() and f.with_suffix(".json").exists():
                continue
            cm = dict(c["cm"])
            cm["Q"] = q
            cm["solver_offset"] = 0
            B, W, P, scale = bm.dictionary(cm)
            X, info = bm.solve(B, P, D, sig, scale, K, SEL["a"], SEL["mu_relative"], live)
            info.update(well=well, run=run, label=label, Q=q.tolist(), qscan_factor=factor,
                        predictor_adjoint_error=adj)
            savemat(f.with_suffix(".mat"), {"X": X})
            bm.writejson(f.with_suffix(".json"), info)
            print(f"{well} {run} {label} Q={np.round(q, 1).tolist()}: iter {info['iterations']} KKT {info['kkt_ratio']:.4f} "
                  f"conv {info['converged']} {info['seconds']:.0f}s (total {time.monotonic() - started:.0f}s)", flush=True)
    print("DONE", well, flush=True)


if __name__ == "__main__":
    main(sys.argv[1])
