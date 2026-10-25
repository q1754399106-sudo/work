"""Pre-run input check: the Python dictionary of SC (common.dictionary, ff95a69) against the frozen
MATLAB prune_dictionary exported by dict_check.m, for H1-H6 at Q0 and Q0x0.7.
Writes analysis/dictionary_check.json. No method output is computed or read.
"""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT.parent / "novelty_evidence_20261025" / "cache" / "deps"))
sys.path.insert(0, str(ROOT / "cache" / "inputs" / "novelty_evidence_20261025" / "scripts"))
import common as cmn  # noqa: E402
import numpy as np  # noqa: E402
from scipy.io import loadmat  # noqa: E402

rows = []
for m in range(1, 7):
    cfg = dict(loadmat(ROOT / "cache/inputs/fd2d_20261025/models" / f"H{m}.mat", simplify_cells=True)["config"])
    for qs in (1.0, 0.7):
        c = dict(cfg)
        c["Q"] = np.asarray(cfg["Q"], float) * qs
        B, W, P, scale = cmn.dictionary(c)
        ref = loadmat(ROOT / "cache/dictcheck" / f"H{m}_q{round(100 * qs):03d}.mat")
        rel = lambda a, b: float(np.linalg.norm(a - b) / np.linalg.norm(b))
        rows.append(dict(model=f"H{m}", q_scale=qs, B=rel(B, ref["B"]), W=rel(W, ref["W"]),
                         scale=rel(scale, np.ravel(ref["scale"]))))
        print(rows[-1], flush=True)
worst = max(max(r["B"], r["W"], r["scale"]) for r in rows)
cmn.writejson(ROOT / "analysis" / "dictionary_check.json", dict(rows=rows, worst_relative_difference=worst,
                                                                pass_1e_8=worst < 1e-8))
print("worst", worst)
