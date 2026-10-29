"""Shared setup: imports of the experiment-13 input snapshot (unchanged code of experiments 7-13) and data paths."""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
E13 = ROOT.parent / "field_baselines_volve_20261029"; IN = E13 / "cache" / "inputs"
E12OUT = ROOT.parent / "baselines_volve_20261029" / "cache" / "out"
sys.path.insert(0, str(ROOT.parent / "novelty_evidence_20261025" / "cache" / "deps"))   # scipy
for p in ("novelty_evidence_20261025/scripts", "field_application_20261028/scripts", "semisynthetic_volve_20261028/scripts",
          "specslope_q_volve_20261028/scripts", "multiwell_tie_volve_20261028/scripts", "srsc_volve_20261028/scripts"):
    sys.path.insert(0, str(IN / p))
import sc11  # noqa: E402
from sc11 import np, lib, semi, sq, bm  # noqa: E402,F401

for mod in (sc11, sq, semi):
    mod.ROOT, mod.IN = E13, IN
sc11.wq.ROOT, sc11.wq.IN, sc11.wq.E9 = E13, IN, IN / "multiwell_tie_volve_20261028"


def datasets():
    return sorted(q for q in E12OUT.glob("T_*.npz") if "__" not in q.name)
