"""Imports of the experiment-13 input snapshot (unchanged code of experiments 7-13), as experiment 14's env14.py."""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
E13 = ROOT.parent / "field_baselines_volve_20261029"; IN = E13 / "cache" / "inputs"
sys.path.insert(0, str(ROOT.parent / "novelty_evidence_20261025" / "cache" / "deps"))   # scipy
for p in ("novelty_evidence_20261025/scripts", "field_application_20261028/scripts", "semisynthetic_volve_20261028/scripts",
          "specslope_q_volve_20261028/scripts"):
    sys.path.insert(0, str(IN / p))
import semi  # noqa: E402
from semi import np, lib  # noqa: E402,F401
import gen  # noqa: E402,F401

semi.ROOT, semi.IN = E13, IN
