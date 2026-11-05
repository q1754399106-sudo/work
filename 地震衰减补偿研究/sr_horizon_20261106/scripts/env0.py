"""Imports for stage 0. The source of earlier experiments is read from Git at a fixed commit, because their
script folders are cleaned locally, and run with __file__ set to its original location so that the cached
inputs of those experiments resolve as usual. Nothing is imported from the working tree of other experiments."""
import subprocess
import sys
import types
from pathlib import Path

EXP = Path(__file__).resolve().parents[1]
RES = EXP.parent
REPO = RES.parent
COMMIT = "2063dd1"
OUT = EXP / "results"
OUT.mkdir(exist_ok=True)


def git_module(name, exp, fname):
    if name in sys.modules:
        return sys.modules[name]
    rel = f"{RES.name}/{exp}/scripts/{fname}"
    src = subprocess.check_output(["git", "show", f"{COMMIT}:{rel}"], cwd=REPO).decode("utf-8")
    m = types.ModuleType(name)
    m.__file__ = str(RES / exp / "scripts" / fname)
    sys.modules[name] = m
    exec(compile(src, m.__file__, "exec"), m.__dict__)
    return m


envf = git_module("envf", "semisynth_fault_q_20261030", "envf.py")
np, lib, semi, gen = envf.np, envf.lib, envf.semi, envf.gen
fq = git_module("fq", "semisynth_fault_q_20261030", "fq.py")
import sq  # noqa: E402  (experiment 10, from the experiment-13 input snapshot on envf's path)

semi.ROOT, semi.IN = envf.E13, envf.IN
TW, DT = semi.TW, semi.DT


def _load(did, _meta={}):
    """As fq.load, with the metadata read from Git (the analysis folder is cleaned locally)."""
    import json
    if not _meta:
        rel = f"{RES.name}/semisynth_fault_q_20261030/analysis/datasets.json"
        _meta.update(json.loads(subprocess.check_output(["git", "show", f"{COMMIT}:{rel}"], cwd=REPO).decode("utf-8")))
    return dict(np.load(RES / "semisynth_fault_q_20261030" / "cache" / "data" / f"{did}.npz")), _meta[did]


fq.load = _load
