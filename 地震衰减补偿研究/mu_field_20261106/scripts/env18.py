"""Imports of experiment 18: the frozen experiment-13 code, read from Git at a fixed commit and run with __file__ set
to its original location, so that experiment 13's cached inputs (cache/inputs) and outputs (cache/out) resolve.
Nothing in experiment 13's folders is written."""
import json
import subprocess
import sys
import types
from pathlib import Path

EXP = Path(__file__).resolve().parents[1]
RES = EXP.parent
REPO = RES.parent
COMMIT = "27a065a"
E13 = "field_baselines_volve_20261029"
CACHE13 = RES / E13 / "cache" / "out"
CACHE = EXP / "cache" / "out"
OUT = EXP / "results"
MU_NEW = 0.25


def git_text(rel):
    return subprocess.check_output(["git", "show", f"{COMMIT}:{RES.name}/{rel}"], cwd=REPO).decode("utf-8")


def git_module(name, fname):
    if name in sys.modules:
        return sys.modules[name]
    m = types.ModuleType(name)
    m.__file__ = str(RES / E13 / "scripts" / fname)
    sys.modules[name] = m
    exec(compile(git_text(f"{E13}/scripts/{fname}"), m.__file__, "exec"), m.__dict__)
    return m


L = git_module("lib13", "lib13.py")
run13 = git_module("run13", "run13.py")
np, semi, sc11, lib, tie = L.np, L.semi, L.sc11, L.lib, L.tie
assert lib.SEL_SC[1] == 10, "the frozen coherence weight of experiment 13 should be 10"


def lines():
    return json.loads(git_text(f"{E13}/analysis/lines.json"))["lines"]


def qscan_factors():
    return {ln: json.loads(git_text(f"{E13}/results/qscan/IL_{ln}.json"))["c"] for ln in L.QSCAN_LINES}


def solve(D, cm, Q, c, p, mu):
    """SC-a with a = 2^-13 and coherence weight mu: the call of lib13.methods_a/b with mu in place of 10."""
    return lib.solve(cm, Q * c, p["wv"], D, p["sigma"], p["K"], p["live"], L.A_TUNED, mu)[0]
