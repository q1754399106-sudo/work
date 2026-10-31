"""Shared code of experiment 17 (PLAN_ZH.md): choose SC-a's coherence weight mu_relative on development data
(the used first and second fault / complex-Q batches), then confirm the choice on the sealed third batch.

Imported unchanged: lib16 (experiment 16; through it lib15 and the experiment-7 to -13 code).
New code: data loading of the third batch, the mu grid, the method sets.
"""
import sys, json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]; BASE = ROOT.parent
sys.path.insert(0, str(BASE / "mu_confirm_20261031" / "scripts"))
import lib16  # noqa: E402
from lib16 import np, lib, semi, sq, sc11, lib12, TW, SEL, A_TUNED, G_IQ, A_ORIG, MU, score  # noqa: E402,F401

DATA3 = BASE / "semisynth_fault_q3_20261031"
META3 = json.loads((DATA3 / "analysis" / "datasets.json").read_text(encoding="utf-8"))
GRID = (0.0, 0.25, 0.5, 1.0, 2.0, 4.0, 10.0)
DEV = [f"Z{k}_c{c}_r{r}" for k in "ABCDEF" for c in (0.75, 1.05) for r in (0, 1)]      # 24, all used before
TEST = [f"Z{k}_c{c}_r{r}" for k in "GHI" for c in (0.75, 1.05) for r in (0, 1)]       # 12, sealed


def load(did):
    if did not in META3:
        return lib16.load(did)
    return {k: v.astype(float) for k, v in np.load(DATA3 / "cache" / "data" / f"{did}.npz").items()}, META3[did]


def estimates(D, cm, Q, T, wins, seed):
    sr = sq.estimate(D, cm, Q, table=T, wins=wins, nboot=0, seed=seed)
    p = sc11.prepare_sc(D, cm); c_qs, qi = sc11.qscan(D, cm, Q, p)
    return float(sr["c"]), float(c_qs), qi, p


def sc(D, cm, Q, c, p, a, mu):
    return lib.solve(cm, Q * c, p["wv"], D, p["sigma"], p["K"], p["live"], a, mu)[0]
