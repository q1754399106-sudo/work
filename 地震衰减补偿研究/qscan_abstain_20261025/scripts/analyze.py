"""Judgments A1-A2 of PLAN_ZH.md and the descriptive statistics, from analysis/gather_scores.csv."""
from pathlib import Path
import hashlib, json
import numpy as np, pandas as pd

ROOT = Path(__file__).resolve().parents[1]
A = ROOT / "analysis"
seal = json.loads((ROOT / "decisions" / "SEALED_SHA256.json").read_text())
changed = [p for p, h in seal["hashes"].items() if p != "scripts/analyze.py" and hashlib.sha256((ROOT / p).read_bytes()).hexdigest() != h]
assert not changed, changed
g = pd.read_csv(A / "gather_scores.csv")
carb = {"K2", "K4", "K6", "K8"}
s = g.groupby(["model", "scene"]).agg(mse_D=("mse_D", "mean"), mse_rule=("mse_rule", "mean"), mse_S2=("mse_S2", "mean"),
                                      mse_B=("mse_B", "mean"), Z1=("Z1", "mean"), shifted=("g", lambda x: 0),
                                      n_shift=("b", "size")).reset_index()
s["n_shift"] = g.assign(sh=g.g != g.b).groupby(["model", "scene"]).sh.sum().to_numpy()
s = s.drop(columns="shifted")
s["rel"] = s.mse_rule / s.mse_D - 1
s["rel_forced_S2"] = s.mse_S2 / s.mse_D - 1
s["rel_rule_vs_B"] = s.mse_rule / s.mse_B - 1
s["carbonate"] = s.model.isin(carb)
s.to_csv(A / "scene_scores.csv", index=False)
prim = s[s.scene.isin(["Q100_14dB", "Q070_14dB"])]
A1 = bool(prim.rel.mean() <= -0.05)
A2 = bool(prim.rel.max() <= 0.05)
sec = {sc: dict(mean_rel=float(d.rel.mean()), worst_rel=float(d.rel.max()), n_shifted_gathers=int(d.n_shift.sum()))
       for sc, d in s[~s.scene.isin(["Q100_14dB", "Q070_14dB"])].groupby("scene")}
split = s.assign(hi=s.Z1 > 72.65656565656568)
out = dict(A1_mean_rel_le_minus_0_05=A1, A2_worst_rel_le_0_05=A2, confirmed=A1 and A2,
           primary=dict(mean_rel=float(prim.rel.mean()), worst_rel=float(prim.rel.max()), n_improved=int((prim.rel < 0).sum()),
                        n=len(prim), n_shifted_gathers=int(prim.n_shift.sum()), mean_rel_forced_S2=float(prim.rel_forced_S2.mean()),
                        worst_rel_forced_S2=float(prim.rel_forced_S2.max())),
           secondary=sec,
           cluster_vs_carbonate_primary=pd.crosstab(prim.Z1 > 72.65656565656568, prim.carbonate).to_dict(),
           Q0_14dB_rule_vs_nominal_B=float(s[s.scene == "Q100_14dB"].rel_rule_vs_B.mean()))
(A / "decisions.json").write_text(json.dumps(out, indent=2, default=str), encoding="utf-8")
print(s.round(4).to_string(index=False))
print(json.dumps(out, indent=2, default=str))
