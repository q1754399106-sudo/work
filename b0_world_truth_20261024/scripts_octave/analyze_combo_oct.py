"""Descriptive (RESUME_ZH.md, not a gate): do the sealed B0 calls, made for FullCalRefit c=1 versus
c=QScan, also predict the combination-level sign, i.e. whether QScan + Standard is worse than
Standard at the public Q? Uses analysis/combo_truth_U*.csv (12 replicates per scene)."""
from pathlib import Path
import json
import numpy as np, pandas as pd

ROOT = Path(__file__).resolve().parents[1]
A = ROOT / "analysis"
c = pd.concat([pd.read_csv(f) for f in sorted(A.glob("combo_truth_U*.csv"))])
t = pd.concat([pd.read_csv(f) for f in sorted(A.glob("truth_U*.csv"))])
k = ["model", "condition", "replicate"]
chk = c.merge(t, on=k, suffixes=("", "_t"))
consistency = float(max((chk.mse_base - chk.mse_base_t).abs().max(), (chk.mse_qscan - chk.mse_qscan_t).abs().max()))
c["dc"] = c.mse_qscan_standard - c.mse_standard
g = c.groupby(["model", "condition"]).agg(mse_standard=("mse_standard", "mean"), mse_qscan_standard=("mse_qscan_standard", "mean"),
                                          dc=("dc", "mean"), worse_fraction=("dc", lambda s: float((s > 0).mean())),
                                          mse_fcr=("mse_base", "mean"), mse_fcr_qscan=("mse_qscan", "mean")).reset_index()
g["rel_combo"] = g.dc / g.mse_standard
g["rel_fcr"] = (g.mse_fcr_qscan - g.mse_fcr) / g.mse_fcr
g["decisive_combo"] = g.rel_combo.abs() >= 0.05
g["same_sign"] = np.sign(g.rel_combo) == np.sign(g.rel_fcr)
g.to_csv(A / "combo_expected_truth.csv", index=False)
s = pd.read_csv(A / "confirm_scored.csv").merge(g[["model", "condition", "rel_combo", "decisive_combo"]], on=["model", "condition"])
s["combo_worse"] = s.rel_combo > 0
calls = {c_: s[c_] for c_ in s.columns if c_.startswith("call_")}
dec = s.decisive_combo
acc = {name.replace("call_", ""): dict(decisive=float((v[dec] == s.combo_worse[dec]).mean()) if dec.any() else None,
                                       all=float((v == s.combo_worse).mean())) for name, v in calls.items()}
out = dict(consistency_with_truth_csv=consistency, scenes=len(g), decisive_scenes_combo=int(g.decisive_combo.sum()),
           decisive_gathers_combo=int(dec.sum()), scenes_combo_worse=int((g.rel_combo > 0).sum()),
           scenes_same_sign_fcr_vs_combo=int(g.same_sign.sum()), accuracy_vs_combo_sign=acc,
           spearman_B0M_vs_combo=float(s.groupby(["model", "condition"]).agg(b=("B0M_rel_delta_mse", "mean"), r=("rel_combo", "first"))
                                       .rank().corr().iloc[0, 1]))
(A / "combo_decisions.json").write_text(json.dumps(out, indent=2), encoding="utf-8")
print(g[["model", "condition", "mse_standard", "mse_qscan_standard", "rel_combo", "worse_fraction", "rel_fcr", "same_sign"]].round(4).to_string(index=False))
print(json.dumps(out, indent=2))
