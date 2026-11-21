"""Frozen rules of PLAN_ZH.md section 3 on results/tradeoff.json and results/s2.json -> results/verdict.json (+ printed tables)."""
import json
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]; R = json.loads((ROOT / "results" / "tradeoff.json").read_text())["rows"]
S2 = json.loads((ROOT / "results" / "s2.json").read_text())
q = lambda x, p: float(np.quantile(x, p))  # noqa: E731
out = {"S1": {}, "desc": {}}
for src in "ABCL":
    rr = [r for r in R if r["src"] == src]
    dg = np.abs([r["gt_mc"] - r["gt_cf"] for r in rr]); dr = np.abs([r["ratio_mc"] - r["ratio_cf"] for r in rr])
    dgH = np.abs([r["gt_mc"] - r["gt_cf_H"] for r in rr])
    ok = bool(np.median(dg) <= 0.01 and np.median(dr) <= 0.05)
    out["S1"][src] = dict(n=len(rr), med_dgt=float(np.median(dg)), p90_dgt=q(dg, 0.9), med_dratio=float(np.median(dr)), p90_dratio=q(dr, 0.9),
                          med_dgt_H_posthoc=float(np.median(dgH)), pass_=ok)
    for sig in (0.447, 0.09):
        s = [r for r in rr if r["sig"] == sig]
        out["desc"][f"{src}_{sig}"] = {k: float(np.median([r[k] for r in s])) for k in
                                       ("sigma_beta", "beta_mean", "beta_sd", "beta_loglin_rms", "lam", "k", "gt_cf", "gt_mc", "gt_mc_se", "ratio_cf", "ratio_mc",
                                        "red_cf_in_mc", "red_mc", "rmse_sat", "rmse_fs", "gtil_at_gamma1", "a0", "h", "f0")}
        out["desc"][f"{src}_{sig}"].update(red_cf_range=[float(min(r["red_cf_in_mc"] for r in s)), float(max(r["red_cf_in_mc"] for r in s))],
                                           hf0_sens_over_first=float(np.median([r["hf0_sens"] / r["hf0_first"] for r in s])))
s1 = all(out["S1"][s]["pass_"] for s in "ABC")
A = [r for r in R if r["src"] == "A" and r["sig"] == 0.447]; redA = float(np.median([r["red_cf_in_mc"] for r in A]))
out["S3"] = dict(median_reduction=redA, range=[float(min(r["red_cf_in_mc"] for r in A)), float(max(r["red_cf_in_mc"] for r in A))], pass_=bool(redA >= 0.30))
out["S2"] = dict(phi_median=S2["phi_median"], ratio=S2["ratio"], ratio_phi0=S2["ratio_phi0"], S2a=S2["S2a"], S2b=S2["S2b"], pass_=S2["S2"])
out["S1_pass"] = s1
out["decision"] = "立项" if (s1 and S2["S2"] and out["S3"]["pass_"]) else "不立项"
(ROOT / "results" / "verdict.json").write_text(json.dumps(out, indent=1, ensure_ascii=False), encoding="utf-8")
print(json.dumps(out, indent=1, ensure_ascii=False))
