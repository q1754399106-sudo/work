"""Experiment 14, descriptive supplements -> analysis/sensitivity.json.

1. Convergence sensitivity (PLAN_ZH.md section 5): per weight, the share of non-converged FISTA solves, and the group
   statistics recomputed on the datasets whose solves at the selected factor and its two grid neighbours converged
   (non-convergence sits almost only at the two smallest factors, 0.40 and 0.46, so requiring all 15 would drop most
   datasets); H3/H4 difference on the datasets where QScan-W(2^-4) converged in that sense.
2. Transfer of the calibration: bias on the development set vs bias on each test group, per weight and for SR.
3. Error by true c in group T for QScan, SR, SR_cal and QScan-W at every weight (raw and calibrated).
"""
import json
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
S = json.loads((ROOT / "analysis" / "summary.json").read_text(encoding="utf-8")); cal = S["calibration"]
test = [json.loads(p.read_text(encoding="utf-8")) for p in sorted((ROOT / "results" / "test").glob("*.json"))]
dev = [json.loads(p.read_text(encoding="utf-8")) for p in sorted((ROOT / "results" / "dev").glob("*.json"))]
lams = [str(l) for l in test[0]["lams"]]; L4 = str(2 ** -4)


def describe(sel, ct):
    e = np.log(sel / ct)
    return dict(n=len(ct), slope=float(np.polyfit(np.log(ct), np.log(sel), 1)[0]) if len(set(ct)) > 1 else None,
                bias_log=float(e.mean()), abs_log=float(np.abs(e).mean()))


def allconv(r, i):
    j = int(np.argmin(r["diag"]["J"][i]))
    return bool(np.all(r["diag"]["converged"][i][max(j - 1, 0):j + 2]))


out = dict(nonconverged_share={}, nonconverged_by_grid={}, sensitivity={}, transfer={}, by_c={})
out["nonconverged_by_grid"] = {l: np.sum([np.logical_not(r["diag"]["converged"][i]) for r in dev + test], axis=0).tolist()
                               for i, l in enumerate(lams)}
for name, rs in (("dev", dev), ("T", [r for r in test if r["group"] == "T"]), ("S", [r for r in test if r["group"] != "T"])):
    out["nonconverged_share"][name] = {l: float(1 - np.mean([np.mean(r["diag"]["converged"][i]) for r in rs]))
                                       for i, l in enumerate(lams)}
    out["nonconverged_share"][name]["all_weights"] = float(np.mean(list(out["nonconverged_share"][name].values())))
    out["nonconverged_share"][name]["max_iterations"] = int(max(np.max(r["diag"]["iterations"]) for r in rs))
for grp in ("T", "S"):
    rs = [r for r in test if (r["group"] == "T") == (grp == "T")]; sens = {}
    for i, l in enumerate(lams):
        keep = [r for r in rs if allconv(r, i)]; ct = np.array([r["c_true"] for r in keep])
        sel = np.array([r["qscan_w"][l] for r in keep])
        sens[l] = dict(raw=describe(sel, ct), cal=describe(sel * np.exp(-cal[l]), ct), ids_dropped=[r["id"] for r in rs if not allconv(r, i)])
    i4 = lams.index(L4); keep = [r for r in rs if allconv(r, i4)]; ct = np.array([r["c_true"] for r in keep])
    w = np.array([r["qscan_w"][L4] for r in keep]) * np.exp(-cal[L4]); sr = np.array([r["c_sr"] for r in keep])
    d = np.abs(np.log(w / ct)) - np.abs(np.log(sr * np.exp(-cal["SR"]) / ct))
    sens["H3_diff_converged"] = dict(n=len(keep), mean=float(d.mean()))
    out["sensitivity"][grp] = sens
    out["transfer"][grp] = {l: dict(dev=cal["dev_describe"][l]["bias_log"], test=S["groups"][grp]["describe"][f"W{float(l)}"]["bias_log"])
                            for l in lams}
    out["transfer"][grp]["SR"] = dict(dev=cal["SR"], test=S["groups"][grp]["describe"]["SR"]["bias_log"])
rs = [r for r in test if r["group"] == "T"]; est = S["groups"]["T"]["estimates"]; ct = np.array([r["c_true"] for r in rs])
for c in sorted(set(ct)):
    m = ct == c
    out["by_c"][str(c)] = {k: dict(mean_sel=float(np.exp(np.mean(np.log(np.array(v)[m])))),
                                   abs_log=float(np.mean(np.abs(np.log(np.array(v)[m] / c))))) for k, v in est.items()}
(ROOT / "analysis" / "sensitivity.json").write_text(json.dumps(out, indent=1, ensure_ascii=False), encoding="utf-8")
print("non-converged share:", json.dumps(out["nonconverged_share"], indent=None))
for grp in ("T", "S"):
    print(grp, "transfer (dev bias -> test bias):", {k: (round(v["dev"], 3), round(v["test"], 3)) for k, v in out["transfer"][grp].items()})
    for l in lams:
        s = out["sensitivity"][grp][l]
        print(f"  {grp} lam={l:<9} conv-only n={s['raw']['n']:2d} raw slope={s['raw']['slope'] or float('nan'):.2f} bias={s['raw']['bias_log']:+.3f} "
              f"|e|={s['raw']['abs_log']:.3f}  cal |e|={s['cal']['abs_log']:.3f}  dropped={s['ids_dropped']}")
    print("  H3 diff on converged:", out["sensitivity"][grp]["H3_diff_converged"])
for c, v in out["by_c"].items():
    print(f"c={c}: " + "  ".join(f"{k}={v[k]['mean_sel']:.2f}" for k in ("QScan", "SR", "W0.5", "W0.25", "W0.0625", "W0.015625")))
