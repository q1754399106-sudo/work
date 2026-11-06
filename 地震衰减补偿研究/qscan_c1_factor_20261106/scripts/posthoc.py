"""Post-hoc checks (not in PLAN_ZH.md): no new scans -> ../analysis/posthoc.json."""
import json, sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).parent))
from analyze import load, num, repo, ROOT, CONDS, F0, DVA

sys.stdout.reconfigure(encoding="utf-8")
T = load(); qts = [50, 70, 100, 140]; seeds = list(range(1201, 1207)); out = {}
NM = dict(vJ="J", v1="E17", vd="chengD", v2="C2")
for sv, crs in (("lasso", ("vJ", "v1", "vd", "v2")), ("irls", ("v1", "vd", "v2"))):
    for cr in crs:
        out[f"{sv}:{NM[cr]}:K_by_Q"] = {str(q): float(np.nanmedian([(num(T[(sv, q, s, "W3_1")][cr]) - num(T[(sv, q, s, "W3_0")][cr])) / DVA for s in seeds])) for q in qts}
    # eq. 17: ratio of Q_sel to the frequency-ratio prediction, crossing picks only
    for c in CONDS:
        out[f"{sv}:E17:ratio:{c}"] = {}
        for q in qts:
            x = [(1 / num(T[(sv, q, s, c)]["v1"])) / (q * T[(sv, q, s, c)]["fw"] / F0) for s in seeds if T[(sv, q, s, c)]["cross"]]
            out[f"{sv}:E17:ratio:{c}"][str(q)] = [round(float(np.median(x)), 3), len(x)] if x else [None, 0]
    # K implied by the frequency ratio, against K measured, both-crossing pairs
    rows = {}
    for q in qts:
        pr, ob = [], []
        for s in seeds:
            a, b = T[(sv, q, s, "W3_0")], T[(sv, q, s, "W3_1")]
            if a["cross"] and b["cross"]:
                pr.append((1 / q) * F0 * (1 / b["fw"] - 1 / a["fw"]) / DVA); ob.append((num(b["v1"]) - num(a["v1"])) / DVA)
        rows[str(q)] = dict(pred=float(np.median(pr)) if pr else None, obs=float(np.median(ob)) if ob else None, n=len(pr))
    out[f"{sv}:E17:K_pred_vs_obs"] = rows
    # how often J / C2 / E17 solutions keep the grid edge, and lasso convergence
    if sv == "lasso":
        out["lasso:converged"] = float(np.mean([np.mean(r["converged"]) for k, r in T.items() if k[0] == "lasso"]))
# fp seen in the test data and the fw of each condition
out["fw_W3"] = {str(q): {c: float(np.median([T[("lasso", q, s, c)]["fw"] for s in seeds])) for c in ("W3_0", "W3_1")} for q in qts}
# repo pipeline: every pick is a first local minimum; how far from balance is it?
rp = {}
for f in sorted((ROOT / "results").glob("repo_*.json")):
    r = json.loads(f.read_text(encoding="utf-8")); e = np.array(r["E1"]); rp[f.stem] = dict(min_abs_E1=float(np.abs(e).min()), v=r["v_e17"], u=r["u"])
out["repo_picks"] = rp
(ROOT / "analysis" / "posthoc.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
for k, v in out.items():
    print(k, json.dumps(v, ensure_ascii=False))
