"""Summary of the experiment-12 development runs and the baseline selection (PLAN_ZH.md section 5).

Rule: inside each family the parameter with the largest mean S = (M1 + M1_lag)/2 over the 20 development datasets;
SCaK = the decoupled SC-a (amplitude on c_SR x k, phase fix to c_SR) at the selected B2 weight a, best k.
Writes analysis/dev_summary.json (all variants) and analysis/dev_selection.json (read by run12.py).
"""
import json
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
recs = [json.loads(p.read_text(encoding="utf-8")) for p in sorted((ROOT / "analysis" / "dev").glob("H_*.json"))]
assert len(recs) == 20, len(recs)
byid = {r["id"]: r for r in recs}
for p in sorted((ROOT / "analysis").glob("dev_ext[23]/H_*.json")):   # second and third extensions, one solve per file
    e = json.loads(p.read_text(encoding="utf-8")); byid[e["id"]]["scores"][e["name"]] = e["score"]
names = list(recs[0]["scores"])
assert all(set(r["scores"]) == set(names) for r in recs)
mean = {m: {s: float(np.mean([r["scores"][m][s] for r in recs])) for s in ("M1", "M1_lag", "S")} for m in names}


def best(prefix):
    fam = [m for m in names if m.startswith(prefix)]
    return max(fam, key=lambda m: mean[m]["S"])


b1, b2, b3, dq, sk = best("IQ_g"), best("SCa_a"), best("VDB_"), best("IQd_"), best("SRSCa_")
a2 = float(b2[5:]); kd = best(f"SCaK_a{a2}_") if any(m.startswith(f"SCaK_a{a2}_") for m in names) else None
sel = dict(B1=dict(name=b1, g_db=int(b1[4:])), B2=dict(name=b2, a=a2),
           SCaK=None if kd is None else dict(name=kd, a=a2, k=float(kd.split("_")[2].split("+")[0])),
           B3=dict(name=b3, win=float(b3.split("_")[1]), eps=float(b3.split("_")[2])),
           IQd=dict(name=dq, k=float(dq.split("_")[1]), g_db=int(dq.split("_g")[1])),
           SRSCa_k=dict(name=sk, k=float(sk.split("_")[1].split("+")[0])))
(ROOT / "analysis" / "dev_summary.json").write_text(json.dumps(dict(n=len(recs), mean=mean), indent=1), encoding="utf-8")
(ROOT / "analysis" / "dev_selection.json").write_text(json.dumps(sel, indent=1), encoding="utf-8")
for m in sorted(names, key=lambda m: -mean[m]["S"]):
    print(f"{m:22s} M1 {mean[m]['M1']:.3f}  M1_lag {mean[m]['M1_lag']:.3f}  S {mean[m]['S']:.3f}")
print(json.dumps(sel, indent=1))
