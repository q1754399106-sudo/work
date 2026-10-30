"""analysis/summary.json and the printed table for PREDICTIONS_ZH.md P1-P3."""
import json
from libc import ROOT, IDS, np

R = {d: json.loads((ROOT / "results" / f"{d}.json").read_text(encoding="utf-8")) for d in IDS}
S = {d: {m: v["S"] for m, v in r["scores"].items()} for d, r in R.items()}
coll = [d for d in IDS if R[d]["fs"] < R[d]["fd"]]
okC = [d for d in IDS if d[1] == "C" and d not in coll]; okAB = [d for d in IDS if d[1] in "AB" and d not in coll]


def diff(ids, a, b):
    x = np.array([S[d][a] - S[d][b] for d in ids]); return dict(mean=float(x.mean()), positive=int((x > 0).sum()), n=len(x),
                                                              values=[round(float(v), 4) for v in x])


out = dict(collapsed=coll, normal_C=okC, normal_AB=okAB,
           reproduce=max(abs(S[d]["SC_0_mu10"] - R[d]["e15_SCa_oracle"]) for d in IDS if "e15_SCa_oracle" in R[d]),
           P1=dict(SC=diff(coll, "SC_1_mu10", "SC_0_mu10"), VDB=diff(coll, "VDB_1", "VDB_0")),
           P2=dict(C=diff(okC, "SC_1_mu0", "SC_1_mu10"), AB=diff(okAB, "SC_1_mu0", "SC_1_mu10"),
                   collapsed=diff(coll, "SC_1_mu0", "SC_1_mu10")),
           P3=dict(before=diff(IDS, "IQ70", "SC_0_mu10"), after=diff(IDS, "IQ70", "SC_1_mu0")),
           extra=dict(VDB1_vs_SC1mu0=diff(IDS, "VDB_1", "SC_1_mu0"), IQ70_vs_VDB1=diff(IDS, "IQ70", "VDB_1"),
                      mu0_orig_band_C=diff(okC, "SC_0_mu0", "SC_0_mu10")),
           mean_S={m: float(np.mean([S[d][m] for d in IDS])) for m in S[IDS[0]]}, S=S)
p1 = out["P1"]; p2 = out["P2"]; p3 = out["P3"]
out["verdicts"] = dict(
    P1="成立" if all(p1[k]["mean"] >= 0.08 and p1[k]["positive"] >= 5 for k in ("SC", "VDB")) else "不成立",
    P2="成立" if p2["C"]["mean"] >= 0.05 and p2["C"]["positive"] == 4 and abs(p2["AB"]["mean"]) < 0.03 else "不成立",
    P3="成立" if p3["after"]["mean"] < 0.5 * p3["before"]["mean"] else "不成立")
(ROOT / "analysis" / "summary.json").write_text(json.dumps(out, indent=1, ensure_ascii=False), encoding="utf-8")
M = list(S[IDS[0]])
print("id".ljust(13), " ".join(m.rjust(9) for m in M))
for d in IDS:
    print(d.ljust(13) + ("*" if d in coll else " "), " ".join(f"{S[d][m]:9.3f}" for m in M))
print({k: v for k, v in out.items() if k not in ("S",)})
