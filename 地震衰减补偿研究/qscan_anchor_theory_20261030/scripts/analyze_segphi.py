"""Analysis of PLAN_SEGPHI_ZH.md: kappa_seg = a0 sum_b C_b (sum_b a / sum_b a^2) / sum_b C_b. -> analysis/segphi.json"""
import json, glob, statistics as st
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
HAT_U = {"3677-32|w2.24": 0.211, "3677-32|w3.00": 0.816, "4508-32|w2.24": 0.255, "4508-32|w3.00": 0.607}
HAT_V = {"B": 0.693, "S": 0.854, "D": 0.592, "D2": 0.624, "M": 0.614}


def curv(h, J):
    return float(2 * np.polyfit(h, J, 2)[0])


def analyse(r):
    h = [-0.5, -0.25, 0.0, 0.25, 0.5]; J0 = r["J0"]
    C = np.array([curv(h, [row[0], row[1], J0, row[2], row[3]]) for row in r["J_seg"]])
    Cu = curv(h, [r["J_uniform"][0], r["J_uniform"][1], J0, r["J_uniform"][2], r["J_uniform"][3]])
    t = np.asarray(r["t"]); a = np.asarray(r["a"]); e = np.asarray(r["seg_edges"]); nb = len(C)
    seg = np.clip(np.searchsorted(e, t, side="right") - 1, 0, nb - 1)
    w = np.array([a[seg == b].sum() / np.sum(a[seg == b] ** 2) for b in range(nb)])
    k = float(r["a0"] * np.sum(C * w) / np.sum(C))
    return dict(kappa_seg=k, share=(C / C.sum()).round(3).tolist(), coupling=float(1 - C.sum() / Cu), n_negative=int((C < 0).sum()),
                C_uniform=Cu, converged=r["all_converged"])


recs = [json.loads(Path(p).read_text()) for p in sorted(glob.glob(str(ROOT / "results22" / "*.json")))]
out = {"U": {}, "V": {}}
for r in recs:
    if r["kind"] == "U":
        key = f"{r['line']}|{r['window']}"; out["U"][key] = dict(kappa_hat=HAT_U[key], **analyse(r))
    else:
        out["V"].setdefault(r["profile"], {})[r["id"]] = analyse(r)
summ = {}
for p, d in out["V"].items():
    ks = [x["kappa_seg"] for x in d.values()]
    summ[p] = dict(n=len(ks), kappa_seg_median=st.median(ks), kappa_seg=sorted(round(x, 3) for x in ks), kappa_hat=HAT_V[p],
                   coupling_median=st.median(x["coupling"] for x in d.values()),
                   share_median=np.median([x["share"] for x in d.values()], axis=0).round(3).tolist())
out["V_summary"] = summ
U = out["U"]
S1 = None if len(U) < 4 else all(abs(x["kappa_seg"] - x["kappa_hat"]) <= 0.10 for x in U.values())
S2 = None if not {"4508-32|w3.00", "3677-32|w3.00"} <= set(U) else U["4508-32|w3.00"]["kappa_seg"] <= U["3677-32|w3.00"]["kappa_seg"] - 0.10
S3 = None if len(summ) < 5 else all(abs(x["kappa_seg_median"] - x["kappa_hat"]) <= 0.10 for x in summ.values())
out["verdicts"] = dict(S1=S1, S2=S2, S3=S3)
(ROOT / "analysis" / "segphi.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
for k, x in sorted(U.items()):
    print(k, "hat %.3f seg %.3f coupling %.2f neg %d conv %s" % (x["kappa_hat"], x["kappa_seg"], x["coupling"], x["n_negative"], x["converged"]))
    print("   share", x["share"])
for p, x in summ.items():
    print(p, "hat %.3f seg median %.3f %s coupling %.2f" % (x["kappa_hat"], x["kappa_seg_median"], x["kappa_seg"], x["coupling_median"]))
    print("   share", x["share_median"])
print(out["verdicts"])
