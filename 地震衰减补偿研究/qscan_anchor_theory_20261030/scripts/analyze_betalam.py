"""Analysis of PLAN_BETALAM_ZH.md: beta(lambda) = beta0 + beta1 lambda, semi-synthetic vs field. -> analysis/betalam.json"""
import json, glob, statistics as st
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]; R = ROOT / "results23"; DV = 1 / 0.7 - 1 / 1.4
LAMS = (0.0625, 0.125, 0.25, 0.5)
recs = {Path(p).stem: json.loads(Path(p).read_text()) for p in glob.glob(str(R / "*.json"))}


def fit(pts):
    pts = [(l, b) for l, b in pts if b is not None]
    if len(pts) < 3:
        return None
    l, b = np.array(pts).T; b1, b0 = np.polyfit(l, b, 1)
    sm = [(x, y) for x, y in pts if x < 0.5]; b0s = float(np.polyfit(*np.array(sm).T, 1)[1]) if len(sm) >= 3 else None
    q = float(np.polyfit(l, b, 2)[0]) if len(pts) >= 4 else None
    return dict(beta0=float(b0), beta1=float(b1), beta0_small=b0s, quad=q, beta_half=dict(pts).get(0.5))


semi, field, n_all, n_edge = {}, {}, 0, 0
for c in (0.7, 1.0, 1.4):
    for s in (9971, 9972, 9973, 9974):
        pts = []
        for lam in LAMS:
            r = recs.get(f"S_c{c}_s{s}_lam{lam}"); n_all += 1
            if r is None or r["edge"]:
                n_edge += 1; pts.append((lam, None)); continue
            pts.append((lam, r["v_sel"] - r["u"]))
        semi[f"c{c}_s{s}"] = dict(points=pts, fit=fit(pts))
for p in sorted(glob.glob(str(ROOT / "results11" / "IL_*__w1.4_a1.0.json"))):
    r5 = json.loads(Path(p).read_text()); ln = r5["inline"]
    va = json.loads(Path(p.replace("a1.0", "a0.7")).read_text()); vb = json.loads(Path(p.replace("a1.0", "a1.4")).read_text())
    k = (va["v_sel"] - vb["v_sel"]) / DV; u = 1 / r5["c_sr"]; pts = []
    for lam in LAMS:
        r = r5 if lam == 0.5 else recs.get(f"F_IL_{ln}_lam{lam}"); n_all += 1
        if r is None or r["edge"]:
            n_edge += 1; pts.append((lam, None)); continue
        pts.append((lam, r["v_sel"] - (1 - k) * u - k))
    field[str(ln)] = dict(kappa=k, u_sr=u, points=pts, fit=fit(pts))


def med(d, key):
    xs = [x["fit"][key] for x in d.values() if x["fit"] and x["fit"][key] is not None]
    return (st.median(xs) if xs else None), len(xs)


summ = {}
for name, d in (("semi", semi), ("field", field)):
    summ[name] = {k: med(d, k) for k in ("beta0", "beta1", "beta0_small", "quad", "beta_half")}
G = summ["field"]["beta_half"][0] - summ["semi"]["beta_half"][0]
P = 0.5 * (summ["field"]["beta1"][0] - summ["semi"]["beta1"][0]); I = summ["field"]["beta0"][0] - summ["semi"]["beta0"][0]
if G < 0.10:
    verdict = "gap < 0.10"
elif P >= 2 * G / 3:
    verdict = "mainly sparsity (lambda-dependent)"
elif I >= 2 * G / 3:
    verdict = "mainly lambda-independent"
else:
    verdict = "both"
out = dict(semi=semi, field=field, summary=summ, G=G, P=P, I=I, verdict=verdict, unreliable_share=n_edge / n_all)
(ROOT / "analysis" / "betalam.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
for name in ("semi", "field"):
    print(name, {k: (None if v[0] is None else round(v[0], 3), v[1]) for k, v in summ[name].items()})
print("G %.3f  P %.3f  I %.3f  ->" % (G, P, I), verdict, " unreliable %.2f" % (n_edge / n_all))
for name, d in (("semi", semi), ("field", field)):
    print(name, "beta(lambda) medians:", [round(st.median([dict(x["points"])[l] for x in d.values() if dict(x["points"])[l] is not None]), 3) for l in LAMS])
