"""Analysis of PLAN_LAMBDA_ZH.md: (beta, kappa_hat) per line, window and lambda on U1320A. -> analysis/lambda.json"""
import json
from pathlib import Path
import numpy as np
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[1]
LAMS = (0.0625, 0.25, 0.5, 1.0, 2.0); LINES = ("3677-32", "4508-32"); VAS = (0.5, 1.0, 2.0)


def select(rec):
    v = np.asarray(rec["v"]); J = np.asarray(rec["J_trace"]).sum(axis=1); jm = int(np.argmin(J))
    if jm in (0, len(v) - 1):
        return None
    a = np.polyfit(v[jm - 1:jm + 2], J[jm - 1:jm + 2], 2)
    return float(np.clip(-a[1] / (2 * a[0]), v[jm - 1], v[jm + 1])) if a[0] > 0 else float(v[jm])


def load(line, w, lam, va):
    p = (ROOT / "results20" / f"{line}__{w}_va{va}.json") if lam == 0.5 else (ROOT / "results21" / f"{line}__{w}_lam{lam}_va{va}.json")
    return json.loads(p.read_text())


out = {}
for line in LINES:
    for w in ("w3.00", "w2.24"):
        rows = []
        for lam in LAMS:
            recs = {va: load(line, w, lam, va) for va in VAS}; s = {va: select(recs[va]) for va in VAS}
            ok = None not in s.values()
            rows.append(dict(lam=lam, v_sel=s, reliable=ok, kappa=((s[2.0] - s[0.5]) / 1.5 if ok else None),
                             beta=(s[1.0] - 1 if ok else None), converged=all(all(recs[va]["converged"]) for va in VAS)))
        good = [r for r in rows if r["reliable"]]
        b = [r["beta"] for r in good]; k = [r["kappa"] for r in good]
        rho = float(spearmanr(b, k)[0]) if len(good) >= 4 else None
        mono = all(b[i] < b[i + 1] for i in range(len(b) - 1)) if len(good) >= 4 else None
        out[f"{line}|{w}"] = dict(rows=rows, n_reliable=len(good), beta_range=(max(b) - min(b) if good else None),
                                  beta_monotone=mono, spearman=rho)

d = {l: out[f"{l}|w3.00"] for l in LINES}
N1 = {l: (None if d[l]["beta_monotone"] is None else bool(d[l]["beta_monotone"] and d[l]["beta_range"] >= 1.0)) for l in LINES}
if not all(N1.values()):
    N2 = N0 = None
else:
    g = [r for r in d["4508-32"]["rows"] if r["reliable"]]; lo = min(g, key=lambda r: r["beta"])
    N2 = bool(lo["kappa"] >= 0.75 and all(d[l]["spearman"] <= -0.8 for l in LINES))
    N0 = bool(all(abs(r["kappa"] - 0.607) <= 0.10 for r in g))
out["verdicts"] = dict(N1=N1, N2=N2, N0=N0)
(ROOT / "analysis" / "lambda.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
for key, r in out.items():
    if key == "verdicts":
        continue
    print(key, "beta range", None if r["beta_range"] is None else round(r["beta_range"], 2), "monotone", r["beta_monotone"],
          "spearman", None if r["spearman"] is None else round(r["spearman"], 2))
    for x in r["rows"]:
        print("   lam %-6s beta %s kappa %s  v %s conv %s" % (x["lam"], None if x["beta"] is None else round(x["beta"], 2),
              None if x["kappa"] is None else round(x["kappa"], 3), {k: None if v is None else round(v, 2) for k, v in x["v_sel"].items()}, x["converged"]))
print(out["verdicts"])
