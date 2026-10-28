"""Development summary of experiment 11 -> analysis/dev_summary.json.

Candidates for the improved pipeline (all with the SR factor c_SR of experiment 10):
  SRSCa          SC-a with c_SR
  SRSCa_k+ph     SC-a with c_SR x k (weaker amplitude compensation), then phase_fix from c_SR x k to c_SR, k = 1.15, 1.3
  QSCa+ph        QScan + SC-a, then phase_fix from the QScan factor to c_SR
Baseline QSCa (current pipeline). Selection: largest mean over datasets of (M1 + M1_lag) / 2.
"""
import sys, json
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import sc11  # noqa: E402
from sc11 import np, semi  # noqa: E402

cm = semi.curtain_cm(); Q = np.asarray(cm["Q"], float); rows = {}
for p in sorted((sc11.ROOT / "analysis" / "dev").glob("*.json")):
    r = json.loads(p.read_text(encoding="utf-8")); did = r["id"]; S = r["scores"]
    z = np.load(sc11.ROOT / "cache" / "dev" / f"{did}.npz"); ideal = np.load(sc11.DEV_OUT / f"{did}.npz")["ideal"].astype(float)
    sc = {"QSCa": S["QSCa"], "SRSCa": S["SCa_SRx1"], "SCa_true": S["SCa_truex1.0"], "IQ30_SR": S["IQ30_SR"], "Input": S["Input"]}
    for k in (1.15, 1.3):
        sc[f"SRSCa_{k}+ph"] = sc11.m1(sc11.phase_fix(z[f"SCa_SRx{k}"].astype(float), cm, Q, r["c_sr"] * k, r["c_sr"]), ideal)
    sc["QSCa+ph"] = sc11.m1(sc11.phase_fix(z["QSCa"].astype(float), cm, Q, r["qscan"], r["c_sr"]), ideal)
    sc["SCa_true_1.33+ph"] = sc11.m1(sc11.phase_fix(z["SCa_truex1.33"].astype(float), cm, Q, r["c_true"] * 1.33, r["c_true"]), ideal)
    rows[did] = dict(c_true=r["c_true"], c_sr=r["c_sr"], qscan=r["qscan"], scores=sc)
    print(did, flush=True)
methods = list(next(iter(rows.values()))["scores"])
summ = {}
for m in methods:
    a = np.array([[v["scores"][m]["M1"], v["scores"][m]["M1_lag"]] for v in rows.values()])
    d = np.array([[v["scores"][m]["M1"] - v["scores"]["QSCa"]["M1"], v["scores"][m]["M1_lag"] - v["scores"]["QSCa"]["M1_lag"]]
                  for v in rows.values()])
    summ[m] = dict(M1=float(a[:, 0].mean()), M1_lag=float(a[:, 1].mean()), score=float(a.mean()),
                   dM1_vs_QSCa=float(d[:, 0].mean()), dM1lag_vs_QSCa=float(d[:, 1].mean()),
                   min_dM1=float(d[:, 0].min()), min_dM1lag=float(d[:, 1].min()))
cands = ["SRSCa", "SRSCa_1.15+ph", "SRSCa_1.3+ph", "QSCa+ph"]
sel = max(cands, key=lambda m: summ[m]["score"])
out = dict(n=len(rows), summary=summ, candidates=cands, selected=sel, rows=rows)
(sc11.ROOT / "analysis" / "dev_summary.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
for m in methods:
    print(f"{m:18s}", {k: round(v, 3) for k, v in summ[m].items()})
print("selected:", sel)
