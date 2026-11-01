"""Score the colour test (BETA_ZH.md step 6) from results16/ -> analysis/beta_color.json."""
import json, glob
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]; RNG = np.random.default_rng(20261108); NB = 2000
R = [json.loads(Path(f).read_text()) for f in sorted(glob.glob(str(ROOT / "results16" / "*.json")))]
L = [r for r in R if r["set"] == "L"]; W = [r for r in R if r["set"] == "W"]


def corr_ci(x, y):
    x, y = np.asarray(x), np.asarray(y); n = len(x); b = []
    for _ in range(NB):
        i = RNG.integers(0, n, n)
        if np.std(x[i]) > 0 and np.std(y[i]) > 0:
            b.append(np.corrcoef(x[i], y[i])[0, 1])
    return float(np.corrcoef(x, y)[0, 1]), [float(np.percentile(b, 5)), float(np.percentile(b, 95))]


out = dict(n_L=len(L), n_W=len(W))
if L:
    bL = [r["beta"] for r in L]; out["L"] = dict(beta_mean=float(np.mean(bL)), beta_sd=float(np.std(bL, ddof=1)),
                                                 kappa_mean=float(np.mean([r["kappa"] for r in L])))
    out["C1_corr"], out["C1_ci"] = corr_ci([r["beta_pred"] for r in L], bL)
    out["C1_slope"] = float(np.polyfit([r["beta_pred"] for r in L], bL, 1)[0])
    out["C3_corr"], out["C3_ci"] = corr_ci([r["kurt_sh"] for r in L], bL)
if W:
    bW = [r["beta"] for r in W]; out["W"] = dict(beta_mean=float(np.mean(bW)), beta_sd=float(np.std(bW, ddof=1)),
                                                 kappa_mean=float(np.mean([r["kappa"] for r in W])))
V = {}
if L:
    V["C1"] = bool(out["C1_ci"][0] > 0); V["C3"] = bool(out["C3_ci"][1] < 0)
if L and W:
    V["C2"] = bool(out["W"]["beta_sd"] <= 0.5 * out["L"]["beta_sd"] and abs(out["W"]["beta_mean"]) <= 0.08)
out["verdicts"] = V; out["rows"] = [{k: r[k] for k in ("job", "beta", "beta_pred", "kurt_sh", "kappa")} for r in R]
(ROOT / "analysis" / "beta_color.json").write_text(json.dumps(out, indent=1))
print({k: v for k, v in out.items() if k != "rows"})
for r in out["rows"]:
    print(r["job"], {k: round(v, 3) for k, v in r.items() if k != "job"})
