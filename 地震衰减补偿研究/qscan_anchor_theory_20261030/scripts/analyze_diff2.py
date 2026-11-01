"""Score the differential test v2 (BETA_ZH.md step 5) from results15/ -> analysis/diff2.json."""
import json, glob
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]; K0 = 0.675; DV = 1 / 0.7 - 1 / 1.4


def quantities(r):
    v1 = r["w1.4_a1.0"]["v_sel"]; a, b = r["w2.6_a0.7"]["v_sel"], r["w2.6_a1.4"]["v_sel"]
    kd = (a - b) / DV; v1d = b + kd * (1 - 1 / 1.4)
    return dict(u_sr=1 / r["c_sr"], v1=v1, kappa_d=kd, u_pair=1 + (v1 - v1d) / (kd - K0))


def slopes(rows):
    """Pooled least-squares slope on Delta u with a separate intercept per base."""
    out = {}
    for q in ("u_sr", "v1", "u_pair"):
        num = den = 0.0
        for base, rs in rows.items():
            du = np.array([r["du"] for r in rs]); y = np.array([quantities(r)[q] for r in rs])
            du, y = du - du.mean(), y - y.mean(); num += float(du @ y); den += float(du @ du)
        out[q] = num / den if den else float("nan")
    return out


if __name__ == "__main__":
    R = {}
    for f in glob.glob(str(ROOT / "results15" / "*.json")):
        r = json.loads(Path(f).read_text()); R.setdefault(r["base"], []).append(r)
    S = {b: v for b, v in R.items() if b.startswith("S") and len(v) == 3}; F = {b: v for b, v in R.items() if b.startswith("F") and len(v) == 3}
    out = dict(n_semi=len(S), n_field=len(F), semi=slopes(S) if S else None, field=slopes(F) if F else None,
               per_base={b: [dict(du=r["du"], u_true=r["u_true"], **{k: round(v, 4) for k, v in quantities(r).items()},
                                  edges=[r[t]["edge"] for t in ("w1.4_a1.0", "w2.6_a0.7", "w2.6_a1.4")])
                             for r in sorted(v, key=lambda x: x["du"])] for b, v in R.items()})
    V = {}
    if S:
        s = out["semi"]; V["V"] = bool(abs(s["v1"] - 0.32) <= 0.12 and s["u_sr"] >= 0.6)
        if V["V"] and F:
            f = out["field"]
            V["D1"] = bool(abs(f["u_sr"] - s["u_sr"]) <= 0.2); V["D2"] = bool(abs(f["v1"] - s["v1"]) <= 0.10)
            V["D3"] = bool(abs(f["u_pair"] - s["u_pair"]) <= 0.3)
    out["verdicts"] = V
    (ROOT / "analysis" / "diff2.json").write_text(json.dumps(out, indent=1))
    print("semi", out["semi"]); print("field", out["field"]); print("verdicts", V)
    for b, rs in out["per_base"].items():
        print(b, [(r["du"], r["u_sr"], r["v1"], r["u_pair"]) for r in rs])
