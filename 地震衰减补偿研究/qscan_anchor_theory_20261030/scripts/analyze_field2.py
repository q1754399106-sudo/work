"""Score the field kappa map (PLAN_FIELD2_ZH.md) from results11/ -> analysis/field2.json (+ printed summary)."""
import json, glob
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]; RNG = np.random.default_rng(20261106); NB = 2000
DV = 1 / 0.7 - 1 / 1.4


def ci_median(x):
    x = np.asarray(x, float); b = [np.median(x[RNG.integers(0, len(x), len(x))]) for _ in range(NB)]
    return [float(np.percentile(b, 5)), float(np.percentile(b, 95))]


if __name__ == "__main__":
    R = {}
    for f in glob.glob(str(ROOT / "results11" / "*.json")):
        r = json.loads(Path(f).read_text(encoding="utf-8")); R.setdefault(r["inline"], {})[r["c_anchor"]] = r
    rows, bad, n = {}, 0, 0
    for ln, s in sorted(R.items()):
        for r in s.values():
            n += 1; bad += bool(r["edge"] or not r["neighbourhood_converged"])
        if len(s) < 3:
            continue
        v = {a: s[a]["v_sel"] for a in s}; k = (v[0.7] - v[1.4]) / DV; c_sr = s[1.0]["c_sr"]; u = 1 / c_sr
        rows[ln] = dict(kappa=k, beta=v[1.0] - (1 - k) * u - k, c_sr=c_sr, c1=1 / v[1.0], ratio_c1_sr=(1 / v[1.0]) / c_sr,
                        monotone=bool(v[0.7] > v[1.0] > v[1.4]), v=v)
    K = [r["kappa"] for r in rows.values()]; B = [r["beta"] for r in rows.values()]
    out = dict(n_lines=len(rows), unreliable_share=bad / max(1, n), kappa_median=float(np.median(K)), kappa_ci=ci_median(K),
               kappa_range=[float(min(K)), float(max(K))], kappa_within_006=sum(abs(k - 0.675) <= 0.06 for k in K),
               monotone=sum(r["monotone"] for r in rows.values()), beta_median=float(np.median(B)), beta_ci=ci_median(B),
               beta_range=[float(min(B)), float(max(B))], beta_positive=sum(b > 0 for b in B),
               corr_beta_csr=float(np.corrcoef(B, [r["c_sr"] for r in rows.values()])[0, 1]),
               corr_beta_inline=float(np.corrcoef(B, list(rows))[0, 1]),
               ratio_c1_sr_median=float(np.median([r["ratio_c1_sr"] for r in rows.values()])),
               lines={str(k): {a: b for a, b in r.items() if a != "v"} for k, r in rows.items()})
    V = {}
    if out["unreliable_share"] <= 0.25:
        V["K1"] = bool(abs(out["kappa_median"] - 0.675) <= 0.12 and out["kappa_within_006"] >= 24)
        V["K2"] = bool(out["monotone"] >= 24)
        V["B1"] = bool(out["beta_positive"] >= 22 and out["beta_ci"][0] > 0.10)
    out["verdicts"] = V
    (ROOT / "analysis" / "field2.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
    print({k: v for k, v in out.items() if k != "lines"})
    for ln, r in out["lines"].items():
        print(ln, {a: (round(b, 3) if isinstance(b, float) else b) for a, b in r.items()})
