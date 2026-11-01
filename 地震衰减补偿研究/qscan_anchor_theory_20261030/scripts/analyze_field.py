"""Score the field kappa test (PLAN_FIELD_ZH.md) from results9/ -> analysis/field.json (+ printed summary)."""
import json, glob, sys
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]; RNG = np.random.default_rng(20261104); NB = 2000
E13 = ROOT.parent / "field_baselines_volve_20261029"
DV = 1 / 0.7 - 1 / 1.4


def load():
    R = {}
    for f in glob.glob(str(ROOT / "results9" / "*.json")):
        r = json.loads(Path(f).read_text(encoding="utf-8")); R.setdefault(r["inline"], {})[(r["window"], r["c_anchor"])] = r
    return R


def sr(ln):
    sys.path.insert(0, str(ROOT / "scripts"))
    import env15  # noqa: F401  (paths of the experiment-13 snapshot)
    sys.path.insert(0, str(E13 / "scripts"))
    import lib13
    return lib13.c_sr(ln)


def ci_median(x):
    x = np.asarray(x, float); b = [np.median(x[RNG.integers(0, len(x), len(x))]) for _ in range(NB)]
    return [float(np.percentile(b, 5)), float(np.percentile(b, 95))]


if __name__ == "__main__":
    R = load(); lines = sorted(R); out = dict(lines={}, n_lines=len(lines))
    e13 = {ln: json.loads((E13 / "results" / "qscan" / f"IL_{ln}.json").read_text())["c"] for ln in lines}
    k1, k2, mono, rep, bad = [], [], 0, 0, {"w1.4": 0, "w2.6": 0}; n = {"w1.4": 0, "w2.6": 0}
    for ln in lines:
        s = R[ln]; v = {k: r["v_sel"] for k, r in s.items()}
        for k, r in s.items():
            n[k[0]] += 1; bad[k[0]] += bool(r["edge"] or not r["neighbourhood_converged"])
        row = dict(v={f"{w}_a{a}": x for (w, a), x in v.items()})
        if ("w1.4", 1.0) in s:
            row["grid_a1"] = s[("w1.4", 1.0)]["v_grid"]; row["c_e13"] = e13[ln]
            rep += abs(1 / row["grid_a1"] - e13[ln]) < 1e-6
        if ("w1.4", 0.7) in v and ("w1.4", 1.4) in v:
            k = (v[("w1.4", 0.7)] - v[("w1.4", 1.4)]) / DV; row["kappa"] = k; k1.append(k)
            if ("w1.4", 1.0) in v:
                mono += v[("w1.4", 0.7)] > v[("w1.4", 1.0)] > v[("w1.4", 1.4)]
                vs = (v[("w1.4", 1.0)] - k) / (1 - k) if k < 1 else float("nan")
                row["v_star"] = vs; row["c_star"] = 1 / vs if vs > 0 else float("nan"); row["c_sr"] = sr(ln)
        if ("w2.6", 0.7) in v and ("w2.6", 1.4) in v:
            row["kappa_deep"] = (v[("w2.6", 0.7)] - v[("w2.6", 1.4)]) / DV; k2.append(row["kappa_deep"])
        out["lines"][str(ln)] = row
    deeper = sum(r.get("kappa_deep", -9) > r.get("kappa", 9) for r in out["lines"].values())
    out.update(kappa_median=float(np.median(k1)) if k1 else None, kappa_ci=ci_median(k1) if k1 else None,
               kappa_deep_median=float(np.median(k2)) if k2 else None, kappa_deep_ci=ci_median(k2) if k2 else None,
               monotone=mono, reproduced=rep, deeper=deeper,
               unreliable={w: bad[w] / max(1, n[w]) for w in bad})
    V = {"F0": rep >= 9}
    if V["F0"]:
        if out["unreliable"]["w1.4"] <= 0.25 and k1:
            V["F1"] = bool(abs(out["kappa_median"] - 0.675) <= 0.12 and mono >= 9)
        if out["unreliable"]["w2.6"] <= 0.25 and k2:
            V["F2"] = bool(abs(out["kappa_deep_median"] - 1.14) <= 0.12 and deeper >= 9)
    out["verdicts"] = V
    (ROOT / "analysis" / "field.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
    print({k: v for k, v in out.items() if k != "lines"})
    for ln, r in out["lines"].items():
        print(ln, {k: (round(x, 3) if isinstance(x, float) else x) for k, x in r.items() if k != "v"})
