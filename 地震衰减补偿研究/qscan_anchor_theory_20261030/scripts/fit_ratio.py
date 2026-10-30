"""Post-hoc test: is the QScan-W bias multiplicative (c_sel = r c_true, i.e. v = B u with no intercept)?

Data (read from git HEAD): experiment 12 sections (exp 14 dev), experiment 14 T and L+N groups (QScan-W at
lambda 1/4, 1/16, 1/64), noise test x0.25 (lambda 1/16; x0.5 picks are identical). Per group and lambda: free fit v = A u + C and
log-space slope, bootstrap over datasets (2000, seed 20261030, 90 %). Across lambda within a group, compare
  add: v = A u + C_lambda   (common slope, lambda only shifts)
  mul: v = B_lambda u       (pure ratio, lambda only scales)
by residual sum of squares (both have one parameter per lambda plus at most one shared).
-> analysis/fit_ratio.json
"""
import json, subprocess
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]; REPO = ROOT.parents[1]; BASE = "地震衰减补偿研究"


def git_json(path):
    return json.loads(subprocess.run(["git", "show", f"HEAD:{BASE}/{path}"], cwd=REPO, capture_output=True, check=True).stdout)


def ls(path):
    out = subprocess.run(["git", "ls-tree", "--name-only", f"HEAD:{BASE}/{path}/"], cwd=REPO, capture_output=True, check=True)
    return [n for n in out.stdout.decode().split() if n.endswith(".json")]


E14 = "qscan_anchor_volve_20261029/results"; E3 = "qscan_anchor_theory_20261030/results3"
LAMS = ["0.25", "0.0625", "0.015625"]
dev = [git_json(f"{E14}/dev/{n}") for n in ls(f"{E14}/dev")]
test = [git_json(f"{E14}/test/{n}") for n in ls(f"{E14}/test")]
noise = [git_json(f"{E3}/{n}") for n in ls(E3) if n.endswith("__qscanw.json")]
rows = lambda rs: [dict(id=r["id"], c_true=r["c_true"], sel={k: v for k, v in r["qscan_w"].items() if k in LAMS}) for r in rs]
groups = {"dev": rows(dev), "T": rows([r for r in test if r["group"] == "T"]), "L+N": rows([r for r in test if r["group"] != "T"]),
          "noise x0.25": [dict(id=r["id"], c_true=r["c_true"], sel={"0.0625": r["c_sel"]}) for r in noise if r["noise_mult"] == 0.25]}


def free(rs, lam):
    u = np.array([1 / r["c_true"] for r in rs]); v = np.array([1 / r["sel"][lam] for r in rs])
    A, C = np.polyfit(u, v, 1); return float(A), float(C), float(np.polyfit(np.log(u), np.log(v), 1)[0])


def joint(rs, lams):
    u = np.array([1 / r["c_true"] for r in rs for _ in lams]); v = np.array([1 / r["sel"][l] for r in rs for l in lams])
    li = np.array([i for _ in rs for i in range(len(lams))]); D = np.eye(len(lams))[li]
    Xa = np.column_stack([u, D]); Xm = D * u[:, None]
    ra = v - Xa @ np.linalg.lstsq(Xa, v, rcond=None)[0]; rm = v - Xm @ np.linalg.lstsq(Xm, v, rcond=None)[0]
    return float(ra @ ra), float(rm @ rm)


def boot(rs, fn):
    rng = np.random.default_rng(20261030); out = []
    for _ in range(2000):
        sub = [rs[i] for i in rng.integers(0, len(rs), len(rs))]
        if len({r["c_true"] for r in sub}) > 1:
            out.append(fn(sub))
    return np.array(out)


ci = lambda x: [float(np.percentile(x, 5)), float(np.percentile(x, 95))]
out = {}
for g, rs in groups.items():
    lams = [l for l in LAMS if l in rs[0]["sel"]]; d = dict(n=len(rs), by_lambda={})
    for l in lams:
        A, C, s = free(rs, l); b = boot(rs, lambda x: free(x, l))
        d["by_lambda"][l] = dict(A=A, C=C, log_slope=s, A_ci=ci(b[:, 0]), C_ci=ci(b[:, 1]), log_slope_ci=ci(b[:, 2]),
                                 ratio=[float(1 / r["sel"][l] * r["c_true"]) for r in rs])
    if len(lams) > 1:
        rss_a, rss_m = joint(rs, lams); b = boot(rs, lambda x: np.subtract(*joint(x, lams)))
        d["joint"] = dict(rss_add=rss_a, rss_mul=rss_m, rss_add_minus_mul_ci=ci(b), frac_add_better=float(np.mean(b < 0)))
    out[g] = d
(ROOT / "analysis").mkdir(exist_ok=True)
(ROOT / "analysis" / "fit_ratio.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
for g, d in out.items():
    print(g, d["n"])
    for l, s in d["by_lambda"].items():
        print(f"  lam {float(l):.4f}: A={s['A']:.2f} {np.round(s['A_ci'], 2)}  C={s['C']:+.2f} {np.round(s['C_ci'], 2)}  "
              f"log={s['log_slope']:.2f} {np.round(s['log_slope_ci'], 2)}")
    if "joint" in d:
        j = d["joint"]; print(f"  joint RSS add={j['rss_add']:.3f} mul={j['rss_mul']:.3f} diff CI {np.round(j['rss_add_minus_mul_ci'], 3)} add better {j['frac_add_better']:.2f}")
