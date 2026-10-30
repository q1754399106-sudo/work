"""Theory figure -> figures/anchor_theory.png (descriptive).

(a) lambda 0.5: 1/c_sel against the anchor 1/c_a for each of the 10 sections (one line per section, colour by true c),
    with the predicted slope kappa = 0.675 drawn through the mean.
(b) slopes in 1/c space, predicted vs measured: anchor coefficient K (anchor experiment), slope on 1/c_true of the
    anchored QScan (anchor experiment and lib.qscan on experiments 12/14), A + K, and the slope of QScan-W.
"""
import json
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
import sys; sys.path.insert(0, str(ROOT.parent / "attenuation_run_20261019" / "dependencies"))
import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "DejaVu Sans"]; plt.rcParams["axes.unicode_minus"] = False
AN = json.loads((ROOT / "analysis" / "anchor.json").read_text(encoding="utf-8"))
FE = json.loads((ROOT / "analysis" / "fit_existing.json").read_text(encoding="utf-8"))
KAPPA = 0.675
recs = [json.loads(p.read_text(encoding="utf-8")) for p in sorted((ROOT / "results").glob("*.json"))]
CT = sorted({r["c_true"] for r in recs}); RAMP = ["#c6dbef", "#9ecae1", "#6baed6", "#3182bd", "#08519c"]

fig, (a, b) = plt.subplots(1, 2, figsize=(11.5, 4.8), gridspec_kw=dict(width_ratios=[1, 1.15]))
for d in AN["datasets"]:
    rs = sorted((r for r in recs if r["id"] == d), key=lambda r: r["c_anchor"]); ct = rs[0]["c_true"]
    a.plot([1 / r["c_anchor"] for r in rs], [1 / r["qscan_a"]["0.5"] for r in rs], color=RAMP[CT.index(ct)], lw=1.6,
           marker="o", ms=5, label=f"真 c = {ct}" if d.endswith("9701") else None)
va = np.array([1 / 1.4, 1.0, 1 / 0.7]); m = np.mean([1 / r["qscan_a"]["0.5"] for r in recs if r["c_anchor"] == 1.0])
a.plot(va, m + KAPPA * (va - 1), color="#c0504d", lw=2, ls="--", label=f"理论斜率 κ = {KAPPA}")
a.set_xlabel("锚点 1/c_a（估子波用的 Q 系数的倒数）"); a.set_ylabel("选出的 1/c")
a.set_title(f"(a) λ = 1/2：选择跟着锚点走（拟合 K = {AN['by_weight']['0.5']['K']:.2f}）", fontsize=11)
a.grid(alpha=0.25, lw=0.6); a.legend(fontsize=8, frameon=False, loc="upper left")

rows = [("锚点系数 K", KAPPA, [("锚点实验 λ=1/2", AN["by_weight"]["0.5"]["K"]), ("锚点实验 λ=1/16", AN["by_weight"]["0.0625"]["K"])]),
        ("锚定 QScan 对真值的斜率 A", 1 - KAPPA, [("锚点实验 λ=1/2", AN["by_weight"]["0.5"]["A"]), ("锚点实验 λ=1/16", AN["by_weight"]["0.0625"]["A"]),
                                         ("实验 12 剖面", FE["dev (exp 12 sections)"]["qscan"]["A"]), ("实验 14 T 组", FE["T"]["qscan"]["A"]),
                                         ("实验 14 S 组", FE["S"]["qscan"]["A"])]),
        ("和 A + K", 1.0, [("锚点实验 λ=1/2", AN["by_weight"]["0.5"]["A_plus_K"]), ("锚点实验 λ=1/16", AN["by_weight"]["0.0625"]["A_plus_K"])]),
        ("QScan-W 斜率", 1.0, [(f"{g} λ={l}", FE[k][f"W{v}"]["A"]) for g, k in (("实验 12", "dev (exp 12 sections)"), ("T", "T"), ("S", "S"))
                              for l, v in (("1/4", "0.25"), ("1/16", "0.0625"))])]
y = 0; ticks = []
for name, pred, obs in rows:
    ys = y + np.arange(len(obs)) * 0.32
    b.scatter([o[1] for o in obs], ys, color="#3182bd", s=28, zorder=3)
    for (lab, val), yy in zip(obs, ys):
        b.text(val + 0.03, yy, lab, fontsize=7, va="center", color="#555555")
    b.plot([pred, pred], [ys[0] - 0.25, ys[-1] + 0.25], color="#c0504d", lw=2)
    ticks.append((ys.mean(), name)); y = ys[-1] + 0.9
b.set_yticks([t[0] for t in ticks]); b.set_yticklabels([t[1] for t in ticks], fontsize=9); b.invert_yaxis()
b.set_xlim(0, 1.75); b.set_xlabel("1/c 空间的系数（红线 = 理论，均匀权重 κ = 0.675）")
b.set_title("(b) 理论值与实测", fontsize=11); b.grid(axis="x", alpha=0.25, lw=0.6)
fig.suptitle("子波锚定理论：1/c_sel = (1 − κ)/c_true + κ/c_a + β", fontsize=12)
fig.tight_layout(); (ROOT / "figures").mkdir(exist_ok=True); fig.savefig(ROOT / "figures" / "anchor_theory.png", dpi=150)
print("figures/anchor_theory.png")
