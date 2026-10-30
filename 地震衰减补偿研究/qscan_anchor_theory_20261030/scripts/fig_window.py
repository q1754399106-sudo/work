"""Anchor-window figure -> figures/window_kappa.png (descriptive).

(a) selected c against true c (geometric mean over the two seeds), anchor c_a = 1.4, lambda 1/2, for the three wavelet
    windows: 1.2-1.6 s (anchor experiment, results/), 1.8-2.2 s and 2.4-2.8 s (results2/).
(b) measured anchor coefficient K (90 % bootstrap interval) against kappa predicted from the model Q, three windows.
"""
import json, subprocess
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]; REPO = ROOT.parents[1]
import sys; sys.path.insert(0, str(ROOT.parent / "attenuation_run_20261019" / "dependencies"))
import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "DejaVu Sans"]; plt.rcParams["axes.unicode_minus"] = False
W = json.loads((ROOT / "analysis" / "window.json").read_text(encoding="utf-8"))
r2 = [json.loads(p.read_text(encoding="utf-8")) for p in sorted((ROOT / "results2").glob("*.json"))]
ids = sorted({r["id"] for r in r2}); CT = sorted({r["c_true"] for r in r2})


def git_json(path):
    return json.loads(subprocess.run(["git", "show", f"HEAD:地震衰减补偿研究/qscan_anchor_theory_20261030/{path}"], cwd=REPO,
                                     capture_output=True, check=True).stdout)


r1 = [git_json(f"results/{d}__a1.4.json") for d in ids]
COL = {"1.4": "#9ecae1", "2.0": "#3182bd", "2.6": "#08306b"}
LAB = {"1.4": "锚窗 1.2–1.6 s（现流程）", "2.0": "锚窗 1.8–2.2 s", "2.6": "锚窗 2.4–2.8 s"}
fig, (a, b) = plt.subplots(1, 2, figsize=(11.5, 4.8))
series = {"1.4": [(r["c_true"], r["qscan_a"]["0.5"]) for r in r1]}
for tc in ("2.0", "2.6"):
    series[tc] = [(r["c_true"], r["c_sel"]) for r in r2 if str(r["window"][2]) == tc and r["c_anchor"] == 1.4]
a.plot(CT, CT, color="#bbbbbb", lw=1, ls="--", label="理想（选出 = 真值）")
for tc, pts in series.items():
    y = [float(np.exp(np.mean(np.log([s for c, s in pts if c == ct])))) for ct in CT]
    a.plot(CT, y, color=COL[tc], lw=2, marker="o", ms=6, label=LAB[tc])
a.set_xscale("log"); a.set_yscale("log"); a.set_xticks(CT); a.set_xticklabels([f"{c:g}" for c in CT])
yt = [0.6, 0.8, 1.0, 1.4, 2.0]; a.set_yticks(yt); a.set_yticklabels([f"{v:g}" for v in yt]); a.minorticks_off()
a.set_xlabel("真 c"); a.set_ylabel("选出的 c（2 个种子几何平均）"); a.set_title("(a) 锚点 c_a = 1.4，λ = 1/2：锚窗越深，斜率越低直至反向", fontsize=10.5)
a.grid(alpha=0.25, lw=0.6); a.legend(fontsize=8, frameon=False, loc="lower right")

pts = [("1.4", 0.675, W["baseline_1p4"]["K"], W["baseline_1p4"]["K_ci"])] + \
      [(tc, W["windows"][tc]["kappa_pred"], W["windows"][tc]["K"], W["windows"][tc]["K_ci"]) for tc in ("2.0", "2.6")]
b.plot([0.5, 1.3], [0.5, 1.3], color="#bbbbbb", lw=1, ls="--", label="实测 = 预测")
for tc, p, k, ci in pts:
    b.errorbar(p, k, yerr=[[k - ci[0]], [ci[1] - k]], fmt="o", ms=8, color=COL[tc], capsize=4, lw=2, label=LAB[tc])
    b.annotate(f"预测 {p:.3f}\n实测 {k:.2f}", (p, k), textcoords="offset points", xytext=(10, -22), fontsize=8, color="#444444")
b.set_xlim(0.5, 1.3); b.set_ylim(0.5, 1.3); b.set_aspect("equal")
b.set_xlabel(r"预测 $\kappa = a_0\langle a\rangle/\langle a^2\rangle$（由模型 Q 计算，无拟合参数）"); b.set_ylabel("实测锚点系数 K（90% 区间）")
b.set_title("(b) 三个锚窗：预测与实测", fontsize=10.5); b.grid(alpha=0.25, lw=0.6); b.legend(fontsize=8, frameon=False, loc="upper left")
fig.suptitle("锚窗深度检验：κ 由公式事先给出，实测与之相差不到 0.02", fontsize=12)
fig.tight_layout(); (ROOT / "figures").mkdir(exist_ok=True); fig.savefig(ROOT / "figures" / "window_kappa.png", dpi=150)
print("figures/window_kappa.png")
