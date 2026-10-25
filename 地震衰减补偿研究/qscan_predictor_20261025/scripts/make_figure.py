"""Descriptive figure (after scoring): each feature against QScan's offset k on the primary scenes,
and the leave-one-model-out rel of the two candidates with the best mean."""
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
A, F = ROOT / "analysis", ROOT / "figures"
plt.rcParams.update({"font.sans-serif": ["Noto Sans CJK SC", "Microsoft YaHei", "DejaVu Sans"], "axes.unicode_minus": False,
                     "font.size": 9, "pdf.fonttype": 42, "savefig.facecolor": "white"})
S = pd.read_csv(A / "scene_features.csv")
P = S[S.scene != "Q100_08dB"]
carb = {"G2", "G4", "G6", "G8", "F2"}
names = {"Z1": "Z1 lasso 支撑数（每道）", "Z2": "Z2 J 的相对曲率", "Z3": "Z3 拟合项占比", "Z4": "Z4 补偿能量比",
         "Z5": "Z5 J 顶点亚网格位置", "Z6": "Z6 数据能量衰减"}
fig, axs = plt.subplots(2, 4, figsize=(13, 6.2), layout="constrained")
for ax, z in zip(axs.flat, names):
    for _, r in P.iterrows():
        c = "#D55E00" if r.model in carb else "#0072B2"
        ax.scatter(r[z], r.k + (0.12 if r.scene == "Q070_14dB" else -0.12), color=c, s=22, marker="o" if r.scene == "Q100_14dB" else "s")
        if r.scene == "Q100_14dB":
            ax.annotate(r.model, (r[z], r.k - 0.12), fontsize=7, xytext=(3, -3), textcoords="offset points")
    rho = P[z].rank().corr(P.k.rank())
    ax.set(xlabel=names[z], ylabel="偏移 k（网格步长）", title=f"Spearman ρ = {rho:+.2f}")
    ax.axhline(0, color="#777777", lw=.7)
R = pd.read_csv(A / "cv_results.csv")
for ax, cand in zip(axs.flat[6:], ["L:Z1", "T:Z2"]):
    d = R[(R.candidate == cand) & R.primary].sort_values("rel")
    ax.barh(np.arange(len(d)), d.rel, color=["#009E73" if v < 0 else "#D55E00" for v in d.rel])
    ax.set_yticks(np.arange(len(d)), [f"{m} {s[:4]}" for m, s in zip(d.held_out, d.scene)], fontsize=6.5)
    ax.axvline(0, color="#777777", lw=.7)
    ax.axvline(-.10, color="#777777", lw=.7, ls="--")
    ax.set(xlabel="留一模型 rel = MSE(修正)/MSE(QScan) − 1", title=f"{cand}：平均 {d.rel.mean():+.3f}，最坏 {d.rel.max():+.3f}")
fig.suptitle("主场景（圆点 Q0/14 dB，方块 Q0×0.7/14 dB；橙色为碳酸盐模型）")
for ext in ["png", "pdf"]:
    fig.savefig(F / f"features_vs_offset.{ext}", dpi=170)
