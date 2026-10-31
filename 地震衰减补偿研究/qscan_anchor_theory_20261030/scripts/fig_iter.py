"""Figures of the iterated re-anchoring test -> figures/iter.png (reads results5/ and analysis/iter.json)."""
import json, glob
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
plt.rcParams.update({"font.sans-serif": ["Noto Sans CJK SC", "WenQuanYi Zen Hei", "SimHei", "DejaVu Sans"],
                     "axes.unicode_minus": False, "font.size": 9})
S = json.loads((ROOT / "analysis" / "iter.json").read_text(encoding="utf-8"))
COL = dict(main="#2a6fbb", mid="#d08a1e", deep="#b8322e")
LAB = dict(main="锚窗 1.2–1.6 s（κ 0.675）", mid="1.8–2.2 s（κ 0.964）", deep="2.4–2.8 s（κ 1.14）")

fig, ax = plt.subplots(1, 3, figsize=(11.5, 3.6))
for arm in ("main", "mid", "deep"):
    if arm not in S:
        continue
    k, first = S[arm]["kappa"], True
    for did, d in S[arm]["per_dataset"].items():
        dv = np.abs(np.diff(d["path"]))
        ax[0].semilogy(np.arange(1, len(dv) + 1), np.maximum(dv, 1e-4), color=COL[arm], alpha=0.35, lw=0.8,
                       label=f"{LAB[arm]}，实测 r = {S[arm]['r']:.2f}" if first else None)
        first = False
    x = np.arange(1, 9); ax[0].semilogy(x, 0.2 * k ** (x - 1), color=COL[arm], lw=2.2, ls="--")
ax[0].set_xlabel("迭代步 k"); ax[0].set_ylabel("|Δv_k| = |v_k − v_{k−1}|"); ax[0].legend(fontsize=7, loc="lower left")
ax[0].set_title("增量按 κ 的速率变化（虚线：理论 κ^k）")

m = S.get("main")
if m:
    u = np.array([d["u"] for d in m["per_dataset"].values()])
    v1 = np.array([d["path"][1] for d in m["per_dataset"].values()])
    vf = np.array([d["path"][-1] for d in m["per_dataset"].values()])
    ex = np.array([d["extrap"] for d in m["per_dataset"].values()])
    uu = np.linspace(0.6, 1.5, 2)
    ax[1].plot(uu, uu, color="k", lw=0.8, label="v = u（无偏）")
    ax[1].scatter(u - 0.012, v1, s=14, color="#888", label=f"第 1 步（现流程）斜率 {m['A_step1']:.2f}")
    ax[1].scatter(u + 0.012, vf, s=14, color=COL["main"], label=f"第 8 步 斜率 {m['A_star']:.2f}")
    ax[1].plot(uu, uu + m["C_star"], color=COL["main"], ls="--", lw=1,
               label=f"u + C*，C* = {m['C_star']:.2f}")
    ax[1].set_xlabel("u = 1/c_true"); ax[1].set_ylabel("选出的 1/c"); ax[1].legend(fontsize=7)
    ax[1].set_title("迭代去掉压缩，留下截距")
    ax[2].scatter(ex, vf, s=14, color=COL["main"]); lim = [min(ex.min(), vf.min()) - 0.05, max(ex.max(), vf.max()) + 0.05]
    ax[2].plot(lim, lim, color="k", lw=0.8)
    ax[2].set_xlabel("一步反解 (v_1 − κ)/(1 − κ)"); ax[2].set_ylabel("迭代 8 步后的 v_8")
    ax[2].set_title(f"迭代 = 一步反解（中位差 {m['extrap_absdev_median']:.3f}）")
fig.tight_layout(); (ROOT / "figures").mkdir(exist_ok=True); fig.savefig(ROOT / "figures" / "iter.png", dpi=150)
print("figures/iter.png")
