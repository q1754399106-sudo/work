"""Descriptive figure (after scoring): rel of the rule and of the forced 2-step shift per model and scene."""
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
plt.rcParams.update({"font.sans-serif": ["Noto Sans CJK SC", "Microsoft YaHei", "DejaVu Sans"], "axes.unicode_minus": False,
                     "font.size": 9, "pdf.fonttype": 42, "savefig.facecolor": "white"})
s = pd.read_csv(ROOT / "analysis" / "scene_scores.csv")
scenes = [("Q100_14dB", "Q0 / 14 dB（主）"), ("Q070_14dB", "Q0×0.7 / 14 dB（主）"), ("Q100_08dB", "Q0 / 8 dB"), ("Q100_20dB", "Q0 / 20 dB")]
fig, axs = plt.subplots(1, 4, figsize=(13, 3.8), layout="constrained", sharey=True)
models = sorted(s.model.unique())
x = np.arange(len(models))
for ax, (tag, title) in zip(axs, scenes):
    d = s[s.scene == tag].set_index("model").loc[models]
    ax.bar(x - .2, d.rel, .4, color="#0072B2", label="规则（Z1>t 才平移 2 步）")
    ax.bar(x + .2, d.rel_forced_S2, .4, color="#BBBBBB", label="无条件平移 2 步")
    ax.axhline(0, color="#555555", lw=.7)
    ax.axhline(.05, color="#D55E00", lw=.8, ls="--")
    ax.set_xticks(x, [m + ("*" if c else "") for m, c in zip(models, d.carbonate)])
    ax.set_title(title)
axs[0].set_ylabel("rel = MSE / MSE(QScan＋标准版) − 1")
axs[0].legend(fontsize=7.5, loc="upper left")
fig.suptitle("新模型 K1–K8（* 为碳酸盐模型；虚线为 A2 门槛 +5%）")
for ext in ["png", "pdf"]:
    fig.savefig(ROOT / "figures" / f"rule_rel.{ext}", dpi=170)
