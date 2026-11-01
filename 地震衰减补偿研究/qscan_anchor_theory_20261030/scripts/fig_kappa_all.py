"""All pre-registered kappa predictions against the measured anchor coefficient -> figures/kappa_all.png.

Measured values are read from analysis/*.json where available; the earlier tests are taken from THEORY_ZH.md
(window.json, trunc.json, anchor.json, noise.json). The deep-lossy profile D is drawn in red (prediction failed);
the exploratory truncated-window points on D are drawn hollow (post hoc, not pre-registered).
"""
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]; AN = ROOT / "analysis"
plt.rcParams.update({"font.sans-serif": ["Noto Sans CJK SC", "WenQuanYi Zen Hei", "DejaVu Sans"],
                     "axes.unicode_minus": False, "font.size": 9})
anc = json.loads((AN / "anchor.json").read_text())["by_weight"]
win = json.loads((AN / "window.json").read_text()); tr = json.loads((AN / "trunc.json").read_text())["anchored"]
prof = json.loads((AN / "profile.json").read_text()); it = json.loads((AN / "iter.json").read_text())


def find(d, key):
    """First value of key in a nested dict (window.json layout differs between windows)."""
    if isinstance(d, dict):
        if key in d:
            return d[key]
        for v in d.values():
            r = find(v, key)
            if r is not None:
                return r
    return None


pts = [("锚点实验 λ=1/2", 0.675, anc["0.5"]["K"], anc["0.5"]["K_ci"], "k"),
       ("锚点实验 λ=1/16", 0.675, anc["0.0625"]["K"], anc["0.0625"]["K_ci"], "k"),
       ("截短分析窗 1.2–2.2 s", 0.804, tr["K"], tr["K_ci"], "k"),
       ("剖面 S（上覆层强吸收）", 0.880, prof["S"]["K"], [prof["S"]["ci_A_K_C"][0][1], prof["S"]["ci_A_K_C"][1][1]], "#2a6fbb"),
       ("剖面 D（深层强吸收）", 0.294, prof["D"]["K"], [prof["D"]["ci_A_K_C"][0][1], prof["D"]["ci_A_K_C"][1][1]], "#b8322e"),
       ("迭代增量比（main）", 0.675, it["main"]["r"], it["main"]["r_ci"], "#d08a1e")]
for tc, lab in (("2.0", "锚窗 1.8–2.2 s"), ("2.6", "锚窗 2.4–2.8 s")):
    w = win["windows"][tc]; pts.append((lab, w["kappa_pred"], w["K"], w["K_ci"], "k"))

fig, ax = plt.subplots(figsize=(5.6, 5.0)); lim = [0.2, 1.25]
ax.plot(lim, lim, color="#999", lw=0.8); ax.fill_between(lim, np.array(lim) - 0.1, np.array(lim) + 0.1, color="#eee")
for name, p, m, ci, col in pts:
    ax.errorbar(p, m, yerr=[[m - ci[0]], [ci[1] - m]], fmt="o", color=col, ms=5, capsize=2)
    ax.annotate(name, (p, m), xytext=(5, -3), textcoords="offset points", fontsize=7, color=col)
fd = AN / "field.json"
if fd.exists():
    F = json.loads(fd.read_text())
    for name, p, key in (("实测 10 条测线，默认锚窗", 0.675, "kappa"), ("实测，深锚窗", 1.14, "kappa_deep")):
        m, ci = F[key + "_median"], F[key + "_ci"]
        ax.errorbar(p, m, yerr=[[m - ci[0]], [ci[1] - m]], fmt="s", color="#2e8b57", ms=6, capsize=2)
        ax.annotate(name, (p, m), xytext=(-95, 4) if p > 1 else (8, -24), textcoords="offset points", fontsize=7, color="#2e8b57")
p2 = AN / "profile2.json"
if p2.exists():
    P2 = json.loads(p2.read_text())
    for prof in ("D2", "M"):
        if prof in P2:
            m, ci = P2[prof]["kappa_median"], P2[prof]["kappa_ci"]
            ax.errorbar(P2[prof]["kappa_uniform"], m, yerr=[[m - ci[0]], [ci[1] - m]], fmt="D", color="#7b4fb0", ms=5, capsize=2)
            ax.annotate(f"剖面 {prof}", (P2[prof]["kappa_uniform"], m), xytext=(5, 3), textcoords="offset points", fontsize=7, color="#7b4fb0")
dt = AN / "dtrunc.json"
if dt.exists():
    for row in json.loads(dt.read_text())["by_t1"]:
        ax.plot(row["kappa_uniform"], row["K"], "o", mfc="none", mec="#b8322e", ms=6)
        ax.annotate(f"D 截至 {row['t1']} s（事后）", (row["kappa_uniform"], row["K"]), xytext=(-78, 5) if row["t1"] < 2.4 else (6, 4),
                    textcoords="offset points", fontsize=6.5, color="#b8322e")
ax.set_xlim(lim); ax.set_ylim(lim); ax.set_xlabel("开跑前由几何算出的 κ = a₀⟨a⟩/⟨a²⟩"); ax.set_ylabel("实测锚点系数 K（90% 区间）")
ax.set_title("κ 的事先预测与实测（灰带：±0.1）")
fig.tight_layout(); fig.savefig(ROOT / "figures" / "kappa_all.png", dpi=160); print("figures/kappa_all.png")
