"""Sliding-wavelet-window field figure -> figures/slide_field.png (descriptive; written after the run, does not affect
the verdicts of PLAN_SLIDE_FIELD_ZH.md). Reads analysis/slide_field.json only.

(a) per-line alpha_a(t0) (thin lines, reliable line-windows only), median with 90 % bootstrap interval, frozen theory,
    alpha_a = 1.
(b) cross-line slope b(t0) of v_sel(c_a = 1) on u_SR with 90 % interval, against 1 - alpha_a(t0).
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
S = json.loads((ROOT / "analysis" / "slide_field.json").read_text(encoding="utf-8"))
T0S = [1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6]; KTH = [S["kappa_theory"][str(t)] for t in T0S]
W = S["windows"]; tw = [t for t in T0S if "kappa_median" in W[str(t)]]

fig, (a, b) = plt.subplots(1, 2, figsize=(11.5, 4.6))
for ln, rec in S["lines"].items():
    pts = [(t, rec[str(t)]["kappa"]) for t in T0S if rec.get(str(t), {}).get("status") == "ok"]
    if pts:
        a.plot(*zip(*pts), color="#c6dbef", lw=0.8, zorder=1)
med = [W[str(t)]["kappa_median"] for t in tw]; lo = [W[str(t)]["kappa_ci"][0] for t in tw]; hi = [W[str(t)]["kappa_ci"][1] for t in tw]
a.fill_between(tw, lo, hi, color="#3182bd", alpha=0.25, lw=0, zorder=2)
a.plot(tw, med, color="#08519c", lw=2, marker="o", ms=6, label="实测中位数（90% 区间）", zorder=3)
a.plot(T0S, KTH, color="#d95f02", lw=1.6, ls="--", marker="s", ms=4, label=r"理论 $\alpha_a$（运行前冻结）", zorder=3)
a.axhline(1, color="#777777", lw=1, ls=":")
a.plot([], [], color="#c6dbef", lw=0.8, label="单条测线")
a.set_xlabel("子波窗中心 $t_0$ (s)"); a.set_ylabel(r"两锚点 $\hat{\alpha}_a$"); a.set_xticks(T0S); a.legend(loc="upper left", fontsize=9)
a.set_title(r"(a) $\hat{\alpha}_a$ 随子波窗深度（26 条测线）", fontsize=11)

bs = [W[str(t)]["slope_v1_usr"] for t in tw]; ci = np.array([W[str(t)]["slope_v1_usr_ci"] for t in tw])
b.errorbar(tw, bs, yerr=[np.array(bs) - ci[:, 0], ci[:, 1] - np.array(bs)], color="#08519c", marker="o", ms=6, lw=1.5,
           capsize=3, label="测线间斜率 b($t_0$)（90% 区间）")
b.plot(T0S, [1 - k for k in KTH], color="#d95f02", lw=1.6, ls="--", marker="s", ms=4, label=r"$1-\alpha_a$（理论）")
b.axhline(0, color="#777777", lw=1, ls=":")
b.set_xlabel("子波窗中心 $t_0$ (s)"); b.set_ylabel("$v_{sel}(c_a=1)$ 对 $u_{SR}$ 的斜率"); b.set_xticks(T0S); b.legend(fontsize=9)
b.set_title("(b) 不用两锚点：以 SR 为真值代理（次要）", fontsize=11)
fig.tight_layout(); (ROOT / "figures").mkdir(exist_ok=True)
fig.savefig(ROOT / "figures" / "slide_field.png", dpi=200); print("figures/slide_field.png")
