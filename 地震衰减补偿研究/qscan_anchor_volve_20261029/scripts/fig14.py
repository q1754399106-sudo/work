"""Experiment 14 summary figure -> figures/selection_T.png (descriptive).

(a) group T: selected factor vs true factor (log axes; geometric mean over the four seeds, thin lines per method;
    QScan-W weights on one sequential hue, dark = small weight).
(b) mean bias log(c_sel / c_true) on the development set vs on test groups T and S, per weight and for SR:
    the offset used for calibration (development) and the bias it meets on new data.
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
S = json.loads((ROOT / "analysis" / "summary.json").read_text(encoding="utf-8"))
X = json.loads((ROOT / "analysis" / "sensitivity.json").read_text(encoding="utf-8"))
LAMS = ["0.5", "0.25", "0.0625", "0.015625"]; LAB = {"0.5": "1/2", "0.25": "1/4", "0.0625": "1/16", "0.015625": "1/64"}
BLUES = ["#9ec5e8", "#5b9bd5", "#2c6fb2", "#123f73"]; GRAY, SRC = "#7a7a7a", "#c0504d"

fig, (a, b) = plt.subplots(1, 2, figsize=(11.5, 4.8))
cs = sorted(float(c) for c in X["by_c"])
a.plot(cs, cs, color="#bbbbbb", lw=1, ls="--", label="理想（选出 = 真值）")
series = [("QScan", "QScan（现流程）", GRAY, "s")] + [(f"W{float(l)}", f"QScan-W，λ = {LAB[l]}", BLUES[i], "o")
                                                    for i, l in enumerate(LAMS)] + [("SR", "SR", SRC, "D")]
for k, lab, col, m in series:
    y = [X["by_c"][str(c)][k]["mean_sel"] for c in cs]
    a.plot(cs, y, color=col, lw=2, marker=m, ms=6, label=f"{lab}  斜率 {S['groups']['T']['describe'][k]['slope']:.2f}")
a.set_xscale("log"); a.set_yscale("log"); a.set_xticks(cs); a.set_xticklabels([f"{c:g}" for c in cs])
yt = [0.4, 0.6, 0.8, 1.0, 1.5, 2.0]; a.set_yticks(yt); a.set_yticklabels([f"{v:g}" for v in yt]); a.minorticks_off()
a.set_xlabel("真 c"); a.set_ylabel("选出的 c（4 个种子几何平均）"); a.set_title("(a) T 组：选出的 c 对真 c", fontsize=11)
a.grid(alpha=0.25, lw=0.6); a.legend(fontsize=8, frameon=False, loc="upper left")

keys = LAMS + ["SR"]; xs = np.arange(len(keys)); wd = 0.26
for j, (name, col) in enumerate((("开发集（定校正量）", "#bdbdbd"), ("检验 T 组", "#d08c2c"), ("检验 S 组", "#7a4a12"))):
    if j == 0:
        v = [X["transfer"]["T"][k]["dev"] for k in keys]
    else:
        v = [X["transfer"]["T" if j == 1 else "S"][k]["test"] for k in keys]
    b.bar(xs + (j - 1) * wd, v, wd - 0.03, color=col, label=name)
b.axhline(0, color="#444444", lw=0.8)
b.set_xticks(xs); b.set_xticklabels([f"QScan-W\nλ = {LAB[k]}" for k in LAMS] + ["SR"])
b.set_ylabel("平均偏差 log(选出 c / 真 c)"); b.set_title("(b) 偏差在开发集与新数据上的差别", fontsize=11)
b.grid(axis="y", alpha=0.25, lw=0.6); b.legend(fontsize=8, frameon=False, loc="lower right")
fig.suptitle("实验 14：每个 c 重估子波后 QScan 的选择不再被压向 1，但校正量搬不到新数据上", fontsize=12)
fig.tight_layout(); (ROOT / "figures").mkdir(exist_ok=True); fig.savefig(ROOT / "figures" / "selection_T.png", dpi=150)
print("figures/selection_T.png")
