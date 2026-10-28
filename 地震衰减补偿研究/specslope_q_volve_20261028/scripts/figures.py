"""Figures of experiment 10 (presentation only). Writes figures/*.png."""
import sys, json
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import sq, run, analyze  # noqa: E402
from sq import np, semi, TW  # noqa: E402
sys.path.insert(0, str(sq.ROOT.parent / "field_application_20261028" / "cache" / "deps"))   # matplotlib
import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib import font_manager  # noqa: E402

font_manager.fontManager.addfont(r"C:\Windows\Fonts\msyh.ttc")
plt.rcParams.update({"font.family": "Microsoft YaHei", "axes.unicode_minus": False, "font.size": 9,
                     "axes.edgecolor": "#8a8a85", "axes.labelcolor": "#3d3d3a", "xtick.color": "#5f5e5a",
                     "ytick.color": "#5f5e5a", "axes.spines.top": False, "axes.spines.right": False})
ROOT = sq.ROOT; FIG = ROOT / "figures"; FIG.mkdir(exist_ok=True)
BLUE, ORANGE, GREEN, YELLOW, PINK, VIOLET = "#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#6250d6"
GRID = "#e8e7e2"; INK = "#8a8a85"


def accuracy():
    H, W = analyze.load("H"), analyze.load("W")
    fig, ax = plt.subplots(figsize=(6.4, 5.2)); rng = np.random.default_rng(0)
    g = np.array([0.45, 1.7]); ax.fill_between(g, g / 1.25, g * 1.25, color=GRID, lw=0); ax.plot(g, g, color=INK, lw=1)
    first = dict(H=True, W=True, Q=True)
    for r in H:
        for i, e in enumerate(r["sr"].values()):
            x = r["c_true"] * np.exp(rng.uniform(-.04, .04))
            ax.errorbar(x, e["c"], yerr=[[e["c"] - e["ci90"][0]], [e["ci90"][1] - e["c"]]], fmt="o", ms=5, color=BLUE,
                        ecolor=BLUE, elinewidth=.8, alpha=.75, mec="white", mew=.8,
                        label="SR（检验集 H，含 90% 区间）" if first["H"] else None); first["H"] = False
        if "qscan" in r:
            ax.plot(r["c_true"] * np.exp(rng.uniform(-.04, .04)), r["qscan"]["factor"], "s", ms=6, color=ORANGE, mec="white",
                    mew=1, label="QScan" if first["Q"] else None); first["Q"] = False
    for r in W:
        e = r["sr"]["n0"]
        ax.plot(r["c_true"] * np.exp(rng.uniform(-.04, .04)), e["c"], "D", ms=5, color=GREEN, mec="white", mew=.8,
                label="SR（真实井段反射系数 W）" if first["W"] else None); first["W"] = False
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlim(*g); ax.set_ylim(0.35, 2.4)
    for a in (ax.xaxis, ax.yaxis):
        a.set_major_formatter(matplotlib.ticker.ScalarFormatter()); a.set_minor_formatter(matplotlib.ticker.NullFormatter())
    ax.set_xticks(run.C_TEST); ax.set_yticks([0.4, 0.55, 0.7, 0.9, 1.2, 1.4, 2.0])
    ax.set_xlabel("真 c（真 Q / 模型 Q）"); ax.set_ylabel("估计的 c"); ax.set_title("估计 Q 系数：SR 与 QScan（灰带为 ±25%）", fontsize=10)
    ax.legend(frameon=False, fontsize=8, loc="upper left"); fig.tight_layout(); fig.savefig(FIG / "accuracy.png", dpi=150); plt.close(fig)


def compensation():
    H = [r for r in analyze.load("H") if "scores" in r]
    ms = ["Input", "IQ30_model", "IQ30_QScan", "IQ30_SR", "SCa_SR", "IQ30_true"]
    lab = ["输入", "反 Q（模型 Q）", "反 Q（QScan）", "反 Q（SR）", "SC-a（SR）", "反 Q（真 Q）"]
    col = [BLUE, INK, ORANGE, VIOLET, PINK, GREEN]
    fig, axs = plt.subplots(1, 2, figsize=(12, 4.2), sharey=True)
    for ax, key, t in zip(axs, ("M1", "M1_lag"), ("零时移", "允许 ±20 ms 整体时移")):
        for j, c in enumerate(run.C_TEST):
            rs = [r for r in H if r["c_true"] == c]
            for i, (m, cl) in enumerate(zip(ms, col)):
                v = [r["scores"][m][key] for r in rs]
                ax.plot(j + (i - 2.5) * 0.13, np.mean(v), "o", ms=7, color=cl, mec="white", mew=1.2,
                        label=lab[i] if j == 0 else None)
                ax.plot([j + (i - 2.5) * 0.13] * 2, [min(v), max(v)], color=cl, lw=1.2)
        ax.set_xticks(range(len(run.C_TEST))); ax.set_xticklabels([f"c = {c}" for c in run.C_TEST])
        ax.set_title(f"M1（{t}；点为 4 个种子平均，竖线为范围）", fontsize=9); ax.grid(axis="y", color=GRID, lw=.6)
        ax.axhline(0, color=INK, lw=.8)
    axs[0].set_ylabel("与理想目标的相关（2.30–2.80 s）"); axs[0].legend(frameon=False, fontsize=8, ncol=2)
    fig.tight_layout(); fig.savefig(FIG / "compensation.png", dpi=150); plt.close(fig)


def real():
    L = json.loads((ROOT / "results" / "real" / "lines.json").read_text(encoding="utf-8"))
    R = json.loads((ROOT / "results" / "real" / "inline10141.json").read_text(encoding="utf-8"))
    fig, axs = plt.subplots(1, 2, figsize=(12, 4.2), gridspec_kw=dict(width_ratios=[1.5, 1]))
    ax = axs[0]; ks = sorted(L, key=int)
    for k in ks:
        e = L[k]; hl = int(k) == semi.INLINE
        ax.errorbar(int(k), e["c"], yerr=[[e["c"] - e["ci90"][0]], [e["ci90"][1] - e["c"]]], fmt="o", ms=6,
                    color=VIOLET if not hl else ORANGE, mec="white", mew=1, elinewidth=1.2, capsize=2)
    med = float(np.median([L[k]["c"] for k in ks if int(k) != semi.INLINE]))
    ax.axhline(med, color=VIOLET, lw=1, ls="--"); ax.text(int(ks[-1]) + 6, med, f"中位数 {med:.2f}", color="#3d3d3a", fontsize=8, va="center")
    ax.axhline(R["qscan"]["factor"], color=ORANGE, lw=1, ls=":")
    ax.text(int(ks[-1]), R["qscan"]["factor"] * 0.97, f"QScan（inline 10141）{R['qscan']['factor']:.2f}", ha="right", va="top",
            color="#3d3d3a", fontsize=8)
    ax.axhline(1, color=INK, lw=.8); ax.set_xlabel("inline"); ax.set_ylabel("SR 估计的 c（含 90% 区间）")
    ax.set_title("Volve 实测：各测线的 Q 系数（橙色为应用测线 10141）", fontsize=9); ax.grid(axis="y", color=GRID, lw=.6)
    ax = axs[1]; f = np.array(R["f"])
    for m, lab, cl in [("Input", "输入", BLUE), ("IQ30_QScan", "反 Q（QScan）", ORANGE), ("IQ30_SR", "反 Q（SR）", VIOLET),
                       ("SCa_SR", "SC-a（SR）", PINK)]:
        ax.semilogy(f, np.maximum(R["spectra"][m]["snr"], 1e-2), color=cl, lw=1.6, label=lab)
    ax.axhline(1, color=INK, lw=.8, ls="--"); ax.set_xlim(0, 110); ax.set_ylim(1e-2, 300)
    ax.set_xlabel("频率 (Hz)"); ax.set_ylabel("相干信噪比（2.0–3.0 s）"); ax.set_title("inline 10141 补偿后", fontsize=9)
    ax.legend(frameon=False, fontsize=8); fig.tight_layout(); fig.savefig(FIG / "real.png", dpi=150); plt.close(fig)
    o = np.load(ROOT / "cache" / "out" / "real_10141.npz"); keys = [("Input", "输入"), ("IQ30_QScan", "反 Q 30 dB（QScan c）"),
                                                                  ("IQ30_SR", "反 Q 30 dB（SR c）"), ("SCa_SR", "SC-a（SR c）")]
    w = semi.win(TW, 2.25, 2.85)
    fig, axs = plt.subplots(1, 4, figsize=(13, 4.4), sharey=True)
    for ax, (k, t) in zip(axs, keys):
        A = o[k][w]; c = np.percentile(np.abs(A), 98)
        ax.imshow(A, aspect="auto", cmap="gray_r", vmin=-c, vmax=c, extent=[2099, 2501, TW[w][-1], TW[w][0]], interpolation="bilinear")
        ax.set_title(t, fontsize=9); ax.set_xlabel("xline")
    axs[0].set_ylabel("双程时间 (s)"); fig.tight_layout(); fig.savefig(FIG / "real_sections.png", dpi=150); plt.close(fig)


if __name__ == "__main__":
    accuracy(); compensation(); real()
