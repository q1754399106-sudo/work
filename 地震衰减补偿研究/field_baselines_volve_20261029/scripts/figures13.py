"""Figures of experiment 13 (presentation only). Writes figures/*.png from analysis/summary.json and cache/out."""
import sys, json
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import lib13 as L  # noqa: E402
from lib13 import np, semi, TW  # noqa: E402
sys.path.insert(0, str(L.ROOT.parent / "field_application_20261028" / "cache" / "deps"))   # matplotlib
import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib import font_manager  # noqa: E402

font_manager.fontManager.addfont(r"C:\Windows\Fonts\msyh.ttc")
plt.rcParams.update({"font.family": "Microsoft YaHei", "axes.unicode_minus": False, "font.size": 9,
                     "axes.edgecolor": "#8a8a85", "axes.labelcolor": "#3d3d3a", "xtick.color": "#5f5e5a",
                     "ytick.color": "#5f5e5a", "axes.spines.top": False, "axes.spines.right": False})
ROOT = L.ROOT; FIG = ROOT / "figures"; FIG.mkdir(exist_ok=True)
BLUE, ORANGE, GREEN, YELLOW, PINK, VIOLET, RED = "#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#6250d6", "#c8463c"
GRID = "#e8e7e2"; INK = "#8a8a85"
S = json.loads((ROOT / "analysis" / "summary.json").read_text(encoding="utf-8"))
STYLE = {"B1_IQ70": ("稳定化反 Q 70 dB", GREEN), "B2_SCa": ("SC-a，a=2^-13（SR）", ORANGE),
         "B3_VDB": ("相位反 Q＋Wiener", PINK), "IQ30": ("稳定化反 Q 30 dB", YELLOW),
         "PROP": ("实验 11 改进方案", BLUE), "QSCa_tuned": ("QScan＋SC-a，a=2^-13", VIOLET),
         "QSCa": ("QScan＋SC-a（原流程）", RED)}


def summary():
    fig, axs = plt.subplots(1, 2, figsize=(11.6, 4.4), gridspec_kw=dict(width_ratios=[1.1, 1]))
    ax = axs[0]; ms = [m for m in STYLE if m in S["vs_input"]]; ms.sort(key=lambda m: S["vs_input"][m]["mean"])
    for y, m in enumerate(ms):
        s = S["vs_input"][m]; col = STYLE[m][1]
        ax.plot(s["ci90"], [y, y], color=col, lw=2, alpha=.55); ax.plot(s["mean"], y, "o", color=col, mec="white", ms=7)
        ax.plot(s["values"], [y + .22] * len(s["values"]), "|", color=col, ms=6, alpha=.7)
    ax.set_yticks(range(len(ms))); ax.set_yticklabels([STYLE[m][0] for m in ms]); ax.axvline(0, color=INK, lw=.8)
    ax.set_xlabel("井震相关 T2 相对不补偿的提高（8 组平均，90% 区间；短竖线为各组）")
    ax.set_title("实测 8 组井", loc="left", fontsize=10); ax.grid(axis="x", color=GRID, lw=.6)
    ax = axs[1]; bk = ["10-25", "25-40", "40-55", "55-70", "70-85"]; xs = np.arange(len(bk))
    for m in ["Input"] + list(STYLE):
        if m not in S["bands_all"]:
            continue
        lab, col = ("不补偿", INK) if m == "Input" else STYLE[m]
        ax.plot(xs, [S["bands_all"][m][b] for b in bk], marker="o", ms=4, color=col, lw=1.2, label=lab)
    ax.set_yscale("log"); ax.axhline(1, color=INK, lw=.7, ls=":"); ax.set_xticks(xs); ax.set_xticklabels([b + " Hz" for b in bk])
    ax.set_ylabel("相干信噪比（2.0–3.0 s，105 条 inline 中位数）"); ax.set_title("各频带信噪比", loc="left", fontsize=10)
    ax.legend(frameon=False, fontsize=7.5); ax.grid(color=GRID, lw=.6)
    fig.tight_layout(); fig.savefig(FIG / "summary.png", dpi=170); plt.close(fig)


def sections():
    """Inline 10141 around the reservoir."""
    z = {**dict(np.load(ROOT / "cache" / "out" / "IL_10141_A.npz")), **dict(np.load(ROOT / "cache" / "out" / "IL_10141_B.npz"))}
    w = semi.win(TW, 2.2, 2.9); ks = [("Input", "不补偿"), ("QSCa", "QScan＋SC-a（原流程）"), ("QSCa_tuned", "QScan＋SC-a，a=2^-13"),
                                     ("B2_SCa", "SC-a，a=2^-13（SR）"), ("B1_IQ70", "反 Q 70 dB"), ("B3_VDB", "相位反 Q＋Wiener")]
    fig, axs = plt.subplots(1, len(ks), figsize=(12.4, 4.6), sharey=True)
    for ax, (k, lab) in zip(axs, ks):
        A = z[k][w].astype(float); v = np.percentile(np.abs(A), 99)
        ax.imshow(A, aspect="auto", cmap="gray_r", vmin=-v, vmax=v, extent=[L.XL0, L.XL1, TW[w][-1], TW[w][0]])
        ax.set_title(lab, loc="left", fontsize=8.5); ax.set_xlabel("xline")
    axs[0].set_ylabel("双程时 (s)"); fig.suptitle("实测 inline 10141", x=0.01, ha="left", fontsize=10)
    fig.tight_layout(); fig.savefig(FIG / "sections.png", dpi=170); plt.close(fig)


def ties():
    """Well-side traces of the best-covered well against the log synthetic."""
    R = json.loads((ROOT / "results" / "tie.json").read_text(encoding="utf-8"))
    n = max(R["wells"], key=lambda k: R["wells"][k]["span"][1] - R["wells"][k]["span"][0]); z = np.load(ROOT / "cache" / "out" / f"tie_{n}.npz")
    t0, t1 = R["wells"][n]["span"]; w = semi.win(TW, t0, t1); ks = ["Input", "QSCa", "QSCa_tuned", "B2_SCa", "B1_IQ70"]
    fig, ax = plt.subplots(figsize=(11.6, 3.6)); r = z["ref"][w] / np.std(z["ref"][w])
    for i, k in enumerate(ks):
        tr = z[k][w] / np.std(z[k][w]); off = -4 * (i + 1)
        ax.plot(TW[w], r + off, color=INK, lw=.8); ax.plot(TW[w], tr + off, color=STYLE.get(k, ("", "#3d3d3a"))[1], lw=1.1)
        ax.text(TW[w][0] - .005, off, f"{STYLE.get(k, ('不补偿',))[0]}  T2 {R['wells'][n]['scores'][k]['T2']:.2f}",
                ha="right", va="center", fontsize=8)
    ax.set_yticks([]); ax.set_xlabel("双程时 (s)"); ax.spines["left"].set_visible(False)
    ax.set_title(f"{n}：井旁道（彩色）与测井合成记录（灰）", loc="left", fontsize=10)
    fig.tight_layout(); fig.savefig(FIG / "ties.png", dpi=170); plt.close(fig)


if __name__ == "__main__":
    summary(); sections(); ties()
