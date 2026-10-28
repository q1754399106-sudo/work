"""Figures of experiment 12 (presentation only). Writes figures/*.png from analysis/summary.json."""
import sys, json
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import lib12 as L  # noqa: E402
from lib12 import np, sc11, semi  # noqa: E402
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
BLUE, ORANGE, GREEN, YELLOW, PINK, VIOLET = "#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#6250d6"
GRID = "#e8e7e2"; INK = "#8a8a85"
S = json.loads((ROOT / "analysis" / "summary.json").read_text(encoding="utf-8"))
SEL = json.loads((ROOT / "analysis" / "dev_selection.json").read_text(encoding="utf-8"))
STYLE = {"SRSCa_1.3+ph": ("改进方案（实验 11 冻结）", BLUE, "o"),
         "SCaK": (f"调好 a 的 SC-a＋解耦（k={SEL['SCaK']['k'] if SEL['SCaK'] else '-'}）", VIOLET, "^"),
         "B2_SCa": (f"SC-a，调好 a（{SEL['B2']['a']:g}）", ORANGE, "s"),
         "B1_IQ": (f"稳定化反 Q（{SEL['B1']['g_db']} dB）", GREEN, "D"),
         "IQd": (f"反 Q＋解耦（k={SEL['IQd']['k']}，{SEL['IQd']['g_db']} dB）", YELLOW, "v"),
         "B3_VDB": ("纯相位反 Q＋时变 Wiener", PINK, "P"),
         "SCa_true": ("SC-a 真 Q（参照）", INK, "x")}


def summary():
    ms = [m for m in STYLE if m in S["mean"]]
    fig, axs = plt.subplots(1, 2, figsize=(11.6, 4.4), gridspec_kw=dict(width_ratios=[1.3, 1]))
    a = axs[0]; byc = S["by_c"]; cs = [float(c) for c in byc]
    for m in ms:
        lab, col, mk = STYLE[m]
        a.plot(cs, [byc[c][m] for c in byc], marker=mk, color=col, lw=1.4 if m in ("SRSCa_1.3+ph", "B1_IQ") else .8,
               ms=6, mec="white" if mk != "x" else col, label=lab)
    a.set_xlabel("真实 Q 系数 c"); a.set_ylabel("综合分 S（各 c 4 个数据集均值）")
    a.set_title("各 c 下的综合分", loc="left", fontsize=10); a.grid(axis="y", color=GRID, lw=.6)
    a.legend(frameon=False, fontsize=7.5, loc="lower right")
    ax = axs[1]; pairs = [k for k in S["comparisons"]]; y = np.arange(len(pairs))[::-1]
    ax.axvline(0, color=INK, lw=.8)
    for yy, k in zip(y, pairs):
        d = S["comparisons"][k]["S"]; col = STYLE.get(k.split(" - ")[0], ("", INK))[1]
        ax.plot([d["lo"], d["hi"]], [yy, yy], color=col, lw=1.8, alpha=.6)
        ax.plot(d["mean"], yy, "o", ms=6, color=col, mec="white")
    lab = {m: STYLE[m][0].split("（")[0] for m in STYLE}; lab["SRSCa"] = "SC-a（a=0.125）"; lab["SRSCa_1.3"] = "改进方案去相位校正"
    ax.set_yticks(y); ax.set_yticklabels([" − ".join(lab.get(p, p) for p in k.split(" - ")) for k in pairs], fontsize=7.5)
    ax.set_xlabel("综合分 S 的配对差（90% 区间）"); ax.set_title("配对比较", loc="left", fontsize=10)
    ax.grid(axis="x", color=GRID, lw=.6)
    fig.tight_layout(); fig.savefig(FIG / "summary.png", dpi=170); plt.close(fig)


def sections():
    """One held-out dataset around the evaluation window."""
    did = "T_c1.05_s9701"; z = np.load(ROOT / "cache" / "out" / f"{did}.npz"); w = semi.win(sc11.TW, 2.2, 2.9)
    r = json.loads((ROOT / "results" / "synth" / f"{did}.json").read_text(encoding="utf-8"))
    ks = [("ideal", "理想目标"), ("SRSCa_1.3+ph", "改进方案"), ("SCaK", "调好 a＋解耦"), ("B1_IQ", "稳定化反 Q"),
          ("B3_VDB", "相位反 Q＋Wiener")]
    ks = [k for k in ks if k[0] in z.files]
    fig, axs = plt.subplots(1, len(ks), figsize=(11.6, 4.4), sharey=True)
    for ax, (k, lab) in zip(axs, ks):
        A = z[k][w][:, 60:141].astype(float); v = np.percentile(np.abs(A), 99)
        ax.imshow(A, aspect="auto", cmap="gray_r", vmin=-v, vmax=v, extent=[60, 140, sc11.TW[w][-1], sc11.TW[w][0]])
        ax.set_title(lab if k == "ideal" else f"{lab}  S {r['scores'][k]['S']:.2f}", loc="left", fontsize=9)
        ax.set_xlabel("道号")
    axs[0].set_ylabel("双程时 (s)")
    fig.suptitle(f"检验数据 {did}（真 c = {r['c_true']}，SR 估 {r['sr']['c']:.2f}）", x=0.01, ha="left", fontsize=10)
    fig.tight_layout(); fig.savefig(FIG / "sections.png", dpi=170); plt.close(fig)


if __name__ == "__main__":
    summary(); sections()
