"""Figures of experiment 9 (presentation only). Writes figures/*.png."""
import sys, json
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import wq, tie  # noqa: E402
from wq import np, semi, TW  # noqa: E402
sys.path.insert(0, str(wq.ROOT.parent / "field_application_20261028" / "cache" / "deps"))   # matplotlib
import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib import font_manager  # noqa: E402

font_manager.fontManager.addfont(r"C:\Windows\Fonts\msyh.ttc")
plt.rcParams.update({"font.family": "Microsoft YaHei", "axes.unicode_minus": False, "font.size": 9,
                     "axes.edgecolor": "#8a8a85", "axes.labelcolor": "#3d3d3a", "xtick.color": "#5f5e5a",
                     "ytick.color": "#5f5e5a", "axes.spines.top": False, "axes.spines.right": False})
ROOT = wq.ROOT; FIG = ROOT / "figures"; FIG.mkdir(exist_ok=True)
BLUE, ORANGE, GREEN, YELLOW, PINK, VIOLET = "#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#6250d6"
GRID = "#e8e7e2"; INK = "#8a8a85"
S = json.loads((ROOT / "analysis" / "summary.json").read_text(encoding="utf-8"))
T = json.loads((ROOT / "results" / "tie.json").read_text(encoding="utf-8"))
PS = json.loads((ROOT / "results" / "posthoc_static.json").read_text(encoding="utf-8"))
LBL = lambda n: n.replace("_", " ")  # noqa: E731


def ties():
    fig, (a, b) = plt.subplots(1, 2, figsize=(10.4, 4.4), gridspec_kw=dict(width_ratios=[1, 1.35]))
    s = S["IQ30_SR_minus_Input"]; g = list(s["by_group"]); y = np.arange(len(g))[::-1]
    a.axvspan(s["ci90"][0], s["ci90"][1], color=GRID, lw=0, label="平均值的 90% 区间")
    a.axvline(0, color=INK, lw=.8); a.axvline(s["mean"], color=BLUE, lw=1.2, ls="--", label=f"平均 {s['mean']:+.2f}")
    seen = {k for k, m in T["groups"].items() if any(n in T["seen"] for n in m)}
    for yy, k in zip(y, g):
        a.plot(s["by_group"][k], yy, "o", ms=8, color=BLUE if k not in seen else "white", mec=BLUE, mew=1.5)
    a.plot([], [], "o", color="white", mec=BLUE, mew=1.5, label="含实验 5–7 用过的井")
    a.set_yticks(y); a.set_yticklabels([" / ".join(LBL(n) for n in T["groups"][k]) for k in g])
    a.set_xlabel("T2(IQ30_SR) − T2(Input)"); a.set_title("H1：8 组井全部为正", loc="left", fontsize=10)
    a.legend(frameon=False, fontsize=8, loc="upper left"); a.set_xlim(-0.05, 0.45)
    W = T["wells"]; names = list(W); x = np.arange(len(names))
    series = [("Input", lambda n: W[n]["scores"]["Input"]["T2"], INK, "o"),
              ("IQ30，c = SR", lambda n: W[n]["scores"]["IQ30_SR"]["T2"], BLUE, "o"),
              ("时不变滤波（事后）", lambda n: PS["wells"][n]["STATIC_SR"], GREEN, "D"),
              ("IQ30，c = 0.77", lambda n: W[n]["scores"]["IQ30_0.77"]["T2"], ORANGE, "s")]
    for i, (lab, fn, col, mk) in enumerate(series):
        b.plot(x + (i - 1.5) * 0.15, [fn(n) for n in names], mk, ms=6, color=col, mec="white", mew=.8, label=lab)
    for n, xx in zip(names, x):
        b.plot([xx - .225, xx - .075], [W[n]["scores"]["Input"]["T2"], W[n]["scores"]["IQ30_SR"]["T2"]], color=GRID, lw=1, zorder=0)
    b.set_xticks(x); b.set_xticklabels([LBL(n) for n in names], rotation=35, ha="right"); b.set_ylabel("井上相关 T2")
    b.set_title("各井：补偿都提高相关，但时不变滤波一样好，c 的取值影响很小", loc="left", fontsize=10)
    b.legend(frameon=False, fontsize=8, ncol=2, loc="upper left"); b.set_ylim(0, 0.8); b.grid(axis="y", color=GRID, lw=.6)
    fig.tight_layout(); fig.savefig(FIG / "ties.png", dpi=170); plt.close(fig)


def traces():
    P = wq.picks(); pick = ["F-4", "F-1_A", "F-11_T2"]
    fig, axs = plt.subplots(1, 3, figsize=(10.4, 5.0), sharey=False)
    for ax, n in zip(axs, pick):
        g = tie.tie_log(n, P); z = np.load(ROOT / "cache" / f"tie_{n}.npz"); w = semi.win(TW, g["t0"] - .02, g["t1"] + .02)
        sc = T["wells"][n]["scores"]
        for k, (lab, col, off) in enumerate((("测井合成记录", INK, 0), ("Input", ORANGE, 1), ("IQ30，c = SR", BLUE, 2))):
            key = {"测井合成记录": "ref", "Input": "Input", "IQ30，c = SR": "IQ30_SR"}[lab]; tr = z[key][w]
            sh = 0 if key == "ref" else sc[key]["shift_ms"] / 1000
            tt = TW[w] - sh
            ax.plot(off + 0.42 * tr / np.max(np.abs(tr)), tt, color=col, lw=1)
            ax.fill_betweenx(tt, off, off + 0.42 * np.clip(tr, 0, None) / np.max(np.abs(tr)), color=col, alpha=.35, lw=0)
            if key != "ref":
                ax.text(off, g["t0"] - .038, f"T2 {sc[key]['T2']:.2f}", ha="center", va="center", fontsize=8, color=col)
        ax.set_xticks([0, 1, 2]); ax.set_xticklabels(["测井合成", "Input", "IQ30"]); ax.set_ylim(g["t1"] + .03, g["t0"] - .05)
        ax.set_title(LBL(n), loc="left", fontsize=10); ax.set_ylabel("双程时 (s)" if n == pick[0] else "")
    fig.suptitle("井旁道与测井合成记录（地震道已按最佳时移对齐）", x=0.01, ha="left", fontsize=10)
    fig.tight_layout(); fig.savefig(FIG / "traces.png", dpi=170); plt.close(fig)


if __name__ == "__main__":
    ties(); traces()
