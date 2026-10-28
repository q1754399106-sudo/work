"""Figures of experiment 11 (presentation only). Writes figures/*.png."""
import sys, json
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import sc11  # noqa: E402
from sc11 import np, TW, semi  # noqa: E402
sys.path.insert(0, str(sc11.ROOT.parent / "field_application_20261028" / "cache" / "deps"))   # matplotlib
import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib import font_manager  # noqa: E402

font_manager.fontManager.addfont(r"C:\Windows\Fonts\msyh.ttc")
plt.rcParams.update({"font.family": "Microsoft YaHei", "axes.unicode_minus": False, "font.size": 9,
                     "axes.edgecolor": "#8a8a85", "axes.labelcolor": "#3d3d3a", "xtick.color": "#5f5e5a",
                     "ytick.color": "#5f5e5a", "axes.spines.top": False, "axes.spines.right": False})
ROOT = sc11.ROOT; FIG = ROOT / "figures"; FIG.mkdir(exist_ok=True)
BLUE, ORANGE, GREEN, YELLOW, PINK, VIOLET = "#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#6250d6"
GRID = "#e8e7e2"; INK = "#8a8a85"
S = json.loads((ROOT / "analysis" / "summary.json").read_text(encoding="utf-8"))
STYLE = {"QSCa": ("QScan＋SC-a（现流程）", ORANGE, "s"), "SRSCa_1.3+ph": ("改进：SR 相位＋保守振幅", BLUE, "o"),
         "SRSCa": ("SR＋SC-a", VIOLET, "^"), "QSCa+ph": ("QScan＋SC-a＋SR 相位", GREEN, "D"),
         "IQ30_SR": ("30 dB 反 Q（SR）", YELLOW, "v"), "SCa_true": ("SC-a 真 Q（参照）", INK, "x")}


def summary():
    fig, axs = plt.subplots(1, 3, figsize=(11.6, 4.2), gridspec_kw=dict(width_ratios=[1.5, 1, 1]))
    a = axs[0]; byc = S["means"]["T_by_c"]; cs = [float(c) for c in byc]
    for m, (lab, col, mk) in STYLE.items():
        a.plot(cs, [byc[str(c) if str(c) in byc else c][m] for c in byc], marker=mk, color=col, lw=1.4 if m in ("QSCa", "SRSCa_1.3+ph") else .8,
               ms=6, mec="white" if mk != "x" else col, label=lab)
    a.set_xlabel("真实 Q 系数 c"); a.set_ylabel("M1（零时移，T 组均值）"); a.set_title("各 c 下的零时移相关", loc="left", fontsize=10)
    a.grid(axis="y", color=GRID, lw=.6); a.legend(frameon=False, fontsize=7.5, loc="lower right")
    for ax, g, ttl in ((axs[1], "T", "T 组（20 个）"), (axs[2], "L", "L 组：地层颜色压力（8 个）")):
        ms = [m for m in STYLE if m != "QSCa"]; y = np.arange(len(ms))[::-1]
        ax.axvline(0, color=INK, lw=.8)
        for yy, m in zip(y, ms):
            for k, dy, mk in (("M1", .12, "o"), ("M1_lag", -.12, "D")):
                d = S["vs_QSCa"][g][m][k]; col = STYLE[m][1]
                ax.plot(d["ci90"], [yy + dy] * 2, color=col, lw=1.6, alpha=.6)
                ax.plot(d["mean"], yy + dy, mk, ms=6, color=col if k == "M1" else "white", mec=col, mew=1.3)
        ax.set_yticks(y); ax.set_yticklabels([STYLE[m][0] for m in ms] if g == "T" else [""] * len(ms))
        ax.set_xlabel("相对 QScan＋SC-a 的差（90% 区间）"); ax.set_title(ttl, loc="left", fontsize=10)
        ax.grid(axis="x", color=GRID, lw=.6)
    axs[2].plot([], [], "o", color=INK, label="M1 零时移"); axs[2].plot([], [], "D", color="white", mec=INK, label="M1 允许时移")
    axs[2].legend(frameon=False, fontsize=7.5, loc="lower right")
    fig.tight_layout(); fig.savefig(FIG / "summary.png", dpi=170); plt.close(fig)


def sections():
    """One T dataset: ideal, QSCa, improved, around the evaluation window."""
    did = "T_c0.75_s9501"; z = np.load(ROOT / "cache" / "out" / f"{did}.npz"); w = semi.win(TW, 2.2, 2.9)
    r = json.loads((ROOT / "results" / "synth" / f"{did}.json").read_text(encoding="utf-8"))
    fig, axs = plt.subplots(1, 4, figsize=(11.6, 4.4), sharey=True)
    for ax, (k, lab) in zip(axs, (("ideal", "理想目标"), ("Input", "输入"), ("QSCa", "QScan＋SC-a"), ("SRSCa_1.3+ph", "改进方案"))):
        A = z[k][w][:, 60:141].astype(float); v = np.percentile(np.abs(A), 99)
        ax.imshow(A, aspect="auto", cmap="gray_r", vmin=-v, vmax=v, extent=[60, 140, TW[w][-1], TW[w][0]])
        t = lab if k == "ideal" else f"{lab}  M1 {r['scores'][k]['M1']:.2f}"
        ax.set_title(t, loc="left", fontsize=9.5); ax.set_xlabel("道号")
    axs[0].set_ylabel("双程时 (s)")
    fig.suptitle(f"检验数据 {did}（真 c = 0.75，SR 估 {r['sr']['c']:.2f}，QScan 选 {r['qscan']:.2f}）", x=0.01, ha="left", fontsize=10)
    fig.tight_layout(); fig.savefig(FIG / "sections.png", dpi=170); plt.close(fig)


def real():
    R = json.loads((ROOT / "results" / "real" / "inline10141.json").read_text(encoding="utf-8")); f = np.array(R["f"])
    z = np.load(ROOT / "cache" / "real_10141.npz"); w = semi.win(TW, 2.2, 2.9)
    fig = plt.figure(figsize=(11.6, 4.6)); gs = fig.add_gridspec(1, 4, width_ratios=[1, 1, 1, 1.3])
    for i, (k, lab) in enumerate((("Input", "输入"), ("QSCa", "QScan＋SC-a"), ("SRSCa_1.3+ph", "改进方案"))):
        ax = fig.add_subplot(gs[i]); A = z[k][w].astype(float); v = np.percentile(np.abs(A), 99)
        ax.imshow(A, aspect="auto", cmap="gray_r", vmin=-v, vmax=v, extent=[0, A.shape[1], TW[w][-1], TW[w][0]])
        ax.set_title(lab, loc="left", fontsize=9.5); ax.set_xlabel("道号"); ax.set_ylabel("双程时 (s)" if i == 0 else "")
    ax = fig.add_subplot(gs[3])
    for k in ("Input", "QSCa", "SRSCa_1.3+ph", "IQ30_SR"):
        lab, col = (("输入", INK) if k == "Input" else STYLE[k][:2])
        ax.plot(f, R["spectra"][k]["snr"], color=col, lw=1.3, label=lab)
    ax.axhline(1, color=INK, lw=.7, ls=":"); ax.set_yscale("log"); ax.set_xlim(0, 100); ax.set_ylim(.05, 50)
    ax.set_xlabel("频率 (Hz)"); ax.set_ylabel("相干信噪比（2.0–3.0 s）"); ax.legend(frameon=False, fontsize=7.5)
    ax.set_title("实测 inline 10141", loc="left", fontsize=9.5); ax.grid(color=GRID, lw=.6)
    fig.tight_layout(); fig.savefig(FIG / "real.png", dpi=170); plt.close(fig)


if __name__ == "__main__":
    summary(); sections(); real()
