"""Figures of experiment 8 (presentation only; no hypothesis depends on them). Writes figures/*.png."""
import sys, json
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import semi  # noqa: E402
from semi import np, lib, TW, C_TRUE  # noqa: E402
import analyze  # noqa: E402
sys.path.insert(0, str(semi.ROOT.parent / "field_application_20261028" / "cache" / "deps"))   # matplotlib
import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib import font_manager  # noqa: E402

font_manager.fontManager.addfont(r"C:\Windows\Fonts\msyh.ttc")
plt.rcParams.update({"font.family": "Microsoft YaHei", "axes.unicode_minus": False, "font.size": 9,
                     "axes.edgecolor": "#8a8a85", "axes.labelcolor": "#3d3d3a", "xtick.color": "#5f5e5a",
                     "ytick.color": "#5f5e5a", "axes.spines.top": False, "axes.spines.right": False})
ROOT = semi.ROOT; FIG = ROOT / "figures"; FIG.mkdir(exist_ok=True)
B = np.load(ROOT / "cache" / "build.npz")
NAMES = {"Input": "输入（未补偿）", "IQ20a": "常规反 Q（20 dB）＋限带", "QIQa": "QScan＋常规反 Q＋限带",
         "SC": "SC（原始 Ricker）", "SCa": "SC-a", "QSCa": "QScan＋SC-a", "SCa_oracle": "SC-a（真 Q）",
         "IQ10": "常规反 Q 10 dB", "IQ20": "常规反 Q 20 dB", "IQ30": "常规反 Q 30 dB"}
SHOW = ["Input", "IQ20a", "QIQa", "IQ30", "SCa", "QSCa"]
COL = dict(zip(SHOW, ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#6250d6"]))


def panel(ax, A, title, clip=None, t=TW, x0=0):
    c = clip if clip is not None else np.percentile(np.abs(A), 98)
    ax.imshow(A, aspect="auto", cmap="gray_r", vmin=-c, vmax=c, extent=[x0 - .5, x0 + A.shape[1] - .5, t[-1], t[0]],
              interpolation="bilinear")
    ax.set_title(title, fontsize=9)


def realism():
    lv = json.loads((ROOT / "analysis" / "build.json").read_text(encoding="utf-8"))["level"]
    fig = plt.figure(figsize=(11, 7.2)); gs = fig.add_gridspec(2, 3, height_ratios=[1.3, 1])
    D = B["D_real"]; X = B["sig_c0.8"] + lv["c0.8"] * B["noise_R0"]
    panel(fig.add_subplot(gs[0, 0]), D, "实测 inline 10141"); panel(fig.add_subplot(gs[0, 1]), X, "半合成（c = 0.8，实测噪声 R0）")
    ax = fig.add_subplot(gs[0, 2]); panel(ax, B["noise_R0"] * lv["c0.8"], "其中的噪声（另一条测线）",
                                          clip=np.percentile(np.abs(X), 98))
    ax = fig.add_subplot(gs[1, :2])
    f, s, _ = lib.coherence_snr(D, TW, 2.0, 3.0); ax.semilogy(f, np.maximum(s, 1e-2), color="#0b0b0b", lw=2, label="实测")
    for c, col in zip(C_TRUE, ["#2a78d6", "#eb6834", "#1baf7a", "#6250d6"]):
        _, s, _ = lib.coherence_snr(B[f"sig_c{c}"] + lv[f"c{c}"] * B["noise_R0"], TW, 2.0, 3.0)
        ax.semilogy(f, np.maximum(s, 1e-2), color=col, lw=1.4, label=f"半合成 c = {c}")
    ax.axhline(1, color="#8a8a85", lw=.8, ls="--"); ax.set_xlim(0, 110); ax.set_ylim(1e-2, 300)
    ax.set_xlabel("频率 (Hz)"); ax.set_ylabel("相干信噪比（2.0–3.0 s）"); ax.legend(frameon=False, ncol=5, fontsize=8)
    ax = fig.add_subplot(gs[1, 2])
    fa, _, Ps = lib.coherence_snr(D, TW, 1.2, 1.6); _, _, Pd = lib.coherence_snr(D, TW, 2.0, 3.0); b = (fa >= 10) & (fa <= 60)
    ax.plot(fa[b], 0.5 * np.log(Pd[b] / Ps[b]), color="#0b0b0b", lw=2, label="实测")
    for c, col in zip(C_TRUE, ["#2a78d6", "#eb6834", "#1baf7a", "#6250d6"]):
        X = B[f"sig_c{c}"] + lv[f"c{c}"] * B["noise_R0"]
        _, _, Ps = lib.coherence_snr(X, TW, 1.2, 1.6); _, _, Pd = lib.coherence_snr(X, TW, 2.0, 3.0)
        ax.plot(fa[b], 0.5 * np.log(Pd[b] / Ps[b]), color=col, lw=1.4, label=f"c = {c}")
    ax.set_xlabel("频率 (Hz)"); ax.set_ylabel("ln(目标层 / 浅层 振幅谱)"); ax.set_title("浅层到目标层的频谱衰减", fontsize=9)
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout(); fig.savefig(FIG / "realism.png", dpi=150); plt.close(fig)


def sections(did):
    o = np.load(ROOT / "cache" / "out" / f"{did}.npz"); keys = ["Input", "IQ20a", "SCa", "QSCa", "ideal"]
    fig, axs = plt.subplots(1, 5, figsize=(14, 5.2), sharey=True)
    for ax, k in zip(axs, keys):
        A = o[k]; panel(ax, A, NAMES.get(k, "理想目标（真反射，真信噪比频带）"))
    axs[0].set_ylabel("双程时间 (s)")
    for ax in axs:
        ax.axvline(semi.CENTRE, color="#e34948", lw=.6, alpha=.5)
    fig.suptitle(f"{did}：中心红线为“井”位置", fontsize=10); fig.tight_layout(); fig.savefig(FIG / f"sections_{did}.png", dpi=150)
    plt.close(fig)
    fig, axs = plt.subplots(1, 5, figsize=(14, 4.2), sharey=True); w = semi.win(TW, 2.25, 2.85)
    for ax, k in zip(axs, keys):
        panel(ax, o[k][w][:, 70:131], NAMES.get(k, "理想目标"), t=TW[w], x0=70)
        ax.axvline(semi.CENTRE, color="#e34948", lw=.6, alpha=.5)
    axs[0].set_ylabel("双程时间 (s)")
    fig.tight_layout(); fig.savefig(FIG / f"zoom_{did}.png", dpi=150); plt.close(fig)


def lag(rows):
    """Post-hoc (descriptive): zero-lag M1 vs M1 with one common shift (+-20 ms), analysis/posthoc_lag.json."""
    L = json.loads((ROOT / "analysis" / "posthoc_lag.json").read_text(encoding="utf-8"))["summary"]
    ms = ["Input", "IQ30", "QSCa", "SCa_oracle"]; col = ["#2a78d6", "#eda100", "#6250d6", "#1baf7a"]
    fig, axs = plt.subplots(1, len(C_TRUE), figsize=(13, 3.8), sharey=True)
    for ax, c in zip(axs, C_TRUE):
        rs = analyze.group(rows, "P", c)
        for i, (m, cl) in enumerate(zip(ms, col)):
            z = np.mean([r["scores"][m]["M1"] for r in rs]); l = L[f"c{c}"][m]["M1_lag"][0]
            ax.plot([i - .15, i + .15], [z, l], color=cl, lw=1.5)
            ax.plot(i - .15, z, "o", ms=7, color="white", mec=cl, mew=2); ax.plot(i + .15, l, "o", ms=7, color=cl, mec="white", mew=1.5)
        ax.set_xticks(range(len(ms))); ax.set_xticklabels([NAMES[m] for m in ms], rotation=30, ha="right", fontsize=8)
        ax.axhline(0, color="#8a8a85", lw=.8); ax.grid(axis="y", color="#e8e7e2", lw=.6)
        ax.set_title(f"真 Q = {c} × 模型 Q", fontsize=9)
    axs[0].set_ylabel("M1（空心：零时移；实心：允许 ±20 ms 整体时移）")
    fig.tight_layout(); fig.savefig(FIG / "posthoc_lag.png", dpi=150); plt.close(fig)


def summary(rows):
    fig, axs = plt.subplots(1, len(C_TRUE), figsize=(13, 4.2), sharey=True)
    for ax, c in zip(axs, C_TRUE):
        rs = analyze.group(rows, "P", c)
        for i, m in enumerate(SHOW):
            v = np.array([r["scores"][m]["M1"] for r in rs]); mm, lo, hi = analyze.ci(v)
            ax.errorbar(i, mm, yerr=[[mm - lo], [hi - mm]], fmt="o", ms=7, color=COL[m], capsize=3, lw=2,
                        mec="white", mew=1.5)
        ax.set_xticks(range(len(SHOW))); ax.set_xticklabels([NAMES[m] for m in SHOW], rotation=40, ha="right", fontsize=8)
        ax.set_title(f"真 Q = {c} × 模型 Q（{len(rs)} 个噪声实现）", fontsize=9); ax.grid(axis="y", color="#e8e7e2", lw=.6)
    axs[0].set_ylabel("M1：与理想目标的相关（2.30–2.80 s）")
    fig.tight_layout(); fig.savefig(FIG / "summary_M1.png", dpi=150); plt.close(fig)
    fig, axs = plt.subplots(1, len(C_TRUE), figsize=(13, 3.8), sharey=True)
    lab = [f"{a}–{b}" for a, b in semi.BANDS]
    for ax, c in zip(axs, C_TRUE):
        rs = analyze.group(rows, "P", c)
        for m in ["Input", "IQ20a", "QIQa", "QSCa"]:
            ax.plot(range(len(lab)), [np.mean([r["scores"][m][f"M2_{a}_{b}"] for r in rs]) for a, b in semi.BANDS],
                    "o-", color=COL[m], lw=2, ms=6, mec="white", label=NAMES[m])
        ax.set_xticks(range(len(lab))); ax.set_xticklabels(lab); ax.set_xlabel("频带 (Hz)")
        ax.set_title(f"真 Q = {c} × 模型 Q", fontsize=9); ax.grid(axis="y", color="#e8e7e2", lw=.6)
    axs[0].set_ylabel("M2：分频带与真响应的相关"); axs[0].legend(frameon=False, fontsize=8)
    fig.tight_layout(); fig.savefig(FIG / "summary_bands.png", dpi=150); plt.close(fig)
    fig, axs = plt.subplots(1, 2, figsize=(10, 4))
    ax = axs[0]
    for c in C_TRUE:
        fac = [r["info"]["qscan"]["factor"] for r in analyze.group(rows, "P", c)]
        ax.scatter(np.full(len(fac), c) * np.exp(np.random.default_rng(0).uniform(-.03, .03, len(fac))), fac,
                   s=40, color="#6250d6", edgecolor="white", lw=1, zorder=3)
    g = np.array([0.5, 2.0]); ax.plot(g, g, color="#8a8a85", lw=1); ax.fill_between(g, g / 1.25, g * 1.25, color="#e8e7e2", lw=0)
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlabel("真 c"); ax.set_ylabel("QScan 选出的系数")
    ax.set_title("QScan 找回 Q（灰带为 ±25%）", fontsize=9)
    ax = axs[1]
    for c, col in zip(C_TRUE, ["#2a78d6", "#eb6834", "#1baf7a", "#6250d6"]):
        rs = analyze.group(rows, "P", c)
        ax.scatter(analyze.diff(rs, "QSCa", "Input"), analyze.diff(rs, "QSCa", "Input", "T2like"), s=36, color=col,
                   edgecolor="white", lw=1, label=f"c = {c}", zorder=3)
    ax.axhline(0, color="#8a8a85", lw=.8); ax.axvline(0, color="#8a8a85", lw=.8)
    ax.set_xlabel("ΔM1（真值评价，QScan＋SC-a − 输入）"); ax.set_ylabel("ΔT2-like（井上相关评价）")
    ax.set_title("井上相关能否看出真实的改进", fontsize=9); ax.legend(frameon=False, fontsize=8)
    fig.tight_layout(); fig.savefig(FIG / "qscan_and_tie.png", dpi=150); plt.close(fig)


if __name__ == "__main__":
    realism()
    rows = analyze.load()
    if rows:
        summary(rows)
        lag(rows)
        for c in C_TRUE:
            if (ROOT / "cache" / "out" / f"P_c{c}_n0.npz").exists():
                sections(f"P_c{c}_n0")
