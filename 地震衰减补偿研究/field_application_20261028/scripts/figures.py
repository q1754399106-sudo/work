"""Figures of experiment 7 (presentation only; no gate depends on them). Writes figures/*.png."""
import sys, json
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import lib  # noqa: E402
from lib import np, bm  # noqa: E402
from scipy.io import loadmat  # noqa: E402
import run  # noqa: E402
sys.path.insert(0, str(lib.ROOT / "cache" / "deps"))
import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib import font_manager  # noqa: E402

font_manager.fontManager.addfont(r"C:\Windows\Fonts\msyh.ttc")
plt.rcParams.update({"font.family": "Microsoft YaHei", "axes.unicode_minus": False, "font.size": 9})
FIG = lib.ROOT / "figures"; FIG.mkdir(exist_ok=True)
NAMES = {"Input": "输入（未补偿）", "Standard": "相干标准版（原始）", "SC": "SC（原始 Ricker 字典）",
         "QSC": "QScan＋SC（原始）", "SCa": "SC-a（实测适配）", "QSCa": "QScan＋SC-a（实测适配）"}
WELLNAME = {"F11": "15/9-F-11 T2", "F15A": "15/9-F-15A"}


def panel(ax, A, tw, title, path=None, x0=0):
    c = np.percentile(np.abs(A), 98)
    ax.imshow(A, aspect="auto", cmap="gray_r", vmin=-c, vmax=c,
              extent=[x0 - .5, x0 + A.shape[1] - .5, tw[-1] + .002, tw[0] - .002], interpolation="bilinear")
    if path is not None:
        ax.plot(path, tw, color="#d62728", lw=1.2, alpha=.8)
    ax.set_title(title, fontsize=9)


def spectrum(A, tw, t0=2.0, t1=3.0):
    w = (tw >= t0 - 1e-9) & (tw <= t1 + 1e-9)
    F = np.abs(np.fft.rfft(A[w] * lib.taper(w.sum())[:, None], 512, axis=0)).mean(axis=1)
    return np.fft.rfftfreq(512, lib.DT), F


def well_figs(well, p, scores, gates):
    c = loadmat(p["curtain"], simplify_cells=True); r = loadmat(p["refs"], simplify_cells=True)
    D = np.asarray(c["D"], float); tw = bm.vec(c["tw"]); idx = bm.vec(c["idx"], int)
    outs = {"Input": D, "SC": loadmat(p["sc"], simplify_cells=True)["X"],
            "SCa": np.load(run.OUT / f"{well}_SCa.npz")["out"], "QSCa": np.load(run.OUT / f"{well}_QSCa.npz")["out"]}
    wv = json.loads((run.AN / f"{well}_wavelets.json").read_text(encoding="utf-8"))
    sc = {s["method"]: s for s in scores if s["well"] == well}
    # Sections.
    fig, axs = plt.subplots(1, 4, figsize=(13, 6.2), sharey=True)
    for ax, (m, A) in zip(axs, outs.items()):
        s = sc[m]
        panel(ax, np.asarray(A, float), tw, f"{NAMES[m]}\nT1 = {s['T1']:.3f}，高频能量比 {s['hf_ratio_60_100Hz_2to3s']:.3g}",
              path=idx - 1)
        ax.set_xlabel("幕帘道号（25 m）")
    axs[0].set_ylabel("双程时 (s)")
    fig.suptitle(f"{WELLNAME[well]}：沿井幕帘剖面（红线为评价道所在的井轨迹）", fontsize=11)
    fig.tight_layout(); fig.savefig(FIG / f"{well}_sections.png", dpi=150); plt.close(fig)
    # Spectra and SNR.
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 4))
    for m, A in outs.items():
        f, F = spectrum(np.asarray(A, float), tw)
        a1.plot(f, 20 * np.log10(F / F.max() + 1e-12), label=NAMES[m])
    a1.axvline(wv["fmax_target"], color="k", ls="--", lw=.8); a1.set_xlim(0, 125); a1.set_ylim(-60, 3)
    a1.set_xlabel("频率 (Hz)"); a1.set_ylabel("2.0–3.0 s 平均振幅谱 (dB，各自归一)"); a1.legend(fontsize=8)
    a1.set_title("补偿后频谱（虚线：目标窗口信噪比 ≥ 1 的上限）")
    f = np.asarray(wv["f"])
    a2.semilogy(f, np.asarray(wv["snr_shallow"]), label="浅层 1.2–1.6 s")
    a2.semilogy(f, np.asarray(wv["snr_target"]), label="目标 2.0–3.0 s")
    a2.axhline(1, color="k", lw=.8); a2.set_xlim(0, 125); a2.set_ylim(1e-3, 1e3)
    a2.set_xlabel("频率 (Hz)"); a2.set_ylabel("相邻道相干信噪比"); a2.legend(fontsize=8)
    a2.set_title(f"输入数据的信噪比（有效频带上限 {wv['fmax_shallow']:.0f} / {wv['fmax_target']:.0f} Hz）")
    fig.tight_layout(); fig.savefig(FIG / f"{well}_spectra.png", dpi=150); plt.close(fig)
    # Tie traces.
    w1 = bm.vec(r["w1"], bool); ref = np.asarray(r["ref1"], float)
    fig, ax = plt.subplots(figsize=(4.8, 7))
    nz = lambda x: (x - x[w1].mean()) / x[w1].std()  # noqa: E731
    ax.plot(nz(ref)[w1], tw[w1], "k", lw=1.4, label="走廊叠加（参照）")
    for k, m in enumerate(["Input", "SC", "SCa", "QSCa"]):
        x = lib.zpfilter(lib.path_trace(np.asarray(outs[m], float), idx), lib.DT, p["fhi"])
        y = np.interp(tw + sc[m]["T1_shift_ms"] / 1000, tw, x)
        ax.plot(nz(y)[w1] + 4 * (k + 1), tw[w1], lw=1, label=f"{NAMES[m]}（{sc[m]['T1']:.2f}）")
    ax.invert_yaxis(); ax.set_xticks([]); ax.set_ylabel("双程时 (s)"); ax.legend(fontsize=7, loc="lower left")
    ax.set_title(f"{WELLNAME[well]}：井旁道与走廊叠加\n（括号内为相关系数 T1）", fontsize=9)
    fig.tight_layout(); fig.savefig(FIG / f"{well}_tie.png", dpi=150); plt.close(fig)


def inline_fig():
    d = np.load(run.OUT / "IL10141_input.npz"); D, tw, xl = d["D"], d["tw"], d["xline"]
    outs = {"Input": D, "SCa": np.load(run.OUT / "IL10141_SCa.npz")["out"], "QSCa": np.load(run.OUT / "IL10141_QSCa.npz")["out"]}
    info = json.loads((run.AN / "inline.json").read_text(encoding="utf-8"))
    fig, axs = plt.subplots(1, 3, figsize=(15, 6.2), sharey=True)
    for ax, (m, A) in zip(axs, outs.items()):
        panel(ax, A, tw, NAMES[m])
        ax.set_xticks(np.linspace(0, A.shape[1] - 1, 5)); ax.set_xticklabels([str(int(v)) for v in np.linspace(xl[0], xl[-1], 5)])
        ax.set_xlabel("xline")
    axs[0].set_ylabel("双程时 (s)")
    fig.suptitle(f"应用测线 inline 10141（201 道，靠近 F-11 T2 与 F-15A）；QScan 因子 {info['qscan']['factor']:.2f}，"
                 f"有效频带上限 {info['fmax_target']:.0f} Hz", fontsize=11)
    fig.tight_layout(); fig.savefig(FIG / "IL10141_sections.png", dpi=150); plt.close(fig)
    # Zoom on the reservoir interval for detail.
    fig, axs = plt.subplots(1, 2, figsize=(12, 5), sharey=True)
    z = (tw >= 2.3) & (tw <= 2.8)
    for ax, m in zip(axs, ["Input", "QSCa"]):
        panel(ax, outs[m][z][:, 60:160], tw[z], NAMES[m] + "（2.3–2.8 s 放大）", x0=60)
        ax.set_xlabel("道号")
    axs[0].set_ylabel("双程时 (s)")
    fig.tight_layout(); fig.savefig(FIG / "IL10141_zoom.png", dpi=150); plt.close(fig)


if __name__ == "__main__":
    scores = json.loads((run.AN / "scores.json").read_text(encoding="utf-8"))
    gates = json.loads((run.AN / "gates.json").read_text(encoding="utf-8"))
    for well, p in run.WELLS.items():
        well_figs(well, p, scores, gates)
    if (run.OUT / "IL10141_QSCa.npz").exists():
        inline_fig()
    print(sorted(x.name for x in FIG.glob("*.png")))
