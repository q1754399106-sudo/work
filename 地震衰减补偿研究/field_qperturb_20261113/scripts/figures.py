"""Figures of the field Q-perturbation test (presentation only). Writes figures/*.png from analysis/summary.json, results/tie.json and cache/out.

usage: figures.py [INLINE]
"""
import sys, json
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT.parent / "field_application_20261028" / "cache" / "deps"))   # matplotlib, if not installed
import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib import font_manager  # noqa: E402

font_manager.fontManager.addfont(r"C:\Windows\Fonts\msyh.ttc")
plt.rcParams.update({"font.family": "Microsoft YaHei", "axes.unicode_minus": False, "font.size": 9, "axes.edgecolor": "#8a8a85", "axes.labelcolor": "#3d3d3a",
                     "xtick.color": "#5f5e5a", "ytick.color": "#5f5e5a", "axes.spines.top": False, "axes.spines.right": False})
FIG = ROOT / "figures"; FIG.mkdir(exist_ok=True); GRID = "#e8e7e2"
S = json.loads((ROOT / "analysis" / "summary.json").read_text(encoding="utf-8"))
STYLE = {"UQW": ("UQ-W", "#2a78d6"), "VDB": ("VDB（相位反 Q＋维纳整形）", "#e87ba4"), "IQ70": ("常规反 Q（70 dB）", "#1baf7a")}
MS = ["0.5", "0.6", "0.7", "0.85", "1.0", "1.2", "1.4", "1.7", "2.0"]; T0, DT, XL0, XSTEP = 1.2, 0.004, 2100, 2


def curves():
    fig, axs = plt.subplots(1, 2, figsize=(11.6, 4.3)); x = [float(m) for m in MS]; ax = axs[0]
    for k, (name, c) in STYLE.items():
        ax.plot(x, [S["group_mean_T2"][k][m] for m in MS], "-o", color=c, lw=2, ms=4.5, label=name)
    ax.axhline(0.220, color="#8a8a85", ls="--", lw=1); ax.text(2.0, 0.226, "不补偿 0.220", ha="right", color="#5f5e5a")
    ax.axvline(1.0, color=GRID, lw=1, zorder=0); ax.set_xscale("log"); ax.set_xticks(x); ax.set_xticklabels(MS); ax.minorticks_off()
    ax.set_xlabel("Q 乘子（<1：Q 给小了，过补偿；>1：Q 给大了，欠补偿）"); ax.set_ylabel("井震相关（8 组井平均）"); ax.set_ylim(0.18, 0.56); ax.grid(axis="y", color=GRID)
    ax.legend(frameon=False, loc="lower left"); ax.set_title("井震相关随 Q 乘子的变化", loc="left")
    ax = axs[1]; d = S["diff_by_m_UQW_minus_VDB"]; mean = np.array([d[m]["mean"] for m in MS]); lo = np.array([d[m]["ci90"][0] for m in MS]); hi = np.array([d[m]["ci90"][1] for m in MS])
    for i, m in enumerate(MS):
        ax.plot([x[i]] * 8, d[m]["values"], ".", color="#b9b8b2", ms=5, zorder=1)
    ax.errorbar(x, mean, yerr=[mean - lo, hi - mean], fmt="o", color="#2a78d6", capsize=3, lw=1.6, ms=5, zorder=2); ax.axhline(0, color="#3d3d3a", lw=1)
    ax.set_xscale("log"); ax.set_xticks(x); ax.set_xticklabels(MS); ax.minorticks_off(); ax.set_xlabel("Q 乘子"); ax.set_ylabel("UQ-W 减 VDB 的井震相关"); ax.grid(axis="y", color=GRID)
    ax.set_title("UQ-W 减 VDB：均值与 90% 区间（灰点是各井组）", loc="left"); fig.tight_layout(); fig.savefig(FIG / "fig1_tie_vs_q.png", dpi=170); plt.close(fig)


def sections(il, ta=2.0, tb=2.8):
    z = np.load(ROOT / "cache" / "out" / f"IL_{il}.npz"); nt, nx = z["Input"].shape; t = T0 + DT * np.arange(nt); k = (t >= ta) & (t <= tb); ext = [XL0, XL0 + XSTEP * (nx - 1), tb, ta]
    fig, axs = plt.subplots(3, 4, figsize=(15.5, 10.5), sharex=True, sharey=True)
    for r, m in enumerate(["0.5", "1.0", "2.0"]):
        for c, key in enumerate(["Input", "UQW", "VDB", "IQ70"]):
            a = z["Input" if key == "Input" else f"{key}_{m}"][k]; a = a / (np.sqrt(np.mean(a ** 2)) + 1e-30); ax = axs[r, c]
            ax.imshow(a, aspect="auto", cmap="gray_r", vmin=-2.5, vmax=2.5, extent=ext, interpolation="bilinear")
            ax.set_title(("输入（不补偿）" if key == "Input" else STYLE[key][0].split("（")[0]) + ("" if key == "Input" else f"，Q 乘子 {m}"), loc="left", fontsize=10)
            if c == 0: ax.set_ylabel("时间 (s)")
            if r == 2: ax.set_xlabel("crossline")
    fig.suptitle(f"Volve inline {il}，{ta:.1f}–{tb:.1f} s；每幅按自身均方根归一，色标 ±2.5", x=0.01, ha="left", fontsize=11); fig.tight_layout(rect=[0, 0, 1, 0.975])
    fig.savefig(FIG / f"fig2_sections_IL{il}.png", dpi=150); plt.close(fig)


def spectra(il, ta=2.0, tb=2.8):
    z = np.load(ROOT / "cache" / "out" / f"IL_{il}.npz"); nt = z["Input"].shape[0]; t = T0 + DT * np.arange(nt); k = (t >= ta) & (t <= tb); n = int(k.sum()); f = np.fft.rfftfreq(512, DT)
    sp = lambda a: np.sqrt(np.mean(np.abs(np.fft.rfft(a[k] * np.hanning(n)[:, None], 512, axis=0)) ** 2, axis=1))
    fig, axs = plt.subplots(1, 3, figsize=(13.5, 4.1), sharey=True); p0 = sp(z["Input"]); p0 = p0 / p0.max()
    for ax, (key, (name, c)) in zip(axs, STYLE.items()):
        ax.plot(f, 20 * np.log10(p0 + 1e-6), color="#8a8a85", lw=1.2, label="输入")
        for m, ls in [("0.5", ":"), ("1.0", "-"), ("2.0", "--")]:
            p = sp(z[f"{key}_{m}"]); ax.plot(f, 20 * np.log10(p / p.max() + 1e-6), color=c, ls=ls, lw=1.8, label=f"Q 乘子 {m}")
        ax.set_xlim(0, 125); ax.set_ylim(-50, 3); ax.grid(color=GRID); ax.set_xlabel("频率 (Hz)"); ax.set_title(name, loc="left"); ax.legend(frameon=False, loc="lower left")
    axs[0].set_ylabel("振幅谱 (dB，各自归一)"); fig.suptitle(f"inline {il}，{ta:.1f}–{tb:.1f} s 的平均振幅谱", x=0.01, ha="left", fontsize=11); fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(FIG / f"fig3_spectra_IL{il}.png", dpi=170); plt.close(fig)


if __name__ == "__main__":
    il = int(sys.argv[1]) if len(sys.argv) > 1 else 10105; curves(); sections(il); spectra(il); print("written", sorted(p.name for p in FIG.glob("*.png")))
