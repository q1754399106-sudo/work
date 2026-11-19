"""Paper figures 1-8 (PAPER_DRAFT_ZH.md). Okabe-Ito colours with redundant line styles / markers, PDF (Type 42 fonts) +
PNG 300 dpi at 174 mm (two-column) or 84 mm width. Reads only saved result files; no data are altered.

usage: make_figs.py [zh|en]
  zh (default): Chinese labels for Oil Geophysical Prospecting -> figures/fig*.pdf|png, figures/MANIFEST.json.
                Fonts: Latin Times New Roman (else Liberation Serif), Chinese SimSun / Songti (else Noto Serif CJK SC),
                maths in the Latin font (STIX only for missing symbols). The first installed font of each list is used, so the same script runs on Windows and Linux.
  en          : the earlier English labels -> figures/en/ (the 2026-11-19 English set was made with the
                scientific-visualization style preset; without it a close fallback style is used).
"""
import sys, re, json, hashlib, logging
from pathlib import Path
import numpy as np
_VENDOR = Path(r"D:/ChatGPT/研究/地震衰减补偿研究/attenuation_run_20261019/dependencies")   # vendored matplotlib 3.11.2 (cp312)
if _VENDOR.exists():
    sys.path.append(str(_VENDOR))
_PRESET = Path(r"C:/Users/qq175/.claude/skills/scientific-visualization/scripts")
if _PRESET.exists():
    sys.path.insert(0, str(_PRESET))
import matplotlib; matplotlib.use("Agg")
import matplotlib.ticker
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.patches import FancyBboxPatch

LANG = sys.argv[1] if len(sys.argv) > 1 else "zh"
assert LANG in ("zh", "en"), LANG
HERE = Path(__file__).resolve().parent; ROOT = HERE.parent; UQ = ROOT.parent / "uq_avo_20261118"
OUT = ROOT / "figures" / ("" if LANG == "zh" else "en"); OUT.mkdir(parents=True, exist_ok=True)
BLUE, VERM, GREEN, PINK, BLACK = "#0072B2", "#D55E00", "#009E73", "#CC79A7", "#000000"
MM = 1 / 25.4; W2, W1 = 174 * MM, 84 * MM; manifest = {}
J = lambda p: json.loads(Path(p).read_text(encoding="utf-8"))

# ------------------------------------------------------------------ labels
EN = dict(
    f1=["Pre-stack gathers\n(no NMO)", "Offset spectral ratio,\nequal-length windows\nalong moveout",
        "$\\hat c=Q_\\mathrm{true}/Q_\\mathrm{model}$\n(50-CMP\nsupergather)", "Path-aware\ncompensation\nwith $\\hat c\\,Q_\\mathrm{model}$",
        "Per-event AVO\n(intercept,\ngradient)"],
    f2_tik="Tikhonov ($\\sigma=0$)", f2_uq="UQ ($\\sigma=0.447$)", f2_x="Numerical $\\Delta(G/A)$ (time-varying operator)",
    f2_y="Closed-form $\\Delta(G/A)$ from $\\Gamma(x)$",
    f3_tik="Tikhonov", f3_oc="Offset-consistent", f3_x="$\\sigma_{\\ln Q}$ (prior width of the Q-level error)", f3_t="{} Hz Ricker",
    f3_y="68th percentile of $|\\Delta(G/A)|$ over the prior\n(median of 3 events)",
    f4_t="{} Hz Ricker", f4_y="Extra $\\Delta(G/A)$ from UQ, $D$\n(UQ minus Tikhonov, Q correct)",
    f5_orig="Compressed far window (original)", f5_eq="Equal-length far window", f5_x0="Ricker peak frequency (Hz)",
    f5_y0="$\\ln(\\hat c/c)$ (noise-free; mean, range of 4 seeds)", f5_x1="Offset (m)",
    f5_y1="Measured / predicted $\\Delta t^*$ per window\n(40 Hz, median)",
    f6_orig="Original\n(M16-M19)", f6_corr="Corrected\n(M20-M23)", f6_y0="$\\ln(\\hat c/c)$, 24 dB, AVO, N = 50 CMPs\n(384 estimates per box)",
    f6_x1="CMPs per supergather, N (log scale)", f6_y1="RMS $\\ln(\\hat c/c)$ (corrected estimator)", f6_leg="{} Hz, {} dB",
    f7=["True Q", "Model Q (sat. damping)", "Model Q (Tikhonov)", "Calibrated $\\hat cQ_\\mathrm{model}$"],
    f7_x="$c=Q_\\mathrm{true}/Q_\\mathrm{model}$ (log scale)", f7_y0="Per-event $|\\Delta(G/A)|$, noisy (log scale)\nmedian of 32 cells",
    f7_y1="Waveform error (log scale)\nmedian of 32 cells",
    f8=["Truth", "True Q", "Model Q", "Calibrated"], f8_t="c = {}, event {} ($\\hat c$ = {:.2f})", f8_y="Amplitude / zero-offset amplitude")
ZH = dict(
    f1=["未动校\n叠前 CMP 道集", "沿时距曲线取窗\n（远道窗与近道等长）\n计算炮检距谱比",
        "估计校正系数 $\\hat c$\n（50 个 CMP\n超道集）", "用 $\\hat c$ 校准后的 Q\n沿射线路径补偿",
        "逐事件 AVO\n（截距、梯度）"],
    f2_tik="一步 Tikhonov（$\\sigma=0$）", f2_uq="期望算子（$\\sigma=0.447$）", f2_x="数值计算的 $\\Delta(G/A)$（时变算子）",
    f2_y="由 $\\Gamma(x)$ 闭式预测的 $\\Delta(G/A)$",
    f3_tik="一步 Tikhonov", f3_oc="炮检距一致算子", f3_x="$\\sigma_{\\ln Q}$（Q 绝对水平误差的先验宽度）", f3_t="{} Hz Ricker 子波",
    f3_y="先验下 $|\\Delta(G/A)|$ 的 68% 分位数\n（3 个事件的中位数）",
    f4_t="{} Hz Ricker 子波", f4_y="期望算子额外产生的 $\\Delta(G/A)$，$D$\n（Q 正确；期望算子减一步 Tikhonov）",
    f5_orig="随时差压缩的远道窗（原估计器）", f5_eq="远道窗与近道等长（修正后）", f5_x0="Ricker 子波主频/Hz",
    f5_y0="$\\ln(\\hat c/c)$\n（无噪声；4 个种子的均值与范围）", f5_x1="炮检距/m",
    f5_y1="逐窗实测与预测 $\\Delta t^*$ 之比\n（40 Hz 子波，中位数）",
    f6_orig="原估计器\n（M16～M19）", f6_corr="修正估计器\n（M20～M23）", f6_y0="$\\ln(\\hat c/c)$（24 dB，有 AVO，N = 50）\n每箱 384 个估计",
    f6_x1="超道集 CMP 数 N（对数坐标）", f6_y1="$\\ln(\\hat c/c)$ 均方根（修正估计器）", f6_leg="{} Hz，{} dB",
    f7=["真值 Q", "模型 Q（频率饱和阻尼）", "模型 Q（一步 Tikhonov）", "用 $\\hat c$ 校准后的 Q"],
    f7_x="$c$（真实 Q 与模型 Q 之比，对数坐标）", f7_y0="有噪逐事件 $|\\Delta(G/A)|$（对数坐标）\n32 个单元的中位数",
    f7_y1="波形相对误差（对数坐标）\n32 个单元的中位数",
    f8=["真值", "真值 Q", "模型 Q", "校准后"], f8_t="$c$ = {}，事件 {}（$\\hat c$ = {:.2f}）", f8_y="振幅与零炮检距振幅之比")
T = ZH if LANG == "zh" else EN
PANEL = "abcdef"; CJK = LATIN = None; SEP = "，" if LANG == "zh" else ", "


def panel(ax, i):
    ax.text(-0.02, 1.02, f"({PANEL[i]})", transform=ax.transAxes, ha="right", va="bottom", fontsize=7)


CJK_RUN = re.compile("[　-鿿＀-￯]+")   # Chinese characters and full-width punctuation


def cjk_math_fix(fig):
    """mathtext lays out text outside $...$ with the label's own font and no glyph fallback (only STIX, no Chinese).
    So in mixed Chinese + maths labels each run of Chinese characters / full-width punctuation outside the dollars is
    put into $\\mathsf{...}$ (mathtext.sf = the Chinese font) at save time, while the label font stays the Latin font:
    Latin letters and digits outside the maths, and the upright maths (mathtext.rm), are then Times New Roman.
    The label strings in ZH are not changed."""
    if LANG != "zh":
        return
    import matplotlib.text as mt
    for t in fig.findobj(mt.Text):
        s = t.get_text()
        if "$" in s and CJK_RUN.search(s):
            parts = s.split("$")
            for i in range(0, len(parts), 2):   # even parts are outside the dollars
                parts[i] = CJK_RUN.sub(lambda m: "$\\mathsf{" + m.group(0) + "}$", parts[i])
            t.set_text("$".join(parts)); t.set_fontfamily(LATIN)
            if "\n" in s:
                t.set_linespacing(1.6)   # mathtext lines sit closer than plain text lines


def save(fig, name, sources, note):
    cjk_math_fix(fig)
    for ext in ("pdf", "png"):
        fig.savefig(OUT / f"{name}.{ext}", dpi=300, facecolor="white")
    manifest[name] = dict(sources=[Path(s).relative_to(ROOT.parent).as_posix() for s in sources],
                          sha256={str(Path(s).name): hashlib.sha256(Path(s).read_bytes()).hexdigest()[:16] for s in sources}, note=note)
    plt.close(fig)


def fig1():
    fig, ax = plt.subplots(figsize=(W2, 34 * MM)); ax.set_axis_off(); ax.set_xlim(0, 10.2); ax.set_ylim(0, 2)
    for txt, x in zip(T["f1"], (0.05, 2.1, 4.15, 6.2, 8.25)):
        ax.add_patch(FancyBboxPatch((x, 0.3), 1.6, 1.4, boxstyle="round,pad=0.05", fc="white", ec=BLACK, lw=0.8)); ax.text(x + 0.8, 1.0, txt, ha="center", va="center", fontsize=6.5, linespacing=1.3)
    for x in (1.7, 3.75, 5.8, 7.85):
        ax.annotate("", xy=(x + 0.36, 1.0), xytext=(x + 0.02, 1.0), arrowprops=dict(arrowstyle="-|>", lw=0.8, color=BLACK))
    save(fig, "fig1_workflow", [], "schematic")


def fig2():
    p = UQ / "results" / "gate1.json"; R = J(p); fig, ax = plt.subplots(figsize=(W1, W1), layout="constrained")
    for sig, col, mk, lab in ((0.0, BLUE, "o", T["f2_tik"]), (0.447, VERM, "^", T["f2_uq"])):
        rr = [r for r in R if r["kind"] == "uq" and r["sig"] == sig and "dga_num" in r and r["model"] in ("M12", "M13", "M14", "M15")]
        ax.scatter([r["dga_num"] for r in rr], [r["dga_closed"] for r in rr], s=8, marker=mk, facecolors="none", edgecolors=col, lw=0.6, label=f"{lab}{SEP}n = {len(rr)}")
    lim = [-2.0, 1.0]; ax.plot(lim, lim, color="0.5", lw=0.6, ls="--", label="1:1"); ax.set(xlim=lim, ylim=lim, xlabel=T["f2_x"], ylabel=T["f2_y"])
    ax.legend(fontsize=6, loc="upper left"); save(fig, "fig2_gamma_check", [p], "M12-M15, both wavelets, 3 lambdas, 3 events, u grid where numeric exists")


def fig3():
    p = UQ / "results" / "fig_sigma.json"; R = J(p); fig, axs = plt.subplots(1, 2, figsize=(W2, 62 * MM), layout="constrained", sharey=True)
    for i, (ax, fp) in enumerate(zip(axs, (60.0, 40.0))):
        for op, col, ls in (("tikhonov", VERM, "-"), ("offset_consistent", BLUE, "--")):
            for M, mk in (("M1", "o"), ("M12", "s")):
                rr = sorted([r for r in R if r["fp"] == fp and r["op"] == op and r["model"] == M], key=lambda r: r["sigma"])
                ax.plot([r["sigma"] for r in rr], [r["q68_median"] for r in rr], color=col, ls=ls, marker=mk, ms=3, lw=1,
                        label=f"{T['f3_tik'] if op == 'tikhonov' else T['f3_oc']}{SEP}{M}")
        ax.set(xlabel=T["f3_x"], title=T["f3_t"].format(int(fp)), xlim=(0, 0.5), ylim=(0, None)); panel(ax, i)
    axs[0].set_ylabel(T["f3_y"]); axs[1].legend(fontsize=6, loc="upper left")
    save(fig, "fig3_avo_vs_sigma", [p], "closed form, lambda 10^-3.5, models used in development; prior truncated at +-3 sigma; q68 because the RMS of the ratio G/A is tail-dominated")


def fig4():
    p1, p2 = UQ / "results" / "lam_scan.json", UQ / "results" / "lam_hi_verdict.json"; L, V = J(p1), J(p2)
    fig, axs = plt.subplots(1, 2, figsize=(W2, 60 * MM), layout="constrained", sharey=True); lg = np.log10(L["grid"])
    for i, (ax, fp) in enumerate(zip(axs, (60, 40))):
        for sig, col, ls, mk in (("0.447", VERM, "-", "o"), ("0.250", BLUE, "--", "s")):
            ax.plot(lg, L["D"][f"sig{sig}_fp{fp}"], color=col, ls=ls, marker=mk, ms=3, lw=1, label=f"$\\sigma={float(sig):.3g}$")
        pos = {}
        for snr in (30, 24, 18, 12):
            pos.setdefault(V["cells"][f"{snr}_aln_{fp}"]["log10lam_median"], []).append(snr)
        for x, ss in pos.items():
            ax.axvline(x, color="0.4", lw=0.8, ls=":"); ax.text(x + 0.04, -0.44, ", ".join(f"{q}" for q in ss) + " dB", rotation=90, fontsize=6, va="bottom", ha="left", color="0.25")
        ax.axhline(0, color="0.6", lw=0.5); ax.set(xlabel="$\\log_{10}\\lambda$", title=T["f4_t"].format(fp), ylim=(-0.45, 0.05)); panel(ax, i)
    axs[0].set_ylabel(T["f4_y"]); axs[0].legend(fontsize=6, loc="lower right")
    axs[0].yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:.1f}".replace("-", "−")))
    save(fig, "fig4_uq_cost_vs_lambda", [p1, p2], "closed-form D (P1 validated); dotted lines = median residual-principle lambda of one-step Tikhonov per SNR (aligned, M12-M15)")


def fig5():
    p = ROOT / "results" / "fig_window.json"; R = J(p); fig, axs = plt.subplots(1, 2, figsize=(W2, 60 * MM), layout="constrained")
    for est, col, ls, mk, lab in (("original", VERM, "-", "o", T["f5_orig"]), ("equal_length", BLUE, "--", "s", T["f5_eq"])):
        fps = sorted({r["fp"] for r in R["bias"]}); v = [[r["ln_err"] for r in R["bias"] if r["fp"] == f and r["est"] == est] for f in fps]
        axs[0].plot(fps, [np.mean(x) for x in v], color=col, ls=ls, marker=mk, ms=3, lw=1, label=lab)
        for f, x in zip(fps, v):
            axs[0].plot([f, f], [min(x), max(x)], color=col, lw=0.8)
        offs = sorted({r["offset"] for r in R["per_offset"]}); m = [np.median([r["y_over_x"] for r in R["per_offset"] if r["est"] == est and r["offset"] == o]) for o in offs]
        axs[1].plot(offs, m, color=col, ls=ls, marker=mk, ms=3, lw=1, label=lab)
    axs[0].axhline(0, color="0.6", lw=0.5); axs[0].set(xlabel=T["f5_x0"], ylabel=T["f5_y0"])
    axs[1].axhline(1, color="0.6", lw=0.5); axs[1].set(xlabel=T["f5_x1"], ylabel=T["f5_y1"])
    for i, ax in enumerate(axs):
        panel(ax, i)
    axs[0].legend(fontsize=6); save(fig, "fig5_window_compression", [p], "used model M1, noise-free, no AVO, c=1")


def fig6():
    pa = [ROOT / "results" / f"a1_{m}.json" for m in ("M16", "M17", "M18", "M19")]; pb = [ROOT / "results" / f"a1b_{m}.json" for m in ("M20", "M21", "M22", "M23")]
    A = [r for q in pa for r in J(q)]; B = [r for q in pb for r in J(q)]
    sel = lambda R, fp, n=50: [np.log(r["c_hat"] / r["c"]) for r in R if r["kind"] == "avo" and r["snr"] == 24 and r["v"] == 1.0 and r["n"] == n and r["fp"] == fp and r["c_hat"] > 0]
    fig, axs = plt.subplots(1, 2, figsize=(W2, 62 * MM), layout="constrained"); data, labs = [], []
    for name, R in ((T["f6_orig"], A), (T["f6_corr"], B)):
        for fp in (40.0, 60.0):
            data.append(sel(R, fp)); labs.append(f"{name}\n{int(fp)} Hz")
    axs[0].boxplot(data, widths=0.5, showfliers=True, flierprops=dict(ms=2)); axs[0].set_xticks(range(1, 5), labs, fontsize=6)
    rng = np.random.default_rng(0)
    for i, d in enumerate(data):
        axs[0].scatter(i + 1 + rng.uniform(-0.15, 0.15, len(d)), d, s=2, color=[VERM, VERM, BLUE, BLUE][i], alpha=0.4, lw=0)
    for y in (-0.1, 0.1):
        axs[0].axhline(y, color="0.6", lw=0.5, ls=":")
    axs[0].axhline(0, color="0.4", lw=0.5); axs[0].set_ylabel(T["f6_y0"])
    for fp, col, ls, mk in ((40.0, BLUE, "-", "o"), (60.0, GREEN, "--", "s")):
        for snr in (24, 30):
            rms = [np.sqrt(np.mean(np.square([np.log(r["c_hat"] / r["c"]) for r in B if r["kind"] == "avo" and r["snr"] == snr and r["v"] == 1.0 and r["n"] == n and r["fp"] == fp and r["c_hat"] > 0]))) for n in (10, 50, 200)]
            axs[1].plot((10, 50, 200), rms, color=col, ls=ls if snr == 24 else ":", marker=mk, ms=3, lw=1, label=T["f6_leg"].format(int(fp), snr))
    axs[1].set(xscale="log", xlabel=T["f6_x1"], ylabel=T["f6_y1"], ylim=(0, None))
    axs[1].axhline(0.1, color="0.6", lw=0.5, ls=":"); axs[1].legend(fontsize=6)
    for i, ax in enumerate(axs):
        panel(ax, i)
    save(fig, "fig6_chat_accuracy", pa + pb, "gate A1 (first attempt) and A1b (second attempt, new models)")


def fig7():
    p = ROOT / "results" / "a2_verdict.json"; D = J(p)["descriptive"]; cs = (0.5, 0.7, 1.0, 1.4, 2.0)
    fig, axs = plt.subplots(1, 2, figsize=(W2, 62 * MM), layout="constrained")
    for m, col, ls, mk, lab in zip(("oracle", "uncal", "tik_uncal", "cal"), (BLACK, VERM, PINK, BLUE), (":", "-", "-.", "--"), ("D", "o", "v", "s"), T["f7"]):
        for ax, q in zip(axs, ("ga", "wave")):
            ax.plot(cs, [D[f"c{c:.1f}"][m][q] for c in cs], color=col, ls=ls, marker=mk, ms=3.5, lw=1, label=lab)
    axs[0].set(xscale="log", yscale="log", xlabel=T["f7_x"], ylabel=T["f7_y0"])
    axs[1].set(xscale="log", yscale="log", xlabel=T["f7_x"], ylabel=T["f7_y1"])
    for i, ax in enumerate(axs):
        ax.set_xticks(cs, [str(c) for c in cs]); ax.minorticks_off(); panel(ax, i)
    axs[0].legend(fontsize=6); save(fig, "fig7_end_to_end", [p], "gate A2, aligned events + thin beds, 24/30 dB, M24-M27, 10 seeds per cell")


def fig8():
    p = ROOT / "results" / "fig_example.json"; R = J(p); fig, axs = plt.subplots(2, 3, figsize=(W2, 100 * MM), layout="constrained", sharex="col")
    for i, c in enumerate((0.5, 2.0)):
        for j in range(3):
            r = [q for q in R if q["c"] == c and q["event"] == j][0]; ax = axs[i, j]; s2 = np.array(r["sin2"])
            for k, col, ls, mk, lab in zip(("truth", "oracle", "uncalibrated", "calibrated"), (BLACK, "0.45", VERM, BLUE), ("-", ":", "-", "--"), (None, "D", "o", "s"), T["f8"]):
                y = np.array(r[k]); ax.plot(s2, y / y[0], color=col, ls=ls, marker=mk, ms=2.5, lw=1, label=lab)
            ax.set_title(T["f8_t"].format(c, j + 1, r["c_hat"]), fontsize=7); panel(ax, 3 * i + j)
            if i == 1:
                ax.set_xlabel("$\\sin^2\\theta$")
            if j == 0:
                ax.set_ylabel(T["f8_y"])
    axs[0, 0].legend(fontsize=6); save(fig, "fig8_example_gather", [p], "M24, 60 Hz, events_aligned, seed 0, 24 dB; one realisation, illustrative; each curve divided by its own zero-offset value (shows the AVO shape; absolute level differs because strong attenuation limits the recoverable band)")


def first_font(names):
    have = {f.name for f in font_manager.fontManager.ttflist}
    return next((n for n in names if n in have), None)


def style():
    """The English set keeps the scientific-visualization preset when it is installed; otherwise (and for zh) plain rcParams."""
    base = {"pdf.fonttype": 42, "font.size": 7, "axes.titlesize": 7, "axes.labelsize": 7, "xtick.labelsize": 6, "ytick.labelsize": 6,
            "axes.spines.top": False, "axes.spines.right": False, "axes.linewidth": 0.8, "legend.frameon": False,
            "xtick.direction": "out", "ytick.direction": "out", "savefig.bbox": "tight", "savefig.pad_inches": 0.02}
    if LANG == "zh":
        latin = first_font(["Times New Roman", "Liberation Serif", "DejaVu Serif"])
        cjk = first_font(["SimSun", "Songti SC", "STSong", "Noto Serif CJK SC", "Source Han Serif SC", "SimHei", "Microsoft YaHei", "WenQuanYi Zen Hei"])
        assert cjk, "no Chinese font found"
        global CJK, LATIN; CJK, LATIN = cjk, latin
        # A label that contains $...$ is laid out entirely by mathtext, whose text outside the dollars uses the label's
        # own font with no per-glyph fallback; cjk_math_fix moves the Chinese runs of such labels into \mathsf (Chinese
        # font), so upright maths and Latin text stay Times; italics (variables) Times italic.
        base.update({"font.family": [latin, cjk], "axes.unicode_minus": True, "mathtext.fontset": "custom",
                     "mathtext.rm": latin, "mathtext.sf": cjk, "mathtext.it": f"{latin}:italic", "mathtext.bf": f"{latin}:bold",
                     "mathtext.fallback": "stix"})
        print("fonts:", latin, "+", cjk)
    return base


if __name__ == "__main__":
    logging.getLogger("matplotlib.font_manager").setLevel(logging.ERROR)
    try:
        from style_presets import style_context
        ctx = style_context("default", palette_name="okabe_ito_on_white") if LANG == "en" else None
    except ImportError:
        ctx = None
    import contextlib
    with (ctx or contextlib.nullcontext()):
        plt.rcParams.update(style())
        for f in (fig1, fig2, fig3, fig4, fig5, fig6, fig7, fig8):
            try:
                f(); print("ok", f.__name__)
            except FileNotFoundError as e:
                print("skip", f.__name__, e)
    (OUT / "MANIFEST.json").write_text(json.dumps(dict(lang=LANG, figures=manifest), indent=1, ensure_ascii=False), encoding="utf-8")
