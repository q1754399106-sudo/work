"""Paper figures 1-8 (PAPER_OUTLINE_ZH.md). Provisional general style (journal not chosen): Okabe-Ito colours with redundant
line styles / markers, English labels, PDF (Type 42 fonts) + PNG 300 dpi at 174 mm (two-column) or 84 mm width.
Reads only saved result files; no data are altered. usage: make_figs.py  -> figures/fig*.pdf|png, figures/MANIFEST.json"""
import sys, json, hashlib
from pathlib import Path
import numpy as np
sys.path.append(r"D:/ChatGPT/研究/地震衰减补偿研究/attenuation_run_20261019/dependencies")  # vendored matplotlib 3.11.2 (cp312)
sys.path.insert(0, r"C:/Users/qq175/.claude/skills/scientific-visualization/scripts")
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
from style_presets import style_context

HERE = Path(__file__).resolve().parent; ROOT = HERE.parent; UQ = ROOT.parent / "uq_avo_20261118"; OUT = ROOT / "figures"; OUT.mkdir(exist_ok=True)
BLUE, VERM, GREEN, PINK, BLACK = "#0072B2", "#D55E00", "#009E73", "#CC79A7", "#000000"
MM = 1 / 25.4; W2, W1 = 174 * MM, 84 * MM; manifest = {}
J = lambda p: json.loads(Path(p).read_text(encoding="utf-8"))


def save(fig, name, sources, note):
    for ext in ("pdf", "png"):
        fig.savefig(OUT / f"{name}.{ext}", dpi=300, facecolor="white")
    manifest[name] = dict(sources=[str(Path(s).relative_to(ROOT.parent)) for s in sources],
                          sha256={str(Path(s).name): hashlib.sha256(Path(s).read_bytes()).hexdigest()[:16] for s in sources}, note=note)
    plt.close(fig)


def fig1():
    fig, ax = plt.subplots(figsize=(W2, 34 * MM)); ax.set_axis_off(); ax.set_xlim(0, 10.2); ax.set_ylim(0, 2)
    boxes = [("Pre-stack gathers\n(no NMO)", 0.05), ("Offset spectral ratio,\nequal-length windows\nalong moveout", 2.1),
             ("$\\hat c=Q_\\mathrm{true}/Q_\\mathrm{model}$\n(50-CMP\nsupergather)", 4.15), ("Path-aware\ncompensation\nwith $\\hat c\\,Q_\\mathrm{model}$", 6.2),
             ("Per-event AVO\n(intercept,\ngradient)", 8.25)]
    for txt, x in boxes:
        ax.add_patch(FancyBboxPatch((x, 0.3), 1.6, 1.4, boxstyle="round,pad=0.05", fc="white", ec=BLACK, lw=0.8)); ax.text(x + 0.8, 1.0, txt, ha="center", va="center", fontsize=6.5)
    for x in (1.7, 3.75, 5.8, 7.85):
        ax.annotate("", xy=(x + 0.36, 1.0), xytext=(x + 0.02, 1.0), arrowprops=dict(arrowstyle="-|>", lw=0.8, color=BLACK))
    save(fig, "fig1_workflow", [], "schematic")


def fig2():
    p = UQ / "results" / "gate1.json"; R = J(p); fig, ax = plt.subplots(figsize=(W1, W1), layout="constrained")
    for sig, col, mk, lab in ((0.0, BLUE, "o", "Tikhonov ($\\sigma=0$)"), (0.447, VERM, "^", "UQ ($\\sigma=0.447$)")):
        rr = [r for r in R if r["kind"] == "uq" and r["sig"] == sig and "dga_num" in r and r["model"] in ("M12", "M13", "M14", "M15")]
        ax.scatter([r["dga_num"] for r in rr], [r["dga_closed"] for r in rr], s=8, marker=mk, facecolors="none", edgecolors=col, lw=0.6, label=f"{lab}, n={len(rr)}")
    lim = [-2.0, 1.0]; ax.plot(lim, lim, color="0.5", lw=0.6, ls="--", label="1:1"); ax.set(xlim=lim, ylim=lim, xlabel="Numerical $\\Delta(G/A)$ (time-varying operator)",
                                                                                          ylabel="Closed-form $\\Delta(G/A)$ from $\\Gamma(x)$")
    ax.legend(fontsize=6, loc="upper left"); save(fig, "fig2_gamma_check", [p], "M12-M15, both wavelets, 3 lambdas, 3 events, u grid where numeric exists")


def fig3():
    p = UQ / "results" / "fig_sigma.json"; R = J(p); fig, axs = plt.subplots(1, 2, figsize=(W2, 62 * MM), layout="constrained", sharey=True)
    for ax, fp in zip(axs, (60.0, 40.0)):
        for op, col, ls in (("tikhonov", VERM, "-"), ("offset_consistent", BLUE, "--")):
            for M, mk in (("M1", "o"), ("M12", "s")):
                rr = sorted([r for r in R if r["fp"] == fp and r["op"] == op and r["model"] == M], key=lambda r: r["sigma"])
                ax.plot([r["sigma"] for r in rr], [r["q68_median"] for r in rr], color=col, ls=ls, marker=mk, ms=3, lw=1,
                        label=f"{'Tikhonov' if op == 'tikhonov' else 'Offset-consistent'}, {M}")
        ax.set(xlabel="$\\sigma_{\\ln Q}$ (prior width of the Q-level error)", title=f"{int(fp)} Hz Ricker", xlim=(0, 0.5), ylim=(0, None))
    axs[0].set_ylabel("68th percentile of $|\\Delta(G/A)|$ over the prior\n(median of 3 events)"); axs[1].legend(fontsize=6, loc="upper left")
    save(fig, "fig3_avo_vs_sigma", [p], "closed form, lambda 10^-3.5, models used in development; prior truncated at +-3 sigma; q68 because the RMS of the ratio G/A is tail-dominated")


def fig4():
    p1, p2 = UQ / "results" / "lam_scan.json", UQ / "results" / "lam_hi_verdict.json"; L, V = J(p1), J(p2)
    fig, axs = plt.subplots(1, 2, figsize=(W2, 60 * MM), layout="constrained", sharey=True); lg = np.log10(L["grid"])
    for ax, fp in zip(axs, (60, 40)):
        for sig, col, ls, mk in (("0.447", VERM, "-", "o"), ("0.250", BLUE, "--", "s")):
            ax.plot(lg, L["D"][f"sig{sig}_fp{fp}"], color=col, ls=ls, marker=mk, ms=3, lw=1, label=f"$\\sigma={float(sig):.3g}$")
        pos = {}
        for snr in (30, 24, 18, 12):
            pos.setdefault(V["cells"][f"{snr}_aln_{fp}"]["log10lam_median"], []).append(snr)
        for x, ss in pos.items():
            ax.axvline(x, color="0.4", lw=0.8, ls=":"); ax.text(x + 0.04, -0.44, ", ".join(f"{q}" for q in ss) + " dB", rotation=90, fontsize=6, va="bottom", ha="left", color="0.25")
        ax.axhline(0, color="0.6", lw=0.5); ax.set(xlabel="$\\log_{10}\\lambda$", title=f"{fp} Hz Ricker", ylim=(-0.45, 0.05))
    axs[0].set_ylabel("Extra $\\Delta(G/A)$ from UQ, $D$\n(UQ minus Tikhonov, Q correct)"); axs[0].legend(fontsize=6, loc="lower right")
    save(fig, "fig4_uq_cost_vs_lambda", [p1, p2], "closed-form D (P1 validated); dotted lines = median residual-principle lambda of one-step Tikhonov per SNR (aligned, M12-M15)")


def fig5():
    p = ROOT / "results" / "fig_window.json"; R = J(p); fig, axs = plt.subplots(1, 2, figsize=(W2, 60 * MM), layout="constrained")
    for est, col, ls, mk, lab in (("original", VERM, "-", "o", "Compressed far window (original)"), ("equal_length", BLUE, "--", "s", "Equal-length far window")):
        fps = sorted({r["fp"] for r in R["bias"]}); v = [[r["ln_err"] for r in R["bias"] if r["fp"] == f and r["est"] == est] for f in fps]
        axs[0].plot(fps, [np.mean(x) for x in v], color=col, ls=ls, marker=mk, ms=3, lw=1, label=lab)
        for f, x in zip(fps, v):
            axs[0].plot([f, f], [min(x), max(x)], color=col, lw=0.8)
        offs = sorted({r["offset"] for r in R["per_offset"]}); m = [np.median([r["y_over_x"] for r in R["per_offset"] if r["est"] == est and r["offset"] == o]) for o in offs]
        axs[1].plot(offs, m, color=col, ls=ls, marker=mk, ms=3, lw=1, label=lab)
    axs[0].axhline(0, color="0.6", lw=0.5); axs[0].set(xlabel="Ricker peak frequency (Hz)", ylabel="$\\ln(\\hat c/c)$ (noise-free; mean, range of 4 seeds)")
    axs[1].axhline(1, color="0.6", lw=0.5); axs[1].set(xlabel="Offset (m)", ylabel="Measured / predicted $\\Delta t^*$ per window\n(40 Hz, median)")
    axs[0].legend(fontsize=6); save(fig, "fig5_window_compression", [p], "used model M1, noise-free, no AVO, c=1")


def fig6():
    pa = [ROOT / "results" / f"a1_{m}.json" for m in ("M16", "M17", "M18", "M19")]; pb = [ROOT / "results" / f"a1b_{m}.json" for m in ("M20", "M21", "M22", "M23")]
    A = [r for q in pa for r in J(q)]; B = [r for q in pb for r in J(q)]
    sel = lambda R, fp, n=50: [np.log(r["c_hat"] / r["c"]) for r in R if r["kind"] == "avo" and r["snr"] == 24 and r["v"] == 1.0 and r["n"] == n and r["fp"] == fp and r["c_hat"] > 0]
    fig, axs = plt.subplots(1, 2, figsize=(W2, 62 * MM), layout="constrained"); data, labs = [], []
    for name, R in (("Original\n(M16-M19)", A), ("Corrected\n(M20-M23)", B)):
        for fp in (40.0, 60.0):
            data.append(sel(R, fp)); labs.append(f"{name}\n{int(fp)} Hz")
    bp = axs[0].boxplot(data, widths=0.5, showfliers=True, flierprops=dict(ms=2)); axs[0].set_xticks(range(1, 5), labs, fontsize=6)
    rng = np.random.default_rng(0)
    for i, d in enumerate(data):
        axs[0].scatter(i + 1 + rng.uniform(-0.15, 0.15, len(d)), d, s=2, color=[VERM, VERM, BLUE, BLUE][i], alpha=0.4, lw=0)
    for y in (-0.1, 0.1):
        axs[0].axhline(y, color="0.6", lw=0.5, ls=":")
    axs[0].axhline(0, color="0.4", lw=0.5); axs[0].set_ylabel("$\\ln(\\hat c/c)$, 24 dB, AVO, N = 50 CMPs\n(384 estimates per box)")
    for fp, col, ls, mk in ((40.0, BLUE, "-", "o"), (60.0, GREEN, "--", "s")):
        for snr, a in ((24, 1.0), (30, 0.5)):
            rms = [np.sqrt(np.mean(np.square([np.log(r["c_hat"] / r["c"]) for r in B if r["kind"] == "avo" and r["snr"] == snr and r["v"] == 1.0 and r["n"] == n and r["fp"] == fp and r["c_hat"] > 0]))) for n in (10, 50, 200)]
            axs[1].plot((10, 50, 200), rms, color=col, ls=ls if snr == 24 else ":", marker=mk, ms=3, lw=1, label=f"{int(fp)} Hz, {snr} dB")
    axs[1].set(xscale="log", xlabel="CMPs per supergather, N (log scale)", ylabel="RMS $\\ln(\\hat c/c)$ (corrected estimator)", ylim=(0, None))
    axs[1].axhline(0.1, color="0.6", lw=0.5, ls=":"); axs[1].legend(fontsize=6); save(fig, "fig6_chat_accuracy", pa + pb, "gate A1 (first attempt) and A1b (second attempt, new models)")


def fig7():
    p = ROOT / "results" / "a2_verdict.json"; D = J(p)["descriptive"]; cs = (0.5, 0.7, 1.0, 1.4, 2.0)
    fig, axs = plt.subplots(1, 2, figsize=(W2, 62 * MM), layout="constrained")
    for m, col, ls, mk, lab in (("oracle", BLACK, ":", "D", "True Q"), ("uncal", VERM, "-", "o", "Model Q (sat. damping)"),
                                ("tik_uncal", PINK, "-.", "v", "Model Q (Tikhonov)"), ("cal", BLUE, "--", "s", "Calibrated $\\hat cQ_\\mathrm{model}$")):
        for ax, q in zip(axs, ("ga", "wave")):
            ax.plot(cs, [D[f"c{c:.1f}"][m][q] for c in cs], color=col, ls=ls, marker=mk, ms=3.5, lw=1, label=lab)
    axs[0].set(xscale="log", yscale="log", xlabel="$c=Q_\\mathrm{true}/Q_\\mathrm{model}$ (log scale)", ylabel="Per-event $|\\Delta(G/A)|$, noisy (log scale)\nmedian of 32 cells")
    axs[1].set(xscale="log", yscale="log", xlabel="$c=Q_\\mathrm{true}/Q_\\mathrm{model}$ (log scale)", ylabel="Waveform error (log scale)\nmedian of 32 cells")
    for ax in axs:
        ax.set_xticks(cs, [str(c) for c in cs]); ax.minorticks_off()
    axs[0].legend(fontsize=6); save(fig, "fig7_end_to_end", [p], "gate A2, aligned events + thin beds, 24/30 dB, M24-M27, 10 seeds per cell")


def fig8():
    p = ROOT / "results" / "fig_example.json"; R = J(p); fig, axs = plt.subplots(2, 3, figsize=(W2, 100 * MM), layout="constrained", sharex="col")
    for i, c in enumerate((0.5, 2.0)):
        for j in range(3):
            r = [q for q in R if q["c"] == c and q["event"] == j][0]; ax = axs[i, j]; s2 = np.array(r["sin2"])
            for k, col, ls, mk, lab in (("truth", BLACK, "-", None, "Truth"), ("oracle", "0.45", ":", "D", "True Q"), ("uncalibrated", VERM, "-", "o", "Model Q"),
                                        ("calibrated", BLUE, "--", "s", "Calibrated")):
                y = np.array(r[k]); ax.plot(s2, y / y[0], color=col, ls=ls, marker=mk, ms=2.5, lw=1, label=lab)
            ax.set_title(f"c = {c}, event {j + 1} ($\\hat c$ = {r['c_hat']:.2f})", fontsize=7)
            if i == 1:
                ax.set_xlabel("$\\sin^2\\theta$")
            if j == 0:
                ax.set_ylabel("Amplitude / zero-offset amplitude")
    axs[0, 0].legend(fontsize=6); save(fig, "fig8_example_gather", [p], "M24, 60 Hz, events_aligned, seed 0, 24 dB; one realisation, illustrative; each curve divided by its own zero-offset value (shows the AVO shape; absolute level differs because strong attenuation limits the recoverable band)")


if __name__ == "__main__":
    with style_context("default", palette_name="okabe_ito_on_white"):
        plt.rcParams.update({"pdf.fonttype": 42, "font.size": 7, "axes.titlesize": 7, "axes.labelsize": 7, "xtick.labelsize": 6, "ytick.labelsize": 6})
        for f in (fig1, fig2, fig3, fig4, fig5, fig6, fig7, fig8):
            try:
                f(); print("ok", f.__name__)
            except FileNotFoundError as e:
                print("skip", f.__name__, e)
    (OUT / "MANIFEST.json").write_text(json.dumps(manifest, indent=1))
