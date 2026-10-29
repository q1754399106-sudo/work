"""Gather figures for experiment 14 (descriptive, not in PLAN_ZH.md).

For every finished test dataset (results/test/<id>.json) the input is rebuilt from its seed exactly as in run14.py,
the ideal target as in run11.py, and the input is compensated with SC-a (sc11.sca, the pipeline of experiments 7-13)
at the Q factor chosen by each method: true c, SR, QScan, QScan-W(2^-4) and QScan-W(2^-4)_cal when
analysis/calibration.json exists. M1 (median trace correlation with the ideal, 2.30-2.80 s) is printed per panel.

usage: gathers14.py [id ...]    default: all finished test datasets without a figure yet
writes figures/gathers/<id>.png and cache/gathers/<id>.npz
"""
import sys, json, time
from env14 import np, semi, sc11, gen, synth, wq, ROOT
from run14 import STRUCT, NOISE, WELLS_L, signal_L

sys.path.insert(0, str(ROOT.parent / "attenuation_run_20261019" / "dependencies"))
import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "DejaVu Sans"]; plt.rcParams["axes.unicode_minus"] = False
LAM = "0.0625"
FIG, CACHE = ROOT / "figures" / "gathers", ROOT / "cache" / "gathers"


def factors(r):
    out = {"真 c": r["c_true"], "SR": r["c_sr"], "QScan": r["qscan"], "QScan-W(λ=1/16)": r["qscan_w"][LAM]}
    cal = ROOT / "analysis" / "calibration.json"
    if cal.exists():
        b = json.loads(cal.read_text(encoding="utf-8"))[LAM]
        out["QScan-W(λ=1/16)_cal"] = r["qscan_w"][LAM] * float(np.exp(-b))
    return out


def build(r, ctx, W):
    cm = ctx["cm"]; Q = np.asarray(cm["Q"], float)
    g = signal_L(ctx, W[r["well"]], r["c_true"], r["seed"]) if r["group"] == "L" else gen.signal(ctx, r["c_true"], r["seed"])
    assert abs(g["level"] - r["level"]) < 1e-12, (r["id"], g["level"], r["level"])
    N0 = r["noise_mult"] * g["level"] * ctx["noise"][0]; D = g["x"] + N0
    if "R" not in g:   # signal_L keeps only x; rebuild R, S, scale with the same seed
        rng = np.random.default_rng(r["seed"]); f = np.fft.rfftfreq(8192, semi.FINE)
        g.update(R=synth.reflectivity(ctx, W[r["well"]], rng), S=semi.source_spectrum(ctx["w_in"], f, 0.0, cm, r["c_true"]),
                 f=f, scale=float(np.sqrt(ctx["Creal"][ctx["band"]].sum() / semi.band_power(g["x"], 1.2, 1.6)[ctx["band"]].sum())))
    f_true = semi.true_band(g["x"], N0)
    ideal = semi.model(g["R"], None, g["S"] * g["scale"] * semi.cut(g["f"], f_true), g["f"], float(cm["fref"]))
    pp = sc11.prepare_sc(D, cm)
    outs = {name: sc11.sca(D, cm, Q, c, pp) for name, c in factors(r).items()}
    return D, ideal, outs


def plot(r, D, ideal, outs, path):
    fac = factors(r); panels = [("理想目标", ideal, None), ("输入（含衰减与噪声）", D, sc11.m1(D, ideal))]
    panels += [(f"SC-a，{k} = {fac[k]:.3f}", o, sc11.m1(o, ideal)) for k, o in outs.items()]
    ncol = 4 if len(panels) > 6 else 3; nrow = -(-len(panels) // ncol)
    t = semi.TW; ext = [0, D.shape[1] - 1, t[-1], t[0]]
    fig, axs = plt.subplots(nrow, ncol, figsize=(4.2 * ncol, 4.6 * nrow), sharex=True, sharey=True, squeeze=False)
    for ax, (title, X, s) in zip(axs.flat, panels):
        v = np.percentile(np.abs(X), 99) or 1.0
        ax.imshow(X, aspect="auto", cmap="gray_r", vmin=-v, vmax=v, extent=ext, interpolation="bilinear")
        ax.axhspan(*semi.WIN_LOG, color="tab:orange", alpha=0.08, lw=0)
        ax.set_title(title if s is None else f"{title}\nM1 = {s['M1']:.3f}", fontsize=10)
    for ax in axs.flat[len(panels):]:
        ax.axis("off")
    for ax in axs[:, 0]:
        ax.set_ylabel("时间 / s")
    for ax in axs[-1]:
        ax.set_xlabel("道号")
    grp = {"T": "T 组（主组）", "L": f"L 组（{r['well']} 长测井）", "N": f"N 组（噪声 ×{r['noise_mult']}）"}[r["group"]]
    fig.suptitle(f"实验 14  {r['id']}  {grp}  真 c = {r['c_true']}  构造线 {STRUCT}，噪声线 {NOISE[0]}"
                 "（各面板按自身 99% 分位裁剪；橙色带为 M1 评分窗 2.30–2.80 s）", fontsize=11)
    fig.tight_layout(); fig.savefig(path, dpi=130); plt.close(fig)


def main(ids):
    FIG.mkdir(parents=True, exist_ok=True); CACHE.mkdir(parents=True, exist_ok=True)
    recs = [json.loads(p.read_text(encoding="utf-8")) for p in sorted((ROOT / "results" / "test").glob("*.json"))]
    recs = [r for r in recs if r["id"] in ids] if ids else [r for r in recs if not (FIG / f"{r['id']}.png").exists()]
    if not recs:
        print("nothing to do"); return
    ctx = gen.prepare(STRUCT, NOISE); P = wq.picks(); W = {w: wq.well(w, P) for w in WELLS_L}
    for r in recs:
        t0 = time.monotonic(); z = CACHE / f"{r['id']}.npz"
        if z.exists():   # re-plot only
            z = dict(np.load(z)); D, ideal = z.pop("Input"), z.pop("ideal")
            outs = {k: z[k] for k in factors(r) if k in z}
            if len(outs) == len(factors(r)):
                plot(r, D, ideal, outs, FIG / f"{r['id']}.png"); print(r["id"], "re-plotted", flush=True); continue
        D, ideal, outs = build(r, ctx, W)
        np.savez_compressed(CACHE / f"{r['id']}.npz", Input=D.astype(np.float32), ideal=ideal.astype(np.float32),
                            **{k: o.astype(np.float32) for k, o in outs.items()})
        plot(r, D, ideal, outs, FIG / f"{r['id']}.png"); print(r["id"], f"{time.monotonic() - t0:.0f}s", flush=True)


if __name__ == "__main__":
    main(sys.argv[1:])
