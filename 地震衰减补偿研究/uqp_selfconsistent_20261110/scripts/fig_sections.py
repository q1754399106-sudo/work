"""Section images for visual inspection (exploratory, NOTE_ZH.md section 12). SR route, semi-synthetic development data.

  compute ID VERSION   sections ideal / input / IQ70 / W0 / SC-a / MW0_K9 -> results/fig_sections/ID_VERSION.npz
                       (SC-a: cal = experiments 12 / 17 setting; nd = notch-tolerant band + selected (2^-9, 2))
  plot ID VERSION MPLDIR   full section and zoom -> figures/ID_VERSION_{full,zoom}.png
Each panel is clipped at its own 99th percentile of |amplitude| (S is a correlation, insensitive to overall scale).
"""
import sys, json
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import uqp, mwdip  # noqa: E402
from uqp import np, lib, semi, ROOT, DT  # noqa: E402
from uqw import uqw, snr_rows, G_IQ, TW, lib12, load  # noqa: E402
from lib17 import score, sc11, A_TUNED  # noqa: E402

OUT = ROOT / "results" / "fig_sections"; FIG = ROOT / "figures"; F_OP = np.fft.rfftfreq(1024, DT); S = semi.gain_sigma2(G_IQ); ANCHOR = 1.4
NAMES = [("ideal", "理想剖面"), ("input", "输入（含衰减和噪声）"), ("IQ70", "IQ70 稳定反 Q"), ("W0", "W0 单道维纳"),
         ("SCa", "SC-a 结构约束稀疏反演"), ("MW0_K9", "多道维纳（9 道）")]


def compute(did, version):
    from pilot_sc_robust import wavelets
    cm = semi.curtain_cm(); Q = np.asarray(cm["Q"], float); arr, meta = load(did)
    D = arr["signal"] + (1.0 if version == "cal" else 4.0) * arr["noise"]
    c = json.loads((ROOT / "results" / "pilot_semi" / f"{did}_{version}.json").read_text(encoding="utf-8"))["sr"]["c_sr"]
    x = lib.absorption(dict(cm, Q=Q * c), TW); dx = x - lib.absorption(dict(cm, Q=Q * c), [ANCHOR])[0]
    g = mwdip.Gather(D, TW, 4); inv = snr_rows(D, TW, F_OP); mod = g.model(x, dx); p = sc11.prepare_sc(D, cm)
    if version == "cal":
        X, info, _ = lib.solve(cm, Q * c, p["wv"], D, p["sigma"], p["K"], p["live"], A_TUNED, 0.25)
    else:
        X, info, _ = lib.solve(cm, Q * c, wavelets(D, cm, False), D, p["sigma"], p["K"], p["live"], 2.0 ** -9, 2.0)
    sec = dict(ideal=arr["ideal"], input=D, IQ70=lib12.iq(D, cm, Q, c, c, G_IQ), W0=uqw(D, cm, Q, c, 0.0, TW, ANCHOR, 0.0, inv_snr=inv),
               SCa=X, MW0_K9=g.apply(x, dx, 0.0, mod, F_OP, cm["fref"], S, 4))
    sc = {k: float(score(v, arr["ideal"])["S"]) for k, v in sec.items()}
    OUT.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(OUT / f"{did}_{version}.npz", **{k: np.asarray(v, np.float32) for k, v in sec.items()},
                        fault_dist=np.asarray(arr["fault_dist"], np.float32))
    (OUT / f"{did}_{version}.json").write_text(json.dumps(dict(id=did, version=version, c_sr=c, ratio=c / meta["c_eff_cum"], S=sc,
                                                               sc_converged=bool(info["converged"])), indent=1), encoding="utf-8")
    print(did, version, {k: round(v, 3) for k, v in sc.items()}, flush=True)


def plot(did, version, mpldir):
    sys.path.insert(0, mpldir)
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    plt.rcParams["font.sans-serif"] = ["Microsoft YaHei"]; plt.rcParams["axes.unicode_minus"] = False
    z = np.load(OUT / f"{did}_{version}.npz"); meta = json.loads((OUT / f"{did}_{version}.json").read_text(encoding="utf-8")); FIG.mkdir(exist_ok=True)
    nt, nx = z["ideal"].shape; noise = "标定噪声" if version == "cal" else "4 倍噪声"
    for tag, (t0, t1, x0, x1) in dict(full=(1.2, 3.0, 0, nx), zoom=(2.1, 2.7, 60, 140)).items():
        r = (TW >= t0) & (TW <= t1); fig, ax = plt.subplots(2, 3, figsize=(15, 9.5), sharex=True, sharey=True)
        for a, (k, name) in zip(ax.ravel(), NAMES):
            v = z[k][r, x0:x1]; clip = np.quantile(np.abs(v), 0.99) + 1e-30
            a.imshow(v, aspect="auto", cmap="gray", vmin=-clip, vmax=clip, extent=(x0, x1, TW[r][-1], TW[r][0]), interpolation="bilinear")
            a.set_title(name if k in ("ideal", "input") and False else (name + ("" if k == "ideal" else "  S = %.3f" % meta["S"][k])), fontsize=12)
        for a in ax[-1]:
            a.set_xlabel("道号")
        for a in ax[:, 0]:
            a.set_ylabel("时间 (s)")
        fig.suptitle("%s，%s，Q 由谱比法估计（估计值/真值 = %.2f）；S 为与理想剖面的相似度（全剖面）" % (did, noise, meta["ratio"]), fontsize=13)
        fig.tight_layout(rect=(0, 0, 1, 0.96)); fig.savefig(FIG / f"{did}_{version}_{tag}.png", dpi=110); plt.close(fig)
    print("saved", did, version)


if __name__ == "__main__":
    if sys.argv[1] == "compute":
        compute(sys.argv[2], sys.argv[3])
    else:
        plot(sys.argv[2], sys.argv[3], sys.argv[4])
