"""Gather figures for W1-W6 and X1-X3 (Moderate noise): figures/overview.png (observed and target per scene)
and figures/<scene>.png (methods with Q correct and with Q -30%). Reads models/, truth/, cache/runner/results
and cache/solved; checks each plotted relative error against analysis/rows.csv. Needs matplotlib in cache/plotdeps.
"""
from pathlib import Path
import csv, sys

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / "cache" / "plotdeps"), str(ROOT.parent / "novelty_evidence_20261025" / "cache" / "deps")]
import numpy as np
import scipy.io as sio
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei"]
plt.rcParams["axes.unicode_minus"] = False
SCENES = ["W1", "W2", "W3", "W4", "W5", "W6", "X1", "X2", "X3"]
NOISE = "Moderate"
OUT = ROOT / "figures"
NAMES = dict(Standard="标准版", SC="SC", FullCalRefit="FullCalRefit", QScanRefit="QScanRefit",
             QScan_Standard="QScan＋标准版", QScan_SC="QScan＋SC")
ROWS = {(r["scene"], r["noise"], int(r["q_error"])): r
        for r in csv.DictReader(open(ROOT / "analysis" / "rows.csv", encoding="utf-8"))}


def load(path):
    return sio.loadmat(path, squeeze_me=True, struct_as_record=False)


def pick(arr, name):
    return next(x for x in np.atleast_1d(arr) if x.name == name)


def scene_data(s):
    m, t = load(ROOT / "models" / f"{s}.mat"), load(ROOT / "truth" / f"{s}_truth.mat")
    obs = np.asarray(pick(m["obs"], NOISE).D, float)
    T = np.asarray(pick(t["tr"], NOISE).target, float)
    return np.asarray(m["config"].t, float), obs, T


def outputs(s, e):
    res = load(ROOT / "cache" / "runner" / "results" / f"{s}_Q{e:+d}_{NOISE}.mat")["gathers"]
    out = {k: np.asarray(getattr(res, k), float) for k in ("Standard", "FullCalRefit", "QScanRefit", "QScan_Standard")}
    for lab, key in (("SC", "SC"), ("QSC", "QScan_SC")):
        out[key] = np.load(ROOT / "cache" / "solved" / f"{s}_{NOISE}_Q{e:+d}_{lab}.npz")["out"]
    return out


def rel(X, T):
    return float(np.linalg.norm(X - T) / np.linalg.norm(T))


def show(ax, D, t, clip, title):
    ax.imshow(D, aspect="auto", cmap="gray_r", vmin=-clip, vmax=clip, interpolation="bilinear",
              extent=[1, D.shape[1], t[-1], t[0]])
    ax.set_title(title, fontsize=9)
    ax.tick_params(labelsize=7)


def clip_of(D):
    return float(np.percentile(np.abs(D), 99.5)) or 1.0


OUT.mkdir(exist_ok=True)
fig, axes = plt.subplots(3, 6, figsize=(16, 10), constrained_layout=True)
for i, s in enumerate(SCENES):
    t, obs, T = scene_data(s)
    a, b = axes[i // 3, 2 * (i % 3)], axes[i // 3, 2 * (i % 3) + 1]
    show(a, obs, t, clip_of(obs), f"{s} 观测（{NOISE}）")
    show(b, T, t, clip_of(T), f"{s} 真值目标")
    a.set_ylabel("时间 (s)", fontsize=8)
for ax in axes[-1]:
    ax.set_xlabel("道号", fontsize=8)
fig.suptitle("V 家族新种子场景：观测道集与真值目标（W 层状模板，X 复杂地质；各面板单独定色标）", fontsize=12)
fig.savefig(OUT / "overview.png", dpi=150)
plt.close(fig)

checked = 0
for s in SCENES:
    t, obs, T = scene_data(s)
    o0, o30 = outputs(s, 0), outputs(s, -30)
    panels = [[("观测", obs, None), ("真值目标", T, None), ("标准版", o0["Standard"], ("Standard", 0)),
               ("SC", o0["SC"], ("SC", 0))],
              [("FullCalRefit", o30["FullCalRefit"], ("FullCalRefit", -30)), ("QScanRefit", o30["QScanRefit"], ("QScanRefit", -30)),
               ("QScan＋标准版", o30["QScan_Standard"], ("QScan_Standard", -30)), ("QScan＋SC", o30["QScan_SC"], ("QScan_SC", -30))]]
    c = clip_of(T)
    fig, axes = plt.subplots(2, 4, figsize=(14, 8), sharex=True, sharey=True, constrained_layout=True)
    for r, row in enumerate(panels):
        for k, (name, D, key) in enumerate(row):
            title = name
            if key:
                v = rel(D, T)
                ref = float(ROWS[(s, NOISE, key[1])][f"rel_{key[0]}"])
                assert abs(v - ref) < 1e-9 * max(1, ref), (s, key, v, ref)
                checked += 1
                title += f"\n相对误差 {v:.3f}"
            show(axes[r, k], D, t, clip_of(obs) if name == "观测" else c, title)
        axes[r, 0].set_ylabel(("Q 正确" if r == 0 else "Q −30%") + "\n时间 (s)", fontsize=9)
    for ax in axes[-1]:
        ax.set_xlabel("道号", fontsize=8)
    fig.suptitle(f"{s}（{NOISE} 噪声）：上排 Q 正确，下排 Q 低估 30%；除观测外共用真值目标的色标", fontsize=12)
    fig.savefig(OUT / f"{s}.png", dpi=150)
    plt.close(fig)
print(f"wrote {len(SCENES) + 1} figures, {checked} errors checked against rows.csv")
