"""figures/overview.png (real line, noise-free signal, data, true Q per structure, c0 = 0.75, r = 0) and
figures/tstar.png (true t* of Z against the model, every structure, c0 = 0.75); as figures.py of the first batch."""
import json
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "field_application_20261028" / "cache" / "deps"))   # matplotlib
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import fq2
from fq2 import fq, np, semi, ROOT
from figures import show
from fq import lib

DATA, FIG = ROOT / "cache" / "data", ROOT / "figures"
TW, NX, C = fq.TW, fq.NX, 0.75


def main():
    FIG.mkdir(exist_ok=True)
    geo = json.loads((ROOT / "analysis" / "geology.json").read_text(encoding="utf-8"))
    meta = json.loads((ROOT / "analysis" / "datasets.json").read_text(encoding="utf-8"))
    keys = list(fq2.STRUCT); fig, axs = plt.subplots(len(keys), 4, figsize=(15, 4.2 * len(keys)), sharey=True)
    for i, k in enumerate(keys):
        real = semi.read_lines([fq2.STRUCT[k]])[fq2.STRUCT[k]]
        did = f"Z{k}_c{C}_r0"; z = np.load(DATA / f"{did}.npz"); m = meta[did]
        show(axs[i, 0], real, f"{k}: real inline {fq2.STRUCT[k]}")
        show(axs[i, 1], z["signal"], f"Z{k} c0={C}: signal (noise-free)")
        show(axs[i, 2], z["D"], f"Z{k} c0={C}: data (noise line {m['noise_line']})")
        im = axs[i, 3].imshow(1 / z["invQ"], aspect="auto", cmap="viridis", vmin=15, vmax=200,
                              extent=[0, NX - 1, TW[-1], TW[0]])
        axs[i, 3].contour(np.arange(NX), TW, z["gas"], [0.5], colors="w", linewidths=1)
        axs[i, 3].set_title(f"true Q (c_eff cum {m['c_eff_cum']:.2f}, diff {m['c_eff_dif']:.2f})", fontsize=9)
        fig.colorbar(im, ax=axs[i, 3], fraction=0.04)
        for f in geo[k]["faults"]:
            for ax in axs[i, 2:]:
                ax.plot(fq.fault_x(f, TW), TW, "r--", lw=0.6)
        for ax in axs[i]:
            ax.set_xlim(0, NX - 1); ax.axvline(fq.CENTRE, color="c", lw=0.5)
        axs[i, 0].set_ylabel("two-way time (s)")
    for ax in axs[-1]:
        ax.set_xlabel("trace (xline 2100-2500, step 2)")
    fig.tight_layout(); fig.savefig(FIG / "overview.png", dpi=110); plt.close(fig)

    cm = semi.curtain_cm(); a = lib.absorption(cm, TW)
    fig, axs = plt.subplots(1, len(keys), figsize=(13, 4.5), sharey=True)
    for ax, k in zip(axs, keys):
        z = np.load(DATA / f"Z{k}_c{C}_r0.npz")
        for j, col in zip((20, 60, 100, 140, 180), plt.cm.plasma(np.linspace(0, 0.85, 5))):
            ax.plot(z["tstar"][:, j] * 1e3, TW, color=col, lw=1, label=f"Z trace {j}")
        ax.plot(a / C * 1e3, TW, "k:", lw=1, label="model / c0")
        ax.set_title(f"structure {k}, c0 = {C}"); ax.set_xlabel("two-way t* (ms)")
    axs[0].invert_yaxis(); axs[0].set_ylabel("two-way time (s)"); axs[-1].legend(fontsize=7)
    fig.tight_layout(); fig.savefig(FIG / "tstar.png", dpi=110); plt.close(fig)


if __name__ == "__main__":
    main()
