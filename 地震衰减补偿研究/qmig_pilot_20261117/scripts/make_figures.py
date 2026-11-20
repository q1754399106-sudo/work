"""Figures (presentation only): model, Kirchhoff images at 14 dB noise, score curves."""
import json
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from common import *

M = model(); ext = [0, M["x"][-1], M["z"][-1], 0]; F = ROOT / "figures"
fig, ax = plt.subplots(1, 2, figsize=(11, 3.4))
for a, (v, t) in zip(ax, ((M["c"], "c0 (m/s)"), (1 / M["qi"], "Q"))):
    im = a.imshow(v.T, extent=ext, aspect="auto", cmap="viridis"); a.set_title(t); a.set_xlabel("x (m)"); a.set_ylabel("z (m)"); plt.colorbar(im, ax=a)
fig.tight_layout(); fig.savefig(F / "model.png", dpi=130); plt.close(fig)

K = np.load(CACHE / "kirch_images.npz"); R = K["REF"]; vm = 0.6 * abs(R[:, 25:]).max()
for tag in ("none", "14dB"):
    fig, ax = plt.subplots(3, 4, figsize=(15, 9), sharex=True, sharey=True)
    panels = [("REF", "reference (no attenuation)"), ("%s_NONE_m1" % tag, "no compensation"),
              ("%s_LP_m1" % tag, "LP, Q x 1"), ("%s_UQ_m1" % tag, "UQ, Q x 1"),
              ("%s_LP_m0.7" % tag, "LP, Q x 0.7"), ("%s_STAB_m0.7" % tag, "STAB, Q x 0.7"), ("%s_WIEN_m0.7" % tag, "WIEN, Q x 0.7"), ("%s_UQ_m0.7" % tag, "UQ, Q x 0.7"),
              ("%s_LP_m2" % tag, "LP, Q x 2"), ("%s_STAB_m2" % tag, "STAB, Q x 2"), ("%s_WIEN_m2" % tag, "WIEN, Q x 2"), ("%s_UQ_m2" % tag, "UQ, Q x 2")]
    for a, (k, t) in zip(ax.ravel(), panels):
        a.imshow(K[k].T, extent=ext, aspect="auto", cmap="gray", vmin=-vm, vmax=vm); a.set_title(t, fontsize=10)
        a.plot([848, 1568, 1568, 848, 848], [600, 600, 1304, 1304, 600], "y-", lw=0.6)
    fig.suptitle("Kirchhoff Q-PSDM, noise %s (same colour scale; box: target zone under the low-Q lens)" % tag)
    fig.tight_layout(); fig.savefig(F / ("kirch_images_%s.png" % tag), dpi=110); plt.close(fig)

r = json.load(open(ROOT / "results" / "kirch.json")); fig, ax = plt.subplots(1, 3, figsize=(14, 4), sharey=True)
for a, tag in zip(ax, ("none", "30dB", "14dB")):
    reals = [k for k in r[tag] if k.startswith("real")]
    for n, st in (("NONE", "k:"), ("LP", "C0o-"), ("STAB", "C1s-"), ("WIEN", "C2^-"), ("UQ", "C3D-"), ("UQ25", "C3D--")):
        a.plot(MS, [np.mean([r[tag][k]["%s_m%g" % (n, m)]["target"]["E"] for k in reals]) for m in MS], st, label=n)
    a.set_xscale("log"); a.set_xticks(MS); a.set_xticklabels([str(m) for m in MS]); a.set_ylim(-0.2, 1.02); a.grid(alpha=0.3)
    a.set_title("noise %s" % tag); a.set_xlabel("Q multiplier of the migration model")
ax[0].set_ylabel("E score, target zone"); ax[0].legend(fontsize=8)
fig.tight_layout(); fig.savefig(F / "kirch_scores.png", dpi=130); plt.close(fig)
