"""Scores of the frequency-domain Q-RTM cross-check (exploratory). Low-pass cutoff tuned at m = 1 per noise level."""
import json
from common import *
from kirch import ZONES, score
from run_kirch import LPS

d = np.load(CACHE / "rtm.npz"); img = d["img"].astype(float); ops = [str(o) for o in d["ops"]]
psig = (abs(np.load(CACHE / "data.npz")["visco"]) ** 2).mean()
post = lambda I: np.gradient(np.gradient(I, axis=1), axis=1)            # removes the low-wavenumber RTM artefact
tap = lambda fc: np.ones(len(FREQ)) if fc is None else 0.5 * (1 + np.cos(np.pi * np.clip((FREQ - fc) / 8.0, 0, 1)))
R = post(img[0, 0].sum(0)); res = {}
def image(name, fc, snr):
    o = ops.index(name); w = tap(fc)[:, None, None]; I = (img[o, 0] * w).sum(0)
    if snr is not None:
        I = I + np.sqrt(psig / 10 ** (snr / 10)) * (img[o, 1] * w).sum(0)
    return post(I)
for tag, snr in (("none", None), ("30dB", 30.0), ("14dB", 14.0)):
    grid = list(LPS) + [None]; sc = [score(image("m1_lam1", fc, snr), R)["target"]["E"] for fc in grid]
    best = grid[int(np.argmax(sc))]; res[tag] = dict(tune=dict(grid=grid, E=sc, best=best))
    print(tag, "LP tuning E(target):", " ".join("%s:%.3f" % z for z in zip(grid, sc)))
    for m in MS:
        a = score(image("m%g_lam0" % m, None, snr), R); b = score(image("m%g_lam1" % m, best, snr), R)
        res[tag]["m%g" % m] = dict(NONE=a, LP=b)
        print("  m=%-4g  no boost corr/E target %.3f/%.3f all %.3f/%.3f | Q-RTM + LP target %.3f/%.3f all %.3f/%.3f" % (
            m, a["target"]["corr"], a["target"]["E"], a["all"]["corr"], a["all"]["E"],
            b["target"]["corr"], b["target"]["E"], b["all"]["corr"], b["all"]["E"]))
json.dump(res, open(ROOT / "results" / "rtm.json", "w"), indent=1)
