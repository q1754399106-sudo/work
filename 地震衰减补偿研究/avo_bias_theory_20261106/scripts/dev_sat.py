"""Development on used models only: frequency-saturated damping variants, own lambda by discrepancy, same data and scores as step 2.
usage: dev_sat.py M [quick] -> results/dev_sat_<M>.json. No rule, descriptive."""
import json, sys, time
from pathlib import Path
import numpy as np
import gate1 as g
from gate1 import C, NX, n, I, a0_ax, atten, avo, bandpass, MODELS, build_data, SIN2, DT, T
import step2 as s2
from step2 import score, GRID
from gate1b import template
import stable, satdamp as sd

OUT = Path(__file__).resolve().parents[1] / "results"
VARIANTS = {"sat.5": (0.5, None), "sat.2": (0.2, None), "sat.05": (0.05, None), "sat.2g6": (0.2, 6.0), "sat.05g6": (0.05, 6.0)}


def run_cell(M, scene, qf, rows):
    Q = MODELS[M] * qf; cols = build_data(M, scene); keep = cols[0][4]; Dat = np.stack([c[2] for c in cols], 2); Tru = np.stack([c[3] for c in cols], 2); sig = np.stack([c[5] for c in cols], 1); meta = [(int(c[0]), int(c[1])) for c in cols]; nc = Dat.shape[2]
    a0, ax = a0_ax(Q); Gd = [atten(a0 + ax[:, i]) for i in range(NX)]; GtG = [G.T @ G for G in Gd]; GtD = [Gd[i].T @ Dat[:, i, :] for i in range(NX)]; target = np.sum(sig ** 2, axis=0) * n; At, Bt, ok = avo(bandpass(Tru), keep)
    ev = scene.startswith("events")
    for v, (rc, gam) in VARIANTS.items():
        ch = [None] * nc; out = np.zeros_like(Dat); last = None
        for lam in GRID:
            fb = sd.fb_of(a0, lam, rc); wn = stable.bound(a0, ax, lam, gam) if gam else None; O = np.zeros_like(Dat)
            for i in range(NX):
                f = fb if wn is None else sd.fb_floor(fb, ax[:, i], wn[:, i]); L = sd.Lmat(ax[:, i], f); O[:, i, :] = np.linalg.solve(GtG[i] + lam * (L.T @ L), GtD[i])
            res = sum(np.sum(((Gd[i] @ O[:, i, :] - Dat[:, i, :]) * keep[:, i, None]) ** 2, axis=0) for i in range(NX)); last = O
            for c in range(nc):
                if ch[c] is None and res[c] <= target[c]:
                    ch[c] = float(lam); out[:, :, c] = O[:, :, c]
            if all(x is not None for x in ch):
                break
        fbk = sum(x is None for x in ch)
        for c in range(nc):
            if ch[c] is None:
                ch[c] = float(GRID[-1]); out[:, :, c] = last[:, :, c]
        s = score(out, Tru, keep, At, Bt, ok)
        for c in range(nc):
            r = dict(M=M, scene=scene, qf=qf, m=v, snr=meta[c][0], seed=meta[c][1], lam=ch[c], fallback=fbk, **{k: float(x[c]) for k, x in s.items()})
            if ev:
                fT = template(Tru[:, :, c], keep); fE = template(out[:, :, c], keep); r["tga"] = float(np.mean(np.abs(fE[:, 1] / fE[:, 0] - fT[:, 1] / fT[:, 0])))
            rows.append(r)


def main(M, quick):
    rows = []; t0 = time.monotonic(); path = OUT / (f"dev_sat_{M}" + ("_q" if quick else "") + ".json")
    scenes = ("events_aligned", "thin_aligned", "dense_aligned") if quick else g.SCENES
    for scene in scenes:
        for qf in ((1.0,) if quick else (0.8, 1.0, 1.2)):
            run_cell(M, scene, qf, rows); print(M, scene, qf, f"{time.monotonic() - t0:.0f}s", flush=True); path.write_text(json.dumps(rows), encoding="utf-8")


if __name__ == "__main__":
    main(sys.argv[1], len(sys.argv) > 2)
