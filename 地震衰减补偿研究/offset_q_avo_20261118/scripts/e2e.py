"""Gate A2 (PLAN_GATE_A2_ZH.md): end-to-end AVO with the offset-spectral-ratio calibration.
Truth Q = c x model Q. For every (model, wavelet, c, SNR, seed): c_hat from a 50-CMP supergather with dense Shuey-AVO
reflectivity (oqa.py, frozen A1b estimator), then AVO gathers (step-2b generator, legacy2b/, unchanged) compensated by
  oracle  saturated damping with the true Q
  uncal   saturated damping with the model Q;   tik_uncal  one-step Tikhonov with the model Q
  cal     saturated damping with model Q x c_hat
all with lambda from the residual principle (step 2b loop, sigma known). Scores: step2b.metrics (per-event G/A 'ga', etc.).
usage: e2e.py run <model> <fp> | dev   (one process per model x wavelet)"""
import sys, json, time
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent; ROOT = HERE.parent
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(HERE / "legacy2b"))
import step2b as s  # noqa: E402
import oqa as o  # noqa: E402

o.EQLEN, o.SMOOTH, o.XWIN = True, 1, False  # frozen A1b estimator
MODELS = {"M24": [42, 70, 55, 95], "M25": [75, 45, 30, 65], "M26": [35, 40, 85, 50], "M27": [95, 60, 45, 35]}
CS = (0.5, 0.7, 1.0, 1.4, 2.0); SNRS = (24, 30); SEEDS = 10; FPS = (60.0, 40.0)
SCENES = ("events_aligned", "thin_aligned", "events_nmo"); EST_BASE = 51_000_000; SUPER = 50


def c_hat(Qm, mnum, c, fp, snr, seed):
    o.CMPS = SUPER; a0m, axm, P = o.rays(Qm)
    rng = np.random.default_rng(EST_BASE + mnum * 100_000 + int(fp) * 1000 + int(round(c * 10)) * 10 + seed)
    A, Bg = o.reflectivity(rng, True); G0 = o.gathers(A, Bg, Qm, c, fp, P); rms = np.array([np.sqrt(np.mean(g ** 2)) for g in G0])
    nrng = np.random.default_rng(rng.integers(1 << 62) + snr); sig = rms / 10 ** (snr / 20)
    ce, _ = o.c_off([g + q * nrng.standard_normal(g.shape) for g, q in zip(G0, sig)], sig, a0m, axm, o.vrms)
    return ce


def columns(Qtrue, mnum, scene, fp, seeds):
    """As step2b.columns with the truth Q passed in; column order (seed, snr)."""
    nmo = scene.endswith("nmo"); D, C0, Tr, Sg, meta = [], [], [], [], []; keep = None
    for seed in seeds:
        rng = np.random.default_rng(8_000_000 + int(fp) * 10000 + mnum * 1000 + s.SCENES.index(scene) * 100 + seed)
        idx, R = s.reflectors(scene, rng); dp, tru = s.synth(idx, R, Qtrue, nmo); rms = np.sqrt(np.sum(dp[:int(2.4 / s.DT)] ** 2, axis=0) / s.n); cl, keep = s.to_t0(dp, nmo)
        for snr in SNRS:
            sg = rms / 10 ** (snr / 20); y, _ = s.to_t0(dp + sg * rng.standard_normal(dp.shape), nmo); D.append(y); C0.append(cl); Tr.append(tru); Sg.append(sg); meta.append((int(snr), int(seed)))
    return np.stack(D, 2), np.stack(C0, 2), np.stack(Tr, 2), np.stack(Sg, 1), meta, keep


def solve(D, Cl, sig, keep, Qop, method):
    """step2b.cell selection loop for one operator Q and one method ('sat' or 'direct'); returns noisy, clean outputs, lambdas."""
    nc = D.shape[2]; a0, ax = s.a0_ax(np.asarray(Qop, float)); RH = np.concatenate([D, Cl], axis=2); n = s.n; I = s.I
    Gd = [s.atten(a0 + ax[:, i]) for i in range(s.NX)]; GtG = [G.T @ G for G in Gd]; GtR = [Gd[i].T @ RH[:, i, :] for i in range(s.NX)]
    target = np.sum(sig ** 2, axis=0) * n; ch = [None] * nc; out = np.zeros_like(RH); O = None
    for lam in s.GRID:
        if method == "direct":
            O = np.stack([np.linalg.solve(GtG[i] + lam * I, GtR[i]) for i in range(s.NX)], axis=1)
        else:
            wn = s.stable.bound(a0, ax, lam, s.GAMMA); O = np.zeros_like(RH)
            for i in range(s.NX):
                L = s.sd.Lmat(ax[:, i], s.sd.fb_floor(np.full(n, s.sd.FF[-1]), ax[:, i], wn[:, i])); O[:, i, :] = np.linalg.solve(GtG[i] + lam * (L.T @ L), GtR[i])
        res = sum(np.sum(((Gd[i] @ O[:, i, :nc] - D[:, i, :]) * keep[:, i, None]) ** 2, axis=0) for i in range(s.NX))
        for c in range(nc):
            if ch[c] is None and res[c] <= target[c]:
                ch[c] = float(lam); out[:, :, c] = O[:, :, c]; out[:, :, nc + c] = O[:, :, nc + c]
        if all(v is not None for v in ch):
            break
    for c in range(nc):
        if ch[c] is None:
            ch[c] = float(s.GRID[-1]); out[:, :, c] = O[:, :, c]; out[:, :, nc + c] = O[:, :, nc + c]
    return out[:, :, :nc], out[:, :, nc:], ch


def run(model, seeds=range(SEEDS), cs=CS, fps=FPS, scenes=SCENES, Qm=None, mnum=None, tag="a2"):
    Qm = np.array(MODELS[model] if Qm is None else Qm, float); mnum = int(model[1:]) if mnum is None else mnum
    rows = []; t0 = time.monotonic(); chats = {}
    for fp in fps:
        s.g.RIC = s.spec(fp); s.th.S2 = s.g.RIC ** 2
        for c in cs:
            for seed in seeds:
                for snr in SNRS:
                    chats[(fp, c, snr, seed)] = c_hat(Qm, mnum, c, fp, snr, seed)
            print(model, fp, c, "c_hat done %.0fs" % (time.monotonic() - t0), flush=True)
            for scene in scenes:
                D, Cl, Tr, sig, meta, keep = columns(c * Qm, mnum, scene, fp, seeds); res = {}
                for name, Qop, meth in (("oracle", c * Qm, "sat"), ("uncal", Qm, "sat"), ("tik_uncal", Qm, "direct")):
                    res[name] = solve(D, Cl, sig, keep, Qop, meth)
                Ocal = np.zeros_like(D); Occal = np.zeros_like(D); lcal = []
                for k, (snr, seed) in enumerate(meta):
                    ce = chats[(fp, c, snr, seed)]; ce = ce if (ce == ce and ce > 0) else 1.0  # failed estimate: fall back to the model
                    a, b, l = solve(D[:, :, k:k + 1], Cl[:, :, k:k + 1], sig[:, k:k + 1], keep, ce * Qm, "sat"); Ocal[:, :, k] = a[:, :, 0]; Occal[:, :, k] = b[:, :, 0]; lcal.append(l[0])
                res["cal"] = (Ocal, Occal, lcal)
                for name, (O, Oc, lam) in res.items():
                    m = s.metrics(O, Oc, Tr, keep, fp)
                    for k, (snr, seed) in enumerate(meta):
                        rows.append(dict(model=model, fp=fp, scene=scene, c=c, snr=snr, seed=seed, method=name, c_hat=chats[(fp, c, snr, seed)],
                                         lam=lam[k], **{q: float(v[k]) for q, v in m.items()}))
                print(model, fp, c, scene, "%.0fs" % (time.monotonic() - t0), flush=True)
                (ROOT / "results" / ("%s_%s_f%d.json" % (tag, model, int(fps[0])) if len(fps) == 1 else "%s_%s.json" % (tag, model))).write_text(json.dumps(rows))


if __name__ == "__main__":
    if sys.argv[1] == "run":
        run(sys.argv[2], fps=(float(sys.argv[3]),))
    elif sys.argv[1] == "dev":  # used model M12 (step 2b), 3 seeds, c 0.5/1/2, 60 Hz, events_aligned only
        SNRS = (24,)
        run("M12", seeds=range(3), cs=(0.5, 1.0, 2.0), fps=(60.0,), scenes=("events_aligned",), Qm=s.NEW["M12"], mnum=12, tag="dev_a2")
