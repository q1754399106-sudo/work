"""PLAN_GATE_A2EB_ZH.md: empirical-Bayes shrinkage of the offset-spectral-ratio calibration towards the model Q.
u = 1/c_hat from the 50-CMP supergather (frozen A1b estimator); s = delete-one-group jackknife (10 groups of 5 CMPs);
u_EB = 1 + max(0, 1 - s^2/(u-1)^2)(u-1), gamma = 1 fixed. Missing s (a jackknife subset failed) -> no shrinkage.
usage: eb.py z            step 0 on used M24-M27 (estimator only) -> results/eb_z.json
       eb.py run <M> <fp> step 1 on new M32-M35                 -> results/a2eb_<M>_f<fp>.json"""
import sys, json, time
from pathlib import Path
from multiprocessing import Pool
import numpy as np
HERE = Path(__file__).resolve().parent; ROOT = HERE.parent
sys.path.insert(0, str(HERE))
import e2e as E  # noqa: E402  (sets the frozen estimator options)
s, o = E.s, E.o

NEW = {"M32": [55, 38, 80, 45], "M33": [85, 50, 65, 32], "M34": [40, 90, 48, 70], "M35": [65, 35, 100, 55]}
CS = (0.5, 0.7, 0.85, 1.0, 1.2, 1.4, 2.0); SEEDS = 6; BASE_NEW = 52_000_000; NG = 10


def c_hat_jk(Qm, rng_seed, c, fp, snr):
    o.CMPS = E.SUPER; a0m, axm, P = o.rays(Qm); rng = np.random.default_rng(rng_seed)
    A, Bg = o.reflectivity(rng, True); G0 = o.gathers(A, Bg, Qm, c, fp, P); rms = np.array([np.sqrt(np.mean(g ** 2)) for g in G0])
    nrng = np.random.default_rng(rng.integers(1 << 62) + snr); sig = rms / 10 ** (snr / 20)
    Gn = [g + q * nrng.standard_normal(g.shape) for g, q in zip(G0, sig)]
    ce, _ = o.c_off(Gn, sig, a0m, axm, o.vrms); uk = []
    for k in range(NG):
        m = np.ones(E.SUPER, bool); m[k * E.SUPER // NG:(k + 1) * E.SUPER // NG] = False
        ck, _ = o.c_off([g[:, m] for g in Gn], sig, a0m, axm, o.vrms); uk.append(1 / ck if (ck == ck and ck > 0) else np.nan)
    uk = np.array(uk); se = float(np.sqrt((NG - 1) / NG * np.sum((uk - uk.mean()) ** 2))) if np.all(np.isfinite(uk)) else float("nan")
    return float(ce), se


def shrink(ce, se):
    if not (ce == ce and ce > 0):
        return 1.0
    if not np.isfinite(se):
        return ce
    d = 1 / ce - 1
    f = max(0.0, 1 - se ** 2 / d ** 2) if d != 0 else 0.0
    return 1 / (1 + f * d)


def z_job(args):
    model, fp = args; Qm = np.array(E.MODELS[model], float); mnum = int(model[1:]); rows = []
    for c in (0.85, 1.0, 1.2):
        for seed in range(10):
            for snr in E.SNRS:
                sd = E.EST_BASE + mnum * 100_000 + int(fp) * 1000 + int(round(c * 10)) * 10 + seed   # e2e.c_hat rule
                ce, se = c_hat_jk(Qm, sd, c, fp, snr)
                rows.append(dict(model=model, fp=fp, c=c, seed=seed, snr=snr, c_hat=ce, se_u=se, z=(1 / ce - 1 / c) / se, c_eb=shrink(ce, se)))
    return rows


def run(model, fp):
    Qm = np.array(NEW[model], float); mnum = int(model[1:]); rows = []; t0 = time.monotonic(); est = {}
    s.g.RIC = s.spec(fp); s.th.S2 = s.g.RIC ** 2; out = ROOT / "results" / ("a2eb_%s_f%d.json" % (model, int(fp)))
    done = set()
    if out.exists():   # resume (added 11-21 after out-of-memory kills): keep only c with all scenes done; seeds depend on c only
        old = json.loads(out.read_text())
        done = {c for c in CS if {r["scene"] for r in old if r["c"] == c} == set(E.SCENES)}
        rows = [r for r in old if r["c"] in done]; print(model, fp, "resume, done c:", sorted(done), flush=True)
    for c in CS:
        if c in done:
            continue
        for seed in range(SEEDS):
            for snr in E.SNRS:
                ce, se = c_hat_jk(Qm, BASE_NEW + mnum * 100_000 + int(fp) * 1000 + int(round(c * 100)) * 10 + seed, c, fp, snr)
                est[(c, snr, seed)] = (ce, se, shrink(ce, se))
        print(model, fp, c, "c_hat done %.0fs" % (time.monotonic() - t0), flush=True)
        for scene in E.SCENES:
            D, Cl, Tr, sig, meta, keep = E.columns(c * Qm, mnum, scene, fp, range(SEEDS)); res = {}
            for name, Qop, meth in (("oracle", c * Qm, "sat"), ("uncal", Qm, "sat"), ("tik_uncal", Qm, "direct")):
                res[name] = E.solve(D, Cl, sig, keep, Qop, meth)
            for name, pick in (("cal", 0), ("cal_eb", 2)):
                O = np.zeros_like(D); Oc = np.zeros_like(D); lam = []
                for k, (snr, seed) in enumerate(meta):
                    ce = est[(c, snr, seed)][pick]; ce = ce if (ce == ce and ce > 0) else 1.0
                    if ce == 1.0:                       # identical to the uncalibrated column (per-column lambda choice)
                        O[:, :, k], Oc[:, :, k] = res["uncal"][0][:, :, k], res["uncal"][1][:, :, k]; lam.append(res["uncal"][2][k]); continue
                    if name == "cal_eb" and ce == est[(c, snr, seed)][0]:
                        O[:, :, k], Oc[:, :, k] = res["cal"][0][:, :, k], res["cal"][1][:, :, k]; lam.append(res["cal"][2][k]); continue
                    a, b, l = E.solve(D[:, :, k:k + 1], Cl[:, :, k:k + 1], sig[:, k:k + 1], keep, ce * Qm, "sat")
                    O[:, :, k], Oc[:, :, k] = a[:, :, 0], b[:, :, 0]; lam.append(l[0])
                res[name] = (O, Oc, lam)
            for name, (O, Oc, lam) in res.items():
                m = s.metrics(O, Oc, Tr, keep, fp)
                for k, (snr, seed) in enumerate(meta):
                    ce, se, cb = est[(c, snr, seed)]
                    rows.append(dict(model=model, fp=fp, scene=scene, c=c, snr=snr, seed=seed, method=name, c_hat=ce, se_u=se, c_eb=cb,
                                     lam=lam[k], **{q: float(v[k]) for q, v in m.items()}))
            print(model, fp, c, scene, "%.0fs" % (time.monotonic() - t0), flush=True)
            out.write_text(json.dumps(rows))


if __name__ == "__main__":
    if sys.argv[1] == "z":
        with Pool(4) as p:
            R = [r for rows in p.map(z_job, [(m, fp) for m in E.MODELS for fp in E.FPS]) for r in rows]
        z = np.array([r["z"] for r in R]); z = z[np.isfinite(z)]
        summ = dict(n=len(R), n_finite=int(len(z)), z_mean=float(np.mean(z)), z_sd=float(np.std(z)),
                    by_c={str(c): dict(z_mean=float(np.mean([r["z"] for r in R if r["c"] == c and np.isfinite(r["z"])])),
                                       z_sd=float(np.std([r["z"] for r in R if r["c"] == c and np.isfinite(r["z"])])),
                                       shrunk_to_1=float(np.mean([r["c_eb"] == 1.0 for r in R if r["c"] == c])),
                                       rms_ln_hard=float(np.sqrt(np.mean([np.log(r["c_hat"] / c) ** 2 for r in R if r["c"] == c]))),
                                       rms_ln_eb=float(np.sqrt(np.mean([np.log(r["c_eb"] / c) ** 2 for r in R if r["c"] == c]))))
                            for c in (0.85, 1.0, 1.2)})
        summ["pass"] = bool(0.7 <= summ["z_sd"] <= 1.5 and abs(summ["z_mean"]) <= 0.5)
        (ROOT / "results" / "eb_z.json").write_text(json.dumps(dict(summary=summ, rows=R), indent=1))
        print(json.dumps(summ, indent=1))
    elif sys.argv[1] == "run":
        run(sys.argv[2], float(sys.argv[3]))
