"""PLAN_ANGLE35_ZH.md: rerun of gate A2 (M24-M27) and gate A3 (M28-M31) with the AVO fit / G/A score restricted to traces whose
incidence angle is <= 35 deg. Data generators, seeds, c_hat estimator, compensation and lambda rule are those of e2e.py / a3.py,
unchanged. Every output is scored twice: with the original mask ('full', must reproduce the stored results) and with the angle
mask ('a35', keys suffixed 35).
Angle: physical incidence at the reflector, theta_j(x) = arcsin(p_j(x) vp[j]) (the Zoeppritz input of legacy2b/gate1.reflectors);
each t0 sample takes the offsets allowed for its nearest event (Voronoi split between event times).
Variant A (main): c_hat unchanged (all offsets).  Variant B (secondary, gate A2 aligned scenes only): c_hat from the spectral
ratio using only offset traces with theta(tc, x) <= 35 deg at the window centre tc; method name 'cal_B'.
usage: angle35.py a2 <M> <fp> | a3 <M> <fp> | mask"""
import sys, json, time
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent; sys.path.insert(0, str(HERE)); sys.path.insert(0, str(HERE / "legacy2b"))
import e2e  # noqa: E402
from e2e import s, o  # noqa: E402

THMAX = 35.0; ROOT = e2e.ROOT; SEEDS2 = range(6)  # gate A2 rerun: seeds 0-5 of the original 10 (cost), paired with the stored rows
g = s.g; EI = np.round(g.C.event_t / g.DT).astype(int)
THETA = np.degrees(np.arcsin(np.stack([[g.GEO[i]["p"][i0] * g.C.vp[j] for i in range(g.NX)] for j, i0 in enumerate(EI)])))  # (3, NX)
mid = 0.5 * (g.C.event_t[1:] + g.C.event_t[:-1]); NEAR = np.searchsorted(mid, g.T)  # nearest event per t0 sample
MASK35 = THETA[NEAR] <= THMAX  # (n, NX)


def scores(O, Oc, Tr, keep, fp):
    a = s.metrics(O, Oc, Tr, keep, fp); b = s.metrics(O, Oc, Tr, keep & MASK35, fp)
    return {**{k: v for k, v in a.items()}, **{k + "35": v for k, v in b.items()}}


# Variant B estimator: o.c_off with one extra condition (angle at the window centre); everything else copied unchanged.
_, _, PRAY = o.rays(np.array(e2e.MODELS["M24"], float))  # ray parameters do not depend on Q
THX = np.degrees(np.arcsin(np.minimum(PRAY * o.vint[:, None], 1.0)))  # (N, NX) angle at t0 with the local interval velocity


def c_off35(G, sig, a0p, axp, vp):
    T, OFF, NX = o.T, o.OFF, o.NX; xs, ys = [], []
    for a, b in o.WIN:
        R0, e = o.power(G[0], a, b); n0 = sig[0] ** 2 * e; P0 = np.maximum(R0 - n0, 1e-12 * R0.max()); tc = 0.5 * (a + b)
        for i in range(1, NX):
            if np.interp(tc, T, THX[:, i]) > THMAX:
                continue
            ta, tm, tb = np.sqrt(np.array([a, tc, b]) ** 2 + (OFF[i] / np.interp([a, tc, b], T, vp)) ** 2)
            if tb > T[-1] - 0.1 or tm / tc > o.MUTE:
                continue
            if o.EQLEN:
                ta, tb = tm - (b - a) / 2, tm + (b - a) / 2
            Rx, ex = o.power(G[i], ta, tb); nx = sig[i] ** 2 * ex; Px = np.maximum(Rx - nx, 1e-12 * Rx.max())
            ok = o.B & (R0 > 3 * n0) & (Rx > 3 * nx); ok &= o.f < (o.f[o.B][np.argmin(ok[o.B])] if not ok[o.B].all() else 1e9)
            if ok.sum() < 12:
                continue
            xs.append(np.pi * o.wmean(axp[:, i], a, b)); ys.append(-np.polyfit(o.f[ok], 0.5 * np.log(Px[ok] / P0[ok]), 1)[0])
    xs, ys = np.array(xs), np.array(ys)
    return (float(np.sum(xs ** 2) / np.sum(xs * ys)) if len(xs) and np.sum(xs * ys) > 0 else float("nan")), len(xs)


def c_hat_both(Qm, mnum, c, fp, snr, seed):
    """e2e.c_hat (same seeds, same supergather), returning the original estimate and the variant-B estimate."""
    o.CMPS = e2e.SUPER; a0m, axm, P = o.rays(Qm)
    rng = np.random.default_rng(e2e.EST_BASE + mnum * 100_000 + int(fp) * 1000 + int(round(c * 10)) * 10 + seed)
    A, Bg = o.reflectivity(rng, True); G0 = o.gathers(A, Bg, Qm, c, fp, P); rms = np.array([np.sqrt(np.mean(x ** 2)) for x in G0])
    nrng = np.random.default_rng(rng.integers(1 << 62) + snr); sig = rms / 10 ** (snr / 20)
    G = [x + q * nrng.standard_normal(x.shape) for x, q in zip(G0, sig)]
    ca, _ = o.c_off(G, sig, a0m, axm, o.vrms); cb, nb = c_off35(G, sig, a0m, axm, o.vrms)
    return ca, cb, nb


def ok_c(ce):
    return ce if (ce == ce and ce > 0) else 1.0


def run_a2(M, fp):
    Qm = np.array(e2e.MODELS[M], float); mnum = int(M[1:]); rows = []; t0 = time.monotonic()
    s.g.RIC = s.spec(fp); s.th.S2 = s.g.RIC ** 2; out = ROOT / "results" / ("angle35_a2_%s_f%d.json" % (M, int(fp)))
    for c in e2e.CS:
        ch = {(snr, seed): c_hat_both(Qm, mnum, c, fp, snr, seed) for seed in SEEDS2 for snr in e2e.SNRS}
        print(M, fp, c, "c_hat done %.0fs" % (time.monotonic() - t0), flush=True)
        for scene in e2e.SCENES:
            D, Cl, Tr, sig, meta, keep = e2e.columns(c * Qm, mnum, scene, fp, SEEDS2); res = {}
            for name, Qop, meth in (("oracle", c * Qm, "sat"), ("uncal", Qm, "sat"), ("tik_uncal", Qm, "direct")):
                res[name] = e2e.solve(D, Cl, sig, keep, Qop, meth)
            variants = (("cal", 0), ("cal_B", 1)) if scene.endswith("aligned") else (("cal", 0),)
            for name, vi in variants:
                Oc1 = np.zeros_like(D); Oc2 = np.zeros_like(D); lam = []
                for k, (snr, seed) in enumerate(meta):
                    ce = ok_c(ch[(snr, seed)][vi]); a, b, l = e2e.solve(D[:, :, k:k + 1], Cl[:, :, k:k + 1], sig[:, k:k + 1], keep, ce * Qm, "sat")
                    Oc1[:, :, k] = a[:, :, 0]; Oc2[:, :, k] = b[:, :, 0]; lam.append(l[0])
                res[name] = (Oc1, Oc2, lam)
            for name, (O, Oc, lam) in res.items():
                m = scores(O, Oc, Tr, keep, fp)
                for k, (snr, seed) in enumerate(meta):
                    cA, cB, nB = ch[(snr, seed)]
                    rows.append(dict(model=M, fp=fp, scene=scene, c=c, snr=snr, seed=seed, method=name, c_hat=cA, c_hat_B=cB, nwin_B=nB,
                                     lam=lam[k], **{q: float(v[k]) for q, v in m.items()}))
            print(M, fp, c, scene, "%.0fs" % (time.monotonic() - t0), flush=True); out.write_text(json.dumps(rows))


def run_a3(M, fp):
    import a3  # sets e2e.EST_BASE = 61 000 000 on import, as in the original gate A3 run
    import bandlim
    Qm = np.array(a3.MODELS[M], float); mnum = int(M[1:]); rows = []; t0 = time.monotonic(); out = ROOT / "results" / ("angle35_a3_%s_f%d.json" % (M, int(fp)))
    s.g.RIC = s.spec(fp); s.th.S2 = s.g.RIC ** 2
    for c in e2e.CS:
        ch = {(snr, seed): e2e.c_hat(Qm, mnum, c, fp, snr, seed) for seed in a3.SEEDS for snr in e2e.SNRS}
        for scene in a3.SCENES:
            D, Cl, Tr, sig, meta, keep = e2e.columns(c * Qm, mnum, scene, fp, a3.SEEDS)
            outs = {m: (np.zeros_like(D), np.zeros_like(D)) for m in ("sat", "direct", "split", "tik0", "band3")}
            for k, (snr, seed) in enumerate(meta):
                Qc = ok_c(ch[(snr, seed)]) * Qm; sl = slice(k, k + 1)
                a, b, _ = e2e.solve(D[:, :, sl], Cl[:, :, sl], sig[:, sl], keep, Qc, "sat"); outs["sat"][0][:, :, k] = a[:, :, 0]; outs["sat"][1][:, :, k] = b[:, :, 0]
                for m in ("direct", "split", "tik0"):
                    a, b = a3.solve_generic(D[:, :, sl], Cl[:, :, sl], sig[:, sl], keep, Qc, m); outs[m][0][:, :, k] = a[:, :, 0]; outs[m][1][:, :, k] = b[:, :, 0]
                outs["band3"][0][:, :, k] = bandlim.apply(D[:, :, k], Qc, sig[:, k], "snr", 3.0); outs["band3"][1][:, :, k] = bandlim.apply(Cl[:, :, k], Qc, sig[:, k], "snr", 3.0)
            for m, (O, Oc) in outs.items():
                q = scores(O, Oc, Tr, keep, fp)
                for k, (snr, seed) in enumerate(meta):
                    rows.append(dict(model=M, fp=fp, scene=scene, c=c, snr=snr, seed=seed, method=m, c_hat=ch[(snr, seed)], **{a: float(v[k]) for a, v in q.items()}))
            print(M, fp, c, scene, "%.0fs" % (time.monotonic() - t0), flush=True); out.write_text(json.dumps(rows))


if __name__ == "__main__":
    if sys.argv[1] == "mask":  # geometry only, no scores
        print("theta (deg) per event x offset\n", np.round(THETA, 1)); print("offsets kept per event", (THETA <= THMAX).sum(1))
        print("variant B: offsets with theta <= 35 at window centres", [int(np.sum(np.interp(0.5 * (a + b), o.T, THX[:, i]) <= THMAX for i in range(1, o.NX))) for a, b in o.WIN])
    elif sys.argv[1] == "a2":
        run_a2(sys.argv[2], float(sys.argv[3]))
    else:
        run_a3(sys.argv[2], float(sys.argv[3]))
