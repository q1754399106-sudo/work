"""PLAN_S2SRC_ZH.md: where does the ~80 % cancellation of the Q sensitivity on gate-A2 NMO gathers come from?
Noise-free data Cl of the events scenes do not depend on the seed, so per (model, fp, c, scene, operator Q, operator kind,
dispersion on/off) the clean output is solved once per lambda of the step-2b grid; cached per-column lambdas are looked up.
Original gate-A2 functions (e2e.solve, step2b.metrics, generator) are called unchanged; dispersion-off swaps in
amplitude-only copies of synth/atten.  usage: s2_source.py table [nproc] | analyse"""
import json, sys, time
from multiprocessing import Pool
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]; OQA = ROOT.parent / "offset_q_avo_20261118"
sys.path.insert(0, str(OQA / "scripts")); sys.path.insert(0, str(OQA / "scripts" / "legacy2b"))
import e2e  # noqa: E402
s = e2e.s
MODELS = e2e.MODELS; FPS = (40.0, 60.0); CS = (0.5, 0.7, 1.4, 2.0); PRIM = (0.7, 1.4); SC = ("events_aligned", "events_nmo")
GRID = s.GRID; TAB = ROOT / "results" / "s2_source_table.json"
_orig = dict(synth=s.synth, atten=s.atten)


def atten_amp(a):
    N = s.C.nfft; f = np.arange(N // 2 + 1) / (N * s.C.dt); H = np.exp(-np.pi * np.outer(f, a)); E = np.fft.rfft(np.eye(s.C.nt), N, axis=0)
    return np.fft.irfft(H * E, N, axis=0)[:s.C.nt]


def synth_amp(idx, R, Q, nmo):
    g = s.g; a0, ax = g.a0_ax(Q); b = -np.pi * g.FS
    d = np.zeros((g.NSYN, g.NX)); tru = np.zeros((g.n, g.NX)); e0 = np.exp(-2j * np.pi * np.outer(g.FS, g.T[idx]))
    for i in range(g.NX):
        tk = g.GEO[i]["tx"][idx] if nmo else g.T[idx]; a = a0[idx] + ax[idx, i]
        X = g.RIC[:, None] * np.exp(np.outer(b, a)) * np.exp(-2j * np.pi * np.outer(g.FS, tk)); d[:, i] = np.fft.irfft(X @ R[:, i], g.NSYN)
        tru[:, i] = np.fft.irfft((g.RIC[:, None] * e0) @ R[:, i], g.NSYN)[:g.n]
    return d, tru


def setdisp(on):
    s.synth = _orig["synth"] if on else synth_amp; s.atten = _orig["atten"] if on else atten_amp
    e2e.s.synth = s.synth; e2e.s.atten = s.atten


def solve_fixed(Cl, keep, Qop, kind, lam):
    """Clean data only: one lambda, no selection (e2e.solve with a one-value grid; residual target irrelevant)."""
    save = s.GRID; s.GRID = np.array([lam])
    try:
        O, _, _ = e2e.solve(Cl, Cl, np.zeros((s.NX, Cl.shape[2])), keep, Qop, kind)
    finally:
        s.GRID = save
    return O


def task(args):
    model, fp, c = args; t0 = time.monotonic(); s.g.RIC = s.spec(fp); s.th.S2 = s.g.RIC ** 2; Qm = np.array(MODELS[model], float); mnum = int(model[1:]); out = []
    for disp in (True, False):
        setdisp(disp); keeps = {}
        data = {}
        for sc in SC:
            _, Cl, Tr, _, _, keep = e2e.columns(c * Qm, mnum, sc, fp, [0]); data[sc] = (Cl[:, :, :1], Tr[:, :, :1]); keeps[sc] = keep
        for sc in SC:
            Cl, Tr = data[sc]
            for meth, Qop in (("oracle", c * Qm), ("uncal", Qm)):
                for kind in (("sat", "direct") if disp else ("sat",)):
                    for lam in GRID:
                        O = solve_fixed(Cl, keeps[sc], Qop, kind, float(lam))
                        mk = {"own": keeps[sc]} if sc == "events_nmo" else {"own": keeps[sc], "nmokeep": keeps["events_nmo"]}
                        for kn, kp in mk.items():
                            m = s.metrics(O, O, Tr, kp, fp)
                            out.append(dict(model=model, fp=fp, c=c, scene=sc, method=meth, kind=kind, disp=disp, lam=float(lam), keep=kn, **{f"gac{j}": float(m[f"gac{j}"][0]) for j in range(3)}))
    setdisp(True); print(model, fp, c, "%.0fs" % (time.monotonic() - t0), flush=True); return out


def table(nproc):
    jobs = [(m, fp, c) for m in MODELS for fp in FPS for c in CS]
    with Pool(nproc) as p:
        rows = sum(p.map(task, jobs, chunksize=1), [])
    TAB.write_text(json.dumps(rows), encoding="utf-8")


def slope(cs, y):
    x = np.log(1 / np.asarray(cs)); x = x - x.mean(); return float((x * (np.asarray(y) - np.mean(y))).sum() / (x ** 2).sum())


def analyse():
    T = json.loads(TAB.read_text()); L = {}
    for r in T:
        L[(r["model"], r["fp"], r["c"], r["scene"], r["method"], r["kind"], r["disp"], r["keep"], round(np.log10(r["lam"]), 3))] = [r[f"gac{j}"] for j in range(3)]
    cache = {}
    for m in MODELS:
        for fp in FPS:
            for r in json.loads((OQA / "results" / f"a2_{m}_f{int(fp)}.json").read_text()):
                if r["method"] in ("uncal", "oracle") and r["scene"] in SC and r["c"] in CS:
                    cache[(m, fp, r["scene"], r["c"], r["snr"], r["seed"], r["method"])] = r
    cols = sorted({k[4:6] for k in cache})
    lg = lambda v: round(float(np.log10(GRID[np.argmin(np.abs(np.log10(GRID) - np.log10(v)))])), 3)  # noqa: E731
    rep = []
    for k, r in cache.items():
        m, fp, sc, c, snr, seed, meth = k; g = L[(m, fp, c, sc, meth, "sat", True, "own", lg(r["lam"]))]; rep += [abs(g[j] - r[f"gac{j}"]) for j in range(3)]

    def phi(policy, kind="sat", disp=True, alkeep="own"):
        """policy: 'cache' | 'common' | float lambda.  Returns phi median and S medians over the 24 cells."""
        cells = []
        for m in MODELS:
            for fp in FPS:
                for j in range(3):
                    S = {}
                    for sc in SC:
                        kp = alkeep if sc == "events_aligned" else "own"; med = {}
                        for c in PRIM:
                            ds = []
                            for snr, seed in (cols if policy in ("cache", "common") else [(0, 0)]):
                                if policy == "cache":
                                    lu, lo = lg(cache[(m, fp, sc, c, snr, seed, "uncal")]["lam"]), lg(cache[(m, fp, sc, c, snr, seed, "oracle")]["lam"])
                                elif policy == "common":
                                    lu = lo = lg(cache[(m, fp, sc, c, snr, seed, "oracle")]["lam"])
                                else:
                                    lu = lo = lg(policy)
                                ds.append(L[(m, fp, c, sc, "uncal", kind, disp, kp, lu)][j] - L[(m, fp, c, sc, "oracle", kind, disp, kp, lo)][j])
                            med[c] = float(np.median(ds))
                        S[sc] = slope(PRIM, [med[c] for c in PRIM])
                    cells.append((S["events_aligned"], S["events_nmo"], 1 - S["events_nmo"] / S["events_aligned"]))
        a = np.array(cells); return dict(phi=float(np.median(a[:, 2])), S_al=float(np.median(a[:, 0])), S_nmo=float(np.median(a[:, 1])),
                                         phi_q25_75=[float(np.quantile(a[:, 2], q)) for q in (0.25, 0.75)])
    res = dict(R0_median_abs_diff=float(np.median(rep)), R0_p95=float(np.quantile(rep, 0.95)),
               V0=phi("cache"), V1=phi("common"), V5=phi("cache", alkeep="nmokeep"),
               V2={f"{l:.0e}": phi(float(l)) for l in GRID}, V3={f"{l:.0e}": phi(float(l), kind="direct") for l in GRID},
               V4={f"{l:.0e}": phi(float(l), disp=False) for l in GRID}, V4_cache=phi("cache", disp=False),
               V5_fixed={f"{l:.0e}": phi(float(l), alkeep="nmokeep") for l in GRID},
               V4V5_fixed={f"{l:.0e}": phi(float(l), disp=False, alkeep="nmokeep") for l in GRID},
               V4V5_cache=phi("cache", disp=False, alkeep="nmokeep"))
    lam_hist = {sc: {meth: np.unique([lg(r["lam"]) for k, r in cache.items() if k[2] == sc and k[6] == meth], return_counts=True) for meth in ("oracle", "uncal")} for sc in SC}
    res["lam_log10_counts"] = {sc: {m: dict(zip(map(str, v[0]), map(int, v[1]))) for m, v in d.items()} for sc, d in lam_hist.items()}
    (ROOT / "results" / "s2_source.json").write_text(json.dumps(res, indent=1), encoding="utf-8"); print(json.dumps(res, indent=1))


if __name__ == "__main__":
    if sys.argv[1] == "table":
        table(int(sys.argv[2]) if len(sys.argv) > 2 else 8)
    else:
        analyse()
