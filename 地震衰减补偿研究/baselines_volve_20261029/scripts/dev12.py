"""Development of experiment 12 on the experiment-10 held-out datasets H (seen; also experiment 11's development set).

usage: dev12.py K N    worker K of N; writes analysis/dev/<id>.json (resumable)
Per dataset (input = noise realisation n0, as in experiments 10 and 11; c_SR = the experiment-10 SR factor):
  IQ_g<g>          Wang stabilised inverse Q, one factor c_SR, maximum gain g dB            (baseline B1)
  IQd_<k>_g<g>     same filter, amplitude on c_SR x k, dispersion on c_SR                   (decoupled IQ)
  VDB_<w>_<eps>    phase-only inverse Q on c_SR, then zero-phase time-varying Wiener shaping (baseline B3)
  SCa_a<a>         SC-a on c_SR with sparsity weight a (0.125 is the frozen SC-a)            (baseline B2)
  SRSCa_<k>+ph     SC-a on c_SR x k, then phase fix to c_SR (k = 1.3 is the experiment-11 method)
SC-a outputs with k = 1, 1.15, 1.3 are the cached experiment-11 development outputs.
"""
import sys, json, time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import lib12 as L  # noqa: E402
from lib12 import np, sc11, semi, lib  # noqa: E402

G_DB = [10, 20, 30, 40, 50, 60]
DEC_K, DEC_G = [1.3, 1.6, 2.0, 2.5], [30, 40, 50]
VDB_WIN, VDB_EPS = [0.3, 0.5], [1e-4, 3e-4, 1e-3, 3e-3, 1e-2, 3e-2, 0.1]
SC_A = [0.25, 0.5, 1.0]
K_NEW = [1.6, 2.0]
C = [0.55, 0.7, 0.9, 1.2, 1.4]
SEEDS = [9101, 9102, 9103, 9104]
OUT = L.ROOT / "analysis" / "dev"


def main(k, n):
    OUT.mkdir(parents=True, exist_ok=True); cm = semi.curtain_cm(); Q = np.asarray(cm["Q"], float)
    jobs = [f"H_c{c}_s{s}" for c in C for s in SEEDS]
    for i, did in enumerate(jobs):
        path = OUT / f"{did}.json"
        if i % n != k or path.exists():
            continue
        t0 = time.monotonic(); z = np.load(L.DEV_E10 / f"{did}.npz"); D = z["Input"].astype(float)
        ideal = z["ideal"].astype(float); zc = np.load(L.DEV_E11 / f"{did}.npz")
        r = json.loads((L.E10 / "results" / "synth" / f"{did}.json").read_text(encoding="utf-8"))
        c, c_sr = r["c_true"], r["sr"]["n0"]["c"]; p = sc11.prepare_sc(D, cm); wv = p["wv"]
        sc = {}
        for g in G_DB:
            sc[f"IQ_g{g}"] = L.score(L.iq(D, cm, Q, c_sr, c_sr, g), ideal)
        for kk in DEC_K:
            for g in DEC_G:
                sc[f"IQd_{kk}_g{g}"] = L.score(L.iq(D, cm, Q, c_sr * kk, c_sr, g), ideal)
        ph = L.iq(D, cm, Q, None, c_sr, 0)
        for w in VDB_WIN:
            for e in VDB_EPS:
                sc[f"VDB_{w}_{e:g}"] = L.score(L.tvw(ph, wv["w_out"], e, win=w), ideal)
        outs = {}
        for a in SC_A:
            outs[f"SCa_a{a}"] = lib.solve(cm, Q * c_sr, wv, D, p["sigma"], p["K"], p["live"], a, lib.SEL_SC[1])[0]
        for kk in K_NEW:
            outs[f"SCa_SRx{kk}"] = sc11.sca(D, cm, Q, c_sr * kk, p)
        sc["SCa_a0.125"] = L.score(zc["SCa_SRx1"].astype(float), ideal)
        for kk, key in ((1.0, "SCa_SRx1"), (1.15, "SCa_SRx1.15"), (1.3, "SCa_SRx1.3")):
            sc[f"SRSCa_{kk}+ph"] = L.score(sc11.phase_fix(zc[key].astype(float), cm, Q, c_sr * kk, c_sr), ideal)
        for kk in K_NEW:
            sc[f"SRSCa_{kk}+ph"] = L.score(sc11.phase_fix(outs[f"SCa_SRx{kk}"], cm, Q, c_sr * kk, c_sr), ideal)
        for a in SC_A:
            sc[f"SCa_a{a}"] = L.score(outs[f"SCa_a{a}"], ideal)
        sc["Input"] = L.score(D, ideal)
        (L.ROOT / "cache" / "dev").mkdir(parents=True, exist_ok=True)
        np.savez_compressed(L.ROOT / "cache" / "dev" / f"{did}.npz", **{m: o.astype(np.float32) for m, o in outs.items()})
        rec = dict(id=did, c_true=c, c_sr=c_sr, scores=sc, seconds=time.monotonic() - t0)
        path.write_text(json.dumps(rec, indent=1), encoding="utf-8"); print(did, f"{rec['seconds']:.0f}s", flush=True)


if __name__ == "__main__":
    main(int(sys.argv[1]), int(sys.argv[2]))
