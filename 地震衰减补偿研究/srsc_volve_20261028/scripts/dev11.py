"""Development of experiment 11 on the experiment-10 held-out datasets H (seen; outputs cached in experiment 10).

usage: dev11.py K N    worker K of N; writes analysis/dev/<id>.json (resumable)
Per dataset (input = noise realisation n0, as in experiment 10): SC-a with c_true x M_TRUE (sensitivity of SC-a to
the Q error), SC-a with c_SR x K_SR (conservative shift of the SR factor), QScan + SC-a (QScan factor of experiment 10).
SC-a with c_SR x 1 is the cached SCa_SR output.
"""
import sys, json, time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import sc11  # noqa: E402
from sc11 import np, semi  # noqa: E402

M_TRUE = [0.75, 1.0, 1.33]
K_SR = [1.15, 1.3]
C = [0.55, 0.7, 0.9, 1.2, 1.4]
SEEDS = [9101, 9102, 9103, 9104]
OUT = sc11.ROOT / "analysis" / "dev"


def main(k, n):
    OUT.mkdir(parents=True, exist_ok=True); cm = semi.curtain_cm(); Q = np.asarray(cm["Q"], float)
    jobs = [f"H_c{c}_s{s}" for c in C for s in SEEDS]
    for i, did in enumerate(jobs):
        path = OUT / f"{did}.json"
        if i % n != k or path.exists():
            continue
        t0 = time.monotonic(); z = np.load(sc11.DEV_OUT / f"{did}.npz"); D = z["Input"].astype(float); ideal = z["ideal"].astype(float)
        r = json.loads((sc11.E10 / "results" / "synth" / f"{did}.json").read_text(encoding="utf-8"))
        c, c_sr, fq = r["c_true"], r["sr"]["n0"]["c"], r["qscan"]["factor"]; p = sc11.prepare_sc(D, cm)
        rec = dict(id=did, c_true=c, c_sr=c_sr, qscan=fq, scores={})
        outs = {"SCa_SRx1": z["SCa_SR"].astype(float)}
        for m in M_TRUE:
            outs[f"SCa_truex{m}"] = sc11.sca(D, cm, Q, c * m, p)
        for kk in K_SR:
            outs[f"SCa_SRx{kk}"] = sc11.sca(D, cm, Q, c_sr * kk, p)
        outs["QSCa"] = sc11.sca(D, cm, Q, fq, p)
        rec["scores"].update({m: sc11.m1(o, ideal) for m, o in outs.items()})
        (sc11.ROOT / "cache" / "dev").mkdir(parents=True, exist_ok=True)
        np.savez_compressed(sc11.ROOT / "cache" / "dev" / f"{did}.npz", **{m: o.astype(np.float32) for m, o in outs.items()})
        for key in ("IQ30_SR", "IQ30_true", "IQ30_model", "IQ30_QScan", "Input"):
            rec["scores"][key] = sc11.m1(z[key].astype(float), ideal)
        rec["seconds"] = time.monotonic() - t0
        path.write_text(json.dumps(rec, indent=1), encoding="utf-8"); print(did, f"{rec['seconds']:.0f}s", flush=True)


if __name__ == "__main__":
    main(int(sys.argv[1]), int(sys.argv[2]))
