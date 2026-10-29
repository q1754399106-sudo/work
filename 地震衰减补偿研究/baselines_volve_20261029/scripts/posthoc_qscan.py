"""Post-hoc addition to experiment 12 (not in PLAN_ZH.md, descriptive only): the user's pipeline QScan + SC-a
on the same held-out sections, so it sits in the same table as the baselines.

For each T dataset (input and ideal target from cache/out/<id>.npz, saved by run12.py as float32):
  QSCa       QScan factor (lib.qscan, unchanged) + SC-a with a = 0.125 (the pipeline of experiments 7-11)
  QSCa_tuned QScan factor + SC-a with the experiment-12 B2 weight (dev_selection.json)
  QSCa+ph    QSCa with the phase fix from the QScan factor to c_SR (experiment 11)
Input is rescored from the float32 copy as a check against results/synth/<id>.json.

usage: posthoc_qscan.py K N    worker K of N; writes results/posthoc_qscan/<id>.json
"""
import sys, json, time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import lib12 as L  # noqa: E402
from lib12 import np, sc11, semi, lib  # noqa: E402

OUT = L.ROOT / "results" / "posthoc_qscan"
SEL = json.loads((L.ROOT / "analysis" / "dev_selection.json").read_text(encoding="utf-8"))


def main(k, n):
    OUT.mkdir(parents=True, exist_ok=True); cm = semi.curtain_cm(); Q = np.asarray(cm["Q"], float)
    for i, p in enumerate(sorted((L.ROOT / "results" / "synth").glob("T_*.json"))):
        did = p.stem; path = OUT / f"{did}.json"
        if i % n != k or path.exists():
            continue
        t0 = time.monotonic(); r = json.loads(p.read_text(encoding="utf-8")); c_sr = r["sr"]["c"]
        z = np.load(L.ROOT / "cache" / "out" / f"{did}.npz"); D = z["Input"].astype(float); ideal = z["ideal"].astype(float)
        pp = sc11.prepare_sc(D, cm); fq, qi = sc11.qscan(D, cm, Q, pp); wv = pp["wv"]
        outs = dict(QSCa=sc11.sca(D, cm, Q, fq, pp))
        outs["QSCa_tuned"] = lib.solve(cm, Q * fq, wv, D, pp["sigma"], pp["K"], pp["live"], SEL["B2"]["a"],
                                       lib.SEL_SC[1])[0]
        outs["QSCa+ph"] = sc11.phase_fix(outs["QSCa"], cm, Q, fq, c_sr)
        rec = dict(id=did, c_true=r["c_true"], c_sr=c_sr, qscan=fq, qscan_boundary=bool(qi["at_boundary"]),
                   scores={m: L.score(o, ideal) for m, o in outs.items()},
                   input_check=dict(float32=L.score(D, ideal), run12=r["scores"]["Input"]), seconds=time.monotonic() - t0)
        np.savez_compressed(L.ROOT / "cache" / "out" / f"{did}__posthoc_qscan.npz",
                            **{m: o.astype(np.float32) for m, o in outs.items()})
        path.write_text(json.dumps(rec, indent=1), encoding="utf-8"); print(did, f"{rec['seconds']:.0f}s", flush=True)


if __name__ == "__main__":
    main(int(sys.argv[1]), int(sys.argv[2]))
