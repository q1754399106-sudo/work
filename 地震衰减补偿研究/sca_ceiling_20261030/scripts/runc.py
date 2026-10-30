"""SC-a ceiling runs (PREDICTIONS_ZH.md). usage: runc.py K N  (worker K of N over the 18 Z datasets, resumable);
       runc.py ids <id> ...  (the listed datasets, resumable; added to re-split the work with fewer processes).
results/<id>.json; no output sections are kept."""
import sys, json, time
from libc import ROOT, IDS, META, E15, np, semi, sc11, load, score, wavelets, methods


def main(k, n, ids=IDS):
    out = ROOT / "results"; out.mkdir(exist_ok=True); cm = semi.curtain_cm(); Q = np.asarray(cm["Q"], float)
    for i, did in enumerate(ids):
        path = out / f"{did}.json"
        if i % n != k or path.exists():
            continue
        t0 = time.monotonic(); arr, meta = load(did); D = arr["D"]
        p = sc11.prepare_sc(D, cm); wv0 = p["wv"]; wv1 = wavelets(D, cm, True)
        outs = methods(D, cm, Q, meta["c_eff_cum"], p, wv0, wv1)
        rec = dict(id=did, c_eff_cum=meta["c_eff_cum"], f_true=meta["f_true"],
                   fs=wv0["fmax_shallow"], fs_robust=wv1["fmax_shallow"], fd=wv0["fmax_target"],
                   scores={m: score(o, arr["ideal"]) for m, o in outs.items()}, seconds=time.monotonic() - t0)
        f15 = E15 / "test" / f"{did}.json"
        if f15.exists():
            rec["e15_SCa_oracle"] = json.loads(f15.read_text(encoding="utf-8"))["scores"]["SCa_oracle"]["S"]
        path.write_text(json.dumps(rec, indent=1), encoding="utf-8"); print(did, f"{rec['seconds']:.0f}s", flush=True)


if __name__ == "__main__":
    if sys.argv[1] == "ids":
        main(0, 1, sys.argv[2:])
    else:
        main(int(sys.argv[1]), int(sys.argv[2]))
