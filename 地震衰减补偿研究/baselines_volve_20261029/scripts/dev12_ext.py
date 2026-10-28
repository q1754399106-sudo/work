"""Grid extension of the experiment-12 development runs (PLAN_ZH.md section 5): the first pass put the optimum of
four families on the grid edge (IQ g = 60 dB, IQd g = 50 dB, VDB eps = 1e-4, SC-a a = 0.125).

usage: dev12_ext.py K N    worker K of N; adds the new variants to analysis/dev/<id>.json (resumable per dataset)
"""
import sys, json, time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import lib12 as L  # noqa: E402
from lib12 import np, sc11, semi, lib  # noqa: E402

G_DB = [70, 80, 100, 120]
DEC_K, DEC_G = [1.15, 1.3, 1.6], [60, 70, 80, 100, 120]
VDB_WIN, VDB_EPS = [0.3, 0.5], [3e-5, 1e-5, 3e-6, 1e-6]
SC_A = [0.0625, 0.03125]
OUT = L.ROOT / "analysis" / "dev"


def main(k, n):
    cm = semi.curtain_cm(); Q = np.asarray(cm["Q"], float)
    for i, path in enumerate(sorted(OUT.glob("H_*.json"))):
        rec = json.loads(path.read_text(encoding="utf-8")); did = rec["id"]
        if i % n != k or rec.get("ext"):
            continue
        t0 = time.monotonic(); z = np.load(L.DEV_E10 / f"{did}.npz"); D = z["Input"].astype(float)
        ideal = z["ideal"].astype(float); c_sr = rec["c_sr"]; p = sc11.prepare_sc(D, cm); wv = p["wv"]; sc = rec["scores"]
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
            sc[f"SCa_a{a}"] = L.score(outs[f"SCa_a{a}"], ideal)
        old = dict(np.load(L.ROOT / "cache" / "dev" / f"{did}.npz"))
        np.savez_compressed(L.ROOT / "cache" / "dev" / f"{did}.npz", **old, **{m: o.astype(np.float32) for m, o in outs.items()})
        rec["ext"] = True; rec["seconds_ext"] = time.monotonic() - t0
        path.write_text(json.dumps(rec, indent=1), encoding="utf-8"); print(did, f"{rec['seconds_ext']:.0f}s", flush=True)


if __name__ == "__main__":
    main(int(sys.argv[1]), int(sys.argv[2]))
