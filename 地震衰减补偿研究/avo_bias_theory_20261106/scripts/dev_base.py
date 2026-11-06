"""Development on used M1 only: export, SC/SC0, score against the four linear methods. usage: dev_base.py export | sc | score"""
import json, sys, time
import numpy as np
import step2b as s
from step2b import NEW, columns, metrics, cell, spec
import gate1 as g
from gate1 import n, NX
import theory as th
import base_common as b

NEW["M1"] = np.array([40.0, 50, 70, 100]); s.SEEDS = 3; s.SNRS = (18, 12); CELLS = [("events_aligned", 1.0), ("dense_aligned", 1.0), ("dense_nmo", 1.0), ("events_aligned", 0.8)]; FP = 60.0; g.RIC = spec(FP); th.S2 = g.RIC ** 2


def name(sc, qf):
    return f"dev_M1_{sc}_{qf:.1f}"


if sys.argv[1] == "export":
    for sc, qf in CELLS:
        D, Cl, Tr, sig, meta, keep = columns("M1", sc, FP); b.write_cell(name(sc, qf), D, NEW["M1"] * qf, FP, sig); print(name(sc, qf), D.shape, flush=True)
elif sys.argv[1] == "sc":
    W_ = None
    for sc, qf in CELLS:
        D, Cl, Tr, sig, meta, keep = columns("M1", sc, FP); W, dicts = b.dictionaries(NEW["M1"] * qf, FP)
        for lab, a, mu in (("SC", 0.125, 10.0), ("SC0", 0.5, 0.0)):
            O = np.zeros_like(D); t0 = time.monotonic(); its = []
            for c in range(D.shape[2]):
                O[:, :, c], info = b.solve_sc(D[:, :, c], sig[:, c], W, dicts, a, mu); its.append((info["it"], info["conv"]))
            (b.EX / f"{name(sc, qf)}_{lab}.bin").write_bytes(np.asfortranarray(O).tobytes(order="F")); print(name(sc, qf), lab, f"{time.monotonic() - t0:.0f}s", its[:3], "all conv", all(x[1] for x in its), flush=True)
else:
    K = ("ga", "eg", "egl", "mg", "cg", "wave")
    for sc, qf in CELLS:
        D, Cl, Tr, sig, meta, keep = columns("M1", sc, FP); rows = []; t0 = time.monotonic(); cell("M1", sc, qf, FP, rows); print("\n", sc, qf, "linear methods", f"{time.monotonic() - t0:.0f}s")
        for m in ("sat", "split", "method", "direct"):
            x = [r for r in rows if r["m"] == m]; print(f"  {m:13s}" + " ".join(f"{k} {np.mean([r[k] for r in x]):.3f}" for k in K))
        for lab in ("FCR", "FrozenSplit", "FrozenDirect", "SC", "SC0"):
            try:
                O = b.read_out(name(sc, qf), lab, D.shape[2])
            except FileNotFoundError:
                print("  missing", lab); continue
            mt = metrics(O, O, Tr, keep, FP); print(f"  {lab:13s}" + " ".join(f"{k} {np.mean(mt[k]):.3f}" for k in K))
