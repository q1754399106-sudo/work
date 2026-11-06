"""Development on used M1 only, second round: Q wrong and wavelet error. usage: dev_base2.py export | sc I | score"""
import sys, time
import numpy as np
import step2b as s
from step2b import NEW, columns, metrics, cell, spec
import gate1 as g
import theory as th
import base_common as b

NEW["M1"] = np.array([40.0, 50, 70, 100]); s.SEEDS = 3; FP = 60.0; g.RIC = spec(FP); th.S2 = g.RIC ** 2
CELLS = [("dense_aligned", 0.8, 1.0), ("dense_aligned", 1.2, 1.0), ("thin_aligned", 1.0, 1.0), ("events_aligned", 1.0, 0.8), ("events_aligned", 1.0, 1.2), ("dense_aligned", 1.0, 0.8), ("dense_aligned", 1.0, 1.2), ("thin_aligned", 1.0, 0.8)]


def name(sc, qf, ff):
    return f"dev2_M1_{sc}_{qf:.1f}_{ff:.1f}"


if sys.argv[1] == "export":
    for sc, qf, ff in CELLS:
        D, Cl, Tr, sig, meta, keep = columns("M1", sc, FP); b.write_cell(name(sc, qf, ff), D, NEW["M1"] * qf, FP * ff, sig)
elif sys.argv[1] == "sc":
    for sc, qf, ff in CELLS[int(sys.argv[2])::2]:
        D, Cl, Tr, sig, meta, keep = columns("M1", sc, FP); W, dicts = b.dictionaries(NEW["M1"] * qf, FP * ff)
        for lab, a, mu in (("SC", 0.125, 10.0), ("SC0", 0.5, 0.0)):
            O = np.zeros_like(D); t0 = time.monotonic(); its = []
            for c in range(D.shape[2]):
                O[:, :, c], info = b.solve_sc(D[:, :, c], sig[:, c], W, dicts, a, mu); its.append(info["conv"])
            (b.EX / f"{name(sc, qf, ff)}_{lab}.bin").write_bytes(np.asfortranarray(O).tobytes(order="F")); print(name(sc, qf, ff), lab, f"{time.monotonic() - t0:.0f}s", "conv", sum(its), "of", len(its), flush=True)
else:
    K = ("ga", "eg", "egl", "mg", "cg", "wave")
    for sc, qf, ff in CELLS:
        D, Cl, Tr, sig, meta, keep = columns("M1", sc, FP); rows = []; print("\n", sc, "Q x", qf, "dict f x", ff)
        if ff == 1.0 or sc:
            cell("M1", sc, qf, FP, rows)
            for m in ("sat", "split", "method", "direct"):
                x = [r for r in rows if r["m"] == m]; print(f"  {m:13s}" + " ".join(f"{k} {np.mean([r[k] for r in x]):.3f}" for k in K))
        for lab in ("FCR", "FrozenSplit", "SC", "SC0"):
            try:
                O = b.read_out(name(sc, qf, ff), lab, D.shape[2])
            except FileNotFoundError:
                print("  missing", lab); continue
            mt = metrics(O, O, Tr, keep, FP); print(f"  {lab:13s}" + " ".join(f"{k} {np.mean(mt[k]):.3f}" for k in K))
