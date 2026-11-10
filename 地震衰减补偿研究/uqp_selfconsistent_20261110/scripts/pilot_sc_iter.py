"""Does the iteration cap depress SC-a at 4 x noise? (exploratory, NOTE_ZH.md section 12). Variant R1_sel of
pilot_sc_robust.py (notch-tolerant shallow band, a 2^-9, mu 2) with benchmark.MAXITER x 4, on 4 of the 16 evaluation
datasets (one per structure family with the default order: ZA, ZC, ZE, ZF at c 0.75). Records S, KKT ratio, iterations.

usage: pilot_sc_iter.py run K N
"""
import sys, json, time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import uqp  # noqa: E402,F401
from uqp import np, lib, semi, ROOT  # noqa: E402
from uqw import load  # noqa: E402
from lib17 import score, sc11  # noqa: E402
from pilot_sc_robust import wavelets  # noqa: E402

OUT = ROOT / "results" / "pilot_sc_iter"; IDS = ["ZA_c0.75_r1", "ZC_c0.75_r1", "ZE_c0.75_r0", "ZF_c0.75_r0"]; FACTOR = 4

if __name__ == "__main__":
    k, n = int(sys.argv[2]), int(sys.argv[3]); OUT.mkdir(parents=True, exist_ok=True); cm = semi.curtain_cm(); Q = np.asarray(cm["Q"], float)
    base = lib.bm.MAXITER; lib.bm.MAXITER = base * FACTOR
    for j, did in enumerate(IDS):
        path = OUT / f"{did}_x{FACTOR}.json"
        if j % n != k or path.exists():
            continue
        t0 = time.monotonic(); arr, meta = load(did); D = arr["signal"] + 4.0 * arr["noise"]; p = sc11.prepare_sc(D, cm)
        c = json.loads((ROOT / "results" / "pilot_semi" / f"{did}_nd.json").read_text(encoding="utf-8"))["sr"]["c_sr"]
        X, info, _ = lib.solve(cm, Q * c, wavelets(D, cm, False), D, p["sigma"], p["K"], p["live"], 2.0 ** -9, 2.0)
        ref = json.loads((ROOT / "results" / "pilot_sc_robust" / f"{did}_R1_sel.json").read_text(encoding="utf-8"))["S"]
        out = dict(id=did, maxiter=base * FACTOR, base_maxiter=base, S=float(score(X, arr["ideal"])["S"]), S_base=ref, iterations=int(info["iterations"]),
                   kkt_ratio=float(info["kkt_ratio"]), converged=bool(info["converged"]), objective=float(info["objective"]), seconds=time.monotonic() - t0)
        path.write_text(json.dumps(out, indent=1), encoding="utf-8"); print(out, flush=True)
