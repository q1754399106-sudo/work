"""SC-a baseline (a 2^-13, mu_relative 0.25: experiments 12 and 17, tuned on this data family) on the 4 x noise version
of the semi-synthetic development set, SR route only (c_sr of pilot_semi.py). Exploratory, NOTE_ZH.md section 9.
usage: pilot_semi_sc.py run K N"""
import sys, json, time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import uqp  # noqa: E402,F401
from uqp import np, lib, semi, ROOT  # noqa: E402
from uqw import load  # noqa: E402
from lib17 import DEV, score, sc11, A_TUNED  # noqa: E402

OUT = ROOT / "results" / "pilot_semi_sc"
if __name__ == "__main__":
    k, n = int(sys.argv[2]), int(sys.argv[3]); OUT.mkdir(parents=True, exist_ok=True); cm = semi.curtain_cm(); Q = np.asarray(cm["Q"], float)
    for j, did in enumerate(DEV):
        path = OUT / f"{did}_nd.json"
        if j % n != k or path.exists():
            continue
        t0 = time.monotonic(); arr, meta = load(did); D = arr["signal"] + 4.0 * arr["noise"]; p = sc11.prepare_sc(D, cm)
        c_sr = json.loads((ROOT / "results" / "pilot_semi" / f"{did}_nd.json").read_text(encoding="utf-8"))["sr"]["c_sr"]
        X, info, _ = lib.solve(cm, Q * c_sr, p["wv"], D, p["sigma"], p["K"], p["live"], A_TUNED, 0.25)
        out = dict(id=did, c_sr=c_sr, SCa=float(score(X, arr["ideal"])["S"]), converged=bool(info["converged"]), seconds=time.monotonic() - t0)
        path.write_text(json.dumps(out, indent=1), encoding="utf-8"); print(out, flush=True)
