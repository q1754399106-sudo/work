"""Stage 0 (PLAN_STAGE0_ZH.md): loss of tuned compensation against a Q multiplier m (c = m x c_eff_cum).

usage: stage0.py smoke          IQ70 / SCa at c_SR on ZA_c0.75_r0 must reproduce experiment 17's dev scores
       stage0.py run K N        worker K of N over (dataset, method, m) jobs; results/stage0/<id>__<method>_m<m>.json
"""
import sys, json, time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]; BASE = ROOT.parent
sys.path.insert(0, str(BASE / "mu_tune_20261031" / "scripts"))
import lib17  # noqa: E402
from lib17 import np, lib, semi, sq, sc11, lib12, SEL, A_TUNED, G_IQ, DEV, load, score  # noqa: E402

MU = 0.25
M_ALL = (0.5, 0.6, 0.7, 0.85, 1.0, 1.2, 1.4, 1.7, 2.0); M_SC = (0.5, 0.7, 1.0, 1.4, 2.0)
OUT = ROOT / "results" / "stage0"


def setup():
    cm = semi.curtain_cm(); return cm, np.asarray(cm["Q"], float)


def compensate(method, D, cm, Q, c, p):
    if method == "IQ70":
        return lib12.iq(D, cm, Q, c, c, G_IQ)
    if method == "VDB":
        return lib12.tvw(lib12.iq(D, cm, Q, None, c, 0), p["wv"]["w_out"], SEL["B3"]["eps"], win=SEL["B3"]["win"])
    return lib.solve(cm, Q * c, p["wv"], D, p["sigma"], p["K"], p["live"], A_TUNED, MU)[0]


def jobs():
    J = []
    for did in DEV:                                     # dataset-major: a partial run covers whole datasets
        J += [(did, "SCa", m) for m in M_SC] + [(did, meth, m) for meth in ("IQ70", "VDB") for m in M_ALL]
    return J


if __name__ == "__main__":
    cm, Q = setup()
    if sys.argv[1] == "smoke":
        did = "ZA_c0.75_r0"; arr, meta = load(did); D = arr["D"]; p = sc11.prepare_sc(D, cm)
        ref = json.loads((BASE / "mu_tune_20261031" / "results" / "dev" / f"{did}.json").read_text(encoding="utf-8"))
        c_sr = ref["c_sr"]; t0 = time.monotonic()
        d_iq = score(compensate("IQ70", D, cm, Q, c_sr, p), arr["ideal"])["S"] - ref["scores"]["IQ70"]["S"]
        d_sc = score(compensate("SCa", D, cm, Q, c_sr, p), arr["ideal"])["S"] - ref["scores"]["SR_mu0.25"]["S"]
        print("smoke dS IQ70 %.2e  SCa %.2e  c_eff_cum %.4f  %.0fs" % (d_iq, d_sc, meta["c_eff_cum"], time.monotonic() - t0))
        sys.exit(0 if abs(d_iq) <= 1e-9 and abs(d_sc) <= 1e-3 else 1)  # PLAN_STAGE0_ZH.md: SCa tolerance revised pre-run
    k, n = int(sys.argv[2]), int(sys.argv[3]); OUT.mkdir(parents=True, exist_ok=True); cache = {}
    for i, (did, meth, m) in enumerate(jobs()):
        path = OUT / f"{did}__{meth}_m{m:g}.json"
        if i % n != k or path.exists():
            continue
        t0 = time.monotonic()
        if did not in cache:
            cache.clear(); arr, meta = load(did); cache[did] = (arr, meta, sc11.prepare_sc(arr["D"], cm))
        arr, meta, p = cache[did]; c = m * meta["c_eff_cum"]
        s = score(compensate(meth, arr["D"], cm, Q, c, p), arr["ideal"])
        rec = dict(id=did, method=meth, m=m, c=c, c_eff_cum=meta["c_eff_cum"], score=s, seconds=time.monotonic() - t0)
        path.write_text(json.dumps(rec, indent=1), encoding="utf-8")
        print(did, meth, m, "S %.4f" % s["S"], f"{rec['seconds']:.0f}s", flush=True)
