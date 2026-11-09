"""Stage 1 (PLAN_STAGE1_ZH.md): UQ-IQ against IQ70 / VDB / SC-a (stage 0) and two deterministic conservative controls
on the 24 development datasets, Q multiplier m (c = m x c_eff_cum).

usage: stage1.py run K N   worker K of N over datasets -> results/stage1/<id>.json (one record list per dataset)
"""
import sys, json, time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from uqiq import ROOT, BASE, np, lib12, semi, G_IQ, load, uqiq  # noqa: E402
from stage0 import DEV, M_ALL, score, sc11, SEL  # noqa: E402

SIG = (0.1225, 0.2235, 0.245, 0.447, 0.894)        # half A, half B, A, B, double B (log scale)
SIG_A, SIG_B = 0.245, 0.447                          # SD of log m over the scenario grids A and B
OUT = ROOT / "results" / "stage1"


def records(did, cm, Q):
    arr, meta = load(did); D, ideal = arr["D"], arr["ideal"]; ce = meta["c_eff_cum"]; R = []
    dev = json.loads((BASE / "mu_tune_20261031" / "results" / "dev" / f"{did}.json").read_text(encoding="utf-8"))
    add = lambda method, m, c, sig, out, **kw: R.append(dict(method=method, m=m, c=c, sig=sig, S=score(out, ideal)["S"], **kw))
    for m in M_ALL:
        c = m * ce; ref = lib12.iq(D, cm, Q, c, c, G_IQ)
        x0 = uqiq(D, cm, Q, c, 0.0)                                                    # P0 equivalence
        add("UQ0", m, c, 0.0, x0, rel_diff=float(np.abs(x0 - ref).max() / np.abs(ref).max()))
        for law in ("gauss", "lognormal"):
            for sig in SIG:
                add(f"UQ_{law}", m, c, sig, uqiq(D, cm, Q, c, sig, law=law))
        for sig in (SIG_A, SIG_B):
            add("UQ_gauss_abs", m, c, sig, uqiq(D, cm, Q, c, sig, law="gauss_abs"))
            c2 = c * np.exp(sig); add("C2", m, c, sig, lib12.iq(D, cm, Q, c2, c2, G_IQ))   # Q taken one sigma larger
        add("C1", m, c, None, lib12.iq(D, cm, Q, c, c, 30))                          # 30 dB stabilisation
    c = dev["c_sr"]; wv = sc11.prepare_sc(D, cm)["wv"]                                 # as stage 0's VDB
    add("SR_IQ70", None, c, 0.0, lib12.iq(D, cm, Q, c, c, G_IQ))
    add("SR_VDB", None, c, None, lib12.tvw(lib12.iq(D, cm, Q, None, c, 0), wv["w_out"], SEL["B3"]["eps"], win=SEL["B3"]["win"]))
    for law in ("gauss", "lognormal"):
        add(f"SR_UQ_{law}", None, c, SIG_A, uqiq(D, cm, Q, c, SIG_A, law=law))
    return R


if __name__ == "__main__" and sys.argv[1] == "run":
    k, n = int(sys.argv[2]), int(sys.argv[3]); OUT.mkdir(parents=True, exist_ok=True)
    cm = semi.curtain_cm(); Q = np.asarray(cm["Q"], float)
    for i, did in enumerate(DEV):
        path = OUT / f"{did}.json"
        if i % n != k or path.exists():
            continue
        t0 = time.monotonic(); R = records(did, cm, Q)
        path.write_text(json.dumps(dict(id=did, records=R, seconds=time.monotonic() - t0), indent=1), encoding="utf-8")
        print(did, len(R), f"{time.monotonic() - t0:.0f}s", flush=True)
