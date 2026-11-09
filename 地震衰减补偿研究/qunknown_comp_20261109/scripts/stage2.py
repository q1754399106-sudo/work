"""Stage 2 (PLAN_STAGE2_ZH.md): confirm UQ-IQ (lognormal, scenario B, sig 0.447) on the sealed fourth batch.

usage: stage2.py run K N   worker K of N over the 12 datasets -> results/stage2/<id>.json
Run through lowprio.py with OMP/OPENBLAS/MKL_NUM_THREADS=1 (low resource use).
"""
import sys, json, time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from uqiq import ROOT, BASE, np, lib12, semi, G_IQ, uqiq  # noqa: E402
from stage0 import M_ALL, score, sc11, SEL  # noqa: E402

DATA4 = BASE / "semisynth_fault_q4_20261109"
META4 = json.loads((DATA4 / "analysis" / "datasets.json").read_text(encoding="utf-8"))
IDS = [f"Z{k}_c{c}_r{r}" for k in "JKL" for c in (0.75, 1.05) for r in (0, 1)]
SIG_B = 0.447
SIG_LOGN = (0.245, 0.447, 0.894)                        # 0.447 decides; 0.245 and 0.894 descriptive
OUT = ROOT / "results" / "stage2"


def load(did):
    return {k: v.astype(float) for k, v in np.load(DATA4 / "cache" / "data" / f"{did}.npz").items()}, META4[did]


def records(did, cm, Q):
    arr, meta = load(did); D, ideal = arr["D"], arr["ideal"]; ce = meta["c_eff_cum"]; R = []
    wv = sc11.prepare_sc(D, cm)["wv"]                                              # as stage 0's VDB
    add = lambda method, m, c, sig, out, **kw: R.append(dict(method=method, m=m, c=c, sig=sig, S=score(out, ideal)["S"], **kw))
    for m in M_ALL:
        c = m * ce; ref = lib12.iq(D, cm, Q, c, c, G_IQ); add("IQ70", m, c, None, ref)
        x0 = uqiq(D, cm, Q, c, 0.0, law="lognormal")
        add("UQ0", m, c, 0.0, x0, rel_diff=float(np.abs(x0 - ref).max() / np.abs(ref).max()))
        add("VDB", m, c, None, lib12.tvw(lib12.iq(D, cm, Q, None, c, 0), wv["w_out"], SEL["B3"]["eps"], win=SEL["B3"]["win"]))
        add("C1", m, c, None, lib12.iq(D, cm, Q, c, c, 30))
        c2 = c * np.exp(SIG_B); add("C2", m, c, SIG_B, lib12.iq(D, cm, Q, c2, c2, G_IQ))
        for sig in SIG_LOGN:
            add("UQ_lognormal", m, c, sig, uqiq(D, cm, Q, c, sig, law="lognormal"))
        add("UQ_gauss", m, c, SIG_B, uqiq(D, cm, Q, c, SIG_B, law="gauss"))         # descriptive
    return R


if __name__ == "__main__" and sys.argv[1] == "run":
    k, n = int(sys.argv[2]), int(sys.argv[3]); OUT.mkdir(parents=True, exist_ok=True)
    cm = semi.curtain_cm(); Q = np.asarray(cm["Q"], float)
    for i, did in enumerate(IDS):
        path = OUT / f"{did}.json"
        if i % n != k or path.exists():
            continue
        t0 = time.monotonic(); R = records(did, cm, Q)
        path.write_text(json.dumps(dict(id=did, records=R, seconds=time.monotonic() - t0), indent=1), encoding="utf-8")
        print(did, len(R), f"{time.monotonic() - t0:.0f}s", flush=True)
