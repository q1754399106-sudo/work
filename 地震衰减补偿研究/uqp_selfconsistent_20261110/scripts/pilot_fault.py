"""Does lateral mixing blur faults? (exploratory, NOTE_ZH.md section 11). Local scores of experiment 15 (lib15.masks
windows): correlation with the ideal section within 3 traces of a fault plane and at least 15 traces away, 2.0-3.0 s.
SR route on the 24 development datasets, both noise versions. Ratio fault / far shows any extra loss at faults.

usage: pilot_fault.py run K N
"""
import sys, json, time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import uqp, mwdip  # noqa: E402
from uqp import np, lib, semi, ROOT, DT  # noqa: E402
from uqw import uqw, snr_rows, G_IQ, TW, lib12, load  # noqa: E402
from lib17 import DEV, score  # noqa: E402

OUT = ROOT / "results" / "pilot_fault"; F_OP = np.fft.rfftfreq(1024, DT); S = semi.gain_sigma2(G_IQ); ANCHOR = 1.4


def local(o, ideal, mk):
    return {k: float(np.corrcoef(o[v], ideal[v])[0, 1]) for k, v in mk.items()}


def run(did, version, cm, Q):
    arr, meta = load(did); D = arr["signal"] + (1.0 if version == "cal" else 4.0) * arr["noise"]; ideal = arr["ideal"]
    c = json.loads((ROOT / "results" / "pilot_semi" / f"{did}_{version}.json").read_text(encoding="utf-8"))["sr"]["c_sr"]
    w = semi.win(TW, 2.0, 3.0)[:, None]; fd_ = arr["fault_dist"]; mk = dict(fault=(fd_ <= 3) & w, far=(fd_ >= 15) & w)
    x = lib.absorption(dict(cm, Q=Q * c), TW); dx = x - lib.absorption(dict(cm, Q=Q * c), [ANCHOR])[0]
    g = mwdip.Gather(D, TW, 4); g0 = mwdip.Gather(D, TW, 4, steer=False); inv = snr_rows(D, TW, F_OP)
    outs = dict(IQ70=lib12.iq(D, cm, Q, c, c, G_IQ), W0=uqw(D, cm, Q, c, 0.0, TW, ANCHOR, 0.0, inv_snr=inv),
                MW0_K1=g.apply(x, dx, 0.0, g.model(x, dx), F_OP, cm["fref"], S, 0),
                MW0_K9=g.apply(x, dx, 0.0, g.model(x, dx), F_OP, cm["fref"], S, 4),
                MW0_K9_nosteer=g0.apply(x, dx, 0.0, g0.model(x, dx), F_OP, cm["fref"], S, 4))
    return dict(id=did, version=version, c_sr=c, n_fault=int(mk["fault"].sum()), n_far=int(mk["far"].sum()),
                scores={k: dict(S=float(score(v, ideal)["S"]), **local(v, ideal, mk)) for k, v in outs.items()})


if __name__ == "__main__":
    k, n = int(sys.argv[2]), int(sys.argv[3]); OUT.mkdir(parents=True, exist_ok=True); cm = semi.curtain_cm(); Q = np.asarray(cm["Q"], float)
    for j, (d, v) in enumerate([(d, v) for v in ("cal", "nd") for d in DEV]):
        path = OUT / f"{d}_{v}.json"
        if j % n != k or path.exists():
            continue
        t0 = time.monotonic(); r = run(d, v, cm, Q); r["seconds"] = time.monotonic() - t0
        path.write_text(json.dumps(r, indent=1), encoding="utf-8")
        print(d, v, {m: {a: round(b, 3) for a, b in s.items()} for m, s in r["scores"].items()}, f"{r['seconds']:.0f}s", flush=True)
