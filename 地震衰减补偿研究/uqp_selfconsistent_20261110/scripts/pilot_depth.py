"""Where does the multichannel gain come from? (exploratory diagnostic, NOTE_ZH.md section 11). SR route, 24 development
datasets, both noise versions. Per 0.4 s depth window: median over traces of the zero-lag trace correlation with the
ideal section (the M1 part of S, restricted to the window) and the pooled correlation of the window (amplitude
sensitive, as the fault scores); plus the rms lateral amplitude error: std over traces of log(rms output / rms ideal).

usage: pilot_depth.py run K N
"""
import sys, json, time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import uqp, mwdip  # noqa: E402
from uqp import np, lib, semi, ROOT, DT  # noqa: E402
from uqw import uqw, snr_rows, G_IQ, TW, lib12, load  # noqa: E402
from lib17 import DEV  # noqa: E402

OUT = ROOT / "results" / "pilot_depth"; F_OP = np.fft.rfftfreq(1024, DT); S = semi.gain_sigma2(G_IQ); ANCHOR = 1.4
WINS = [(1.2, 1.6), (1.6, 2.0), (2.0, 2.4), (2.4, 2.8), (2.6, 3.0)]


def metrics(o, ideal):
    out = {}
    for a, b in WINS:
        w = semi.win(TW, a, b); O, I = o[w], ideal[w]; live = np.abs(I).sum(axis=0) > 0
        tc = [np.corrcoef(O[:, j], I[:, j])[0, 1] for j in np.where(live)[0]]
        amp = np.log(np.sqrt((O[:, live] ** 2).mean(axis=0)) / np.sqrt((I[:, live] ** 2).mean(axis=0)))
        out[f"{a:.1f}-{b:.1f}"] = dict(trace=float(np.median(tc)), pooled=float(np.corrcoef(O[:, live].ravel(), I[:, live].ravel())[0, 1]),
                                       amp_sd=float(np.std(amp)))
    return out


def run(did, version, cm, Q):
    arr, meta = load(did); D = arr["signal"] + (1.0 if version == "cal" else 4.0) * arr["noise"]; ideal = arr["ideal"]
    c = json.loads((ROOT / "results" / "pilot_semi" / f"{did}_{version}.json").read_text(encoding="utf-8"))["sr"]["c_sr"]
    x = lib.absorption(dict(cm, Q=Q * c), TW); dx = x - lib.absorption(dict(cm, Q=Q * c), [ANCHOR])[0]
    g = mwdip.Gather(D, TW, 4); inv = snr_rows(D, TW, F_OP); mod = g.model(x, dx)
    outs = dict(W0=uqw(D, cm, Q, c, 0.0, TW, ANCHOR, 0.0, inv_snr=inv), MW0_K1=g.apply(x, dx, 0.0, mod, F_OP, cm["fref"], S, 0),
                MW0_K9=g.apply(x, dx, 0.0, mod, F_OP, cm["fref"], S, 4))
    return dict(id=did, version=version, c_sr=c, scores={k: metrics(v, ideal) for k, v in outs.items()})


if __name__ == "__main__":
    k, n = int(sys.argv[2]), int(sys.argv[3]); OUT.mkdir(parents=True, exist_ok=True); cm = semi.curtain_cm(); Q = np.asarray(cm["Q"], float)
    for j, (d, v) in enumerate([(d, v) for v in ("cal", "nd") for d in DEV]):
        path = OUT / f"{d}_{v}.json"
        if j % n != k or path.exists():
            continue
        t0 = time.monotonic(); r = run(d, v, cm, Q); r["seconds"] = time.monotonic() - t0
        path.write_text(json.dumps(r, indent=1), encoding="utf-8"); print(d, v, f"{r['seconds']:.0f}s", flush=True)
