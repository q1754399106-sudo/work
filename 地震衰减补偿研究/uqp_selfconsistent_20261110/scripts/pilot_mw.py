"""Exploratory pilot of UQ-MW (mw.py) on the used FD data H1-H6 (no pre-registered rules, NOTE_ZH.md).

At every Q multiplier m (assumed Q = m x true Q): K = 1, 5, 9 traces x {plug-in sig 0, prior sig 0.447}, signal
model fitted at the assumed Q. Same noise seeds as pilot_fd.py, so UQW / W0 / ORA of that run are the references
(ORA: single-trace ceiling with the true Q, the true noise power and the elastic target's signal power).

usage: pilot_mw.py run K N LEVEL
"""
import sys, json, time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import uqp, mw  # noqa: E402
from uqp import np, lib, semi, ROOT, DT, SIG0  # noqa: E402
import fdcheck as fd  # noqa: E402
from uqw import G_IQ  # noqa: E402
from stage0 import M_ALL  # noqa: E402

OUT = ROOT / "results" / "pilot_mw"; F_OP = np.fft.rfftfreq(1024, DT); S = semi.gain_sigma2(G_IQ)
LEVELS = {"14dB": (14.0, 20261113000), "30dB": (30.0, 20261121000)}; MW = (0, 2, 4)        # K = 1, 5, 9


def run(name, level, r):
    D0, T, cm, tw, fp = fd.load_fd(name); db, base = LEVELS[level]; seed = base + 1000 * (fd.MODELS.index(name) + 1) + r
    D = D0 + fd.noise(D0, db, seed); x1 = lib.absorption(cm, tw); f, cen, P, C = mw.lag_spectra(D, tw); rec = []
    for m in M_ALL:
        x = x1 / m; ls0, q, N = mw.signal_model(f, cen, P, C, tw, x, x); d = dict(m=m, q_med=float(np.median(q[:, (f >= 10) & (f <= 60)])))
        for k in MW:
            for tag, sig in (("MW0", 0.0), ("UQMW", SIG0)):
                X = mw.uqmw(D, tw, x, x, sig, ls0, q, N, f, cen, F_OP, cm["fref"], S, k)
                d[f"{tag}_K{2 * k + 1}"] = float(fd.score_fd(X, T, tw)[0])
        rec.append(d)
    return dict(model=name, level=level, r=r, seed=seed, records=rec)


if __name__ == "__main__":
    k, n, level = int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]; OUT.mkdir(parents=True, exist_ok=True)
    jobs = [(h, r) for h in fd.MODELS for r in fd.R]
    for j, (h, r) in enumerate(jobs):
        path = OUT / f"{h}_{level}_r{r}.json"
        if j % n != k or path.exists():
            continue
        t0 = time.monotonic(); out = run(h, level, r); out["seconds"] = time.monotonic() - t0
        path.write_text(json.dumps(out, indent=1), encoding="utf-8"); q = out["records"][4]
        print(h, level, r, " ".join("%s %.3f" % (key, v) for key, v in q.items() if key != "m"), f"{out['seconds']:.0f}s", flush=True)
