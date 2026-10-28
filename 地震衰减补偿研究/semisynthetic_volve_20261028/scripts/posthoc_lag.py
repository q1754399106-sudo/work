"""Post-hoc, descriptive (not in PLAN_ZH.md): M1 allowing one common time shift (+-20 ms, 1 ms) per output,
to separate the bulk dispersion delay (fref convention) from waveform improvement. Writes analysis/posthoc_lag.json."""
import sys, json
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import semi  # noqa: E402
from semi import np, TW, C_TRUE  # noqa: E402
import analyze  # noqa: E402

ROOT = semi.ROOT; w1 = semi.win(TW, *semi.WIN_LOG)


def shift(A, s):
    return np.stack([np.interp(TW + s, TW, a) for a in A.T], axis=1)


def lagged(out, ideal):
    best = (-np.inf, 0.0)
    for s in np.arange(-20, 21) / 1000:
        m = float(np.median(semi.colcorr(shift(out, s)[w1], ideal[w1])))
        best = max(best, (m, s * 1000))
    return best


res = {}
for c in C_TRUE:
    for i in range(10):
        did = f"P_c{c}_n{i}"; o = np.load(ROOT / "cache" / "out" / f"{did}.npz")
        res[did] = {m: lagged(o[m], o["ideal"]) for m in analyze.METHODS}
        print(did, {m: (round(v[0], 3), v[1]) for m, v in res[did].items()}, flush=True)
summ = {}
for c in C_TRUE:
    ids = [f"P_c{c}_n{i}" for i in range(10)]
    summ[f"c{c}"] = {m: dict(M1_lag=analyze.ci([res[d][m][0] for d in ids]),
                             shift_ms=float(np.median([res[d][m][1] for d in ids]))) for m in analyze.METHODS}
    summ[f"c{c}"]["QSCa_minus_Input"] = analyze.ci([res[d]["QSCa"][0] - res[d]["Input"][0] for d in ids])
    summ[f"c{c}"]["QSCa_minus_IQ30"] = analyze.ci([res[d]["QSCa"][0] - res[d]["IQ30"][0] for d in ids])
(ROOT / "analysis" / "posthoc_lag.json").write_text(json.dumps(dict(per_dataset=res, summary=summ), indent=1), encoding="utf-8")
for c, s in summ.items():
    print(c, {m: (round(v["M1_lag"][0], 3), v["shift_ms"]) for m, v in s.items() if m in analyze.METHODS},
          "QSCa-In", np.round(s["QSCa_minus_Input"], 3), "QSCa-IQ30", np.round(s["QSCa_minus_IQ30"], 3))
