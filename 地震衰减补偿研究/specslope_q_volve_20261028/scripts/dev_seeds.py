"""Development (before the plan): SR on new reflectivity seeds 7001-7006 with the experiment-8 structure
(inline 10141), experiment-8 noise lines 10041/10061/10081 and the experiment-8 c values. Measures the bias
spread over reflectivity realizations and the coverage of the 90% interval. Writes analysis/dev_seeds.json."""
import sys, json, time
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import sq, gen  # noqa: E402
from sq import np  # noqa: E402

t0 = time.monotonic(); ctx = gen.prepare(10141, [10041, 10061, 10081]); cm = ctx["cm"]; Q = np.asarray(cm["Q"], float)
f = np.fft.rfftfreq(sq.NFFT, sq.DT); wins = sq.windows(); T = sq.mw_table(cm, Q, f, wins)
res = {}
for seed in range(7001, 7007):
    for c in [0.64, 0.8, 1.0, 1.58]:
        g = gen.signal(ctx, c, seed); x, lev = g["x"], g["level"]
        for i, N in enumerate(ctx["noise"]):
            e = sq.estimate(x + lev * N, cm, Q, table=T, wins=wins, nboot=50, seed=seed)
            res[f"s{seed}_c{c}_n{i}"] = dict(c_true=c, seed=seed, level=lev, **e)
        r = [res[f"s{seed}_c{c}_n{i}"] for i in range(3)]
        print(seed, c, "lev", round(lev, 2), "c_hat", [round(v["c"], 3) for v in r], "ci", [np.round(v["ci90"], 2).tolist() for v in r][:1],
              round(time.monotonic() - t0), "s", flush=True)
(sq.ROOT / "analysis" / "dev_seeds.json").write_text(json.dumps(res, indent=1), encoding="utf-8")
v = list(res.values()); lr = np.array([np.log(r["c"] / r["c_true"]) for r in v])
cov = np.mean([r["ci90"][0] <= r["c_true"] <= r["ci90"][1] for r in v])
print("log ratio mean", lr.mean().round(3), "sd", lr.std().round(3), "within 25%", np.mean(np.abs(lr) <= np.log(1.25)).round(3), "coverage", cov.round(3))
for c in [0.64, 0.8, 1.0, 1.58]:
    l = [np.log(r["c"] / c) for r in v if r["c_true"] == c]; print(c, np.round(np.mean(l), 3), np.round(np.std(l), 3))
