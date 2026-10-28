"""Development (seen data): SR on the experiment-8 build (c = 0.64, 0.8, 1.0, 1.58; noises R0-R9 and G0-G4 at
the calibrated levels, and noise-free) and on the real inline 10141, primary multi-window estimator and the
two-window variant. Writes analysis/dev.json."""
import sys, json
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import sq  # noqa: E402
from sq import np  # noqa: E402

B = np.load(json.loads((sq.ROOT / "cache" / "SNAPSHOT.json").read_text(encoding="utf-8"))["local"]["build8.npz"]["path"])
lv = json.loads((sq.IN / "semisynthetic_volve_20261028" / "analysis" / "build.json").read_text(encoding="utf-8"))["level"]
cm = sq.semi.curtain_cm(); Q = np.asarray(cm["Q"], float)
f = np.fft.rfftfreq(sq.NFFT, sq.DT); wins = sq.windows(); T = sq.mw_table(cm, Q, f, wins); T2 = sq.decay_table(cm, Q, f)
res = {}
for c in sq.semi.C_TRUE:
    for nk in ["clean"] + [f"R{i}" for i in range(10)] + [f"G{i}" for i in range(5)]:
        X = B[f"sig_c{c}"] + (0 if nk == "clean" else lv[f"c{c}"] * B[f"noise_{nk}"])
        e = sq.estimate(X, cm, Q, table=T, wins=wins, nboot=100 if nk in ("R0", "R1") else 0)
        e2 = sq.estimate2(X, cm, Q, table=T2); e2.pop("_curve")
        res[f"c{c}_{nk}"] = dict(c_true=c, mw=e, two=e2)
    r = [res[f"c{c}_R{i}"]["mw"]["c"] for i in range(10)]; g = [res[f"c{c}_G{i}"]["mw"]["c"] for i in range(5)]
    print(c, "clean", round(res[f"c{c}_clean"]["mw"]["c"], 3), "R", np.round([np.median(r), min(r), max(r)], 3),
          "G", np.round(np.median(g), 3), "two-window R median", round(float(np.median([res[f"c{c}_R{i}"]["two"]["c"] for i in range(10)])), 3),
          "se_jack", round(res[f"c{c}_R0"]["mw"]["se_jack_log"], 3), "se_boot", round(res[f"c{c}_R0"]["mw"]["se_boot_log"], 3),
          "ci", np.round(res[f"c{c}_R0"]["mw"]["ci90"], 3), "fhi", res[f"c{c}_R0"]["mw"]["f_hi"], flush=True)
e = sq.estimate(B["D_real"], cm, Q, table=T, wins=wins, nboot=200); res["real10141"] = e
print("real", {k: v for k, v in e.items()}, flush=True)
(sq.ROOT / "analysis" / "dev.json").write_text(json.dumps(res, indent=1), encoding="utf-8")
