"""Experiment 10 decisions E1-E4 (PLAN_ZH.md) and the descriptive groups. Writes analysis/summary.json."""
import sys, json
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import sq, run  # noqa: E402
from sq import np, bm  # noqa: E402
from scipy import stats  # noqa: E402

ROOT = sq.ROOT; AN = ROOT / "analysis"; L25 = np.log(1.25)


def ci(d):
    d = np.asarray(d, float); m = d.mean(); h = stats.t.ppf(0.95, len(d) - 1) * d.std(ddof=1) / np.sqrt(len(d))
    return [float(m), float(m - h), float(m + h)]


def load(group):
    return [json.loads(p.read_text(encoding="utf-8")) for p in sorted((ROOT / "results" / "synth").glob(f"{group}_*.json"))]


def accuracy(rows):
    est = [(r["c_true"], e) for r in rows for e in r["sr"].values()]
    lr = np.array([np.log(e["c"] / c) for c, e in est])
    cover = [e["ci90"][0] <= c <= e["ci90"][1] for c, e in est]
    per_c = {}
    for c in sorted({c for c, _ in est}):
        x = np.array([np.log(e["c"] / cc) for cc, e in est if cc == c])
        per_c[str(c)] = dict(n=len(x), mean_log_ratio=float(x.mean()), median_abs_log=float(np.median(np.abs(x))),
                             within25=float(np.mean(np.abs(x) <= L25)))
    return dict(n=len(lr), within25=float(np.mean(np.abs(lr) <= L25)), mean_log_ratio=float(lr.mean()),
                sd_log_ratio=float(lr.std(ddof=1)), coverage90=float(np.mean(cover)),
                median_ci_width_log=float(np.median([np.log(e["ci90"][1] / e["ci90"][0]) for _, e in est])), per_c=per_c)


def main():
    H, W = load("H"), load("W"); out = dict(n_H=len(H), n_W=len(W))
    out["H_accuracy"] = accuracy(H); out["W_accuracy"] = accuracy(W)
    heavy = [r for r in H if "scores" in r]
    d = [abs(np.log(r["sr"]["n0"]["c"] / r["c_true"])) - abs(np.log(r["qscan"]["factor"] / r["c_true"])) for r in heavy]
    out["E2"] = dict(n=len(d), abs_log_SR_minus_QScan=ci(d),
                     SR_abs_log_mean=float(np.mean([abs(np.log(r["sr"]["n0"]["c"] / r["c_true"])) for r in heavy])),
                     QScan_abs_log_mean=float(np.mean([abs(np.log(r["qscan"]["factor"] / r["c_true"])) for r in heavy])),
                     qscan_factors={r["id"]: r["qscan"]["factor"] for r in heavy})
    comp = {}
    for key in ("M1", "M1_lag"):
        ms = list(heavy[0]["scores"])
        comp[key] = dict(means={m: float(np.mean([r["scores"][m][key] for r in heavy])) for m in ms},
                         SR_minus_QScan=ci([r["scores"]["IQ30_SR"][key] - r["scores"]["IQ30_QScan"][key] for r in heavy]),
                         SR_minus_model=ci([r["scores"]["IQ30_SR"][key] - r["scores"]["IQ30_model"][key] for r in heavy]),
                         SR_minus_true=ci([r["scores"]["IQ30_SR"][key] - r["scores"]["IQ30_true"][key] for r in heavy]),
                         SCaSR_minus_IQ30SR=ci([r["scores"]["SCa_SR"][key] - r["scores"]["IQ30_SR"][key] for r in heavy]))
    out["compensation"] = comp
    a = out["H_accuracy"]; e2 = out["E2"]["abs_log_SR_minus_QScan"]; e4 = comp["M1"]["SR_minus_QScan"]
    out["verdict"] = dict(
        E1=bool(a["within25"] >= 0.8),
        E2=bool(e2[2] < 0),
        E3=bool(a["coverage90"] >= 0.75),
        E4="SR better" if e4[1] > 0 else "QScan better" if e4[2] < 0 else "no difference")
    bm.writejson(AN / "summary.json", out)
    print(json.dumps(out["verdict"]), flush=True)
    print("H", {k: (round(v, 3) if isinstance(v, float) else v) for k, v in a.items() if k != "per_c"})
    print("W", {k: (round(v, 3) if isinstance(v, float) else v) for k, v in out["W_accuracy"].items() if k != "per_c"})
    print("E2", np.round(e2, 3), "M1", {m: round(v, 3) for m, v in comp["M1"]["means"].items()})


if __name__ == "__main__":
    main()
