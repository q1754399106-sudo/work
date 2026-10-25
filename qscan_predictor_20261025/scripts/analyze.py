"""Development scoring of PLAN_ZH.md: reproduction check, leave-one-model-out cross-validation of
the shift rules C, L and T on features Z1-Z6, judgment P1 and descriptive statistics."""
from pathlib import Path
import itertools, json
import numpy as np, pandas as pd

ROOT = Path(__file__).resolve().parents[1]
IN, A = ROOT / "cache" / "inputs", ROOT / "analysis"
FEATS = ["Z1", "Z2", "Z3", "Z4", "Z5", "Z6"]
PRIMARY = ["Q100_14dB", "Q070_14dB"]
KS = range(-4, 5)
FSCENE = {("F1", "Q100_14dB"): "S01_F1_Q100_14dB", ("F1", "Q100_08dB"): "S03_F1_Q100_08dB", ("F1", "Q070_14dB"): "S05_F1_Q070_14dB",
          ("F2", "Q100_14dB"): "S02_F2_Q100_14dB", ("F2", "Q100_08dB"): "S04_F2_Q100_08dB", ("F2", "Q070_14dB"): "S06_F2_Q070_14dB"}


def curves():
    g = pd.read_csv(IN / "qscan_offset_20261025/analysis/offset_curves.csv")
    f = pd.read_csv(IN / "combo_bias_variance_20261024/analysis/sweep_curves.csv")
    f = f[(f.method == "Standard") & f.scene.isin(FSCENE.values())].copy()
    inv = {v: k for k, v in FSCENE.items()}
    f["model"] = [inv[s][0] for s in f.scene]
    f["scene"] = [inv[s][1] for s in f.scene]
    c = pd.concat([g[["model", "scene", "index", "mse", "times_selected"]], f[["model", "scene", "index", "mse", "times_selected"]]])
    out = {}
    for (m, s), d in c.groupby(["model", "scene"]):
        d = d.sort_values("index")
        grid = d[d["index"] >= 2]
        mse = grid.mse.to_numpy()                       # grid index 1..13 -> mse[g-1]
        sel = grid.times_selected.to_numpy()
        out[(m, s)] = dict(mse=mse, mode=int(np.argmax(sel)) + 1, best=int(np.argmin(mse)) + 1, counts=sel)
    return out


def rel(cv, key, k):
    c = cv[key]
    g = int(np.clip(c["mode"] - k, 1, 13))
    return c["mse"][g - 1] / c["mse"][c["mode"] - 1] - 1


def fit(rule, train, cv):
    """train: DataFrame rows (model, scene, Z, k) of training primary scenes."""
    keys = list(zip(train.model, train.scene))
    if rule == "C":
        best = min(KS, key=lambda k0: (np.mean([rel(cv, q, k0) for q in keys]), abs(k0)))
        return lambda z: best, dict(k0=best)
    z, k = train.Z.to_numpy(), train.k.to_numpy()
    if rule == "L":
        b, a = np.polyfit(z, k, 1)
        return (lambda zz: int(np.clip(np.round(a + b * zz), -4, 4))), dict(a=a, b=b)
    u = np.unique(z)
    ts = (u[:-1] + u[1:]) / 2
    cand = []
    for t, k1, k2 in itertools.product(ts, KS, KS):
        r = np.mean([rel(cv, q, k1 if zz <= t else k2) for q, zz in zip(keys, z)])
        cand.append((r, abs(k1) + abs(k2), t, k1, k2))
    r, _, t, k1, k2 = min(cand)
    return (lambda zz: k1 if zz <= t else k2), dict(t=t, k1=k1, k2=k2)


def main():
    cv = curves()
    F = pd.concat([pd.read_csv(f) for f in sorted(A.glob("features_*.csv"))])
    # reproduction: Octave selections versus archived counts (G) or per-replicate factors (F)
    rep = pd.read_csv(IN / "combo_bias_variance_20261024/analysis/replicates.csv").drop_duplicates(["scene", "replicate"])
    repro = []
    for (m, s), d in F.groupby(["model", "scene"]):
        if m.startswith("G"):
            oct_counts = np.bincount(d.best, minlength=14)[1:]
            ok = bool((oct_counts == cv[(m, s)]["counts"]).all())
        else:
            arch = rep[rep.scene == FSCENE[(m, s)]].set_index("replicate").qscan_factor
            ok = bool(all(abs(arch[r] - f) < 1e-12 for r, f in zip(d.replicate, d.qscan_factor)))
        repro.append(dict(model=m, scene=s, identical=ok, octave_mode=int(d.best.mode().iloc[0]), archived_mode=cv[(m, s)]["mode"]))
    repro = pd.DataFrame(repro)
    S = F.groupby(["model", "scene"])[FEATS].mean().reset_index()
    S["k"] = [cv[(m, s)]["mode"] - cv[(m, s)]["best"] for m, s in zip(S.model, S.scene)]
    S.to_csv(A / "scene_features.csv", index=False)
    models = sorted(S.model.unique())
    results, rows = {}, []
    for feat, rule in [("-", "C")] + [(f, r) for f in FEATS for r in ["L", "T"]]:
        name = "C" if rule == "C" else f"{rule}:{feat}"
        for held in models:
            tr = S[(S.model != held) & S.scene.isin(PRIMARY)]
            tr = tr.assign(Z=tr[feat] if feat != "-" else 0.0)
            predict, params = fit(rule, tr, cv)
            for _, row in S[S.model == held].iterrows():
                kh = predict(row[feat] if feat != "-" else 0.0)
                rows.append(dict(candidate=name, held_out=held, scene=row.scene, k_true=row.k, k_hat=kh,
                                 rel=rel(cv, (held, row.scene), kh), primary=row.scene in PRIMARY))
    R = pd.DataFrame(rows)
    R.to_csv(A / "cv_results.csv", index=False)
    summ = []
    for name, d in R.groupby("candidate", sort=False):
        p = d[d.primary]
        summ.append(dict(candidate=name, mean_rel=p.rel.mean(), worst_rel=p.rel.max(), n_improved=int((p.rel < 0).sum()),
                         n=len(p), mean_rel_08dB=d[~d.primary].rel.mean(), worst_rel_08dB=d[~d.primary].rel.max(),
                         P1=bool(p.rel.mean() <= -0.10 and p.rel.max() <= 0.10 and (p.rel < 0).sum() >= 14)))
    summ = pd.DataFrame(summ)
    summ.to_csv(A / "cv_summary.csv", index=False)
    passing = summ[summ.P1].copy()
    selected, final = None, None
    if len(passing):
        passing["isT"] = passing.candidate.str.startswith("T")
        selected = passing.sort_values(["mean_rel", "isT", "worst_rel"], ascending=[True, False, True]).iloc[0].candidate
        rule, feat = selected.split(":")
        tr = S[S.scene.isin(PRIMARY)].assign(Z=S[feat])
        _, final = fit(rule, tr, cv)
        final = dict(rule=rule, feature=feat, **{k: float(v) for k, v in final.items()})
    prim = S[S.scene.isin(PRIMARY)]
    spear = {f: float(prim[f].rank().corr(prim.k.rank())) for f in FEATS}
    a, b = S[S.scene == "Q100_14dB"].set_index("model"), S[S.scene == "Q070_14dB"].set_index("model")
    inv = {f: float(np.median(np.abs(np.log(b[f] / a[f])) if f in ["Z1", "Z2", "Z3", "Z4"] else np.abs(b[f] - a[f]))) for f in FEATS}
    out = dict(reproduction_all_identical=bool(repro.identical.all()), reproduction=repro.to_dict("records"),
               P1_passed=bool(len(passing)), selected=selected, final_rule=final, summary=summ.to_dict("records"),
               spearman_primary=spear, invariance_median=inv)
    (A / "decisions.json").write_text(json.dumps(out, indent=2, ensure_ascii=False, default=float), encoding="utf-8")
    print(summ.to_string(index=False))
    print(json.dumps({k: out[k] for k in ["reproduction_all_identical", "P1_passed", "selected", "final_rule", "spearman_primary", "invariance_median"]},
                     indent=2, default=float))


if __name__ == "__main__":
    main()
