"""Analysis of PLAN_WLEN_ZH.md: realisation spread of beta against wavelet-window length at fixed centre.
-> analysis/wlen.json"""
import json, statistics as st
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]; DV = 1 / 0.7 - 1 / 1.4
WINS = ("W04_14", "W04_16", "W08_16", "W04_18", "W12_18")
OLD = list(range(10001, 10013)); NEW = list(range(10013, 10025)); SEEDS = OLD + NEW
FIELD = [9985, 10065, 10141, 10209, 10289, 10369]
KTH = {"W04_14": 0.675, "W04_16": 0.771, "W08_16": 0.771, "W04_18": 0.868, "W12_18": 0.868}
bsr = json.loads((ROOT / "analysis" / "betasr.json").read_text())["rows"]


def load(path):
    return json.loads(path.read_text(encoding="utf-8"))


def semi_scans(s, w):
    if s in OLD and w == "W04_14":
        return load(ROOT / "results16" / f"L{s}.json")["scans"]
    if s in OLD and w == "W04_16":
        return load(ROOT / "results25" / f"S_L{s}__w1.6.json")["scans"]
    return load(ROOT / "results27" / f"S_L{s}__{w}.json")["scans"]


def field_scans(ln, w):
    if w == "W04_14":
        r = {ca: load(ROOT / "results11" / f"IL_{ln}__w1.4_a{ca}.json") for ca in ("0.7", "1.0", "1.4")}
        return {ca: dict(v_sel=x["v_sel"], edge=x["edge"]) for ca, x in r.items()}, 1 / r["1.0"]["c_sr"]
    if w == "W04_16":
        sc = load(ROOT / "results25" / f"F_{ln}__w1.6.json")["scans"]
    else:
        sc = load(ROOT / "results27" / f"F_{ln}__{w}.json")["scans"]
    return sc, 1 / load(ROOT / "results11" / f"IL_{ln}__w1.4_a1.0.json")["c_sr"]


def kb(v, u_sr):
    if any(x["edge"] for x in v.values()):
        return None
    k = (v["0.7"]["v_sel"] - v["1.4"]["v_sel"]) / DV
    return dict(kappa=k, beta_true=v["1.0"]["v_sel"] - 1.0, beta_sr=v["1.0"]["v_sel"] - (1 - k) * u_sr - k)


c_sr = {s: bsr[f"L{s}"]["c_sr"] for s in OLD}
c_sr.update({s: load(ROOT / "results27" / f"S_L{s}__W04_14.json")["c_sr"] for s in NEW})
S = {w: {s: kb(semi_scans(s, w), 1 / c_sr[s]) for s in SEEDS} for w in WINS}
F = {}
for w in WINS:
    F[w] = {}
    for ln in FIELD:
        v, u = field_scans(ln, w); F[w][ln] = kb(v, u)

rng = np.random.default_rng(20261109); NB = 2000


def pair(wl, ws, key):
    ok = [s for s in SEEDS if S[wl][s] and S[ws][s]]
    a = np.array([S[wl][s][key] for s in ok]); b = np.array([S[ws][s][key] for s in ok])
    return ok, a, b


def ratio_boot(key):
    oA, aA, bA = pair("W08_16", "W04_16", key); oB, aB, bB = pair("W12_18", "W04_18", key)
    sd = lambda x: float(np.std(x, ddof=1))
    RA, RB = sd(aA) / sd(bA), sd(aB) / sd(bB); R = (RA * RB) ** 0.5; bs = []
    for _ in range(NB):
        iA = rng.integers(0, len(aA), len(aA)); iB = rng.integers(0, len(aB), len(aB))
        bs.append((sd(aA[iA]) / sd(bA[iA]) * sd(aB[iB]) / sd(bB[iB])) ** 0.5)
    lo, hi = np.quantile(bs, [0.05, 0.95])
    return dict(R_A=RA, R_B=RB, R=R, R_lo=float(lo), R_hi=float(hi), n_A=len(oA), n_B=len(oB),
                contains_0639=bool(lo <= 0.639 <= hi), contains_0875=bool(lo <= 0.875 <= hi))


out = dict(windows={}, field={})
for w in WINS:
    ok = [x for x in S[w].values() if x]; bt = np.array([x["beta_true"] for x in ok]); bsrv = np.array([x["beta_sr"] for x in ok])
    km = float(np.median([x["kappa"] for x in ok])); sdt = float(np.std(bt, ddof=1))
    fb = [x["beta_sr"] for x in F[w].values() if x]
    out["windows"][w] = dict(n=len(ok), kappa_median=km, kappa_theory=KTH[w], sd_beta_true=sdt, sd_beta_sr=float(np.std(bsrv, ddof=1)),
                             median_beta_true=float(np.median(bt)), median_beta_sr=float(np.median(bsrv)), floor=sdt / (1 - km),
                             field_n=len(fb), field_median_beta_sr=st.median(fb) if fb else None,
                             field_kappa_median=st.median(x["kappa"] for x in F[w].values() if x) if fb else None,
                             gap=(st.median(fb) - float(np.median(bsrv))) if fb else None)
f0 = out["windows"]["W04_14"]["floor"]
for w in WINS:
    out["windows"][w]["floor_ratio"] = out["windows"][w]["floor"] / f0
out["true"] = ratio_boot("beta_true"); out["sr"] = ratio_boot("beta_sr")
R = out["true"]
L1 = "supported" if R["R_hi"] < 1.0 else ("not supported" if R["R"] >= 0.95 else "undecided")
L2 = None if L1 != "supported" else ("mainly sampling" if R["R"] <= 0.70 else ("partly sampling" if R["R"] < 0.90 else "sampling minor"))
L3 = dict(W08_16_ge_0p9=out["windows"]["W08_16"]["floor_ratio"] >= 0.9, W12_18_gt_1=out["windows"]["W12_18"]["floor_ratio"] > 1)
out["verdicts"] = dict(L1=L1, L2=L2, L3_prediction_holds=L3)
out["per_seed"] = {w: {str(s): S[w][s] for s in SEEDS} for w in WINS}
out["per_line"] = {w: {str(ln): F[w][ln] for ln in FIELD} for w in WINS}
(ROOT / "analysis" / "wlen.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
for w in WINS:
    x = out["windows"][w]
    print(w, "n %d kappa %.3f (th %.3f) SD_true %.3f SD_sr %.3f med_true %.3f med_sr %.3f floor %.3f (x%.2f) | field med %s gap %s" % (
        x["n"], x["kappa_median"], x["kappa_theory"], x["sd_beta_true"], x["sd_beta_sr"], x["median_beta_true"], x["median_beta_sr"],
        x["floor"], x["floor_ratio"], None if x["field_median_beta_sr"] is None else round(x["field_median_beta_sr"], 3),
        None if x["gap"] is None else round(x["gap"], 3)))
print("beta_true ratios", {k: (round(v, 3) if isinstance(v, float) else v) for k, v in out["true"].items()})
print("beta_sr ratios", {k: (round(v, 3) if isinstance(v, float) else v) for k, v in out["sr"].items()})
print(out["verdicts"])
