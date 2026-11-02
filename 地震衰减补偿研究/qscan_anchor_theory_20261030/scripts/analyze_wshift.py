"""Analysis of PLAN_WSHIFT_ZH.md: kappa_hat and beta_SR per wavelet window, field (6 lines) vs semi (12 L). -> analysis/wshift.json"""
import json, statistics as st
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]; DV = 1 / 0.7 - 1 / 1.4
FIELD = [9985, 10065, 10141, 10209, 10289, 10369]; SEMI = [f"L{s}" for s in range(10001, 10013)]
bsr = json.loads((ROOT / "analysis" / "betasr.json").read_text())["rows"]


def kb(v, u_sr):
    if any(x["edge"] for x in v.values()):
        return None
    k = (v["0.7"]["v_sel"] - v["1.4"]["v_sel"]) / DV; return dict(kappa=k, beta_sr=v["1.0"]["v_sel"] - (1 - k) * u_sr - k)


out = {}
for w in ("w1.4", "w1.5", "w1.6", "w1.7"):
    F, S = {}, {}
    for ln in FIELD:
        r11 = {ca: json.loads((ROOT / "results11" / f"IL_{ln}__w1.4_a{ca}.json").read_text()) for ca in ("0.7", "1.0", "1.4")}
        u = 1 / r11["1.0"]["c_sr"]
        v = {ca: dict(v_sel=r11[ca]["v_sel"], edge=r11[ca]["edge"]) for ca in r11} if w == "w1.4" else \
            json.loads((ROOT / "results25" / f"F_{ln}__{w}.json").read_text())["scans"]
        F[ln] = kb(v, u)
    for j in SEMI:
        u = 1 / bsr[j]["c_sr"]
        v = json.loads((ROOT / "results16" / f"{j}.json").read_text())["scans"] if w == "w1.4" else \
            json.loads((ROOT / "results25" / f"S_{j}__{w}.json").read_text())["scans"]
        S[j] = kb(v, u)
    fb = [x["beta_sr"] for x in F.values() if x]; sb = [x["beta_sr"] for x in S.values() if x]
    out[w] = dict(field=F, semi=S, field_median=st.median(fb), field_n=len(fb), semi_median=st.median(sb), semi_max=max(sb),
                  semi_min=min(sb), semi_n=len(sb), gap=st.median(fb) - st.median(sb),
                  field_kappa=st.median(x["kappa"] for x in F.values() if x), semi_kappa=st.median(x["kappa"] for x in S.values() if x))
sh = [out[w] for w in ("w1.5", "w1.6", "w1.7")]
W1 = sum(x["field_median"] <= 0.20 for x in sh) >= 2
W2 = all(x["field_median"] >= 0.25 and x["field_median"] > x["semi_max"] for x in sh)
out["verdicts"] = dict(W1=W1, W2=W2, undecided=not (W1 or W2))
(ROOT / "analysis" / "wshift.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
for w in ("w1.4", "w1.5", "w1.6", "w1.7"):
    x = out[w]; print(w, "field beta_SR %.3f (n %d, kappa %.3f)  semi median %.3f [%.3f, %.3f] (n %d, kappa %.3f)  gap %.3f" % (
        x["field_median"], x["field_n"], x["field_kappa"], x["semi_median"], x["semi_min"], x["semi_max"], x["semi_n"], x["semi_kappa"], x["gap"]))
print(out["verdicts"])
