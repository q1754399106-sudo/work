"""Analysis of PLAN_PHI2_ZH.md: per-sample curvature C_t of the split objective, kappa_curv = a0 sum(C/a) / sum(C).

Fit: 5 grid points centred on the total-J grid minimum (sorted by v), least-squares quadratic, C_t = 2 x leading
coefficient; fewer than two evaluated points on either side -> unreliable. Block version: 10 equal-length blocks.
Comparators: U1320A two-anchor kappa (analysis/u1320a.json); Volve S/D (results7) and D2/M (results10) median
two-anchor kappa over the 12 datasets; base: two-anchor kappa of this run (c_a 0.7 vs 1.4). -> analysis/phi2.json
"""
import json, glob, statistics as st
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]; R = ROOT / "results19"; DV = 1 / 0.7 - 1 / 1.4


def kcurv(r):
    v = np.asarray(r["v"]); J = np.asarray(r["J"]); jm = int(np.argmin(J))
    if jm < 2 or jm > len(v) - 3:
        return None
    Jt = np.asarray(r["Jres"]) + np.asarray(r["Jpen"]); sl = slice(jm - 2, jm + 3)
    C = 2 * np.polyfit(v[sl], Jt[sl], 2)[0]; a = np.asarray(r["a"]); t = np.asarray(r["t"]); a0 = r["a0"]; m = a > 0
    k = float(a0 * np.sum(C[m] / a[m]) / np.sum(C[m]))
    edges = np.linspace(t[0], t[-1] + 1e-9, 11); b = np.clip(np.searchsorted(edges, t, side="right") - 1, 0, 9)
    Cb = np.array([C[b == i].sum() for i in range(10)]); ab = np.array([a[b == i].mean() for i in range(10)])
    kb = float(a0 * np.sum(Cb / ab) / np.sum(Cb))
    return dict(kappa_curv=k, kappa_block=kb, share=(Cb / Cb.sum()).round(4).tolist(),
                block_start=edges[:-1].round(3).tolist(), negative_share=float(-C[C < 0].sum() / np.abs(C).sum()),
                total_curv=float(C.sum()), conv_nbhd=bool(all(r["converged"][jm - 1:jm + 2])))


recs = {Path(p).stem: json.loads(Path(p).read_text()) for p in sorted(glob.glob(str(R / "*.json")))}
out = {"U": {}, "V": {}, "repro": {}}

# stage 1
u18 = json.loads((ROOT / "analysis" / "u1320a.json").read_text())
for line in ("3677-32", "4508-32"):
    for w in ("w2.24", "w3.00"):
        row = {"kappa_hat": u18[f"{line}|{w}"]["kappa"]}
        for va in (0.5, 1.0, 2.0):
            r = recs.get(f"U_{line}__{w}_va{va}")
            if r is None:
                continue
            row[str(va)] = kcurv(r)
            old = json.loads((ROOT / "results18" / f"{line}__{w}_va{va}.json").read_text())["v_sel"]
            out["repro"][f"U_{line}__{w}_va{va}"] = round(r["v_sel"] - old, 4)
        out["U"][f"{line}|{w}"] = row

# stage 2
def two_anchor(dirn, prof):
    ks = []
    for f in sorted(glob.glob(str(ROOT / dirn / f"{prof}_*__a0.7.json"))):
        a = json.loads(Path(f).read_text()); b = json.loads(Path(f.replace("a0.7", "a1.4")).read_text())
        ks.append((a["v_sel"] - b["v_sel"]) / DV)
    return ks


comp = {"S": two_anchor("results7", "S"), "D": two_anchor("results7", "D"), "D2": two_anchor("results10", "D2"),
        "M": two_anchor("results10", "M")}
kb = []
for k, r in recs.items():
    if k.startswith("V_B_") and k.endswith("__a0.7"):
        r2 = recs.get(k.replace("a0.7", "a1.4"))
        if r2:
            kb.append((r["v_sel"] - r2["v_sel"]) / DV)
comp["B"] = kb
for prof, src in (("B", None), ("S", "results7"), ("D", "results7"), ("D2", "results10"), ("M", "results10")):
    per = {}
    for k, r in recs.items():
        if k.startswith(f"V_{prof}_") and k.endswith("__a1.0"):
            per[r["id"]] = kcurv(r)
            if src:
                old = json.loads((ROOT / src / f"{prof}_{r['id']}__a1.0.json").read_text())["v_sel"]
                out["repro"][k] = round(r["v_sel"] - old, 4)
    other = {}
    if prof == "B":
        for ca in ("0.7", "1.4"):
            other[ca] = [kcurv(r)["kappa_curv"] for k, r in recs.items()
                         if k.startswith("V_B_") and k.endswith(f"__a{ca}") and kcurv(r)]
    ok = [x["kappa_curv"] for x in per.values() if x]
    out["V"][prof] = dict(n=len(per), n_reliable=len(ok), kappa_curv_median=(st.median(ok) if ok else None),
                          kappa_block_median=(st.median([x["kappa_block"] for x in per.values() if x]) if ok else None),
                          kappa_curv=sorted(round(x, 3) for x in ok), two_anchor_median=(st.median(comp[prof]) if comp[prof] else None),
                          n_two_anchor=len(comp[prof]), other_anchor_medians={c: st.median(v) for c, v in other.items() if v},
                          share_median=(np.median([x["share"] for x in per.values() if x], axis=0).round(3).tolist() if ok else None),
                          per_dataset=per)

# verdicts
U = out["U"]
def kc(key):
    x = U.get(key, {}).get("1.0"); return x["kappa_curv"] if x else None

p1 = [(kc(k) is not None, kc(k) is not None and abs(kc(k) - U[k]["kappa_hat"]) <= 0.10) for k in U]
P1 = None if len(p1) < 4 or not all(a for a, _ in p1) else all(b for _, b in p1)
a, b = kc("4508-32|w3.00"), kc("3677-32|w3.00"); P2 = None if a is None or b is None else a <= b - 0.10
p3 = []
for prof, row in out["V"].items():
    if row["n_reliable"] < max(1, row["n"] / 2) or row["two_anchor_median"] is None:
        p3.append(None)
    else:
        p3.append(abs(row["kappa_curv_median"] - row["two_anchor_median"]) <= 0.10)
P3 = None if (len(p3) < 5 or any(x is None for x in p3)) else all(p3)
out["verdicts"] = dict(P1=P1, P2=P2, P3=P3)
(ROOT / "analysis" / "phi2.json").write_text(json.dumps(out, indent=1), encoding="utf-8")

for k, row in U.items():
    print(k, "hat %.3f" % row["kappa_hat"], " ".join(f"va{va}: {row[va]['kappa_curv']:.3f}/{row[va]['kappa_block']:.3f}"
                                                  for va in ("0.5", "1.0", "2.0") if row.get(va)))
    if row.get("1.0"):
        print("   share", row["1.0"]["share"], "neg %.3f" % row["1.0"]["negative_share"])
for prof, row in out["V"].items():
    print(prof, {k: row[k] for k in ("n", "n_reliable", "kappa_curv_median", "kappa_block_median", "two_anchor_median",
                                    "other_anchor_medians")})
    print("   share", row["share_median"])
print("repro max |dv|", max((abs(x) for x in out["repro"].values()), default=None))
print(out["verdicts"])
