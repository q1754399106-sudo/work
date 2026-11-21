"""PLAN_ANGLE35_ZH.md section 4-5: reproduction check, E1-E5 / F1-F4 with the full mask (paired original) and the <= 35 deg mask
(variant A; variant B = cal_B, gate A2 aligned scenes), 68 % quantiles, c_hat accuracy, and the section-5 judgement."""
import json
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[1]; RES = ROOT / "results"; CS = (0.5, 0.7, 1.0, 1.4, 2.0); ALN = ("events_aligned", "thin_aligned")
G2 = ["M24", "M25", "M26", "M27"]; G3 = ["M28", "M29", "M30", "M31"]; FP = (60, 40)
load = lambda pat, ms: [r for m in ms for fp in FP for r in json.loads((RES / (pat % (m, fp))).read_text())]
N2, N3 = load("angle35_a2_%s_f%d.json", G2), load("angle35_a3_%s_f%d.json", G3)
O2, O3 = [r for r in load("a2_%s_f%d.json", G2) if r["seed"] < 6], load("a3_%s_f%d.json", G3)
out = {}

# reproduction: full-mask scores of the rerun vs stored rows
key = lambda r: (r["model"], r["fp"], r["scene"], r["c"], r["snr"], r["seed"], r["method"])
for name, new, old in (("a2", N2, O2), ("a3", N3, O3)):
    d = {key(r): r for r in new}; m = [(d[key(r)], r) for r in old if key(r) in d]
    out["repro_" + name] = dict(matched=len(m), stored=len(old), max_abs_ga=float(max(abs(a["ga"] - b["ga"]) for a, b in m)),
                                max_abs_wave=float(max(abs(a["wave"] - b["wave"]) for a, b in m)),
                                max_rel_chat=float(max(abs(a["c_hat"] / b["c_hat"] - 1) for a, b in m if b["c_hat"] == b["c_hat"])))
    out["repro_" + name]["ok"] = bool(m and len(m) == len(old) and out["repro_" + name]["max_abs_ga"] <= 1e-4 and out["repro_" + name]["max_rel_chat"] <= 1e-9)


def cells(R, method, q):
    d = {}
    for r in R:
        if r["method"] == method:
            d.setdefault((r["model"], r["fp"], r["scene"], r["c"], r["snr"]), []).append(r[q])
    return {k: float(np.mean(v)) for k, v in d.items()}


def gate2(R, sfx, cal="cal"):
    ga = {m: cells(R, m, "ga" + sfx) for m in ("oracle", "uncal", "tik_uncal", cal)}; wv = {m: cells(R, m, "wave" + sfx) for m in ("uncal", cal)}
    keys = [k for k in ga[cal] if k[2] in ALN]; o = {}
    e1 = []
    for g in sorted({(k[0], k[1], k[2], k[4]) for k in keys}):
        ks = [(g[0], g[1], g[2], c, g[3]) for c in CS]
        e1.append(np.mean([ga[cal][k] for k in ks]) / min(np.mean([ga["uncal"][k] for k in ks]), np.mean([ga["tik_uncal"][k] for k in ks])))
    o["E1_frac"] = float(np.mean(np.array(e1) <= 0.5)); o["E1_median_ratio"] = float(np.median(e1)); o["E1_pass"] = o["E1_frac"] >= 0.9
    o["E2"] = float(np.median([ga[cal][k] - ga["oracle"][k] for k in keys])); o["E2_pass"] = o["E2"] <= 0.10
    k1 = [k for k in keys if k[3] == 1.0]
    o["E3"] = float(np.median([ga[cal][k] - ga["uncal"][k] for k in k1])); o["E3_pass"] = o["E3"] <= 0.05
    o["E4"] = float(np.median([wv[cal][k] - wv["uncal"][k] for k in k1])); o["E4_pass"] = o["E4"] <= 0.02
    kn = [k for k in ga[cal] if k[2] == "events_nmo" and k[3] != 1.0]
    if kn:
        o["E5"] = float(np.mean([ga[cal][k] <= ga["uncal"][k] for k in kn])); o["E5_pass"] = o["E5"] >= 0.8
    o["by_c"] = {str(c): {m: float(np.median([ga[m][k] for k in keys if k[3] == c])) for m in ga} for c in CS}
    for c in CS:
        b = o["by_c"][str(c)]; b["reduction"] = 1 - b[cal] / min(b["uncal"], b["tik_uncal"])
    ncs = [c for c in CS if c != 1.0]; o["reduction_c_ne_1"] = [o["by_c"][str(c)]["reduction"] for c in ncs]
    o["q68"] = {m: float(np.percentile([ga[m][k] for k in keys], 68)) for m in ga}
    o["median"] = {m: float(np.median([ga[m][k] for k in keys])) for m in ga}
    return o


def gate3(R, sfx):
    S = cells(R, "sat", "ga" + sfx); o = {"n_cells": len(S)}; ops = ("sat", "direct", "split", "band3", "tik0"); C = {m: cells(R, m, "ga" + sfx) for m in ops}
    for nm, m in (("F1", "split"), ("F2", "band3")):
        d = np.array([S[k] - C[m][k] for k in S]); q25, q75 = np.percentile(d, [25, 75])
        o[nm] = dict(median_diff=float(np.median(d)), q25=float(q25), q75=float(q75)); o[nm + "_pass"] = bool(abs(np.median(d)) <= 0.02 and q25 >= -0.03 and q75 <= 0.03)
    r = np.array([C["tik0"][k] / S[k] for k in S]); o["F3"] = dict(median_ratio=float(np.median(r)), frac_ge2=float(np.mean(r >= 2)), min_ratio=float(r.min()))
    o["F3_pass"] = bool(np.median(r) >= 5 and np.mean(r >= 2) >= 0.9)
    r = np.array([C["direct"][k] / S[k] for k in S]); o["F4"] = dict(median_ratio=float(np.median(r))); o["F4_pass"] = bool(np.median(r) >= 1.2)
    o["median"] = {m: float(np.median(list(C[m].values()))) for m in ops}; o["q68"] = {m: float(np.percentile(list(C[m].values()), 68)) for m in ops}
    o["q68_ratio_to_sat"] = {m: o["q68"][m] / o["q68"]["sat"] for m in ops}
    o["median_cell_ratio_to_sat"] = {m: float(np.median([C[m][k] / S[k] for k in S])) for m in ops}
    o["by_c"] = {str(c): {m: float(np.median([v for k, v in C[m].items() if k[3] == c])) for m in ops} for c in CS}
    w = {m: cells(R, m, "wave" + sfx) for m in ops}; o["wave_median"] = {m: float(np.median(list(w[m].values()))) for m in ops}
    return o


def chat_acc(R, col):
    rr = [r for r in R if r["method"] == "cal" and r["scene"] == "events_aligned"]; o = {}
    for c in CS:
        v = [r[col] for r in rr if r["c"] == c]; ok = [x for x in v if x == x and x > 0]; e = np.log(np.array(ok) / c)
        o[str(c)] = dict(rms=float(np.sqrt(np.mean(e ** 2))), mean=float(np.mean(e)), fail=len(v) - len(ok))
    e = np.concatenate([np.log(np.array([r[col] for r in rr if r["c"] == c and r[col] == r[col] and r[col] > 0]) / c) for c in CS])
    o["all"] = dict(rms=float(np.sqrt(np.mean(e ** 2))), mean=float(np.mean(e))); return o


out["a2_full"] = gate2(N2, ""); out["a2_A35"] = gate2(N2, "35"); out["a2_B35"] = gate2(N2, "35", "cal_B"); out["a2_B_full"] = gate2(N2, "", "cal_B")
out["a3_full"] = gate3(N3, ""); out["a3_A35"] = gate3(N3, "35")
out["chat_A"] = chat_acc(N2, "c_hat"); out["chat_B"] = chat_acc(N2, "c_hat_B")
out["stored_a2_10seeds"] = {k: v for k, v in json.loads((RES / "a2_verdict.json").read_text()).items() if k != "descriptive"}


def judge(base, new, stat, passk, absolute=False):
    b, n = stat(base), stat(new); flip = bool(base[passk]) != bool(new[passk])
    ch = abs(n - b) if absolute else abs(n - b) / abs(b); lim = 0.01 if absolute else 0.20
    return dict(orig=b, new=n, orig_pass=bool(base[passk]), new_pass=bool(new[passk]), change=float(ch),
                verdict="结论变了" if flip else ("与角度无关" if ch < lim else "方向不变但幅度变"))


def judge_all(b2, n2, b3=None, n3=None):
    j = {"E1": judge(b2, n2, lambda o: o["E1_median_ratio"], "E1_pass"), "E2": judge(b2, n2, lambda o: o["E2"], "E2_pass"),
         "E3": judge(b2, n2, lambda o: o["E3"], "E3_pass"), "E4": judge(b2, n2, lambda o: o["E4"], "E4_pass")}
    if "E5" in n2:
        j["E5"] = judge(b2, n2, lambda o: o["E5"], "E5_pass")
    rb, rn = np.mean(b2["reduction_c_ne_1"]), np.mean(n2["reduction_c_ne_1"]); dirn = all(x > 0 for x in n2["reduction_c_ne_1"])
    j["reduction"] = dict(orig=float(rb), new=float(rn), change=float(abs(rn - rb) / rb), verdict=("结论变了" if not dirn else ("与角度无关" if abs(rn - rb) / rb < 0.2 else "方向不变但幅度变")))
    if b3:
        j["F1"] = judge(b3, n3, lambda o: o["F1"]["median_diff"], "F1_pass", True); j["F2"] = judge(b3, n3, lambda o: o["F2"]["median_diff"], "F2_pass", True)
        j["F3"] = judge(b3, n3, lambda o: o["F3"]["median_ratio"], "F3_pass"); j["F4"] = judge(b3, n3, lambda o: o["F4"]["median_ratio"], "F4_pass")
    return j


out["judge_A"] = judge_all(out["a2_full"], out["a2_A35"], out["a3_full"], out["a3_A35"])
out["judge_B"] = judge_all(out["a2_full"], out["a2_B35"])
(RES / "angle35_verdict.json").write_text(json.dumps(out, indent=1, ensure_ascii=False))
print(json.dumps({k: out[k] for k in ("repro_a2", "repro_a3", "judge_A", "judge_B")}, indent=1, ensure_ascii=False))
