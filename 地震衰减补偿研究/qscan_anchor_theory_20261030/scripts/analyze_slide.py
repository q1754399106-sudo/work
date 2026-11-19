"""Score the sliding-wavelet-window field test (PLAN_SLIDE_FIELD_ZH.md) from results28/ -> analysis/slide_field.json.

usage: analyze_slide.py [results_dir]    (default results28; a different directory is only for testing the scoring code)
Pure numpy/json; does not need the field data. Theory values KTH are frozen in the plan and re-derived here from the
curtain model (field_vsp_volve_f11_20261026/inputs/posthoc_curtain_meta.json) as a consistency check.
"""
import sys, json, glob
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]; NB = 2000; SEED = 20261120
T0S = (1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6); CAS = (0.7, 1.0, 1.4); DV = 1 / 0.7 - 1 / 1.4
KTH = {1.4: 0.675, 1.6: 0.771, 1.8: 0.868, 2.0: 0.964, 2.2: 1.028, 2.4: 1.084, 2.6: 1.140}   # frozen (plan, table 1)
T_CROSS_TH = 2.1


def kappa_theory():
    m = json.loads((ROOT.parent / "field_vsp_volve_f11_20261026" / "inputs" / "posthoc_curtain_meta.json").read_text(encoding="utf-8"))
    q, v, h = (np.asarray(m[k], float) for k in ("Q_VSP", "vp", "h"))
    tw = np.round(np.arange(300, 751) * 0.004, 6); b = np.r_[0, np.cumsum(2 * h / v[:3]), np.inf]

    def a(t):
        t = np.atleast_1d(t); return np.sum(np.maximum(0, np.minimum(t[:, None], b[1:]) - b[:-1]) / q, axis=1)
    A = a(tw)
    return {t0: float(a(t0)[0] * A.sum() / (A ** 2).sum()) for t0 in T0S}


def boot_ci(x, stat, rng):
    x = np.asarray(x, float); b = [stat(x[rng.integers(0, len(x), len(x))]) for _ in range(NB)]
    return [float(np.percentile(b, 5)), float(np.percentile(b, 95))]


def slope(x, y):
    return float(np.polyfit(x, y, 1)[0])


def crossing(t, k):
    """First t0 where the median curve reaches 1 (linear interpolation between neighbouring windows); None if never."""
    for i in range(1, len(t)):
        if k[i - 1] < 1 <= k[i]:
            return float(t[i - 1] + (1 - k[i - 1]) * (t[i] - t[i - 1]) / (k[i] - k[i - 1]))
    return None


def spearman(x, y):
    rx = np.argsort(np.argsort(x)); ry = np.argsort(np.argsort(y)); return float(np.corrcoef(rx, ry)[0, 1])


def reproduction(R):
    """S0: new scans against the same geometry scanned before (results11 t0 1.4, results25 t0 1.6, results9 t0 2.6)."""
    pairs = []
    for (ln, t0, ca), r in R.items():
        old = None
        if t0 == 1.4:
            p = ROOT / "results11" / f"IL_{ln}__w1.4_a{ca}.json"
            old = json.loads(p.read_text(encoding="utf-8"))["v_sel"] if p.exists() else None
        elif t0 == 1.6:
            p = ROOT / "results25" / f"F_{ln}__w1.6.json"
            old = json.loads(p.read_text(encoding="utf-8"))["scans"][str(ca)]["v_sel"] if p.exists() else None
        elif t0 == 2.6:
            p = ROOT / "results9" / f"IL_{ln}__w2.6_a{ca}.json"
            old = json.loads(p.read_text(encoding="utf-8"))["v_sel"] if p.exists() else None
        if old is not None:
            pairs.append(dict(inline=ln, t0=t0, c_anchor=ca, new=r["v_sel"], old=old, diff=r["v_sel"] - old))
    d = np.abs([p["diff"] for p in pairs]) if pairs else np.array([np.inf])
    return dict(n_pairs=len(pairs), share_within_0002=float(np.mean(d <= 0.002)), max_abs_diff=float(d.max()), pairs=pairs)


if __name__ == "__main__":
    rdir = ROOT / (sys.argv[1] if len(sys.argv) > 1 else "results28"); rng = np.random.default_rng(SEED)
    kth_model = kappa_theory()
    assert all(abs(kth_model[t] - KTH[t]) < 5e-4 for t in T0S), kth_model
    R = {}
    for f in glob.glob(str(rdir / "*.json")):
        r = json.loads(Path(f).read_text(encoding="utf-8")); R[(r["inline"], round(r["t0"], 1), r["c_anchor"])] = r
    lines = sorted({k[0] for k in R}); csr = {k[0]: r["c_sr"] for k, r in R.items()}
    out = dict(results_dir=rdir.name, n_scans=len(R), kappa_theory=KTH, windows={}, lines={}, verdicts={})

    # per line and window
    K = {}       # (ln, t0) -> kappa_hat for reliable line-windows
    for ln in lines:
        out["lines"][str(ln)] = {}
        for t0 in T0S:
            s = {ca: R.get((ln, t0, ca)) for ca in CAS}
            if any(x is None for x in s.values()):
                out["lines"][str(ln)][str(t0)] = dict(status="missing"); continue
            ok = all(not x["edge"] and x["neighbourhood_converged"] for x in s.values())
            v = {ca: s[ca]["v_sel"] for ca in CAS}; k = (v[0.7] - v[1.4]) / DV
            rec = dict(status="ok" if ok else "unreliable", kappa=k, v=v, monotone=bool(v[0.7] > v[1.0] > v[1.4]),
                       beta_sr=v[1.0] - (1 - k) / csr[ln] - k, c1=1 / v[1.0])
            out["lines"][str(ln)][str(t0)] = rec
            if ok:
                K[(ln, t0)] = k

    # per window
    med = {}
    for t0 in T0S:
        recs = [out["lines"][str(ln)].get(str(t0), {}) for ln in lines]
        n_all = sum(r.get("status") in ("ok", "unreliable") for r in recs); n_bad = sum(r.get("status") == "unreliable" for r in recs)
        ks = [K[(ln, t0)] for ln in lines if (ln, t0) in K]
        w = dict(n_lines=n_all, n_unreliable=n_bad, determinable=bool(n_all > 0 and n_bad <= n_all / 4 and len(ks) > 0))
        if ks:
            ok_l = [ln for ln in lines if (ln, t0) in K]
            u = np.array([1 / csr[ln] for ln in ok_l]); v1 = np.array([out["lines"][str(ln)][str(t0)]["v"][1.0] for ln in ok_l])
            xy = np.c_[u, v1]
            w.update(kappa_median=float(np.median(ks)), kappa_ci=boot_ci(ks, np.median, rng), kappa_sd=float(np.std(ks, ddof=1)),
                     diff_theory=float(np.median(ks) - KTH[t0]), n_below1=int(sum(k < 1 for k in ks)), n_above1=int(sum(k > 1 for k in ks)),
                     n_ok=len(ks), monotone_anchor=int(sum(out["lines"][str(ln)][str(t0)]["monotone"] for ln in ok_l)),
                     beta_sr_median=float(np.median([out["lines"][str(ln)][str(t0)]["beta_sr"] for ln in ok_l])),
                     c1_median=float(np.median(1 / v1)),
                     slope_v1_usr=slope(u, v1),
                     slope_v1_usr_ci=boot_ci(np.arange(len(u)), lambda i: slope(xy[i.astype(int), 0], xy[i.astype(int), 1]), rng),
                     slope_theory=1 - KTH[t0])
            med[t0] = w["kappa_median"]
        out["windows"][str(t0)] = w

    W = out["windows"]; det = [t0 for t0 in T0S if W[str(t0)]["determinable"]]
    out["n_windows_determinable"] = len(det)
    V = out["verdicts"]

    # S0 reproduction
    rep = reproduction(R); out["reproduction"] = rep
    V["S0"] = bool(rep["n_pairs"] > 0 and rep["share_within_0002"] >= 0.95 and rep["max_abs_diff"] <= 0.01)

    if len(det) < 5:
        for s in ("S1", "S2", "S3", "S4"):
            V[s] = "undeterminable"
    else:
        # S1 tracking: every determinable window within +-0.10; slope of medians on theory in [0.85, 1.20]
        d = [W[str(t)]["diff_theory"] for t in det]; b = slope([KTH[t] for t in det], [med[t] for t in det])
        out["S1_slope_medians_on_theory"] = b; out["S1_max_abs_diff"] = float(np.max(np.abs(d)))
        V["S1"] = bool(np.max(np.abs(d)) <= 0.10 and 0.85 <= b <= 1.20)
        # S2 crossing: median curve crosses 1 in [1.9, 2.3] s; >= 90% of reliable lines below 1 at t0 <= 1.8, above 1 at t0 >= 2.4
        tc = crossing(det, [med[t] for t in det]); out["S2_crossing_t0"] = tc
        sh = [t for t in det if t <= 1.8]; dp = [t for t in det if t >= 2.4]
        frac_sh = [W[str(t)]["n_below1"] / W[str(t)]["n_ok"] for t in sh]; frac_dp = [W[str(t)]["n_above1"] / W[str(t)]["n_ok"] for t in dp]
        out["S2_share_below1_shallow"] = frac_sh; out["S2_share_above1_deep"] = frac_dp
        V["S2"] = bool(tc is not None and 1.9 <= tc <= 2.3 and sh and dp and min(frac_sh) >= 0.9 and min(frac_dp) >= 0.9)
        # S3 per-line monotonicity in t0: Spearman >= 0.9 on lines with >= 6 reliable windows; >= 90% of such lines
        rho = {}
        for ln in lines:
            ts = [t for t in det if (ln, t) in K]
            if len(ts) >= 6:
                rho[str(ln)] = spearman(ts, [K[(ln, t)] for t in ts])
        out["S3_spearman"] = rho; n3 = len(rho); g3 = sum(r >= 0.9 for r in rho.values())
        V["S3"] = bool(n3 >= 20 and g3 >= 0.9 * n3)
        # S4 (secondary): cross-line slope of v_sel(c_a = 1) on u_SR falls with t0 and changes sign
        bs = [W[str(t)]["slope_v1_usr"] for t in det]; out["S4_spearman_slope_t0"] = spearman(det, bs)
        V["S4"] = bool(out["S4_spearman_slope_t0"] <= -0.8 and det[0] <= 1.6 and det[-1] >= 2.4 and bs[0] > 0 and bs[-1] < 0)

    out["notes"] = "S0 invalidates S1-S4 when false; S4 is secondary (SR used as proxy for the true u)."
    (ROOT / "analysis" / ("slide_field.json" if rdir.name == "results28" else f"slide_field_{rdir.name}.json")).write_text(
        json.dumps(out, indent=1, default=float), encoding="utf-8")
    print("scans", len(R), "lines", len(lines), "windows determinable", len(det))
    print("t0   th     median  [90%]           sd     n_ok bad below1 above1 slope_v1~uSR (th)    beta_sr  c1")
    for t0 in T0S:
        w = W[str(t0)]
        if "kappa_median" in w:
            print(f"{t0:.1f}  {KTH[t0]:.3f}  {w['kappa_median']:.3f}  [{w['kappa_ci'][0]:.3f}, {w['kappa_ci'][1]:.3f}]  "
                  f"{w['kappa_sd']:.3f}  {w['n_ok']:4d} {w['n_unreliable']:3d} {w['n_below1']:6d} {w['n_above1']:6d}  "
                  f"{w['slope_v1_usr']:+.3f} ({w['slope_theory']:+.3f})  {w['beta_sr_median']:+.3f}  {w['c1_median']:.3f}")
    print("reproduction", {k: v for k, v in rep.items() if k != "pairs"})
    print("verdicts", V, "crossing", out.get("S2_crossing_t0"))
