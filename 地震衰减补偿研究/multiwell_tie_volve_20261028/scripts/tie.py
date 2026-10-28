"""Multi-well tie of experiment 9 (PLAN_ZH.md): does compensation raise the well-tie correlation on real data?

Per well: the AC/DEN reflectivity (5-point running medians, gaps up to GAP_TIE bridged) in two-way time from the
TVD-based sonic stretched to the picks (wq.well), kept between the first and last pick that carries a TWT (anchored
span). Tie reference: RC * Ricker(fpeak of the experiment-6 curtain model), as the T2-like tie of
experiment 8. Seismic: at every 4 ms sample, the mean of the MIG_RAW trace nearest to the hole and its four
neighbours (inline +-2, xline +-2), each trace compensated as a whole (1.2-3.0 s) before averaging. Score: Pearson
correlation over the anchored span, best of +-20 ms shifts (lib.tie).

usage: tie.py groups    wells, spans and well groups (no seismic)
       tie.py run       read traces, compensate, tie; writes results/tie.json
"""
import sys, json
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import wq  # noqa: E402
from wq import np, semi, lib, TW  # noqa: E402

TIE_WELLS = ["F-1", "F-1_A", "F-11_A", "F-11_T2", "F-12", "F-14", "F-15", "F-15_A", "F-15_B", "F-4", "F-5"]
SEEN = ["F-11_T2", "F-15_A"]             # tied in experiments 5-7
GAP_TIE = 5.0                            # m
MIN_SPAN = 0.15                          # s, anchored span kept
GROUP_DIST = 250.0                       # m, single-linkage distance of span midpoints
G_DB = 30
OUT = wq.ROOT / "results"


def tie_log(name, P):
    """AC/DEN RC on the time axis of wq.well (sonic over the whole AC log, stretched to the picks): the AC&DEN piece
    (gaps up to GAP_TIE bridged) with the longest overlap with the anchored MD range, cut to the anchored span."""
    W = wq.well(name, P); am = np.array([a[0] for a in W["anchors"]]); at = np.array([a[1] for a in W["anchors"]])
    A = np.loadtxt(wq.E9 / "inventory" / "logs" / f"15_9-{name}.csv", delimiter=",", skiprows=1)
    md, ac, den = A.T; ok = np.isfinite(ac) & np.isfinite(den); md, ac, den = md[ok], ac[ok], den[ok]
    pieces = np.split(np.arange(len(md)), np.where(np.diff(md) > GAP_TIE)[0] + 1)
    ov = [max(0.0, min(md[p[-1]], am.max()) - max(md[p[0]], am.min())) for p in pieces]
    s = pieces[int(np.argmax(ov))]; md, ac, den = md[s], wq.running_median(ac[s]), wq.running_median(den[s])
    t = np.interp(md, W["md"], W["t"])                   # W["md"] is increasing along W["t"]
    lo, hi = max(t[0], at.min()), min(t[-1], at.max())
    tt = np.round(np.arange(np.ceil(lo / wq.FINE), np.floor(hi / wq.FINE) + 1) * wq.FINE, 6)
    z = np.log(np.interp(tt, t, 0.3048 / (ac * 1e-6) * den))
    return dict(name=name, t=tt, rc=0.5 * np.diff(z), E=np.interp(tt, W["t"], W["E"]), N=np.interp(tt, W["t"], W["N"]),
                t0=float(tt[0]), t1=float(tt[-1]), n_anchors=int(((am >= md[0]) & (am <= md[-1])).sum()))


def groups(logs):
    """Single-linkage clusters of span midpoints closer than GROUP_DIST; names sorted, groups by first member."""
    names = sorted(logs); mid = {n: (np.mean(logs[n]["E"]), np.mean(logs[n]["N"])) for n in names}
    lab = {n: n for n in names}
    for a in names:
        for b in names:
            if a < b and np.hypot(mid[a][0] - mid[b][0], mid[a][1] - mid[b][1]) < GROUP_DIST:
                la, lb = lab[a], lab[b]
                for n in names:
                    if lab[n] == lb:
                        lab[n] = la
    out = {}
    for n in names:
        out.setdefault(lab[n], []).append(n)
    return {g[0]: g for g in sorted(out.values())}, mid


def cube():
    snap = json.loads((wq.ROOT / "cache" / "SNAPSHOT.json").read_text(encoding="utf-8"))
    d = np.memmap(snap["cube"], dtype=np.uint8, mode="r"); ns = 1126; tb = 240 + 4 * ns; n = (len(d) - 3600) // tb
    h = np.asarray(d[3600:3600 + n * tb]).reshape(n, tb)
    sc = h[:, 70:72].copy().view(">i2").ravel().astype(float); sc = np.where(sc < 0, -1 / sc, np.where(sc == 0, 1, sc))
    x = h[:, 180:184].copy().view(">i4").ravel() * sc; y = h[:, 184:188].copy().view(">i4").ravel() * sc
    il = h[:, 188:192].copy().view(">i4").ravel(); xl = h[:, 192:196].copy().view(">i4").ravel()
    return h, x, y, il, xl


def path(log, x, y, il, xl, pad=0.024):
    """For each TW sample of the span +-pad (positions clamped to the span): cube index of the nearest trace and
    its four neighbours."""
    w = semi.win(TW, log["t0"] - pad, log["t1"] + pad); key = {(a, b): i for i, (a, b) in enumerate(zip(il, xl))}; rows = []
    for t in TW[w]:
        e, n = np.interp(t, log["t"], log["E"]), np.interp(t, log["t"], log["N"]); j = int(np.argmin((x - e) ** 2 + (y - n) ** 2))
        rows.append([key.get((il[j] + di, xl[j] + dx), j) for di, dx in ((0, 0), (-2, 0), (2, 0), (0, -2), (0, 2))])
    return w, np.array(rows)


def sr_factor(inline):
    """Experiment-10 SR factor interpolated in inline number (results/real/lines.json, 346b562)."""
    L = json.loads((wq.IN / "specslope_q_volve_20261028" / "results" / "real" / "lines.json").read_text(encoding="utf-8"))
    k = sorted((int(a), v["c"]) for a, v in L.items() if a != "10141")
    return float(np.exp(np.interp(inline, [a for a, _ in k], [np.log(c) for _, c in k])))


def show_groups():
    P = wq.picks(); logs = {n: tie_log(n, P) for n in TIE_WELLS}
    for n, g in logs.items():
        print(n, round(g["t0"], 3), round(g["t1"], 3), round(g["t1"] - g["t0"], 3), "anchors", g["n_anchors"])
    G, mid = groups({n: g for n, g in logs.items() if g["t1"] - g["t0"] >= MIN_SPAN})
    print(G); print({n: (round(m[0]), round(m[1])) for n, m in mid.items()})


def run():
    P = wq.picks(); logs = {n: tie_log(n, P) for n in TIE_WELLS}
    logs = {n: g for n, g in logs.items() if g["t1"] - g["t0"] >= MIN_SPAN}; G, mid = groups(logs)
    h, x, y, il, xl = cube(); cm = semi.curtain_cm(); Q = np.asarray(cm["Q"], float); fp = float(cm["fpeak"])
    wt = semi.win(np.arange(1126) * wq.DT, semi.T0, semi.T1); rows = {}
    for n, g in logs.items():
        w, idx = path(g, x, y, il, xl); uniq = np.unique(idx)
        D = np.stack([semi.ibm(h[j, 240:].copy().view(">u4")) for j in uniq], axis=1)[wt]
        jm = int(np.argmin((x - mid[n][0]) ** 2 + (y - mid[n][1]) ** 2)); c_sr = sr_factor(int(il[jm]))
        outs = {"Input": D}
        for k, c in (("IQ30_SR", c_sr), ("IQ30_0.77", 0.77), ("IQ30_1.0", 1.0)):
            outs[k] = semi.inverse_q(D, cm, Q * c, G_DB)
        pos = {j: i for i, j in enumerate(uniq)}
        rf = np.zeros(len(semi.TF)); tm = np.round(0.5 * (g["t"][1:] + g["t"][:-1]) - 0.0005, 6)
        i0 = int(round((tm[0] - semi.TF[0]) / wq.FINE)); rf[i0:i0 + len(g["rc"])] = g["rc"]
        ref = semi.ricker_ref(rf, fp); sc, trs = {}, {}; ws = semi.win(TW, g["t0"], g["t1"])
        for k, A in outs.items():
            tr = np.zeros(len(TW)); tr[w] = [A[ii, [pos[j] for j in r]].mean() for ii, r in zip(np.where(w)[0], idx)]
            trs[k] = tr; sc[k] = dict(zip(("T2", "shift_ms"), lib.tie(tr, ref, TW, ws)))
        for k in outs:
            if k != "Input":
                sc[k]["boot90_vs_input"] = lib.tie_boot(trs[k], trs["Input"], ref, TW, ws, sc[k]["shift_ms"],
                                                        sc["Input"]["shift_ms"])
        (wq.ROOT / "cache").mkdir(exist_ok=True)
        np.savez_compressed(wq.ROOT / "cache" / f"tie_{n}.npz", ref=ref, **{k: v for k, v in trs.items()})
        rows[n] = dict(span=[g["t0"], g["t1"]], inline_mid=int(il[jm]), c_sr=c_sr, n_traces=int(len(uniq)), scores=sc)
        print(n, "done", flush=True)                     # no scores in the log
    OUT.mkdir(exist_ok=True)
    bm_json = dict(groups=G, seen=SEEN, fpeak=fp, wells=rows)
    (OUT / "tie.json").write_text(json.dumps(bm_json, indent=1), encoding="utf-8")


if __name__ == "__main__":
    {"groups": show_groups, "run": run}[sys.argv[1]]()
