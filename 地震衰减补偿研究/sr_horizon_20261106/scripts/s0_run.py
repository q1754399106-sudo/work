"""Stage 0 of the along-horizon SR (STAGE0_ZH.md). usage: s0_run.py debug | geometry | precision | trend [s c0 r] | colour
Results go to results/*.json. 'debug' only checks that everything runs on ZA_c1.05_r0 and prints no estimates."""
import json
import sys
import time

from env0 import np, lib, semi, sq, fq, OUT, RES, REPO, COMMIT, git_module, subprocess
import hsr


def save(name, obj):
    (OUT / name).write_text(json.dumps(obj, indent=1), encoding="utf-8")


def model():
    cm = semi.curtain_cm(); return cm, np.asarray(cm["Q"], float)


def ids():
    return [f"{v}A_c{c}_r{r}" for v in "ZY" for c in fq.C0 for r in (0, 1)]


def one(D, cm, Q, truth):
    sr = sq.estimate(D, cm, Q)
    h = hsr.estimate(D, cm, Q, sr["f_hi"])
    return dict(truth=truth, sr=sr["c"], f_hi=sr["f_hi"], hsr=h["c"], rms=h.get("rms"), curvature=h.get("curvature"),
                at_grid_edge=h.get("at_grid_edge"), windows=len(h["used"]), used=h["used"],
                log_err_sr=float(np.log(sr["c"] / truth)), log_err_hsr=float(np.log(h["c"] / truth)))


def debug():
    cm, Q = model(); arr, meta = fq.load("ZA_c1.05_r0")
    t0 = time.monotonic(); r = one(arr["D"].astype(float), cm, Q, meta["c_eff_dif"])
    print("ran in", round(time.monotonic() - t0, 1), "s; windows used", r["windows"], "edge", r["at_grid_edge"], flush=True)
    b = hsr.tilt(np.random.default_rng(0).standard_normal((len(fq.TF), 3)), 0.0099)
    print("tilt ok", b.shape, flush=True)


def geometry():
    L = json.loads(subprocess.check_output(["git", "show", f"{COMMIT}:{RES.name}/specslope_q_volve_20261028/results/real/lines.json"],
                                           cwd=REPO).decode("utf-8"))
    lines = sorted(int(k) for k in L) + [10141]
    lines = sorted(set(lines)); data = semi.read_lines(lines); out = {}
    for line, D in data.items():
        tau = hsr.structure(D); row = []
        for a, b in sq.windows():
            tb, valid = hsr.window_shifts(tau, a, b)
            w = semi.win(semi.TW, a, b); X = hsr.flatten(D, tb)[w]
            x = X - X.mean(axis=0, keepdims=True)
            cc = np.sum(x[:, :-1] * x[:, 1:], axis=0) / np.sqrt(np.sum(x[:, :-1] ** 2, axis=0) * np.sum(x[:, 1:] ** 2, axis=0))
            row.append(dict(win=[a, b], spread=hsr.spread(tb, valid), valid=float(valid.mean()),
                            pair_corr=float(np.median(cc)), tb_min=float(tb.min()), tb_max=float(tb.max())))
        out[str(line)] = row; print(line, [round(r["spread"], 3) for r in row], flush=True)
    S = np.array([r["spread"] for v in out.values() for r in v])
    res = dict(lines=out, n=int(S.size), frac_ge_0p1=float(np.mean(S >= 0.1)), quantiles={str(q): float(np.quantile(S, q)) for q in (0.1, 0.25, 0.5, 0.75, 0.9)},
               g1_pass=bool(np.mean(S >= 0.1) >= 0.5),
               per_window={f"{a}-{b}": float(np.median([v[k]["spread"] for v in out.values()])) for k, (a, b) in enumerate(sq.windows())})
    save("geometry.json", res); print({k: v for k, v in res.items() if k != "lines"}, flush=True)


def precision():
    cm, Q = model(); out = {}
    for did in ids():
        t0 = time.monotonic(); arr, meta = fq.load(did); truth = meta["c_eff_dif"]
        noisy = one(arr["D"].astype(float), cm, Q, truth); clean = one(arr["signal"].astype(float), cm, Q, truth)
        out[did] = dict(noisy=noisy, clean=clean, c_eff_cum=meta["c_eff_cum"])
        print(did, f"{time.monotonic() - t0:.0f}s", flush=True)
        save("precision.json", out)
    e = np.array([v["noisy"]["log_err_hsr"] for v in out.values()])
    sm = {k: dict(sd_hsr=float(np.std([v[k]["log_err_hsr"] for v in out.values()], ddof=1)),
                  mean_hsr=float(np.mean([v[k]["log_err_hsr"] for v in out.values()])),
                  rms_hsr=float(np.sqrt(np.mean([v[k]["log_err_hsr"] ** 2 for v in out.values()]))),
                  sd_sr=float(np.std([v[k]["log_err_sr"] for v in out.values()], ddof=1)),
                  mean_sr=float(np.mean([v[k]["log_err_sr"] for v in out.values()])),
                  rms_sr=float(np.sqrt(np.mean([v[k]["log_err_sr"] ** 2 for v in out.values()])))) for k in ("noisy", "clean")}
    sm["g2_pass"] = bool(np.std(e, ddof=1) <= 0.20); sm["n"] = int(e.size)
    save("precision_summary.json", sm); print(json.dumps(sm, indent=1), flush=True)


def trend(args):
    """Datasets with a colour trend, built with the unmodified generator after tilting its base reflectivity."""
    grid = [(float(s), c, r) for s in (0.0033, -0.0033, 0.0099, -0.0099) for c in (0.75, 1.05) for r in (0, 1)]
    if args:
        grid = [(float(args[0]), float(args[1]), int(args[2]))]
    cm, Q = model(); ctx = fq.prepare("A"); orig = fq.base_series; p = OUT / "trend.json"
    out = json.loads(p.read_text(encoding="utf-8")) if p.exists() else {}
    for s, c0, r in grid:
        key = f"s{s:+.4f}_c{c0}_r{r}"
        if key in out:
            continue
        t0 = time.monotonic(); fq.base_series = lambda rng: hsr.tilt(orig(rng), s)
        built = fq.build(ctx, c0, r); row = {}
        for var, o in built.items():
            m, a = o["meta"], o["arrays"]; truth = m["c_eff_dif"]
            row[var] = dict(noisy=one(a["D"], cm, Q, truth), clean=one(a["signal"], cm, Q, truth), level=m["level"], f_true=m["f_true"])
        out[key] = row; save("trend.json", out); print(key, f"{time.monotonic() - t0:.0f}s", flush=True)


def colour():
    """Slope of ln amplitude spectrum over 10-60 Hz in 0.2 s windows along each well, and its trend with time."""
    wells = git_module("wells", "specslope_q_volve_20261028", "wells.py")
    segs = wells.all_segments(); res = {}
    f = np.fft.rfftfreq(512, 0.001); m = (f >= 10) & (f <= 60); L = 200; H = 100; tp = lib.taper(L)

    def slopes(rc, t0):
        ts, ss = [], []
        for st in range(0, len(rc) - L + 1, H):
            P = np.abs(np.fft.rfft(rc[st:st + L] * tp, 512)) ** 2
            ss.append(float(np.polyfit(f[m], 0.5 * np.log(lib.smooth(P, 9)[m] + 1e-30), 1)[0])); ts.append(st * 0.001 if t0 is None else t0 + st * 0.001)
        return np.array(ts), np.array(ss)

    pooled = []
    for w, sg in segs.items():
        rows = []
        for rc in sg:
            if len(rc) < 400:
                continue
            t, s = slopes(rc, None); b = float(np.polyfit(t, s, 1)[0]); rows.append(dict(seconds=len(rc) * 0.001, trend=b, sd=float(np.std(s, ddof=1)), n=int(len(s))))
            pooled.append((len(rc) * 0.001, b))
        if rows:
            res[w] = rows
    P = np.array(pooled) if pooled else np.zeros((0, 2))
    summary = dict(wells=res, n_segments=int(len(P)), median_abs_trend=float(np.median(np.abs(P[:, 1]))) if len(P) else None,
                   length_weighted_mean=float(np.average(P[:, 1], weights=P[:, 0])) if len(P) else None,
                   sd_trend=float(np.std(P[:, 1], ddof=1)) if len(P) > 1 else None, note="time axis is the time within a segment; only the slope with time is used")
    save("colour.json", summary); print(json.dumps({k: v for k, v in summary.items() if k != "wells"}, indent=1), flush=True); print(json.dumps(res), flush=True)


if __name__ == "__main__":
    mode = sys.argv[1]
    if mode == "debug":
        debug()
    elif mode == "geometry":
        geometry()
    elif mode == "precision":
        precision()
    elif mode == "trend":
        trend(sys.argv[2:])
    elif mode == "colour":
        colour()
