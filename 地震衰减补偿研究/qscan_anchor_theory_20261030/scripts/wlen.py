"""Wavelet-window length (PLAN_WLEN_ZH.md): does a longer wavelet window shrink the realisation spread of beta?

usage: wlen.py [--smoke]    worker; jobs claimed with results27/<job>.claim, resumable.
Windows (t1, t2, tc): W04_14 1.2-1.6/1.4 (default), W04_16 1.4-1.8/1.6, W08_16 1.2-2.0/1.6, W04_18 1.6-2.0/1.8,
W12_18 1.2-2.4/1.8. Semi: beta_color.py L recipe (structure 10227, noise line 10149 x 1, c = 1.0); old seeds
10001-10012 run the three new windows (W04_14 in results16, W04_16 in results25), new seeds 10013-10024 run all
five and get c_SR (sq.estimate, as betasr.py). Field: the six wshift.py lines, three new windows.
Scans: c_a 0.7 / 1.0 / 1.4, lambda 1/2, qprofile.scan with anchor2.wavelets_at. -> results27/<job>.json
--smoke: L10001 with wavelets_at(1.2, 1.6, 1.4) must reproduce results16/L10001.json (|dv| <= 1e-6).
"""
import os, sys, json, time
from env15 import np, lib, semi, sc11, bm, gen, sq, ROOT
from anchor2 import wavelets_at
import qprofile as qp
from field_kappa2 import L13

OUT = ROOT / "results27"; bm.MAXITER = 20000; qp.LAM = 0.5
WINS = {"W04_14": (1.2, 1.6, 1.4), "W04_16": (1.4, 1.8, 1.6), "W08_16": (1.2, 2.0, 1.6),
        "W04_18": (1.6, 2.0, 1.8), "W12_18": (1.2, 2.4, 1.8)}
NEW = ("W08_16", "W04_18", "W12_18"); CAS = (0.7, 1.0, 1.4)
OLD_SEEDS = list(range(10001, 10013)); NEW_SEEDS = list(range(10013, 10025))
FIELD = [9985, 10065, 10141, 10209, 10289, 10369]


def scans(D, cm, Q, win):
    pp = sc11.prepare_sc(D, cm); out = {}
    for ca in CAS:
        wv = wavelets_at(D, semi.TW, dict(cm, Q=Q * ca), *WINS[win]); r = qp.scan(D, cm, Q, pp, lambda cc: wv)
        out[str(ca)] = dict(v_sel=r["v_sel"], edge=r["edge"], ok=r["neighbourhood_converged"])
    return out


def jobs():
    # seed-major so that a partial run covers whole realisations; field last
    J = [(f"S_L{s}__{w}", ("s", s, w)) for s in OLD_SEEDS for w in NEW]
    J += [(f"S_L{s}__{w}", ("s", s, w)) for s in NEW_SEEDS for w in WINS]
    J += [(f"F_{ln}__{w}", ("f", ln, w)) for ln in FIELD for w in NEW]
    return J


_ctx = None


def semi_data(seed):
    global _ctx
    _ctx = _ctx or gen.prepare(10227, [10149, 10273, 10355])
    g = gen.signal(_ctx, 1.0, seed); return g["x"] + g["level"] * _ctx["noise"][0], _ctx["cm"]


if __name__ == "__main__":
    OUT.mkdir(exist_ok=True)
    if "--smoke" in sys.argv:
        D, cm = semi_data(10001); Q = np.asarray(cm["Q"], float); t0 = time.monotonic()
        new = scans(D, cm, Q, "W04_14"); old = json.loads((ROOT / "results16" / "L10001.json").read_text())["scans"]
        d = {k: new[k]["v_sel"] - old[k]["v_sel"] for k in new}
        print("smoke dv", d, "max |dv| %.2e" % max(abs(x) for x in d.values()), f"{time.monotonic() - t0:.0f}s", flush=True)
        sys.exit(0 if max(abs(x) for x in d.values()) <= 1e-6 else 1)
    cmf = semi.curtain_cm()
    for job, spec in jobs():
        path = OUT / f"{job}.json"
        if path.exists():
            continue
        try:
            os.close(os.open(OUT / f"{job}.claim", os.O_CREAT | os.O_EXCL | os.O_WRONLY))
        except FileExistsError:
            continue
        t0 = time.monotonic()
        if spec[0] == "f":
            D = L13.read_inline(spec[1]); cm = cmf; rec = dict(kind="field", inline=spec[1], window=spec[2])
        else:
            D, cm = semi_data(spec[1]); rec = dict(kind="semi", seed=spec[1], window=spec[2])
            if spec[2] == "W04_14" and spec[1] in NEW_SEEDS:
                rec["c_sr"] = float(sq.estimate(D, cm, np.asarray(cm["Q"], float))["c"])
        rec["scans"] = scans(D, cm, np.asarray(cm["Q"], float), spec[2]); rec["seconds"] = time.monotonic() - t0
        path.write_text(json.dumps(rec, indent=1), encoding="utf-8"); (OUT / f"{job}.claim").unlink()
        print(job, {k: round(v["v_sel"], 3) for k, v in rec["scans"].items()}, f"{rec['seconds']:.0f}s", flush=True)
