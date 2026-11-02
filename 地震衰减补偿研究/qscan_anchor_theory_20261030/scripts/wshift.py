"""Wavelet-window shift (PLAN_WSHIFT_ZH.md): anchor scans with the wavelet window at 1.3-1.7, 1.4-1.8, 1.5-1.9 s.

usage: wshift.py    worker; jobs claimed with results25/<job>.claim, resumable.
Field: 6 of the 26 SR inlines (field_kappa2.LINES[::5]), read as field_kappa.py. Semi: the 12 L realizations of
beta_color.py (rebuilt as betasr.py). c_a = 0.7 / 1.0 / 1.4, lambda 1/2, qprofile.scan. -> results25/<job>.json
"""
import os, json, time
from env15 import np, lib, semi, sc11, bm, gen, ROOT
from anchor2 import wavelets_at
import qprofile as qp
from field_kappa2 import LINES, L13
import beta_color as bc

OUT = ROOT / "results25"; bm.MAXITER = 20000; qp.LAM = 0.5
WINS = {"w1.5": (1.3, 1.7, 1.5), "w1.6": (1.4, 1.8, 1.6), "w1.7": (1.5, 1.9, 1.7)}; CAS = (0.7, 1.0, 1.4)
FIELD = LINES[::5]; SEMI = [j for j in bc.JOBS if j[0] == "L"]


def scans(D, cm, Q, win):
    pp = sc11.prepare_sc(D, cm); out = {}
    for ca in CAS:
        wv = wavelets_at(D, semi.TW, dict(cm, Q=Q * ca), *WINS[win]); r = qp.scan(D, cm, Q, pp, lambda cc: wv)
        out[str(ca)] = dict(v_sel=r["v_sel"], edge=r["edge"], ok=r["neighbourhood_converged"])
    return out


def jobs():
    return [(f"F_{ln}__{w}", ("f", ln, w)) for w in WINS for ln in FIELD] + [(f"S_{j}__{w}", ("s", j, w)) for w in WINS for j in SEMI]


if __name__ == "__main__":
    OUT.mkdir(exist_ok=True); cm = semi.curtain_cm(); Q = np.asarray(cm["Q"], float); ctx = None
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
            D = L13.read_inline(spec[1]); rec = dict(kind="field", inline=spec[1], window=spec[2])
            cmu = cm
        else:
            ctx = ctx or gen.prepare(10227, [10149, 10273, 10355]); seed = int(spec[1][1:])
            g = gen.signal(ctx, 1.0, seed); D = g["x"] + g["level"] * ctx["noise"][0]; cmu = ctx["cm"]
            rec = dict(kind="semi", job=spec[1], window=spec[2])
        rec["scans"] = scans(D, cmu, np.asarray(cmu["Q"], float), spec[2]); rec["seconds"] = time.monotonic() - t0
        path.write_text(json.dumps(rec, indent=1), encoding="utf-8"); (OUT / f"{job}.claim").unlink()
        print(job, {k: round(v["v_sel"], 3) for k, v in rec["scans"].items()}, f"{rec['seconds']:.0f}s", flush=True)
