"""Field kappa map (PLAN_FIELD2_ZH.md): two-anchor kappa and the implied offset beta on the 26 SR inlines.

usage: field_kappa2.py    worker; jobs claimed with results11/<job>.claim, resumable.
Lines: the 26 inlines on which experiment 10 measured the SR factor directly (results/real/lines.json of
specslope_q_volve_20261028, 9985-10369 every 16). Read, model Q, lambda and scan exactly as field_kappa.py;
default wavelet window only, c_a = 0.7 / 1.0 / 1.4.
-> results11/IL_<inline>__w1.4_a<c_a>.json
"""
import os, json
from field_kappa import run, L13, OUT as _OUT9, np, semi, bm, ROOT

OUT = ROOT / "results11"; bm.MAXITER = 20000
SRL = json.loads((L13.IN / "specslope_q_volve_20261028" / "results" / "real" / "lines.json").read_text(encoding="utf-8"))
LINES = sorted(int(k) for k in SRL)


def jobs():
    return [(ln, "w1.4", a) for ln in LINES for a in (0.7, 1.0, 1.4)]


if __name__ == "__main__":
    OUT.mkdir(exist_ok=True); cm = semi.curtain_cm(); Q = np.asarray(cm["Q"], float)
    for ln, w, ca in jobs():
        job = f"IL_{ln}__{w}_a{ca}"; path = OUT / f"{job}.json"
        if path.exists():
            continue
        try:
            os.close(os.open(OUT / f"{job}.claim", os.O_CREAT | os.O_EXCL | os.O_WRONLY))
        except FileExistsError:
            continue
        rec = run(ln, w, ca, cm, Q); rec["c_sr"] = float(SRL[str(ln)]["c"])
        path.write_text(json.dumps(rec, indent=1), encoding="utf-8"); (OUT / f"{job}.claim").unlink()
        print(job, f"v {rec['v_sel']:.3f}", "edge" if rec["edge"] else "", f"{rec['seconds']:.0f}s", flush=True)
