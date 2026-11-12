"""Field Q-perturbation test (PLAN_ZH.md): UQ-W, VDB, IQ70 on the experiment-13 Volve inlines with the baseline Q
(curtain model x SR factor of the inline) multiplied by m in M_ALL.

Data access, well geometry and the SR factor: experiment 13 (lib13, unchanged). Operators: confirmW.methods (unchanged;
UQ-W sig 0.447, gamma 0, anchor 1.4 s; VDB and IQ70 with the experiment-12 parameters).
Per inline: cache/out/IL_<L>.npz (Input and <method>_<m>), results/lines/IL_<L>.json (bands, self-correlation with m = 1
in 2.3-2.8 s). Logs print no scores.

usage: run.py smoke      one inline, m = 1: timing and equality with the archived experiment-13 outputs; saves nothing
       run.py run K N    worker K of N over analysis lines of experiment 13 (resumable)
"""
import sys, json, time, subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]; BASE = ROOT.parent; E13 = BASE / "field_baselines_volve_20261029"
sys.path.insert(0, str(E13 / "scripts"))
import lib13 as L  # noqa: E402
from lib13 import np, semi, sc11, lib, TW, DT  # noqa: E402
sys.path.insert(0, str(BASE / "qunknown_comp_20261109" / "scripts"))
import confirmW as CW  # noqa: E402

METHODS = ("UQW", "VDB", "IQ70"); M_ALL = CW.M_ALL; OUT = ROOT / "results" / "lines"; CACHE = ROOT / "cache" / "out"; F_OP = np.fft.rfftfreq(1024, DT)
WS = (TW >= 2.3) & (TW <= 2.8)


def lines():
    txt = subprocess.run(["git", "show", "HEAD:./analysis/lines.json"], cwd=E13, capture_output=True, check=True).stdout.decode("utf-8")
    return json.loads(txt)["lines"]


def compensate(D, cm, Q, c, ms):
    w_out = sc11.prepare_sc(D, cm)["wv"]["w_out"]; inv = CW.snr_rows(D, TW, F_OP); out = {"Input": D}
    for m in ms:
        r = CW.methods(D, cm, Q, m * c, TW, 1.4, w_out, inv)
        for k in METHODS:
            out[f"{k}_{m}"] = r[k]
    return out


def corr(a, b):
    a, b = a[WS].ravel(), b[WS].ravel(); return float(a @ b / np.sqrt((a @ a) * (b @ b)))


if __name__ == "__main__":
    cm = semi.curtain_cm(); Q = np.asarray(cm["Q"], float); mode = sys.argv[1]
    if mode == "smoke":
        ln = 10141; t0 = time.monotonic(); D = L.read_inline(ln); o = compensate(D, cm, Q, L.c_sr(ln), (1.0,)); dt = time.monotonic() - t0
        A = np.load(E13 / "cache" / "out" / f"IL_{ln}_A.npz")
        for mine, theirs in (("Input", "Input"), ("VDB_1.0", "B3_VDB"), ("IQ70_1.0", "B1_IQ70")):
            print(mine, "max |difference| / max |archived| = %.2e" % (np.abs(o[mine] - A[theirs]).max() / np.abs(A[theirs]).max()))
        print("shape", D.shape, "dead traces", int((np.abs(D).sum(0) == 0).sum()), "finite", bool(np.isfinite(o["UQW_1.0"]).all()), "seconds for one multiplier %.0f" % dt)
    else:
        k, n = int(sys.argv[2]), int(sys.argv[3]); OUT.mkdir(parents=True, exist_ok=True); CACHE.mkdir(parents=True, exist_ok=True)
        for i, ln in enumerate(lines()):
            path = OUT / f"IL_{ln}.json"
            if i % n != k or path.exists():
                continue
            t0 = time.monotonic(); D = L.read_inline(ln); c = L.c_sr(ln); o = compensate(D, cm, Q, c, M_ALL)
            rec = dict(inline=ln, c_sr=c, bands={a: L.bands(x) for a, x in o.items()},
                       self_corr={f"{a}_{m}": corr(o[f"{a}_{m}"], o[f"{a}_1.0"]) for a in METHODS for m in M_ALL}, seconds=time.monotonic() - t0)
            np.savez_compressed(CACHE / f"IL_{ln}.npz", **{a: x.astype(np.float32) for a, x in o.items()})
            path.write_text(json.dumps(rec, indent=1), encoding="utf-8"); print(ln, f"{rec['seconds']:.0f}s", flush=True)
