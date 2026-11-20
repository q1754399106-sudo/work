"""Field test of UQW-S (PLAN_ZH.md): UQW-S, UQ-W, VDB (frozen target), VDBin (VDB with UQW-S's robust target), HWin and IQ70
on the experiment-13 Volve inlines, baseline Q (curtain model x SR factor of the inline) multiplied by m in M_ALL.
field_qperturb_20261113/scripts/run.py with the method set of uqw_snr_20261114/scripts/confirmS.py (operators unchanged:
uqws.shape, uqw.uqw through confirmW, lib12.tvw, iq). Per inline: cache/out/IL_<L>.npz, results/lines/IL_<L>.json (bands,
self-correlation with m = 1 in 2.3-2.8 s, fs_r). Logs print no scores.
usage: run.py smoke      one inline, m = 1: timing, finiteness, UQW / VDB / IQ70 equal to the operators of the qperturb test
       run.py run K N    worker K of N over the analysis lines of experiment 13 (resumable)
"""
import sys, json, time, subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]; BASE = ROOT.parent; E13 = BASE / "field_baselines_volve_20261029"
sys.path.insert(0, str(E13 / "scripts"))
import lib13 as L  # noqa: E402
from lib13 import np, semi, sc11, lib, TW, DT  # noqa: E402
sys.path.insert(0, str(BASE / "uqw_snr_20261114" / "scripts"))
import confirmS as CS  # noqa: E402
import uqws as S  # noqa: E402
sys.path.insert(0, str(BASE / "qunknown_comp_20261109" / "scripts"))
import confirmW as CW  # noqa: E402

METHODS = ("UQWS", "UQW", "VDB", "VDBin", "HWin", "IQ70"); M_ALL = CS.MS; OUT = ROOT / "results" / "lines"; CACHE = ROOT / "cache" / "out"
F_OP = np.fft.rfftfreq(1024, DT); WS = (TW >= 2.3) & (TW <= 2.8)


def lines():
    txt = subprocess.run(["git", "show", "HEAD:./analysis/lines.json"], cwd=E13, capture_output=True, check=True).stdout.decode("utf-8")
    return json.loads(txt)["lines"]


def compensate(D, cm, Q, c, ms):
    w_frozen = sc11.prepare_sc(D, cm)["wv"]["w_out"]; inv = CW.snr_rows(D, TW, F_OP); w_t, fs_r = S.target(D, cm); out = {"Input": D}
    for m in ms:
        r = CS.methods(D, cm, Q, m * c, TW, 1.4, w_t, fs_r, w_frozen, inv)
        for k in METHODS:
            out[f"{k}_{m}"] = r[k]
    return out, fs_r, (w_frozen, inv)


def corr(a, b):
    a, b = a[WS].ravel(), b[WS].ravel(); return float(a @ b / np.sqrt((a @ a) * (b @ b)))


if __name__ == "__main__":
    cm = semi.curtain_cm(); Q = np.asarray(cm["Q"], float); mode = sys.argv[1]
    if mode == "smoke":
        ln = 10141; t0 = time.monotonic(); D = L.read_inline(ln); c = L.c_sr(ln); o, fs_r, (w, inv) = compensate(D, cm, Q, c, (1.0,)); dt = time.monotonic() - t0
        q = CW.methods(D, cm, Q, c, TW, 1.4, w, inv)                         # the operators of field_qperturb_20261113 (run.compensate)
        for a in ("UQW", "VDB", "IQ70"):
            print(a, "vs qperturb operator max |d| / max = %.1e" % (np.abs(o[f"{a}_1.0"] - q[a]).max() / np.abs(q[a]).max()))
        print("fs_r %.0f Hz; UQWS == HWin: %s; finite %s; %.0f s for one multiplier" % (fs_r, bool(np.array_equal(o["UQWS_1.0"], o["HWin_1.0"])),
              all(bool(np.isfinite(x).all()) for x in o.values()), dt))
    else:
        k, n = int(sys.argv[2]), int(sys.argv[3]); OUT.mkdir(parents=True, exist_ok=True); CACHE.mkdir(parents=True, exist_ok=True)
        for i, ln in enumerate(lines()):
            path = OUT / f"IL_{ln}.json"
            if i % n != k or path.exists():
                continue
            t0 = time.monotonic(); D = L.read_inline(ln); c = L.c_sr(ln); o, fs_r, _ = compensate(D, cm, Q, c, M_ALL)
            rec = dict(inline=ln, c_sr=c, fs_r=fs_r, bands={a: L.bands(x) for a, x in o.items()},
                       self_corr={f"{a}_{m}": corr(o[f"{a}_{m}"], o[f"{a}_1.0"]) for a in METHODS for m in M_ALL}, seconds=time.monotonic() - t0)
            np.savez_compressed(CACHE / f"IL_{ln}.npz", **{a: x.astype(np.float32) for a, x in o.items()})
            path.write_text(json.dumps(rec), encoding="utf-8"); print(ln, f"{rec['seconds']:.0f}s", flush=True)
