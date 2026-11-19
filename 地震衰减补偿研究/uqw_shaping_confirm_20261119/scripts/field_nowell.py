"""Field no-well description of PLAN_CONFIRM_SHAPING_ZH.md section 3.1: the 12 never-used inlines of lines6.json
(field_nowell). No well, no truth; descriptive only.

usage: field_nowell.py seal      read each inline (experiment-13 reader lib13.read_inline: xline 2100-2500 every 2nd, 1.2-3.0 s)
                                 and its SR factor (lib13.c_sr), save cache/field/IL_<L>.npz, write SEALED_FIELD.txt
                                 (SHA-256 of each saved file); no compensation
       field_nowell.py run K N   worker K of N: all methods of run_confirm.compensate for the 9 multipliers (Q = model Q x m x c_SR);
                                 per method and m: adjacent-trace correlation 2.3-2.8 s (median), coherent SNR in 15 Hz bands
                                 2.0-3.0 s and the 5-95% energy band (lib13.bands), share of energy above 80 Hz;
                                 for the input: 10-55 Hz curvature of ln(A 2.3-2.8 s / A 1.2-1.6 s) (specform definition)
                                 -> results/field_nowell/IL_<L>.json. Logs print progress only.
"""
import sys, json, time, hashlib
from pathlib import Path

HERE = Path(__file__).resolve().parent; ROOT = HERE.parent; BASE = ROOT.parent
LINES = json.loads((ROOT / "lines6.json").read_text(encoding="utf-8"))["field_nowell"]
CACHE = ROOT / "cache" / "field"; SEAL = ROOT / "SEALED_FIELD.txt"; OUT = ROOT / "results" / "field_nowell"


def seal():
    sys.path.insert(0, str(BASE / "field_baselines_volve_20261029" / "scripts"))
    import lib13 as L  # noqa: E402
    CACHE.mkdir(parents=True, exist_ok=True); rows = []
    for ln in LINES:
        D = L.read_inline(ln); p = CACHE / f"IL_{ln}.npz"; np_ = L.np
        np_.savez(p, D=D, c_sr=np_.array(L.c_sr(ln)), dead=np_.array(int((np_.abs(D).sum(0) == 0).sum())))
        rows.append(f"{hashlib.sha256(p.read_bytes()).hexdigest()}  cache/field/IL_{ln}.npz"); print(ln, "read", D.shape, flush=True)
    SEAL.write_text("\n".join(rows) + "\n", encoding="utf-8")


def run(k, n):
    import run_confirm as RC  # noqa: E402
    from run_confirm import np, CW, TW, DT, TA, probe, semi  # noqa: E402
    sys.path.insert(0, str(BASE / "semisynth_fault_q6_20261119" / "scripts"))
    import fq6  # noqa: E402
    table = {r.split("  ")[1]: r.split("  ")[0] for r in SEAL.read_text(encoding="utf-8").split("\n") if r.strip()}
    cm = semi.curtain_cm(); Q = np.asarray(cm["Q"], float); OUT.mkdir(parents=True, exist_ok=True)

    def amp(D, a, c):
        kk = (TW >= a) & (TW <= c); live = np.abs(D).sum(0) > 0
        return np.sqrt(TA.lib.smooth(np.mean(np.abs(np.fft.rfft(D[kk][:, live] * np.hanning(kk.sum())[:, None], fq6.NF, axis=0)) ** 2, axis=1), 5))
    for i, ln in enumerate(LINES):
        path = OUT / f"IL_{ln}.json"
        if i % n != k or path.exists():
            continue
        p = CACHE / f"IL_{ln}.npz"
        if hashlib.sha256(p.read_bytes()).hexdigest() != table[f"cache/field/IL_{ln}.npz"]:
            raise SystemExit(f"seal mismatch {ln}")
        t0 = time.monotonic(); z = np.load(p); D = z["D"].astype(float); c = float(z["c_sr"])
        w_out = CW.sc11.prepare_sc(D, cm)["wv"]["w_out"]; inv = CW.snr_rows(D, TW, RC.F1024)
        desc = lambda X: dict(latcoh=TA.lateral_coherence(X), bands=TA.L.bands(X), hf=probe.hf(X))  # noqa: E731
        rec = dict(inline=ln, c_sr=c, dead=int(z["dead"]), input=desc(D),
                   input_curvature=fq6.curvature(np.log(amp(D, 2.3, 2.8) / amp(D, 1.2, 1.6))), methods={})
        for m in RC.M_ALL:
            for meth, X in RC.compensate(D, cm, Q, m * c, inv, w_out).items():
                rec["methods"][f"{meth}_{m}"] = dict(desc(X), finite=bool(np.isfinite(X).all()))
        path.write_text(json.dumps(rec), encoding="utf-8"); print(ln, "%.0fs" % (time.monotonic() - t0), flush=True)


if __name__ == "__main__":
    sys.path.insert(0, str(HERE))
    if sys.argv[1] == "seal":
        seal()
    else:
        run(int(sys.argv[2]), int(sys.argv[3]))
