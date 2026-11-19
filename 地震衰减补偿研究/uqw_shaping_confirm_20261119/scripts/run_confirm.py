"""Compensation and scoring of PLAN_CONFIRM_SHAPING_ZH.md on the sixth semi-synthetic batch.

Cases: 24 base datasets (fq6.IDS) x spectral variant K0 / TV / ST x noise state cal (D = signal' + noise) / nd (signal' + 4 noise)
= 144; each with the 9 Q multipliers M_ALL (Q = curtain-model Q x m x c_eff_cum).
Solver side: observed() returns D only; compensate() receives D, cm, Q, c, the SNR estimate and w_out and nothing else.
Truth (ideal, R) is loaded by score() after compensate() has returned. The source text of compensate() and shaping.py is
checked for truth array names at start (logged). Seal: every dataset file must match SEALED_SHA256.txt before use.
Logs print progress and timing only, never scores. Output: results/cases/<case>.json (resumable).

usage: run_confirm.py smoke DEV_ID   used development set (e.g. ZA_c0.75_r0 of batch 1): code paths, finiteness, timing,
                                     HAW equal to dev6.anchor_shape; no score is computed or saved
       run_confirm.py run K N        worker K of N
"""
import sys, json, time, hashlib, inspect
from pathlib import Path

HERE = Path(__file__).resolve().parent; ROOT = HERE.parent; BASE = ROOT.parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(BASE / "semisynth_fault_q6_20261119" / "scripts"))
sys.path.insert(0, str(BASE / "uqw_improve_20261113" / "scripts"))
import dev6  # noqa: E402  (anchor_shape = HAW)
import probe  # noqa: E402  (t2, hf)
from probe import CW, np, semi, TW, DT  # noqa: E402


def _import_tie_ablation():
    """I2's tie_ablation imports field_qperturb's run.py under the module name "run", which sc11 (already imported, via
    confirmW) also uses for experiment 10's run.py; lib13 also re-points ROOT / IN of shared modules. Import it with its own
    "run", then restore the module table, the path and those attributes, so the solver modules are exactly as confirmW set them."""
    keep_run = sys.modules.pop("run", None); saved = {}
    for name, mod in list(sys.modules.items()):
        for a in ("ROOT", "IN", "E9", "E10"):
            if hasattr(mod, a) and not name.startswith("_") and name not in ("__main__",):
                saved[(name, a)] = getattr(mod, a)
    p = str(BASE / "field_qperturb_20261113" / "scripts"); sys.path.insert(0, p)
    import tie_ablation  # noqa: E402
    sys.path.remove(p); sys.modules.pop("run", None)
    if keep_run is not None:
        sys.modules["run"] = keep_run
    for (name, a), v in saved.items():
        setattr(sys.modules[name], a, v)
    return tie_ablation


TA = _import_tie_ablation()   # ti_shape of I2 (UQW_TIw) and lateral_coherence, unchanged
import shaping  # noqa: E402
import fq6  # noqa: E402

M_ALL = (0.5, 0.6, 0.7, 0.85, 1.0, 1.2, 1.4, 1.7, 2.0); SIG, GAMMA, ANCHOR = 0.447, 0.0, 1.4; NOISE_X = 4.0
VARS = {"K0": "signal", "TV": "signal_TV", "ST": "signal_ST"}; STATES = {"cal": 1.0, "nd": NOISE_X}
CASES = [f"{d}__{v}__{s}" for d in fq6.IDS for v in VARS for s in STATES]
F1024 = np.fft.rfftfreq(1024, DT); EPS, WIN = CW.SEL["B3"]["eps"], CW.SEL["B3"]["win"]
OUT = ROOT / "results" / "cases"; FORBIDDEN = ('"ideal"', '"R"', '"signal', '"noise"', '"tstar"', '"invQ"', '"src"', "c_eff")


def sealed():
    rows = (fq6.ROOT / "SEALED_SHA256.txt").read_text(encoding="utf-8").split("\n")
    return {r.split("  ")[1]: r.split("  ")[0] for r in rows if r.strip()}


def check_seal(did, table):
    p = fq6.ROOT / "cache" / "data" / f"{did}.npz"
    if hashlib.sha256(p.read_bytes()).hexdigest() != table[f"cache/data/{did}.npz"]:
        raise SystemExit(f"seal mismatch {did}")


def observed(case):
    """D of a case and the experimental Q factor c_eff_cum (sets the multiplier m; not passed to any solver)."""
    did, var, st = case.split("__"); a, meta = fq6.load(did)
    return (a[VARS[var]].astype(float) + STATES[st] * a["noise"].astype(float)), float(meta["c_eff_cum"])


def compensate(D, cm, Q, c, inv, w_out):
    """All methods of plan section 2.3 for one Q factor c (= m x experimental factor)."""
    U = CW.uqw(D, cm, Q, c, SIG, TW, ANCHOR, GAMMA, inv_snr=inv)
    PH = CW.fd.iq_t(D, cm, Q, None, c, 0, TW); T = shaping.target(U, TW, DT)
    return dict(VDB=CW.lib12.tvw(PH, w_out, EPS, win=WIN), IQ70=CW.fd.iq_t(D, cm, Q, c, c, CW.G_IQ, TW), UQW=U,
                UQW_TI=shaping.ti_shape(U, T, TW, DT), PH_TI=shaping.ti_shape(PH, T, TW, DT),
                HAW=dev6.anchor_shape(U, TW, inv), UQW_TIw=TA.ti_shape(U, w_out, EPS, 1.2, 3.0))


def truth_path_check():
    src = inspect.getsource(compensate) + inspect.getsource(shaping)
    bad = [t for t in FORBIDDEN if t in src]; print("truth-path check (compensate + shaping.py):", "clean" if not bad else bad, flush=True)
    if bad:
        raise SystemExit("truth names in solver code")


def mismatched(R, seed):
    """T2-mis reference reflectivity (descriptive, C2): per trace a constant shift U(-8, 8) ms (1 ms grid) and 30% of an
    independent white series of the same rms: sqrt(1 - 0.09) R(t - s) + 0.3 e."""
    rng = np.random.default_rng(seed); out = np.empty_like(R)
    for j in range(R.shape[1]):
        r = np.roll(R[:, j], int(rng.integers(-8, 9))); e = rng.standard_normal(len(r)) * np.sqrt(np.mean(r ** 2))
        out[:, j] = np.sqrt(1 - 0.09) * r + 0.3 * e
    return out


def score(case, outs, fp):
    """Scores of all outputs of one case (truth loaded here, after compensation)."""
    did = case.split("__")[0]; a, _ = fq6.load(did); ideal, R = a["ideal"].astype(float), a["R"].astype(float)
    Rm = mismatched(R, 20261119); rec = []
    for (meth, m), X in outs.items():
        rec.append(dict(method=meth, m=m, T2=probe.t2(X, R, fp), S=CW.score(X, ideal)["S"], T2mis=probe.t2(X, Rm, fp),
                        latcoh=TA.lateral_coherence(X), hf=probe.hf(X)))
    return rec


def run_case(case, cm, Q, fp):
    D, ce = observed(case); w_out = CW.sc11.prepare_sc(D, cm)["wv"]["w_out"]; inv = CW.snr_rows(D, TW, F1024)
    outs = {("Input", 0.0): D}
    for m in M_ALL:
        for k, X in compensate(D, cm, Q, m * ce, inv, w_out).items():
            outs[(k, m)] = X
    bad = [k for k, X in outs.items() if not np.isfinite(X).all()]
    did, var, st = case.split("__")
    return dict(case=case, did=did, structure=did[1], variant=var, state=st, nonfinite=[list(k) for k in bad],
                rows=score(case, outs, fp))


if __name__ == "__main__":
    cm = semi.curtain_cm(); Q = np.asarray(cm["Q"], float); fp = float(cm["fpeak"]); mode = sys.argv[1]; truth_path_check()
    if mode == "smoke":
        a, meta = probe.load(sys.argv[2]); D = a["D"].astype(float); t0 = time.monotonic()
        w_out = CW.sc11.prepare_sc(D, cm)["wv"]["w_out"]; inv = CW.snr_rows(D, TW, F1024); t1 = time.monotonic()
        o = compensate(D, cm, Q, meta["c_eff_cum"], inv, w_out); t2 = time.monotonic()
        haw_ref = dev6.anchor_shape(CW.uqw(D, cm, Q, meta["c_eff_cum"], SIG, TW, ANCHOR, GAMMA, inv_snr=inv), TW, inv)
        print("smoke", sys.argv[2], "finite", {k: bool(np.isfinite(v).all()) for k, v in o.items()})
        print("HAW vs dev6.anchor_shape max |diff| = %.1e" % np.abs(o["HAW"] - haw_ref).max())
        ref = CW.methods(D, cm, Q, meta["c_eff_cum"], TW, ANCHOR, w_out, inv)
        for k in ("UQW", "VDB", "IQ70"):
            print("%s vs confirmW.methods max |diff| / max = %.1e" % (k, np.abs(o[k] - ref[k]).max() / np.abs(ref[k]).max()))
        print("sc11.run is", getattr(CW.sc11.run, "__file__", "?"))
        print("seconds: prepare %.0f, one multiplier %.0f -> per case about %.0f" % (t1 - t0, t2 - t1, (t1 - t0) + 9 * (t2 - t1)))
    else:
        k, n = int(sys.argv[2]), int(sys.argv[3]); OUT.mkdir(parents=True, exist_ok=True); table = sealed(); checked = set()
        for i, case in enumerate(CASES):
            path = OUT / f"{case}.json"
            if i % n != k or path.exists():
                continue
            did = case.split("__")[0]
            if did not in checked:
                check_seal(did, table); checked.add(did)
            t0 = time.monotonic(); rec = run_case(case, cm, Q, fp)
            path.write_text(json.dumps(rec), encoding="utf-8")
            print(case, "%.0fs" % (time.monotonic() - t0), "nonfinite", len(rec["nonfinite"]), flush=True)
