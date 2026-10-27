"""Experiment 7 runs (PLAN_ZH.md): SC-a and QScan+SC-a on the F-11 T2 and F-15A curtains, scoring, gates,
and the descriptive application to inline 10141. Truth-free; wells are used only for scoring.

usage: run.py wells | run.py inline   (resumable: solved outputs in results/ are reused)
"""
import sys, json, struct, time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import lib  # noqa: E402  (sets up the snapshot solver and NumPy)
from lib import np, bm  # noqa: E402
from scipy.io import loadmat  # noqa: E402

ROOT = lib.ROOT
IN = ROOT / "cache" / "inputs"
E6, E5 = IN / "field_vsp_volve_f11_20261026", IN / "field_vsp_volve_20261025"
WELLS = {
    "F11": dict(curtain=E6 / "inputs/posthoc_curtain.mat", refs=E6 / "analysis/posthoc_references.mat", fhi=65,
                standard=E6 / "results_posthoc/V0.mat", sc=E6 / "results_sc_posthoc/F11_V0_SC.mat",
                qsc=E6 / "results_sc_posthoc/F11_V0_QSC.mat"),
    "F15A": dict(curtain=E5 / "inputs/curtain.mat", refs=E5 / "analysis/references.mat", fhi=70,
                 standard=E5 / "results/V0.mat", sc=E6 / "results_sc/F15A_V0_SC.mat",
                 qsc=E6 / "results_sc/F15A_V0_QSC.mat"),
}
RESOLUTION = 0.064   # |T1(Input, W1) - T1(Input, W3)| on F-11 (experiment 6), PLAN_ZH.md gate A2
OUT = ROOT / "results"; AN = ROOT / "analysis"


def setup(D):
    live = np.abs(D).sum(axis=0) > 0
    K, _, adj = bm.predictor(D, live)
    return live, K, adj


def solve_pair(tag, cm, D, sigma, live, K):
    """SC-a at the model Q and QScan+SC-a; cached in results/<tag>_*.npz and .json."""
    wv = lib.wavelets(D, cm["tw"], cm)
    res = {}
    f = OUT / f"{tag}_SCa.npz"
    if not f.exists():
        X, info, _ = lib.solve(cm, cm["Q"], wv, D, sigma, K, live, *lib.SEL_SC)
        np.savez_compressed(f, out=X); bm.writejson(f.with_suffix(".json"), info)
    res["SCa"] = np.load(f)["out"]
    fq = OUT / f"{tag}_qscan.json"
    if not fq.exists():
        c, qi = lib.qscan(cm, cm["Q"], wv, D, sigma, K, live)
        bm.writejson(fq, dict(factor=c, **qi))
    qs = json.loads(fq.read_text(encoding="utf-8"))
    f = OUT / f"{tag}_QSCa.npz"
    if not f.exists():
        X, info, _ = lib.solve(cm, np.asarray(cm["Q"]) * qs["factor"], wv, D, sigma, K, live, *lib.SEL_SC)
        np.savez_compressed(f, out=X); bm.writejson(f.with_suffix(".json"), info)
    res["QSCa"] = np.load(f)["out"]
    wjson = {k: (v.tolist() if isinstance(v, np.ndarray) else v) for k, v in wv.items()}
    bm.writejson(AN / f"{tag}_wavelets.json", wjson)
    return res, qs, wv


def curtain_cm(c):
    cm = dict(c["cm"]); cm["tw"] = bm.vec(c["tw"]); cm["Q"] = bm.vec(cm["Q"]); cm["solver_offset"] = 0
    return cm


def wells():
    rows, gates = [], {}
    for well, p in WELLS.items():
        t0 = time.monotonic()
        c = loadmat(p["curtain"], simplify_cells=True); r = loadmat(p["refs"], simplify_cells=True)
        D = np.asarray(c["D"], float); sigma = bm.vec(c["sigma"]); tw = bm.vec(c["tw"]); idx = bm.vec(c["idx"], int)
        cm = curtain_cm(c)
        live, K, _ = setup(D)
        new, qs, wv = solve_pair(well, cm, D, sigma, live, K)
        outs = {"Input": D, "Standard": loadmat(p["standard"], simplify_cells=True)["G"]["Standard"],
                "SC": loadmat(p["sc"], simplify_cells=True)["X"], "QSC": loadmat(p["qsc"], simplify_cells=True)["X"], **new}
        w1 = bm.vec(r["w1"], bool); w2 = bm.vec(r["w2"], bool)
        xin = lib.path_trace(D, idx); t1in, l1in = lib.tie(lib.zpfilter(xin, lib.DT, p["fhi"]), r["ref1"], tw, w1)
        t2in, l2in = lib.tie(xin, r["ref2"], tw, w2)
        sc = {}
        for m, A in outs.items():
            A = np.asarray(A, float); x = lib.path_trace(A, idx); xf = lib.zpfilter(x, lib.DT, p["fhi"])
            t1, l1 = lib.tie(xf, r["ref1"], tw, w1); t2, l2 = lib.tie(x, r["ref2"], tw, w2)
            ci1 = lib.tie_boot(xf, lib.zpfilter(xin, lib.DT, p["fhi"]), r["ref1"], tw, w1, l1, l1in) if m != "Input" else [0, 0]
            ci2 = lib.tie_boot(x, xin, r["ref2"], tw, w2, l2, l2in) if m != "Input" else [0, 0]
            hf, b05, b95 = lib.spectral(A, D, tw)
            row = dict(well=well, method=m, T1=t1, T1_shift_ms=l1, T1_minus_input_90CI=ci1, T2=t2, T2_shift_ms=l2,
                       T2_minus_input_90CI=ci2, hf_ratio_60_100Hz_2to3s=hf, band05_Hz=b05, band95_Hz=b95)
            if "ref1w3" in r:
                row["T1w3"] = lib.tie(xf, r["ref1w3"], tw, w1)[0]
            rows.append(row); sc[m] = row
        g = {}
        for m in ("SCa", "QSCa"):
            g[m] = dict(A1_noise=sc[m]["hf_ratio_60_100Hz_2to3s"] <= 10,
                        A2_tie=sc[m]["T1"] >= sc["Input"]["T1"] - RESOLUTION,
                        A3_broader=sc[m]["band95_Hz"] >= sc["Input"]["band95_Hz"] + 5)
            g[m]["all"] = all(g[m].values())
        g["qscan"] = qs; g["f_peak"] = wv["f_peak"]; g["fmax_shallow"] = wv["fmax_shallow"]; g["fmax_target"] = wv["fmax_target"]
        g["Q_model"] = cm["Q"].tolist(); g["seconds"] = time.monotonic() - t0
        gates[well] = g
        print(well, json.dumps({m: g[m] for m in ("SCa", "QSCa")}), f"qscan {qs['factor']:.3f} boundary {qs['at_boundary']}",
              f"fmax {wv['fmax_shallow']:.1f}/{wv['fmax_target']:.1f} Hz", flush=True)
    for m in ("SCa", "QSCa"):
        gates[f"verdict_{m}"] = "usable on both wells" if all(gates[w][m]["all"] for w in WELLS) else "not shown"
    bm.writejson(AN / "scores.json", rows); bm.writejson(AN / "gates.json", gates)
    print(json.dumps({k: v for k, v in gates.items() if k.startswith("verdict")}), flush=True)


def ibm(u):
    u = u.astype(np.uint64); s = 1 - 2 * ((u >> 31) & 1).astype(float)
    e = ((u >> 24) & 127).astype(float) - 64; f = (u & 0xFFFFFF).astype(float) / 16777216
    return s * f * 16.0 ** e


def haar_sigma(x):
    d = (x[1::2][: len(x) // 2] - x[0::2][: len(x) // 2]) / np.sqrt(2)
    return float(np.median(np.abs(d)) / 0.6745)


def inline():
    """Inline 10141, xlines 2100-2500 (25 m), 1.2-3.0 s; Q model of F-11 (post-hoc), wavelets and QScan from the line."""
    snap = json.loads((ROOT / "cache" / "SNAPSHOT.json").read_text(encoding="utf-8"))
    d = Path(snap["cube"]).read_bytes(); ns = 1126; tb = 240 + 4 * ns; n = (len(d) - 3600) // tb
    hdr = np.frombuffer(d[3600:], dtype=np.uint8)[: n * tb].reshape(n, tb)
    il = hdr[:, 188:192].copy().view(">i4").ravel(); xl = hdr[:, 192:196].copy().view(">i4").ravel()
    sel = np.where((il == 10141) & (xl >= 2100) & (xl <= 2500))[0]; sel = sel[np.argsort(xl[sel])]
    X = np.stack([ibm(hdr[j, 240:].copy().view(">u4")) for j in sel], axis=1)
    t = np.arange(ns) * lib.DT; w = (t >= 1.2 - 1e-9) & (t <= 3.0 + 1e-9); D = X[w]; tw = t[w]
    c = loadmat(WELLS["F11"]["curtain"], simplify_cells=True); cm = curtain_cm(c); cm["tw"] = tw; cm["t"] = tw
    sigma = np.array([haar_sigma(D[:, j]) for j in range(D.shape[1])])
    live, K, _ = setup(D)
    new, qs, wv = solve_pair("IL10141", cm, D, sigma, live, K)
    hf = {m: lib.spectral(A, D, tw) for m, A in new.items()}
    np.savez_compressed(OUT / "IL10141_input.npz", D=D, tw=tw, xline=xl[sel])
    bm.writejson(AN / "inline.json", dict(traces=int(D.shape[1]), xline=[int(xl[sel][0]), int(xl[sel][-1])],
                                          qscan=qs, f_peak=wv["f_peak"], fmax_shallow=wv["fmax_shallow"],
                                          fmax_target=wv["fmax_target"], spectral={m: v for m, v in hf.items()}))
    print("inline", D.shape, f"qscan {qs['factor']:.3f}", {m: [round(v, 2) for v in hf[m]] for m in hf}, flush=True)


if __name__ == "__main__":
    OUT.mkdir(exist_ok=True); AN.mkdir(exist_ok=True)
    {"wells": wells, "inline": inline}[sys.argv[1]]()
