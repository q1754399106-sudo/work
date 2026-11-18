"""Post hoc mechanism diagnosis (PLAN_ABLATION_ZH.md; the verdicts R1-R5 of PLAN_ZH.md are unchanged): why does UQ-W lose to
VDB on the Volve well ties? Time-invariant versus time-variant spectral shaping of the archived UQ-W outputs, and two tie
references (Ricker as before, and the section's own zero-phase anchor-window wavelet w_out, which is VDB's shaping target).

No compensation is recomputed except the phase-only Kolsky correction PH (input of VDB), and only on the inlines the well
paths touch. Reads cache/out/IL_<L>.npz of run.py. Writes results/ablation_tie.json and results/ablation_desc.json;
prints no scores (analyze_ablation.py gives the verdict).

usage: tie_ablation.py smoke   one inline, m = 1: VDB rebuilt from PH equals the archived VDB; timing; saves nothing
       tie_ablation.py run
"""
import sys, json, time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import run as R  # noqa: E402
from run import L, np, semi, lib, TW, DT, CW  # noqa: E402
import run13  # noqa: E402
from lib13 import tie  # noqa: E402

M_ALL = R.M_ALL; EPS = CW.SEL["B3"]["eps"]; WIN = CW.SEL["B3"]["win"]
TI_FULL = (1.2, 3.0); TI_DEEP = (2.0, 3.0)          # analysis windows of the two time-invariant filters (fixed before any tie)
WS = semi.win(TW, 2.3, 2.8)


def ti_shape(X, w_out, eps, t0, t1, nfft=1024, nsm=21):
    """Time-invariant version of lib12.tvw: ONE zero-phase filter |W_out| A / (A^2 + eps max A^2) per section, A the
    section-mean amplitude spectrum of the live traces in [t0, t1] (smoothed over nsm bins ~ the 5 Hz of tvw), applied to
    the whole trace. Same target wavelet and eps as VDB."""
    n = X.shape[0]; live = np.abs(X).sum(axis=0) > 0; w = semi.win(TW, t0, t1)
    Wo = np.abs(np.fft.rfft(w_out, nfft)); Wo = Wo / Wo.max()
    A = np.sqrt(lib.smooth(np.mean(np.abs(np.fft.rfft(X[w][:, live], nfft, axis=0)) ** 2, axis=1), nsm))
    G = Wo * A / (A ** 2 + eps * np.max(A ** 2))
    return np.fft.irfft(np.fft.rfft(X, nfft, axis=0) * G[:, None], nfft, axis=0)[:n]


def variants(ln, cm, Q, ms):
    """All outputs of one inline: archived Input/UQW/VDB/IQ70 plus PH, VDB_re, UQW_TI, UQW_TId, UQW_TV, PH_TI."""
    z = {k: v.astype(float) for k, v in np.load(R.CACHE / f"IL_{ln}.npz").items()}; D = z["Input"]
    w_out = lib.wavelets(D, TW, cm)["w_out"]; c = L.c_sr(ln); out = {"Input": D}
    for m in ms:
        U, PH = z[f"UQW_{m}"], CW.fd.iq_t(D, cm, Q, None, m * c, 0, TW)
        out.update({f"UQW_{m}": U, f"VDB_{m}": z[f"VDB_{m}"], f"IQ70_{m}": z[f"IQ70_{m}"], f"PH_{m}": PH,
                    f"VDB_re_{m}": CW.lib12.tvw(PH, w_out, EPS, win=WIN),
                    f"UQW_TI_{m}": ti_shape(U, w_out, EPS, *TI_FULL), f"UQW_TId_{m}": ti_shape(U, w_out, EPS, *TI_DEEP),
                    f"UQW_TV_{m}": CW.lib12.tvw(U, w_out, EPS, win=WIN), f"PH_TI_{m}": ti_shape(PH, w_out, EPS, *TI_FULL)})
    return {k: v.astype(np.float32) for k, v in out.items()}, w_out


def wavelet_ref(rf, w_out):
    """Log reflectivity (1 ms) convolved with w_out (4 ms, 201 samples, centred) resampled to 1 ms; sampled at TW."""
    k = np.arange(-400, 401) * semi.FINE; wf = np.interp(k, np.arange(-100, 101) * DT, w_out)
    return np.interp(TW, semi.TF, np.convolve(rf, wf, mode="same"))


def lateral_coherence(X):
    """Median correlation of adjacent live traces in 2.3-2.8 s (descriptive, no well)."""
    a = X[WS]; live = np.abs(X).sum(axis=0) > 0; c = []
    for j in range(a.shape[1] - 1):
        if live[j] and live[j + 1]:
            u, v = a[:, j] - a[:, j].mean(), a[:, j + 1] - a[:, j + 1].mean(); c.append(u @ v / np.sqrt((u @ u) * (v @ v) + 1e-300))
    return float(np.median(c))


def main(mode):
    logs = run13.well_logs(); G, _ = tie.groups(logs); h, x, y, il, xl = L.cube(); cm = semi.curtain_cm(); Q = np.asarray(cm["Q"], float)
    fp = float(cm["fpeak"]); paths = {n: tie.path(g, x, y, il, xl) for n, g in logs.items()}
    if mode == "smoke":
        ln = int(il[paths[sorted(logs)[0]][1][0, 0]]); t0 = time.monotonic(); o, _ = variants(ln, cm, Q, (1.0,))
        print("inline", ln, "VDB_re vs archived VDB: max |diff| / max |VDB| = %.2e" % (np.abs(o["VDB_re_1.0"] - o["VDB_1.0"]).max() / np.abs(o["VDB_1.0"]).max()),
              "finite", all(np.isfinite(v).all() for v in o.values()), "seconds %.0f" % (time.monotonic() - t0)); return
    sec, wo = {}, {}

    def section(ln):
        if ln not in sec:
            sec[ln], wo[ln] = variants(ln, cm, Q, M_ALL); print("inline", ln, "prepared", flush=True)
        return sec[ln]
    rows = {}
    for n, g in logs.items():
        w, idx = paths[n]; rf = np.zeros(len(semi.TF)); tm = np.round(0.5 * (g["t"][1:] + g["t"][:-1]) - 0.0005, 6)
        i0 = int(round((tm[0] - semi.TF[0]) / L.wq.FINE)); rf[i0:i0 + len(g["rc"])] = g["rc"]
        ln0 = int(il[idx[0, 0]]); keys = list(section(ln0)); ws = semi.win(TW, g["t0"], g["t1"])
        refs = {"ricker": semi.ricker_ref(rf, fp), "wout": wavelet_ref(rf, wo[ln0])}; sc = {r: {} for r in refs}
        for k in keys:
            tr = np.zeros(len(TW))
            for ii, r in zip(np.where(w)[0], idx):
                tr[ii] = np.mean([section(int(il[j]))[k][ii, (xl[j] - L.XL0) // L.XSTEP] for j in r])
            for rn, ref in refs.items():
                sc[rn][k] = lib.tie(tr, ref, TW, ws)[0]
        rows[n] = sc; print(n, "done", flush=True)
    desc = {"lateral_coherence_median": {k: float(np.median([lateral_coherence(s[k]) for s in sec.values()])) for k in next(iter(sec.values()))},
            "lines": sorted(sec), "vdb_rebuild_max_rel": max(float(np.abs(s[f"VDB_re_{m}"] - s[f"VDB_{m}"]).max() / np.abs(s[f"VDB_{m}"]).max()) for s in sec.values() for m in M_ALL)}
    (R.ROOT / "results" / "ablation_tie.json").write_text(json.dumps(dict(groups=G, wells=rows, ti_windows=[TI_FULL, TI_DEEP]), indent=1), encoding="utf-8")
    (R.ROOT / "results" / "ablation_desc.json").write_text(json.dumps(desc, indent=1), encoding="utf-8")
    print("lines", len(sec), "VDB rebuild max rel diff %.2e" % desc["vdb_rebuild_max_rel"])


if __name__ == "__main__":
    main(sys.argv[1])
