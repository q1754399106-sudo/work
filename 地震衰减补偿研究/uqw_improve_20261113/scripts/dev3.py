"""Step 3, third probe (development data only): how much of the shaping does UQ-W need?

cache: UQ-W outputs of the 24 development sets (variants O and F of probe2.py, 5 Q multipliers), w_out, ideal and R are
       written to CACHE_DIR (outside the repository), so that shaping variants can be compared without recomputing UQ-W.
eval:  shaping variants applied to the cached UQ-W outputs; S and T2 as probe.py -> results/dev3/<tag>.json
         tvw   eps, win, nsm      lib12.tvw (the VDB shaping; UQWS of probe2.py is eps, win of SEL["B3"], nsm 11)
         lam   lambda             filter of tvw raised to the power lambda inside the band (partial shaping)
         poly  degree             only the polynomial trend (in f, band LO-HI) of ln(W_out / A) is corrected in each window;
                                  outside the band the gain is tapered to zero (removes the out-of-band gain of UQ-W)
         lp                       no shaping, only the band taper of "poly" (reference for the out-of-band defect)

usage: dev3.py cache CACHE_DIR ID [ID ...]
       dev3.py eval CACHE_DIR TAG KIND [PARAMS ...]
"""
import sys, json, time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import probe2 as P2  # noqa: E402
from probe2 import np, semi, score, sc11, TW, DT, snr_rows, CW, load, lib, ROOT, t2, lib12, SEL, tvfilt, CURVE, MS  # noqa: E402

LO, HI, NF = 8.0, 72.0, 512


def taper(f):
    """1 inside LO-HI, cosine to 0 over 4 Hz below LO and 12 Hz above HI."""
    g = np.ones_like(f); a = (f < LO); g[a] = 0.5 * (1 + np.cos(np.pi * np.clip((LO - f[a]) / 4.0, 0, 1)))
    b = (f > HI); g[b] = 0.5 * (1 + np.cos(np.pi * np.clip((f[b] - HI) / 12.0, 0, 1))); return g


def shape(X, w_out, kind, par, win=0.5, hop=0.1, nsm=11):
    n = X.shape[0]; L = int(round(win / DT)); step = int(round(hop / DT)); h = np.hanning(L); f = np.fft.rfftfreq(NF, DT)
    live = np.abs(X).sum(axis=0) > 0; Wo = np.abs(np.fft.rfft(w_out, NF)); Wo = Wo / Wo.max(); band = (f >= LO) & (f <= HI); tp = taper(f)
    out = np.zeros_like(X); ws = np.zeros(n); u = (f - 0.5 * (LO + HI)) / (0.5 * (HI - LO))
    for s0 in range(-L // 2, n - L // 2, step):
        i0, i1 = max(s0, 0), min(s0 + L, n); hw = h[i0 - s0:i1 - s0]; Fx = np.fft.rfft(X[i0:i1] * hw[:, None], NF, axis=0)
        A = np.sqrt(lib.smooth(np.mean(np.abs(Fx[:, live]) ** 2, axis=1), nsm)); A = A / A[band].max()
        if kind == "lp":
            G = tp
        elif kind == "lam":
            G0 = Wo * A / (A ** 2 + SEL["B3"]["eps"] * np.max(A ** 2)); G = np.where(band, G0 ** par[0], G0)
        elif kind == "poly":
            d = np.log(np.maximum(Wo[band], 1e-6) / np.maximum(A[band], 1e-6)); wgt = Wo[band]                   # weight: where the target has energy
            c = np.polynomial.legendre.legfit(u[band], d, int(par[0]), w=wgt); c[0] = 0.0                       # level is not corrected
            G = np.exp(np.polynomial.legendre.legval(np.clip(u, -1, 1), c)) * tp
        out[i0:i1] += np.fft.irfft(Fx * G[:, None], NF, axis=0)[:i1 - i0]; ws[i0:i1] += hw
    return out / np.maximum(ws, 1e-3 * ws.max())[:, None]


if __name__ == "__main__" and sys.argv[1] == "cache":
    cdir = Path(sys.argv[2]); cdir.mkdir(parents=True, exist_ok=True); cm = semi.curtain_cm(); Q = np.asarray(cm["Q"], float)
    r = np.asarray(json.loads(CURVE.read_text(encoding="utf-8"))["r"]); F = np.fft.rfftfreq(1024, DT)
    for did in sys.argv[3:]:
        if (cdir / f"{did}.npz").exists():
            continue
        t0 = time.monotonic(); arr, meta = load(did); d = dict(ideal=arr["ideal"].astype(np.float32), R=arr["R"].astype(np.float32))
        for var, D in (("O", arr["D"]), ("F", tvfilt(arr["signal"], r) + arr["noise"])):
            d["w_" + var] = sc11.prepare_sc(D, cm)["wv"]["w_out"]; inv = snr_rows(D, TW, F)
            d["U_" + var] = np.stack([CW.uqw(D, cm, Q, m * meta["c_eff_cum"], CW.SIG, TW, 1.4, CW.GAMMA, inv_snr=inv) for m in MS]).astype(np.float32)
        np.savez(cdir / f"{did}.npz", **d); print(did, "%.0fs" % (time.monotonic() - t0), flush=True)

if __name__ == "__main__" and sys.argv[1] == "eval":
    cdir = Path(sys.argv[2]); tag, kind = sys.argv[3], sys.argv[4]; par = [float(x) for x in sys.argv[5:]]; fp = float(semi.curtain_cm()["fpeak"]); rec = []
    for p in sorted(cdir.glob("*.npz")):
        z = np.load(p); ideal, R_ = z["ideal"].astype(float), z["R"].astype(float)
        for var in "OF":
            for i, m in enumerate(MS):
                U = z["U_" + var][i].astype(float); w = z["w_" + var]
                X = U if kind == "none" else lib12.tvw(U, w, par[0], win=par[1], nsm=int(par[2])) if kind == "tvw" else shape(U, w, kind, par)
                rec.append(dict(id=p.stem, var=var, m=m, S=score(X, ideal)["S"], T2=t2(X, R_, fp)))
    out = ROOT / "results" / "dev3"; out.mkdir(parents=True, exist_ok=True); (out / f"{tag}.json").write_text(json.dumps(rec), encoding="utf-8")
    for var in "OF":
        for k in ("T2", "S"):
            v = np.array([[q[k] for q in rec if q["var"] == var and q["m"] == m] for m in MS])
            print("%-10s %s %-2s " % (tag, var, k) + " ".join("%.3f" % x for x in v.mean(1)) + " | mean %.3f worst %.3f  n %d" % (v.mean(), v.mean(1).min(), v.shape[1]), flush=True)
