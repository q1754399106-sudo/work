"""SID-UQ step 0 (PLAN_ZH.md section 3): physical-premise diagnostics. No compensation, no well-tie scoring.

Z2a  per-line quadratic coefficient of the field log spectral ratio (2.3-2.8 s over 1.2-1.6 s), 10-55 Hz.
Z2b  constrained fit of the mean ratio, 10-72 Hz:  a - pi f D - eta f^2/(f^2+fc^2),  D, eta >= 0, fc grid 2-200 Hz.
Z2c  exact 1-D normal-incidence transmission (all internal multiples, flux-normalised, 1 ms TWT layers) through the
     well-log reflectivity of F-12 and F-14 between 1.4 and 2.55 s; two-way ln|T| differences vs the field excess drop.
D3   reflectivity colour of the 2.3-2.8 s log windows vs the field residual.
D4   adjacent-trace coherent SNR per frequency in the deep and shallow windows.

The Volve cache is read (read-only) from the main checkout, as the worktree has no cache.
usage: step0.py
"""
import sys, json, subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]; REPO = ROOT.parents[1]
MAIN = Path(r"D:/ChatGPT/研究/地震衰减补偿研究/field_baselines_volve_20261029")
sys.path.insert(0, str(MAIN / "scripts"))
import lib13 as L  # noqa: E402
from lib13 import np, lib, semi, wq, tie, TW, DT  # noqa: E402

NF = 512; F = np.fft.rfftfreq(NF, DT); FIT = (F >= 10) & (F <= 55); FULL = (F >= 10) & (F <= 72); HF = (F >= 36) & (F <= 72)
FS = np.arange(8, 73, 4); T_SH, T_DP = (1.2, 1.6), (2.3, 2.8); PATH = (1.4, 2.55)


def lines():
    txt = subprocess.run(["git", "show", "HEAD:地震衰减补偿研究/field_baselines_volve_20261029/analysis/lines.json"],
                         cwd=REPO, capture_output=True, check=True).stdout.decode("utf-8")
    return json.loads(txt)["lines"]


def spec(D, a, b):
    """Complex spectra (NF/2+1, live traces) of the hann-tapered window a-b s."""
    k = (TW >= a) & (TW <= b); live = np.abs(D).sum(0) > 0
    return np.fft.rfft(D[k][:, live] * np.hanning(k.sum())[:, None], NF, axis=0)


def amp(X):
    return np.sqrt(lib.smooth(np.mean(np.abs(X) ** 2, axis=1), 5))


def coh_snr(X):
    """Adjacent-trace magnitude-squared coherence pooled along the line, SNR = C/(1-C)."""
    a, b = X[:, :-1], X[:, 1:]
    num = np.abs(lib.smooth(np.sum(a * np.conj(b), axis=1).real, 5) + 1j * lib.smooth(np.sum(a * np.conj(b), axis=1).imag, 5)) ** 2
    den = lib.smooth(np.sum(np.abs(a) ** 2, axis=1), 5) * lib.smooth(np.sum(np.abs(b) ** 2, axis=1), 5)
    C = np.clip(num / den, 0, 0.9999); return C / (1 - C)


def at(v, f):
    return float(v[int(np.argmin(np.abs(F - f)))])


def fit_sid(m):
    """min ||m - (a - pi f D - eta phi)|| over D, eta >= 0, fc on a grid (a free); returns best dict."""
    f = F[FULL]; y = m[FULL]; best = None
    for fc in np.exp(np.linspace(np.log(2), np.log(200), 400)):
        phi = f ** 2 / (f ** 2 + fc ** 2); A = np.column_stack([np.ones_like(f), -np.pi * f, -phi])
        for act in ((1, 1), (1, 0), (0, 1), (0, 0)):         # active-set NNLS on (D, eta), a unconstrained
            cols = [0] + [i + 1 for i, s in enumerate(act) if s]; c, *_ = np.linalg.lstsq(A[:, cols], y, rcond=None)
            if np.any(c[1:] < 0):
                continue
            full = np.zeros(3); full[cols] = c; r = y - A @ full; rss = float(r @ r)
            if best is None or rss < best["rss"]:
                best = dict(fc=float(fc), a=float(full[0]), D=float(full[1]), eta=float(full[2]), rss=rss, rms=float(np.sqrt(rss / len(y))))
    p1 = np.polyfit(f, y, 1); r1 = y - np.polyval(p1, f); best["rms_line"] = float(np.sqrt(np.mean(r1 ** 2))); best["D_line"] = float(-p1[0] / np.pi)
    return best


def transmission(rc, freqs):
    """|T(f)| of a Goupillaud stack (one-way layer time = 0.5 ms for a 1 ms TWT grid), flux-normalised interfaces."""
    tau = 0.0005; out = []
    for fr in freqs:
        z = np.exp(1j * 2 * np.pi * fr * tau); M = np.eye(2, dtype=complex)
        for r in rc:
            M = M @ (np.array([[1, r], [r, 1]]) / np.sqrt(1 - r * r)) @ np.diag([z, 1 / z])
        out.append(1 / abs(M[0, 0]))
    return np.array(out)


if __name__ == "__main__":
    cm = semi.curtain_cm(); x = lib.absorption(cm, list(PATH)); Da = float(x[1] - x[0]); res = dict(Da_curtain=Da)
    print("curtain Delta a (1.4 -> 2.55 s) %.4f s" % Da, flush=True)

    # ---------- field spectra ----------
    LN = lines(); R = []; C2 = []; snr_d = []; snr_s = []
    for i, ln in enumerate(LN):
        D = L.read_inline(ln); Xs, Xd = spec(D, *T_SH), spec(D, *T_DP); r = np.log(amp(Xd) / amp(Xs)); R.append(r)
        C2.append(float(np.polyfit(F[FIT], r[FIT], 2)[0])); snr_d.append(coh_snr(Xd)); snr_s.append(coh_snr(Xs))
        if i % 15 == 0:
            print("  line %d (%d/%d) c2 %+.2e" % (ln, i + 1, len(LN), C2[-1]), flush=True)
    R = np.array(R); C2 = np.array(C2); m = R.mean(0); snr_d = np.median(snr_d, 0); snr_s = np.median(snr_s, 0)
    frac_pos = float(np.mean(C2 > 0))
    print("\nZ2a: lines %d, quadratic c2 (10-55 Hz) median %+.2e, IQR [%+.2e, %+.2e], fraction c2>0 %.3f -> %s" % (
        len(C2), np.median(C2), *np.percentile(C2, [25, 75]), frac_pos, "REJECT (eta_hat<0)" if frac_pos >= 0.8 else "not rejected"))
    c2m = float(np.polyfit(F[FIT], m[FIT], 2)[0]); print("     mean-ratio c2 %+.2e (curvature over 10-55 Hz %+.2f)" % (c2m, c2m * (45 / 2) ** 2))
    res["Z2a"] = dict(n=len(C2), c2_median=float(np.median(C2)), c2_iqr=[float(v) for v in np.percentile(C2, [25, 75])], frac_pos=frac_pos,
                      c2_mean_ratio=c2m, reject=frac_pos >= 0.8)

    b = fit_sid(m); u_hat = b["D"] / Da
    print("Z2b: SID fit 10-72 Hz: D %.4f s (u_hat %.2f), eta %.2f, fc %.1f Hz, rms %.3f (straight line %.3f, D_line %.4f) -> %s" % (
        b["D"], u_hat, b["eta"], b["fc"], b["rms"], b["rms_line"], b["D_line"], "REJECT (u_hat<0.5)" if u_hat < 0.5 else "not rejected"))
    b["u_hat"] = u_hat; b["reject"] = bool(u_hat < 0.5); res["Z2b"] = b

    p_hf = np.polyfit(F[HF], m[HF], 1); D_hf = max(0.0, -p_hf[0] / np.pi)
    E_field = (at(m, 10) - at(m, 35)) - np.pi * 25 * D_hf
    print("     high-frequency (36-72 Hz) slope D_HF %.4f s (u %.2f); field drop 10->35 Hz %.2f, excess E_field %.2f" % (
        D_hf, D_hf / Da, at(m, 10) - at(m, 35), E_field))
    print("     ln ratio (re 10 Hz) at " + " ".join("%d" % f for f in FS)); print("     " + " ".join("%+.2f" % (at(m, f) - at(m, 10)) for f in FS))
    res["field"] = dict(D_HF=D_hf, E_field=E_field, drop_10_35=at(m, 10) - at(m, 35), ln_ratio_re10={int(f): at(m, f) - at(m, 10) for f in FS})

    # ---------- wells ----------
    P = wq.picks(); fr = np.arange(2, 81, 1.0); res["Z2c"] = {}; ok_wells = 0; reject_c = True
    for n in ("F-12", "F-14"):
        W = wq.well(n, P); tm = 0.5 * (W["t"][1:] + W["t"][:-1]); k = (tm >= PATH[0]) & (tm <= PATH[1]); rc = W["rc"][k]
        T = transmission(rc, fr); lt = 2 * np.log(T)     # two-way
        g = lambda f: float(lt[int(np.argmin(np.abs(fr - f)))])  # noqa: E731
        E_well = g(10) - g(35); hi = g(35) - g(72)
        cov = (float(tm[k].min()), float(tm[k].max()))
        print("Z2c: %s path %.3f-%.3f s (%d layers, sum R^2 %.3f): two-way ln|T| drop 10->35 Hz %.3f, 35->72 Hz %.3f, at 10/35/72 Hz %.3f/%.3f/%.3f;"
              " E_well/E_field %.3f" % (n, *cov, len(rc), float((rc ** 2).sum()), E_well, hi, g(10), g(35), g(72), E_well / E_field))
        res["Z2c"][n] = dict(coverage=cov, layers=int(len(rc)), sumR2=float((rc ** 2).sum()), E_well=E_well, drop_35_72=hi,
                             ln_T2={int(f): g(f) for f in (5, 10, 20, 35, 50, 72)}, ratio=E_well / E_field)
        reject_c &= E_well < 0.25 * E_field
    res["Z2c"]["reject"] = bool(reject_c); print("     Z2c ->", "REJECT (log-predicted scattering < 0.25 of field excess)" if reject_c else "not rejected")

    # ---------- D3: reflectivity colour in the deep window ----------
    Rf = []
    for n in tie.TIE_WELLS:
        try:
            W = wq.well(n, P)
        except Exception as e:  # noqa: BLE001
            print("  D3 skip", n, e); continue
        tm = 0.5 * (W["t"][1:] + W["t"][:-1]); k = (tm >= T_DP[0]) & (tm <= T_DP[1])
        if k.sum() < 300:
            continue
        rc = W["rc"][k]; ff = np.fft.rfftfreq(1024, 0.001); A = np.abs(np.fft.rfft(rc * np.hanning(len(rc)), 1024))
        Rf.append(np.log(np.sqrt(lib.smooth(np.interp(F, ff, A ** 2), 5))))
    if Rf:
        rf = np.mean(Rf, 0); y1 = m[FULL] - np.polyval(np.polyfit(F[FULL], m[FULL], 1), F[FULL]); y2 = rf[FULL] - np.polyval(np.polyfit(F[FULL], rf[FULL], 1), F[FULL])
        cc = float(np.corrcoef(y1, y2)[0, 1])
        print("D3: %d wells with >=0.3 s in 2.3-2.8 s; residual correlation (10-72 Hz) field ratio vs deep-window log reflectivity colour %.2f" % (len(Rf), cc))
        print("     log-RC ln amp minus line at " + " ".join("%d" % f for f in FS)); print("     " + " ".join("%+.2f" % (at(rf, f) - np.polyval(np.polyfit(F[FULL], rf[FULL], 1), f)) for f in FS))
        res["D3"] = dict(n_wells=len(Rf), resid_corr=cc)

    # ---------- D4 ----------
    print("D4: median adjacent-trace coherent SNR, deep / shallow window, at " + " ".join("%d" % f for f in FS))
    print("     deep    " + " ".join("%.1f" % at(snr_d, f) for f in FS)); print("     shallow " + " ".join("%.1f" % at(snr_s, f) for f in FS))
    res["D4"] = dict(deep={int(f): at(snr_d, f) for f in FS}, shallow={int(f): at(snr_s, f) for f in FS})

    z1 = True; stop = z1 and (res["Z2a"]["reject"] or res["Z2b"]["reject"] or res["Z2c"]["reject"])
    res["stop"] = bool(stop); print("\nSTOP RULE (Z1 and any of Z2a/b/c):", stop)
    (ROOT / "results" / "step0.json").write_text(json.dumps(res, indent=1), encoding="utf-8")
