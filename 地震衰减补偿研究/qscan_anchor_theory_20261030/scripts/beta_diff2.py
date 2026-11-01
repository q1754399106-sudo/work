"""Differential test v2 (BETA_ZH.md step 5): add a known attenuation by re-synthesis in the reflectivity domain.

usage: beta_diff2.py    worker; jobs claimed with results15/<job>.claim, resumable.
Base sections on the wide window 1.0-3.2 s:
  S<seed>  semi-synthetic validation: experiment-14 recipe (structure 10227, noise line 10149 x 1 inside 1.2-3.0 s),
           c0 = 1.0, seeds 9995-9997, forward modelled on the wide window (semi.model with output times 1.0-3.2 s);
  F<line>  real inlines 10141, 10193, 10241, 10305.
Pipeline: u_hat = 1 / SR(1.2-3.0 s); reflectivity R = stabilised inverse Q at u_hat (70 dB) then division by |w_in|
inside the coherent band (beta_field_refl.reflectivity); re-synthesis with the generator's forward model at
u_hat + Delta_u, Delta_u = 0, 0.3, 0.6 (no new noise: R carries the original noise). Truth of a re-synthesised
section: u_true + Delta_u (u_true = 1/c0 for S, unknown for F).
Estimators: SR; fixed-anchor QScan default window c_a = 1; deep window (2.4-2.8 s) c_a = 0.7 and 1.4 (lambda 1/2).
-> results15/<base>_d<Delta_u>.json
"""
import os, json, time
from env15 import np, lib, semi, sc11, bm, gen, ROOT, sq
from anchor2 import wavelets_at
from qprofile import scan
import beta_field_refl as bf

OUT = ROOT / "results15"; bm.MAXITER = 20000
BASES = ["S9995", "S9996", "S9997", "F10141", "F10193", "F10241", "F10305"]; DUS = (0.0, 0.3, 0.6)
TW, TWIDE = semi.TW, bf.TWIDE; INNER = (TWIDE >= 1.2 - 1e-9) & (TWIDE <= 3.0 + 1e-9)


def model_wide(R, tstar, S, f, fref, N=8192):
    """semi.model with output times TWIDE."""
    aa = np.clip((120 - f) / 20, 0, 1); aa = 0.5 - 0.5 * np.cos(np.pi * aa)
    b = -np.pi * f + 2j * f * np.log(np.maximum(f, np.finfo(float).tiny) / fref); b[0] = 0
    off = np.round((TWIDE[:, None] - semi.TF[None, :]) / semi.FINE).astype(int) % N; M = np.empty((len(TWIDE), len(semi.TF)))
    for k in range(len(semi.TF)):
        H = S * aa * np.exp(b * tstar[k]); H[-1] = H[-1].real; M[:, k] = np.fft.irfft(H, N)[off[:, k]]
    return M @ R


def base_section(name, ctx):
    cm = ctx["cm"]
    if name.startswith("F"):
        return bf.read_wide(int(name[1:])), None
    seed = int(name[1:]); g = gen.signal(ctx, 1.0, seed); f = np.fft.rfftfreq(8192, semi.FINE)
    x = model_wide(g["R"], semi.tstar_true(cm, 1.0), g["S"], f, float(cm["fref"])) * g["scale"]
    assert np.allclose(x[INNER], g["x"], atol=1e-6 * np.abs(g["x"]).max()), "wide model differs from gen.signal"
    x[INNER] += g["level"] * ctx["noise"][0]; return x, 1.0


def reflectivity(Dw, c, cm, Q, w_in):
    X = bf.inverse_q(Dw, cm, Q, c); D = Dw[INNER]
    f4, snr, _ = lib.coherence_snr(D, TW, 2.0, 3.0); _, _, P = lib.coherence_snr(D, TW, 1.2, 1.6)
    fpk = float(f4[np.argmax(np.sqrt(P) * (f4 >= 5))]); hi = float(min(lib.snr_band(f4, snr, fpk), 80.0))
    N = 2048; f = np.fft.rfftfreq(N, semi.DT); W = np.abs(np.fft.rfft(np.roll(np.pad(w_in, (0, N - len(w_in))), -100), N))
    m = bf.band_mask(f, 6.0, hi); inv = np.where(W > 1e-3 * W.max(), m / np.maximum(W, 1e-12), 0) * W.max()
    R4 = np.fft.irfft(np.fft.rfft(X, N, axis=0) * inv[:, None], N, axis=0)[:len(TWIDE)]
    live = np.abs(Dw).sum(axis=0) > 0
    R1 = np.fft.irfft(np.fft.rfft(R4 * live[None, :], axis=0), len(TWIDE) * 4, axis=0)[:len(semi.TF)] * 4
    return R1, hi


def resynth(R, u, ctx):
    cm = ctx["cm"]; f = np.fft.rfftfreq(8192, semi.FINE); c = 1.0 / u
    S = semi.source_spectrum(ctx["w_in"], f, 0.0, cm, c); return semi.model(R, semi.tstar_true(cm, c), S, f, float(cm["fref"]))


def estimators(D, cm, Q):
    out = dict(c_sr=float(sq.estimate(D, cm, Q)["c"])); pp = sc11.prepare_sc(D, cm)
    for tag, ca, win in (("w1.4_a1.0", 1.0, None), ("w2.6_a0.7", 0.7, (2.4, 2.8, 2.6)), ("w2.6_a1.4", 1.4, (2.4, 2.8, 2.6))):
        acm = dict(cm, Q=Q * ca); wv = lib.wavelets(D, TW, acm) if win is None else wavelets_at(D, TW, acm, *win)
        r = scan(D, cm, Q, pp, lambda cc: wv); out[tag] = dict(v_sel=r["v_sel"], edge=r["edge"], ok=r["neighbourhood_converged"])
    return out


if __name__ == "__main__":
    OUT.mkdir(exist_ok=True); ctx = gen.prepare(10227, [10149, 10273, 10355]); cm = ctx["cm"]; Q = np.asarray(cm["Q"], float)
    for name in BASES:
        for du in DUS:
            job = f"{name}_d{du}"; path = OUT / f"{job}.json"
            if path.exists():
                continue
            try:
                os.close(os.open(OUT / f"{job}.claim", os.O_CREAT | os.O_EXCL | os.O_WRONLY))
            except FileExistsError:
                continue
            t0 = time.monotonic(); Dw, c0 = base_section(name, ctx); live = np.abs(Dw[INNER]).sum(axis=0) > 0
            u_hat = 1.0 / float(sq.estimate(Dw[INNER], cm, Q)["c"]); R, hi = reflectivity(Dw, 1.0 / u_hat, cm, Q, ctx["w_in"])
            D = resynth(R, u_hat + du, ctx) * live[None, :]
            rec = dict(base=name, du=du, u_hat=u_hat, u_true=(1.0 / c0 + du) if c0 else None, band_hi=hi, **estimators(D, cm, Q))
            if du == 0:
                rec["corr_resynth_base"] = float(np.median([np.corrcoef(D[:, j], Dw[INNER][:, j])[0, 1] for j in np.where(live)[0][::10]]))
            rec["seconds"] = time.monotonic() - t0
            path.write_text(json.dumps(rec, indent=1), encoding="utf-8"); (OUT / f"{job}.claim").unlink()
            print(job, f"u_hat {u_hat:.3f} c_sr {rec['c_sr']:.3f}", {k: round(v["v_sel"], 3) for k, v in rec.items() if isinstance(v, dict)},
                  f"{rec['seconds']:.0f}s", flush=True)
