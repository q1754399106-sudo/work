"""Exploratory pilot (no frozen rules): can the Q scale u = 1/c be read from dispersion alone?

Phase-only Kolsky correction exp(-i theta(f) [a(t) u + g]) with theta(f) = 2 f ln(f/fref) (the imaginary part of
b(f) in semi.model), a(t) = model t* at c = 1 on the 1.2-3.0 s observation; g absorbs the time-invariant dispersion
(the unidentifiable absolute level, like the wavelet anchor), psi a constant phase rotation (source phase).
Criterion: kurtosis of the corrected, band-limited (8-60 Hz), envelope-normalised data, maximised over (psi, g) for
each u. The amplitude spectrum is never touched, so reflectivity colour cannot enter.
Data: experiment-14 recipe through gen (structure 10227, noise line 10149 x 1), c 0.7 / 1.0 / 1.4, seeds 10001-10004
(already used; exploratory), plus noise-free versions. -> ../results/pilot_kurt.json
"""
import sys, json, time
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[2] / "qscan_anchor_theory_20261030" / "scripts"))
from env15 import np, lib, semi, gen  # noqa: E402

DT = semi.DT; NF = 1024
U = np.round(np.arange(0.0, 2.01, 0.125), 4); G = np.linspace(-0.02, 0.02, 9); PSI = np.deg2rad(np.arange(0, 180, 15))


def corrector(cm):
    f = np.fft.rfftfreq(NF, DT); fref = float(cm["fref"])
    th = 2 * f * np.log(np.maximum(f, 1e-9) / fref); th[0] = 0
    band = np.clip((f - 6) / 4, 0, 1) * np.clip((62 - f) / 4, 0, 1)            # 8-60 Hz cosine-free trapezoid
    a = lib.absorption(cm, semi.TW)
    E = np.exp(2j * np.pi * f[None, :] * (np.arange(len(semi.TW)) * DT)[:, None])  # (nt, nf) synthesis kernel
    return f, th, band, a, E


def kurt_curve(D, cm):
    f, th, band, a, E = corrector(cm); nt = D.shape[0]
    X = np.fft.rfft(D, NF, axis=0) * band[:, None]                               # (nf, nx)
    out = np.zeros(len(U))
    for i, u in enumerate(U):
        best = -np.inf
        for g in G:
            P = np.exp(-1j * th[None, :] * (a[:, None] * u + g))                    # (nt, nf)
            Y = (E * P) @ X                                                          # analytic-ish signal, (nt, nx)
            Y = 2 * Y / NF
            env = np.sqrt(lib.smooth(np.mean(np.abs(Y) ** 2, axis=1), 51))[:, None]  # 0.2 s lateral-mean envelope
            for p in PSI:
                y = np.real(Y * np.exp(-1j * p)) / env
                k = float(np.mean(y ** 4) / np.mean(y ** 2) ** 2)
                best = max(best, k)
        out[i] = best
    return out


def pick(curve):
    j = int(np.argmax(curve))
    if 0 < j < len(U) - 1:
        c = np.polyfit(U[j - 1:j + 2], curve[j - 1:j + 2], 2)
        if c[0] < 0:
            return float(np.clip(-c[1] / (2 * c[0]), U[j - 1], U[j + 1])), j
    return float(U[j]), j


if __name__ == "__main__":
    t0 = time.monotonic(); ctx = gen.prepare(10227, [10149, 10273, 10355]); cm = ctx["cm"]; rows = []
    for seed in (10001, 10002, 10003, 10004):
        for c in (0.7, 1.0, 1.4):
            g = gen.signal(ctx, c, seed)
            for kind, D in (("noisy", g["x"] + g["level"] * ctx["noise"][0]), ("clean", g["x"])):
                cur = kurt_curve(D, cm); u_sel, j = pick(cur)
                rows.append(dict(seed=seed, c=c, u_true=1 / c, kind=kind, u_sel=u_sel, edge=j in (0, len(U) - 1), curve=cur.tolist()))
                print(seed, c, kind, "u_true %.3f u_sel %.3f" % (1 / c, u_sel), "edge" if j in (0, len(U) - 1) else "",
                      f"{time.monotonic() - t0:.0f}s", flush=True)
    out = HERE.parent / "results"; out.mkdir(exist_ok=True)
    for kind in ("noisy", "clean"):
        r = [x for x in rows if x["kind"] == kind]; ut = np.array([x["u_true"] for x in r]); us = np.array([x["u_sel"] for x in r])
        A = np.polyfit(ut, us, 1); print(kind, "slope %.3f intercept %.3f |log err| median %.3f" % (A[0], A[1], np.median(np.abs(np.log(us / ut)))))
    (out / "pilot_kurt.json").write_text(json.dumps(dict(U=U.tolist(), rows=rows), indent=1), encoding="utf-8")
