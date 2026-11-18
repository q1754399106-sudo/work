"""Gate 1 of direction 1 (PLAN_GATE1_ZH.md): does per-trace UQ compensation create a fake AVO gradient?

Aligned (ideal NMO) gather of core.C (offsets 0-1500 m, events 0.38/0.65/0.90 s), Q model M12-M15, surface anchor.
Path t* of event j at offset x: a(x) = a0 + ax (core.ray). Truth: a_t = u a_m (u multiplies every 1/Q).
Operator per row (aligned time t, offset x), prior u' = exp(sig e), e ~ N(0,1), lognormal (UQ-W law, 41-pt Gauss-Hermite):
    H = conj(E[G]) / (E|G|^2 + rho),   G(a) = exp(b a),  b = -pi f + 2i f ln(f/fref)       (UQ-W, gamma = 0)
sig = 0 gives one-step Tikhonov G*/(|G|^2 + rho). rho = lam is the noise term (N/S); constant over the gather.
Closed-form gain (local stationarity):  Gamma_j(x) = int S^2 Re[G(u a_m) H] df / int S^2 df.
Numeric check: time-varying operator rows applied to isolated unit events, read by a zero-phase matched filter.
dGA = B'/A' - B/A, where (A', B') fit co_j(x) Gamma_j(x) on sin^2 theta_j(x) and co_j are core.reference Shuey coefficients.
OC variant (descriptive): per row time, per frequency, equalise the expected transfer T = |E G|^2/(E|G|^2 + rho) over offsets
to its minimum T*: H = T* conj(E G) / |E G|^2.

usage: gate1.py smoke | run
"""
import sys, json, time
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from core import C, ray, reference  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
N = C.nfft; F = np.fft.rfftfreq(N, C.dt); NT = C.nt
B = -np.pi * F + 2j * F * np.log(np.maximum(F, 1e-300) / C.fref); B[0] = 0
GH_E, GH_W = np.polynomial.hermite_e.hermegauss(41); GH_W = GH_W / GH_W.sum()
MODELS = {"M12": [55, 40, 100, 70], "M13": [28, 65, 40, 90], "M14": [70, 28, 35, 50], "M15": [40, 120, 60, 30],
          "M1": [40, 50, 70, 100]}
GATE_MODELS = ["M12", "M13", "M14", "M15"]
FPS = [60.0, 40.0]; LAMS = [1e-4, 1e-3, 1e-2]; SIG = 0.447
U_GRID = [0.5, 0.7, 0.8, 1.0, 1.2, 1.4, 2.0]
XS = C.offset; NX = len(XS)


def spec(fp):
    return F ** 2 / fp ** 2 * np.exp(-F ** 2 / fp ** 2)


def G(a):
    """a: (...,) -> (..., nf)"""
    return np.exp(np.multiply.outer(np.asarray(a, float), B))


def moments(a_m, sig):
    """E[G], E|G|^2 over u' = exp(sig e); a_m (...,) -> (..., nf) each."""
    if sig == 0:
        g = G(a_m); return g, np.abs(g) ** 2
    us = np.exp(sig * GH_E); EG = 0; EG2 = 0
    for u, w in zip(us, GH_W):
        g = G(np.asarray(a_m) * u); EG = EG + w * g; EG2 = EG2 + w * np.abs(g) ** 2
    return EG, EG2


def op_uq(a_m, sig, lam):
    EG, EG2 = moments(a_m, sig); return np.conj(EG) / (EG2 + lam)


def _oc(EG, EG2, lam):
    """moments with offsets on axis 0: equalise T = |EG|^2/(E|G|^2 + lam) over axis 0 to its minimum."""
    m2 = np.abs(EG) ** 2; Ts = (m2 / (EG2 + lam)).min(axis=0, keepdims=True); ok = m2 > 1e-280
    return np.where(ok, Ts * np.conj(EG) / np.where(ok, m2, 1.0), 0.0)


def op_oc(a_m_x, sig, lam):
    """a_m_x: (NX, ...) path t* over offsets for the same aligned time(s); returns (NX, ..., nf)."""
    EG, EG2 = moments(a_m_x, sig); return _oc(EG, EG2, lam)


def paths(Q):
    """Event path t* a[j, x] and sin^2, Shuey coefficients co[j, x]."""
    a = np.zeros((3, NX)); s2 = np.zeros((3, NX)); co = np.zeros((3, NX))
    for ix, x in enumerate(XS):
        a0, ax, _ = ray(C.event_t, x, np.asarray(Q, float)); a[:, ix] = a0 + ax
        _, co[:, ix], s2[:, ix] = reference(x)
    return a, s2, co


def row_paths(Q):
    """a_m(t_i, x) for every sample time and offset: (NX, NT)."""
    out = np.zeros((NX, NT))
    for ix, x in enumerate(XS):
        a0, ax, _ = ray(C.t, x, np.asarray(Q, float)); out[ix] = a0 + ax
    out[:, 0] = out[:, 1]
    return out


def dga(gam, s2, co):
    """gam, s2, co: (NX,) for one event -> B'/A' - B/A."""
    X = np.c_[np.ones(NX), s2]
    A, Bc = np.linalg.lstsq(X, co, rcond=None)[0]; A2, B2 = np.linalg.lstsq(X, co * gam, rcond=None)[0]
    return B2 / A2 - Bc / A
def lng_slope(gam, s2):
    return np.polyfit(s2, np.log(np.maximum(gam, 1e-300)), 1)[0]


def gamma_closed(a_ev, u, H, S2):
    """a_ev (NX,), H (NX, nf) -> Gamma (NX,)"""
    net = G(u * a_ev) * H
    return (S2[None, :] * net.real).sum(1) / S2.sum()


def build_M(Hrows_x):
    """Hrows_x (NX, NT, nf) -> (NX, NT, NT) matrices; output_i = sum_k h_i(i - k) d_k (uqw.uqw convention)."""
    lag = np.arange(NT); idx = (np.arange(NT)[:, None] - lag[None, :]) % N; Ms = np.empty((NX, NT, NT))
    for ix in range(NX):
        h = Hrows_x[ix].copy(); h[:, 0] = h[:, 0].real; h[:, -1] = h[:, -1].real
        Ms[ix] = np.take_along_axis(np.fft.irfft(h, N, axis=1), idx, axis=1)
    return Ms


def numeric_gamma(Ms, a_ev, u, S):
    """Ms (NX, NT, NT), a_ev (3, NX) -> Gamma (3, NX) by matched filter on isolated unit events."""
    out = np.zeros((3, NX)); ph = np.exp(-2j * np.pi * np.outer(C.event_t, F))
    w = np.fft.irfft(S[None, :] * ph, N, axis=1)[:, :NT]
    for ix in range(NX):
        d = np.zeros(NT)
        for j in range(3):
            d += np.fft.irfft(S * G(u * a_ev[j, ix]) * ph[j], N)[:NT]
        y = Ms[ix] @ d
        for j in range(3):
            win = np.abs(C.t - C.event_t[j]) <= 0.12  # isolate event j (events are >= 0.25 s apart)
            out[j, ix] = (y[win] * w[j, win]).sum() / (w[j, win] ** 2).sum()
    return out


def event_ops(a_ev, kind, sig, lam):
    """Operators at the event points: (3, NX, nf)."""
    if kind == "uq":
        return op_uq(a_ev, sig, lam)
    return np.stack([op_oc(a_ev[j][:, None], sig, lam)[:, 0] for j in range(3)])


def smoke():
    a_ev, s2, co = paths(MODELS["M1"])
    lam = 1e-3; H0 = op_uq(a_ev, 0.0, lam); g = G(a_ev)
    print("sig=0 vs Tikhonov rel %.1e" % (np.abs(H0 - np.conj(g) / (np.abs(g) ** 2 + lam)).max() / np.abs(H0).max()))
    Hs = op_uq(a_ev, 1e-6, lam)
    print("sig=1e-6 vs sig=0 rel %.1e" % (np.abs(Hs - H0).max() / np.abs(H0).max()))
    EG, EG2 = moments(np.array([0.01]), 0.3); print("E|G|^2 >= |EG|^2:", bool((EG2 + 1e-300 >= np.abs(EG) ** 2 * (1 - 1e-12)).all()))
    # dga of a constant gain is zero
    print("dga(const) %.1e" % abs(dga(np.full(NX, 0.7), s2[0], co[0])))
    # numeric pipeline with identity operator returns ~1
    S = spec(60.0); Hid = np.ones((NX, NT, len(F)), complex)
    gn = numeric_gamma(build_M(Hid), a_ev, 0.0, S); print("identity numeric gamma range %.4f-%.4f" % (gn.min(), gn.max()))
    # OC with constant path over offsets equals UQ
    a1 = np.repeat(a_ev[1:2, :1], NX, axis=1)
    d = np.abs(op_oc(a1[0][:, None], 0.447, lam)[:, 0] - op_uq(a1[0], 0.447, lam)).max()
    print("OC == UQ for equal paths: %.1e" % d)


def run():
    t0 = time.monotonic(); res = []
    for mname, Q in MODELS.items():
        a_ev, s2, co = paths(Q); am = row_paths(Q)
        for sig in (0.0, SIG, 0.25):
            EGr, EG2r = moments(am, sig); EGe, EG2e = moments(a_ev, sig)
            for lam in LAMS:
                for kind in ("uq", "oc"):
                    if kind == "oc" and sig == 0.25:
                        continue
                    if kind == "uq":
                        Hrows = np.conj(EGr) / (EG2r + lam); Hev = np.conj(EGe) / (EG2e + lam)
                    else:
                        Hrows = _oc(EGr, EG2r, lam); Hev = np.swapaxes(_oc(np.swapaxes(EGe, 0, 1), np.swapaxes(EG2e, 0, 1), lam), 0, 1)
                    need = [u for u in U_GRID if u == 1.0 or (mname in GATE_MODELS and sig == SIG)]
                    Ms = build_M(Hrows) if need else None; del Hrows
                    for fp in FPS:
                        S = spec(fp); S2 = S ** 2
                        for u in U_GRID:
                            gc = np.stack([gamma_closed(a_ev[j], u, Hev[j], S2) for j in range(3)])
                            gn = numeric_gamma(Ms, a_ev, u, S) if u in need else None
                            for j in range(3):
                                pf = np.polyfit(s2[j], co[j], 1)
                                r = dict(model=mname, fp=fp, lam=lam, kind=kind, sig=sig, u=u, ev=j,
                                         gam_closed=gc[j].tolist(), dga_closed=dga(gc[j], s2[j], co[j]),
                                         slope_closed=lng_slope(gc[j], s2[j]), BA=float(pf[0] / pf[1]))
                                if gn is not None:
                                    r.update(gam_num=gn[j].tolist(), dga_num=dga(gn[j], s2[j], co[j]))
                                res.append(r)
                    print("%s sig %g lam %g %s done %.0fs" % (mname, sig, lam, kind, time.monotonic() - t0), flush=True)
            del EGr, EG2r
    (ROOT / "results").mkdir(exist_ok=True)
    (ROOT / "results" / "gate1.json").write_text(json.dumps(res))
    print("wrote %d rows %.0fs" % (len(res), time.monotonic() - t0))


if __name__ == "__main__":
    {"smoke": smoke, "run": run}[sys.argv[1]]()
