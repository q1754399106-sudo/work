"""Closed-form prediction of the amplitude-versus-offset bias from Tikhonov absorption compensation (THEORY_ZH.md).
Local stationary approximation at an isolated event: compensation multiplies the event's spectrum by
    net(f, x) = H_true(a_t) conj(H_m) / (|H_m|^2 + damping),
and a matched filter with the true wavelet reads amplitude r(x) Gamma(x), Gamma = sum S^2 Re net / sum S^2."""
import numpy as np
from gate1 import C, FS, RIC, GEO, SIN2, NX, a0_ax

B = -np.pi * FS + 2j * FS * np.log(np.maximum(FS, 1e-300) / C.fref); B[0] = 0
S2 = RIC ** 2


def H(a):
    return np.exp(B * a)


def gamma(method, at0, atx, am0, amx, lam, lx=1e-3, fc=None):
    """Gain of a method for one event and one offset. at0, atx: true a0, ax; am0, amx: model a0, ax."""
    Ht = H(at0 + atx)
    if method == "split":
        Hx, H0 = H(amx), H(am0)
        net = Ht * np.conj(Hx) / (np.abs(Hx) ** 2 + lx) * np.conj(H0) / (np.abs(H0) ** 2 + lam)
    else:
        Hd = H(am0 + amx); w = 1.0 if method == "direct" else np.exp(-2 * np.pi * fc * amx)
        net = Ht * np.conj(Hd) / (np.abs(Hd) ** 2 + lam * w)
    return float(np.sum(S2 * net.real) / np.sum(S2))


def fc_star(am0, lam, spec=None):
    """Centroid of the regularisation transition band at zero offset: uses the model Q and the wavelet spectrum only, no truth."""
    rho = 1 / (1 + lam * np.exp(2 * np.pi * FS * am0)); wt = (S2 if spec is None else spec) * rho * (1 - rho)
    return float(np.sum(wt * FS) / np.sum(wt))


def kappa(j, Qm):
    """First-order coefficient: ax ~ kappa sin^2(theta) at event j (small-angle ray; v_loc as in SIN2)."""
    i0 = int(round(C.event_t[j] / C.dt)); tau = GEO[0]["tau"][i0]; vloc = np.sqrt(SIN2[i0, 1]) / GEO[1]["p"][i0] if GEO[1]["p"][i0] > 0 else C.vp[j]
    return float(0.5 * np.sum(tau * C.vp ** 2 / Qm) / vloc ** 2)


def first_order(j, Qm, lam):
    """Delta(G/A) = -2 pi kappa <f (1 - rho)>_{S^2 rho}, Q correct, one-step ridge, linear in ax."""
    i0 = int(round(C.event_t[j] / C.dt)); a0, _ = a0_ax(Qm); rho = 1 / (1 + lam * np.exp(2 * np.pi * FS * a0[i0])); wt = S2 * rho
    return float(-2 * np.pi * kappa(j, Qm) * np.sum(wt * (1 - rho) * FS) / np.sum(wt))


def predict_ga(method, j, Qt, Qm, lam, r, lx=1e-3, fc=None, keep=None):
    """Predicted G/A of the output: amplitude r(x) Gamma(x) fitted against sin^2 over the kept offsets, and of truth."""
    i0 = int(round(C.event_t[j] / C.dt)); a0t, axt = a0_ax(Qt); a0m, axm = a0_ax(Qm); ok = np.ones(NX, bool) if keep is None else keep[i0]
    g = np.array([gamma(method, a0t[i0], axt[i0, i], a0m[i0], axm[i0, i], lam, lx, fc) for i in range(NX)]); x2 = SIN2[i0]
    pe = np.polyfit(x2[ok], (r * g)[ok], 1); pt = np.polyfit(x2[ok], r[ok], 1)
    return pe[0] / pe[1], pt[0] / pt[1]
