"""UQ-A: AVO-consistent UQ operator (design for gate 2, see PLAN_GATE2_ZH.md).
H(f; t, x) = H0(f; t) * R(f; t, x)
  H0 = conj(E G(u a0)) / (E|G(u a0)|^2 + rho)                      zero-offset UQ-W (gamma 0)
  R  = E[G(u a0)] / E[G(u a)],  a = a0 + ax                         ratio of prior-expected responses
Then E_u[H G(u a)] = H0 E_u[G(u a0)] for every x: the prior-expected transfer is offset-invariant at every frequency,
so the expected AVO bias is zero. sig = 0 gives H0 / G(ax) (zero-offset Tikhonov + exact offset correction, the split idea).
Noise safeguard (same principle and gamma as the saturated damping): |R| stops growing above f_n(t, x), the largest f_n with
sum_f |H0 R_sat|^2 <= gamma^2 sum_f |H0|^2; R_sat(f) = |R(min(f, f_n))| exp(i arg R(f))."""
import numpy as np
from gate1 import F, G, moments

GAMMA = 6.0


def ratio(a0, a, sig):
    E0, _ = moments(a0, sig); E1, _ = moments(a, sig)
    ok = np.abs(E1) > 1e-280
    return np.where(ok, E0 / np.where(ok, E1, 1.0), 0.0), E0


def saturate(H0, R, gamma=GAMMA):
    """H0, R: (..., nf). Returns R_sat and the index of f_n (per leading element)."""
    P0 = np.sum(np.abs(H0) ** 2, axis=-1, keepdims=True); A = np.abs(R); ph = np.exp(1j * np.angle(R))
    # cumulative noise gain if saturating at index k: sum_{f<=k} |H0 R|^2 + |R_k|^2 sum_{f>k} |H0|^2
    h2 = np.abs(H0) ** 2; c1 = np.cumsum(h2 * A ** 2, axis=-1); tail = P0 - np.cumsum(h2, axis=-1)
    gain = c1 + A ** 2 * tail; ok = gain <= gamma ** 2 * P0
    # largest k with ok over a monotone prefix: first failure minus one
    bad = ~ok; first_bad = np.where(bad.any(-1), bad.argmax(-1), A.shape[-1] - 1 + 1) - 1
    k = np.maximum(first_bad, 0)
    Ak = np.take_along_axis(A, k[..., None], axis=-1); fidx = np.arange(A.shape[-1])
    Asat = np.where(fidx <= k[..., None], A, Ak)
    return Asat * ph, k


def op_ua(a0, a, sig, rho, sat=True):
    R, E0 = ratio(a0, a, sig); _, E02 = moments(a0, sig); H0 = np.conj(E0) / (E02 + rho)
    if sat:
        R, _ = saturate(H0, R)
    return H0 * R
