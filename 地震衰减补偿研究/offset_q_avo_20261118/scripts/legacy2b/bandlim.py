"""My implementation of the idea in Ratnett, Baldock, Masoomzadeh (2019, SEG): inverse Q without a gain limit up to a time-varying frequency limit, then a high cut there.
The limit is one per time window, the same for every offset (they estimate it on the stack). Their abstract gives no formula, so window sizes, taper and smoothing are my choices.
Limit rules: 'peak' = peak frequency of the gather-mean spectrum (their wording); 'snr' = highest frequency at which the gather-mean spectrum exceeds k times the expected noise
amplitude (their Figure 3 ties the limit to signal-to-noise); the 'snr' rule is mine. Their filter is amplitude-only on NMO-corrected data; here the dispersion phase is also removed,
since the forward model has it. Windowed Fourier overlap-add."""
import numpy as np
from gate1 import C, n, NX, DT, T, a0_ax

L = 50; HOP = 12; NF = 1024; F = np.fft.rfftfreq(NF, DT); BF = -np.pi * F + 2j * F * np.log(np.maximum(F, 1e-300) / C.fref); BF[0] = 0; WIN = np.hanning(L + 1)[:-1]; CT = np.arange(-L // 2, n + L // 2, HOP)
NORM = np.zeros(n)
for c in CT:
    k = np.arange(c - L // 2, c + L // 2); ok = (k >= 0) & (k < n); NORM[k[ok]] += WIN[ok]
NORM = np.maximum(NORM, 1e-6); FIN = 150.0


def limits(D, sig, rule="snr", k_=2.0, floor=15.0):
    Ls = 100; w = np.hanning(Ls); out = np.zeros(len(CT)); b = (F >= 5) & (F <= FIN); fb = F[b]; nmean = float(np.mean(sig)) * np.sqrt(np.sum(w ** 2)) * np.sqrt(np.pi) / 2
    for q, c in enumerate(CT):
        k = np.arange(c - Ls // 2, c + Ls // 2); ok = (k >= 0) & (k < n); X = np.zeros((Ls, D.shape[1])); X[ok] = D[k[ok]]; S = np.mean(np.abs(np.fft.rfft(X * w[:, None], NF, axis=0)), axis=1)[b]
        if rule == "peak":
            out[q] = fb[np.argmax(S)]
        else:
            S = np.convolve(S, np.ones(9) / 9, mode="same"); g = np.where(S >= k_ * nmean)[0]; out[q] = fb[g[-1]] if len(g) else floor
    return np.maximum([np.median(out[max(0, q - 2):q + 3]) for q in range(len(CT))], floor)


def apply(D, Q, sig, rule="snr", k_=2.0, nocut=False):
    """D (n, NX), sig (NX) -> compensated (n, NX)."""
    a0, ax = a0_ax(Q); a = a0[:, None] + ax; lim = limits(D, sig, rule, k_); O = np.zeros_like(D); ii = np.arange(L + 200) - 100
    for q, c in enumerate(CT):
        k = np.arange(c - L // 2, c + L // 2); ok = (k >= 0) & (k < n)
        if not ok.any():
            continue
        cc = int(np.clip(c, 0, n - 1)); X = np.zeros((L, NX)); X[ok] = D[k[ok]]; sp = np.fft.rfft(X * WIN[:, None], NF, axis=0); cut = 1.0 if nocut else np.clip((lim[q] + 10 - F) / 10, 0, 1)[:, None]
        y = np.fft.irfft(sp * np.exp(-np.outer(BF, a[cc])) * cut, NF, axis=0); kk = c - L // 2 + ii; o2 = (kk >= 0) & (kk < n); O[kk[o2]] += y[ii % NF][o2]
    return O / NORM[:, None]
