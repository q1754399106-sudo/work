"""
Minimal test: does Equivariant Imaging (EI) help seismic Q compensation?
Linear reconstructor x_hat = W y (numpy, hand-written gradients).

Compared:
  Tikhonov  : x = (A^T A + mu I)^-1 A^T y, mu tuned with ground truth (oracle-tuned classical baseline)
  MC        : self-supervised, data consistency only        ||A W y - y||^2
  MC+EI     : self-supervised, + shift equivariance          ||W A T_g W y - T_g W y||^2
  Supervised: optimal linear (Wiener) estimator from labels (upper bound)
"""
import numpy as np, json, sys

rng = np.random.default_rng(0)
N, dt, Q, fp = 256, 0.004, 50.0, 25.0
SNR_DB = float(sys.argv[1]) if len(sys.argv) > 1 else 30.0
t = np.arange(N) * dt

# ---------- attenuation operator (Kolsky constant-Q, amplitude + dispersion) ----------
def atten_matrix(Q, N, dt):
    M = 2 * N                                   # pad to avoid wrap-around
    f = np.fft.rfftfreq(M, dt); f[0] = 1e-6
    fh = 1 / (2 * dt)
    A = np.zeros((N, N))
    for j in range(N):
        tau = j * dt
        H = np.exp(-np.pi * f * tau / Q) * np.exp(-1j * 2 * np.pi * f * tau * (1 - np.log(f / fh) / (np.pi * Q)))
        A[:, j] = np.fft.irfft(H, n=M)[:N]
    return A
A = atten_matrix(Q, N, dt)

# ---------- data: sparse reflectivity * Ricker = stationary trace x ----------
tw = (np.arange(-40, 41)) * dt
ricker = (1 - 2 * (np.pi * fp * tw) ** 2) * np.exp(-(np.pi * fp * tw) ** 2)
def make(n):
    r = (rng.random((N, n)) < 0.08) * rng.standard_normal((N, n))
    x = np.stack([np.convolve(r[:, i], ricker, mode="same") for i in range(n)], 1)
    y0 = A @ x
    sig = np.sqrt((y0 ** 2).mean()) * 10 ** (-SNR_DB / 20)
    return x, y0 + sig * rng.standard_normal(y0.shape)
Xtr, Ytr = make(4000)
Xte, Yte = make(500)

# ---------- metrics ----------
deep = t > 0.6
fz = np.fft.rfftfreq(deep.sum(), dt); band = (fz >= 30) & (fz <= 60)
def snr(x, xh, m=slice(None)): return 10 * np.log10((x[m] ** 2).sum() / ((x[m] - xh[m]) ** 2).sum())
def band_ratio(x, xh):  # recovered / true amplitude in 30-60 Hz, deep window
    S = lambda z: np.abs(np.fft.rfft(z[deep], axis=0)).mean(1)
    return float(S(xh)[band].mean() / S(x)[band].mean())
def evaluate(name, Xh):
    return dict(method=name, snr_all=round(snr(Xte, Xh), 2), snr_deep=round(snr(Xte, Xh, deep), 2),
                band_30_60Hz_deep=round(band_ratio(Xte, Xh), 3))

results = []
# Tikhonov, oracle-tuned mu
best = None
for mu in np.logspace(-6, 0, 25):
    Xh = np.linalg.solve(A.T @ A + mu * np.eye(N), A.T @ Yte)
    s = snr(Xte, Xh)
    if best is None or s > best[0]: best = (s, mu, Xh)
results.append(evaluate(f"Tikhonov(mu={best[1]:.1e})", best[2]))

# Supervised linear (Wiener) estimator
W_sup = Xtr @ Ytr.T @ np.linalg.inv(Ytr @ Ytr.T + 1e-6 * np.eye(N))
results.append(evaluate("Supervised (upper bound)", W_sup @ Yte))

# ---------- self-supervised training ----------
def shift_cols(Z, g):           # non-circular shift down by g (zero fill)
    out = np.zeros_like(Z)
    if g >= 0: out[g:] = Z[:N - g]
    else:      out[:N + g] = Z[-g:]
    return out
def shift_cols_T(Z, g): return shift_cols(Z, -g)   # adjoint (up to truncation)

def train(lam_ei, iters=int(sys.argv[2]) if len(sys.argv)>2 else 3000, bs=128, lr=2e-3, gmax=80, seed=1):
    r = np.random.default_rng(seed)
    W = np.zeros((N, N))
    m = np.zeros_like(W); v = np.zeros_like(W); b1, b2 = 0.9, 0.999
    for it in range(1, iters + 1):
        Y = Ytr[:, r.integers(0, Ytr.shape[1], bs)]
        Xh = W @ Y
        R1 = A @ Xh - Y
        G = 2 * A.T @ R1 @ Y.T / bs                         # MC gradient
        if lam_ei > 0:
            for gi in np.array_split(np.arange(bs), 4):     # 4 random shifts per batch
                g = int(r.integers(-gmax, gmax + 1))
                Yg, Xg = Y[:, gi], Xh[:, gi]
                Z = shift_cols(Xg, g)                       # T_g x_hat
                AZ = A @ Z
                R2 = W @ AZ - Z                             # f(A T_g x_hat) - T_g x_hat
                Gw = 2 * R2 @ AZ.T                          # direct dependence on W
                dZ = 2 * (W @ A - np.eye(N)).T @ R2         # via Z = T_g W y
                Gw += shift_cols_T(dZ, g) @ Yg.T
                G += lam_ei * Gw / bs
        m = b1 * m + (1 - b1) * G; v = b2 * v + (1 - b2) * G ** 2
        W -= lr * (m / (1 - b1 ** it)) / (np.sqrt(v / (1 - b2 ** it)) + 1e-8)
    return W

W_mc = train(0.0);  results.append(evaluate("MC only (self-sup.)", W_mc @ Yte))
W_ei = train(1.0);  results.append(evaluate("MC + EI (self-sup.)", W_ei @ Yte))

print(json.dumps(dict(SNR_dB=SNR_DB, Q=Q, results=results), indent=1, ensure_ascii=False))
np.savez(f"/home/claude/ei/out_snr{int(SNR_DB)}.npz", Xte=Xte[:, :8], Yte=Yte[:, :8],
         tik=best[2][:, :8], mc=(W_mc @ Yte)[:, :8], ei=(W_ei @ Yte)[:, :8], sup=(W_sup @ Yte)[:, :8],
         spec_true=np.abs(np.fft.rfft(Xte[deep], axis=0)).mean(1),
         spec_tik=np.abs(np.fft.rfft(best[2][deep], axis=0)).mean(1),
         spec_mc=np.abs(np.fft.rfft((W_mc @ Yte)[deep], axis=0)).mean(1),
         spec_ei=np.abs(np.fft.rfft((W_ei @ Yte)[deep], axis=0)).mean(1),
         spec_sup=np.abs(np.fft.rfft((W_sup @ Yte)[deep], axis=0)).mean(1),
         spec_att=np.abs(np.fft.rfft(Yte[deep], axis=0)).mean(1), fz=fz)
