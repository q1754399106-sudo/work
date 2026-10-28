"""Self-tests of semi.py on synthetic inputs only (no field data)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import semi  # noqa: E402
from semi import np, lib, TW, TF  # noqa: E402

rng = np.random.default_rng(1)
cm = semi.curtain_cm()
n = len(TW)

# 1. 1 ms modelling agrees with the 4 ms solver dictionary for reflectors on the 4 ms grid (layered t*).
w = np.exp(-0.5 * (np.arange(-100, 101) * semi.DT / 0.006) ** 2); w /= w.max()      # smooth zero-phase wavelet
B, P, scale = lib.dictionary(cm, w, w)
r4 = np.zeros(n); r4[rng.choice(np.arange(20, n - 20), 40, replace=False)] = rng.normal(size=40)
d_dict = (B * scale) @ r4
R = np.zeros((len(TF), 1)); R[np.round((TW - TF[0]) / semi.FINE).astype(int), 0] = r4
f = np.fft.rfftfreq(8192, semi.FINE)
S = semi.source_spectrum(w, f)
ts = np.interp(TF, TW, lib.absorption(cm, TW)); ts[TF < TW[0]] = lib.absorption(cm, TF[TF < TW[0]])
d_mod = semi.model(R, ts, S, f, float(cm["fref"]))[:, 0]
c1 = np.corrcoef(d_dict[20:-20], d_mod[20:-20])[0, 1]
amp = np.linalg.norm(d_mod) / np.linalg.norm(d_dict)
print(f"1 modelling vs dictionary: corr {c1:.5f}, amplitude ratio {amp:.4f}")
assert c1 > 0.999 and abs(amp - 1) < 0.02

# 2. Stabilised inverse Q: maximum gain as requested; with a high limit it undoes the attenuation in band.
for g in (10, 20, 30):
    s = semi.gain_sigma2(g); b = np.logspace(-6, 0, 20001); lam = (b + s) / (b * b + s)
    assert abs(20 * np.log10(lam.max()) - g) < 1e-3, g
D0 = semi.model(R, None, S, f, float(cm["fref"]))[:, 0]
X = semi.inverse_q(d_mod[:, None], cm, cm["Q"], 80)[:, 0]
band = lambda x: semi.bandpass(x[:, None], 10, 60)[:, 0]  # noqa: E731
c2 = np.corrcoef(band(X)[60:-20], band(D0)[60:-20])[0, 1]
print(f"2 inverse Q (80 dB) vs unattenuated, 10-60 Hz: corr {c2:.4f}")
assert c2 > 0.98

# 3. Residual noise: unit variance on white noise; small leakage of a dipping band-limited event set.
E = semi.residual_noise(rng.standard_normal((n, 201)), np.zeros((n, 200)))
print(f"3 residual of white noise: variance {E.var():.3f}")
assert abs(E.var() - 1) < 0.03
base = semi.bandpass(rng.standard_normal((n, 1)), 5, 60)[:, 0]
Sd = np.stack([semi.sinc_shift(base, np.arange(n) - 0.6 * j) for j in range(201)], axis=1)
L = semi.local_lags(Sd)
Es = semi.residual_noise(Sd, L)
leak = np.mean(Es[40:-40, 5:-5] ** 2) / np.mean(Sd[40:-40, 5:-5] ** 2)
print(f"3 dip 0.6 samples/trace: estimated lag {np.median(L):.3f}, leakage power {leak:.4f}")
assert abs(np.median(L) - 0.6) < 0.05 and leak < 0.02

# 4. Gaussian noise keeps the window rms of the reference noise.
Eref = np.concatenate([rng.standard_normal((n // 2, 201)), 3 * rng.standard_normal((n - n // 2, 201))])
G = semi.gaussian_noise(Eref, rng)
for a, b in [(1.2, 1.5), (2.7, 3.0)]:
    m = semi.win(TW, a, b); ratio = np.sqrt(np.mean(G[m] ** 2) / np.mean(Eref[m] ** 2))
    print(f"4 gaussian noise rms ratio {a}-{b} s: {ratio:.3f}")
    assert abs(ratio - 1) < 0.1
print("all self-tests passed")
