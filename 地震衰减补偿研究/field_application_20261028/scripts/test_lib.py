"""Synthetic self-test of lib.py (no field data): band-limited zero-phase wavelet, layered Q, white noise.

Checks that the SNR band and wavelet estimate are sensible and that SC-a recovers the band-limited target
better than the input does. Prints the numbers; asserts only gross failures.
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import lib
from lib import np, bm

rng = np.random.default_rng(1)
n, nx = 451, 40
tw = 1.2 + np.arange(n) * lib.DT
cm = dict(dt=lib.DT, nfft=2048, fref=100, t=tw, tw=tw, vp=np.array([1800., 2200., 3200., 3200.]),
          h=np.array([630., 1440., 1240.]), Q=np.array([80., 80., 130., 130.]), solver_offset=0)
f = np.fft.rfftfreq(512, lib.DT)
amp = np.exp(-((f - 30) / 18) ** 2) * (f > 4) * (f < 75)          # true source spectrum, high-cut 75 Hz
w = np.fft.irfft(amp, 512); w = np.roll(w, 100)[:201]; w /= np.abs(w).max()
r = np.zeros((n, nx)); ev = rng.choice(np.arange(20, n - 20), 40, replace=False)
for e in ev:
    r[e + np.round(np.linspace(0, rng.integers(-3, 4), nx)).astype(int), np.arange(nx)] = rng.normal()
B, P, scale = lib.dictionary(cm, w, w)
clean = (B * scale) @ r
D = clean + 0.05 * np.std(clean) * rng.normal(size=clean.shape)
wv = lib.wavelets(D, tw, cm)
print("f_peak %.1f  fmax shallow %.1f  target %.1f" % (wv["f_peak"], wv["fmax_shallow"], wv["fmax_target"]))
target = lib.toeplitz_wavelet(wv["w_out"], n) @ r
live = np.ones(nx, bool); K, _, _ = bm.predictor(D, live)
sigma = np.full(nx, 0.05 * np.std(clean))
out, info, _ = lib.solve(cm, cm["Q"], wv, D, sigma, K, live, *lib.SEL_SC)
cc = lambda a, b: np.corrcoef(a.ravel(), b.ravel())[0, 1]
print("corr(input, target) %.3f  corr(SC-a, target) %.3f  converged %s" % (cc(D, target), cc(out, target), info["converged"]))
c, qi = lib.qscan(cm, cm["Q"], wv, D, sigma, K, live)
print("qscan factor %.3f boundary %s" % (c, qi["at_boundary"]))
assert 15 < wv["f_peak"] < 45 and cc(out, target) > cc(D, target)
