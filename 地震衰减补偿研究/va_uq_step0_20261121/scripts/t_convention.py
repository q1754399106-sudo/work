"""1-D toy: does the velocity/phase-reference convention decide whether Q compensation looks useful,
and does an anchored, velocity-referenced operator (AUQ) remove the field-type loss?

Trace in two-way time (1 ms). Kjartansson attenuation with dispersion referenced to f_d:
  D(f) = sum_j r_j W(f) exp[(-pi f + i k_d) t*_j] exp(-i 2 pi f tau_j) + N(f),   k_d = 2 f ln(f/f_d)
tau_j is the arrival time at f_d. Two conventions for the same physics:
  'syn'   f_d = 100 Hz (stage 1-3 RTM: migration velocity = 100 Hz phase velocity, elastic reference)
  'field' f_d = 25 Hz  (tomographic velocity calibrated in band: events at 25 Hz already at true time)
  'field+P' as 'field' plus a time-invariant whitening of the whole trace (decon-like processing)
Reference: r * W (source wavelet, attenuation free) on the same time axis.

Operators (nonstationary, per output time tau; Q multiplier m: t*_c = t*/m; Bayesian expected Wiener
with lognormal u, sigma 0.447, oracle noise power; s = 70 dB floor):
  R0    identity
  PH    phase only, f_op = 100;  PH25 phase only, f_op = 25
  UQ100 expected Wiener on absolute t*, f_op = 100       (= current UQ-RTM convention)
  UQ25  expected Wiener on absolute t*, f_op = 25        (post-hoc 'A')
  AUQ   expected Wiener on dt* = t* - t*_a below an anchor at 0.8 s, f_op = 25, signal power of the anchor
Scores in 1.2-1.7 s: E1 = (c0 + max_|L|<=4 ms c_L)/2 (stage-2 form), TIE = max_|L|<=20 ms c_L, and
TIEbox = TIE after box-filtering output and reference to R0's own -10 dB band (stage-4 field form).
"""
import json
import sys
import numpy as np

OUT = sys.argv[1] if len(sys.argv) > 1 else 't_convention.json'
dt, N = 1e-3, 2048
t = np.arange(N) * dt
f = np.fft.rfftfreq(N, dt)
Q, FP = 60.0, 25.0
band = (f >= 3) & (f <= 90)
fs = np.where(f > 0, f, 1e-6)
W = (2 / np.sqrt(np.pi)) * fs ** 2 / FP ** 3 * np.exp(-(fs / FP) ** 2)        # Ricker amplitude spectrum
W = W / W.max() * band
TA = 0.8
win = (t >= 1.2) & (t <= 1.7)
MS = [0.5, 0.7, 1.0, 1.4, 2.0]
gh_x, gh_w = np.polynomial.hermite_e.hermegauss(21)
gh_w = gh_w / gh_w.sum()
SIG = 0.447


def k_of(fref):
    return 2 * fs * np.log(fs / fref)


def synth(rng, fd, snr_db, whiten):
    tau = np.arange(0.3, 1.95, dt)
    r = rng.standard_normal(tau.size) * (rng.random(tau.size) < 0.25)       # dense-ish reflectivity
    tstar = tau / Q
    kd = k_of(fd)
    D = np.zeros(f.size, complex)
    for rj, tj, sj in zip(r, tau, tstar):
        D += rj * W * np.exp((-np.pi * fs + 1j * kd) * sj - 2j * np.pi * fs * tj)
    d = np.fft.irfft(D, N)
    ref = np.zeros(N); ref[np.round(tau / dt).astype(int)] = r
    ref = np.fft.irfft(np.fft.rfft(ref) * W, N)
    sp = np.mean(d[(t > 0.3) & (t < 0.9)] ** 2)                             # noise relative to shallow power
    d = d + rng.standard_normal(N) * np.sqrt(sp * 10 ** (-snr_db / 10))
    Dn = np.fft.rfft(d)
    if whiten:                                                              # decon-like whitening, whole trace
        a = np.abs(Dn); a = np.convolve(a, np.ones(31) / 31, 'same')
        P = np.where(band, a.max() / np.maximum(a, 1e-3 * a.max()), 0)
        P = P / np.sqrt(np.mean((P * W)[band] ** 2) / np.mean(W[band] ** 2))
        Dn = Dn * P
    noise_pow = np.mean(np.abs(np.fft.rfft(rng.standard_normal(N) * np.sqrt(sp * 10 ** (-snr_db / 10)))) ** 2)
    sig_pow = float(np.sum(r ** 2))                                          # E|R(f)|^2 of the reflectivity
    return Dn, ref, noise_pow, sig_pow


IB = np.where(band)[0]
FB = fs[IB]
EB = np.exp(2j * np.pi * np.outer(t, FB))                                   # (N, nb)


def apply(Dn, tstar_of_t, fop, m, S, n):
    """y(tau) = Re sum_f K(tau,f) D(f) e^{i2pi f tau}; K = E[G*] S / (E|G|^2 S + n + s S), band bins only."""
    k = 2 * FB * np.log(FB / fop)
    xs = tstar_of_t / m                                                     # model t* (or dt*) per output time
    u = np.exp(SIG * gh_x)[:, None, None]
    G = np.exp((-np.pi * FB + 1j * k)[None, None, :] * (u * xs[None, :, None]))   # (nu, N, nb)
    EG = np.tensordot(gh_w, G, axes=1)
    EG2 = np.tensordot(gh_w, np.abs(G) ** 2, axes=1)
    Sv = S[IB]
    K = np.conj(EG) * Sv[None, :] / (EG2 * Sv[None, :] + n + 1e-7 * Sv.max())
    return np.real((K * Dn[IB][None, :] * EB).sum(axis=1)) * 2 / N


def phase_only(Dn, tstar_of_t, fop, m):
    k = 2 * FB * np.log(FB / fop)
    K = np.exp(-1j * np.outer(tstar_of_t / m, k))
    return np.real((K * Dn[IB][None, :] * EB).sum(axis=1)) * 2 / N


def cc(a, b):
    a = a - a.mean(); b = b - b.mean(); return float(a @ b / np.sqrt((a @ a) * (b @ b)))


def r0_band(y0):
    """-10 dB band of the uncompensated trace in the scoring window (as fld_eval.m band())."""
    seg = y0[win] * np.hanning(win.sum())
    S = np.abs(np.fft.rfft(seg, 1024)); fr = np.fft.rfftfreq(1024, dt)
    i = int(np.argmax(S)); thr = S[i] * 10 ** (-0.5)
    lo = np.where(S[:i] < thr)[0]; hi = np.where(S[i:] < thr)[0]
    return fr[lo[-1]] if lo.size else fr[0], fr[i + hi[0]] if hi.size else fr[-1]


def boxpass(x, lo, hi):
    X = np.fft.rfft(x); return np.fft.irfft(X * ((f >= lo) & (f <= hi)), N)


def score(y, ref, box=None):
    idx = np.where(win)[0]
    c = {L: cc(y[idx + L], ref[idx]) for L in range(-20, 21)}
    out = [0.5 * (c[0] + max(c[L] for L in range(-4, 5))), max(c.values())]
    if box is not None:
        yb, rb = boxpass(y, *box), boxpass(ref, *box)
        out.append(max(cc(yb[idx + L], rb[idx]) for L in range(-20, 21)))
    return out


def run(conv, snr_db, nrep=6, seed=7):
    fd = 100.0 if conv == 'syn' else 25.0
    whiten = conv == 'field+P'
    res = {}
    rng = np.random.default_rng(seed)
    tstar = t / Q
    dts = np.maximum(t - TA, 0) / Q
    for rep in range(nrep):
        Dn, ref, npow, spow = synth(rng, fd, snr_db, whiten)
        # signal power per frequency at the source side: |W|^2 * (reflectivity power); noise power per bin
        S0 = spow * W ** 2
        Sa = S0 * np.exp(-2 * np.pi * fs * TA / Q)                         # anchor signal power (attenuated)
        meths = {'R0': None}
        for m in MS:
            out = {}
            out['R0'] = np.fft.irfft(Dn * band, N)
            # phase only, f_op = 100
            out['PH'] = phase_only(Dn, tstar, 100.0, m)
            out['PH25'] = phase_only(Dn, tstar, 25.0, m)
            out['UQ100'] = apply(Dn, tstar, 100.0, m, S0, npow)
            out['UQ25'] = apply(Dn, tstar, 25.0, m, S0, npow)
            out['AUQ'] = apply(Dn, dts, 25.0, m, Sa, npow)
            box = r0_band(out['R0'])
            for key, y in out.items():
                res.setdefault(key, {}).setdefault(m, []).append(score(y, ref, box))
    summ = {key: {str(m): [float(np.mean([v[i] for v in vals])) for i in range(3)]
                  for m, vals in d.items()} for key, d in res.items()}
    return summ


if __name__ == '__main__':
    allres = {}
    for conv in ['syn', 'field', 'field+P']:
        for snr in [30.0, 14.0]:
            key = f'{conv}_snr{int(snr)}'
            allres[key] = run(conv, snr)
            print(f'--- {key}   (E1 / TIE / TIEbox)   m = ' + ' '.join(f'{m:g}' for m in MS))
            for meth in ['R0', 'PH', 'PH25', 'UQ100', 'UQ25', 'AUQ']:
                row = allres[key][meth]
                print(f'   {meth:6s} ' + '  '.join('/'.join(f'{v:.2f}' for v in row[str(m)]) for m in MS))
    json.dump(allres, open(OUT, 'w'), indent=1)
