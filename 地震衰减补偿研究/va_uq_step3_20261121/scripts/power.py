"""Detection power of the residual Delta t* metric on the Volve line (PLAN_ZH.md section 3).

Quantity: residual Delta t* between the anchor window a (2200-2500 m TVDSS) and the reservoir window r
(2600-3200 m) at F-5, from the log-spectral slope of the two extracted wavelets:
   y(f) = ln|W_r(f)| - ln|W_a(f)| = c0 - pi * dt_res * f,   dt_res = -slope / pi.
Targets (Q = 120, stage-4 nominal; op_pred.json): dt_true(a->r) and 0.3 * dt_true.
   'see any compensation benefit'  : SE <= dt_true / 2.8          (R0 residual != 0 at 80 % power, 5 % two-sided)
   'resolve a 30 % Q error'        : SE <= 0.3 * dt_true / 2.8

Error models (all per window, independent between windows):
 1. noise-limited, reflectivity known exactly at K independent columns (best case well calibration):
      var(ln W_w) per Rayleigh bin (1/T_w) = 1 / (2 SNR_w(f) K)   ->  closed form (lower bound) + Monte Carlo
 2. tie-limited: same with SNR replaced by rho^2/(1-rho^2), rho = in-band F-5 tie of R0 (0.63, stage 4)
 3. statistical (no well): one reflectivity realisation per window (layered geology: lateral averaging
    does not remove it); colour from block bootstrap of the F-5 log, Monte Carlo of the slope.
SNR_w(f) = SNR15 * exp(-2 pi (f - 15) t*_w)   (attenuated signal over white noise, flat designatured source),
t*_a = t*(anchor-window centre), t*_r = t*_a + dt_true (Q = 120).

usage: python -I power.py <uq_rtm_20261117/cache/fld> <op_pred.json> <out.json>"""
import sys, os, json
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import wls_slope, sine_tapers

fld, opj, out = sys.argv[1:4]
OP = json.load(open(opj))['field']
dt_true = OP['dtstar_true_by_Q']['Q120']
ts_a = OP['tstar_Q120_a_s']
Ta = OP['twt_a_s'][1] - OP['twt_a_s'][0]
Tr = OP['twt_r_s'][1] - OP['twt_r_s'][0]
tgt_any, tgt_30 = dt_true / 2.8, 0.3 * dt_true / 2.8
BANDS = [(6, 23), (6, 40), (6, 55)]
SNRS = [10, 15, 20, 30]
rng = np.random.default_rng(20261121)
R = dict(dt_true_Q120_s=dt_true, tstar_a_s=ts_a, T_a_s=Ta, T_r_s=Tr, target_any_s=tgt_any, target_30_s=tgt_30,
         slope_30pct_Np_per_Hz=0.3 * np.pi * dt_true)


def snr_f(f, snr15_db, ts):
    return 10 ** (snr15_db / 10) * np.exp(-2 * np.pi * (f - 15) * ts)


def se_closed(band, snr15_db, K, flat=False, snr_override=None):
    """WLS slope SE of the difference (lower bound): var = sum_w 1 / (T_w int w_w(f) (f - fbar_w)^2 df)."""
    f = np.linspace(band[0], band[1], 2001)
    df = f[1] - f[0]
    v = 0.0
    for T, ts in ((Ta, ts_a), (Tr, ts_a + dt_true)):
        s = np.full_like(f, snr_override) if snr_override is not None else (
            np.full_like(f, 10 ** (snr15_db / 10)) if flat else snr_f(f, snr15_db, ts))
        w = 2 * s * K
        fb = np.sum(w * f) / np.sum(w)
        v += 1 / (T * np.sum(w * (f - fb) ** 2) * df)
    return np.sqrt(v) / np.pi                                              # seconds of Delta t*


# ---- 1. closed form, noise-limited ----
C = {}
for b in BANDS:
    for s in SNRS:
        se1 = se_closed(b, s, 1)
        C['%d-%d Hz, SNR15 %d dB' % (b[0], b[1], s)] = dict(
            se_K1_s=se1, se_K1_flatSNR_s=se_closed(b, s, 1, flat=True),
            ratio_to_target30=se1 / tgt_30, ratio_to_target_any=se1 / tgt_any,
            K_needed_30=(se1 / tgt_30) ** 2, K_needed_any=(se1 / tgt_any) ** 2)
R['closed_form_noise_limited'] = C
# ---- 2. tie-limited ----
rho = 0.63
st = rho ** 2 / (1 - rho ** 2)
R['tie_limited'] = {'rho': rho, 'snr_eff': st, **{'%d-%d Hz' % b: dict(se_K1_s=se_closed(b, 0, 1, snr_override=st),
                                                                      ratio_to_target30=se_closed(b, 0, 1, snr_override=st) / tgt_30)
                                                for b in BANDS}}

# ---- Monte Carlo of the estimator (time domain, 1 ms, sine multitapers) ----
DT = 0.001
NF = 2048
fax = np.fft.rfftfreq(NF, DT)


def wavelet_spec(f, fmax, ts):
    """zero-phase source flat in [4, fmax] with 2 Hz cosine tails (stage-4 style), attenuated by exp(-pi f t*)."""
    w = np.ones_like(f)
    a = f < 6
    w[a] = 0.5 * (1 - np.cos(np.pi * np.clip(f[a] - 4, 0, 2) / 2))
    a = f > fmax - 2
    w[a] = 0.5 * (1 + np.cos(np.pi * np.clip(f[a] - fmax + 2, 0, 2) / 2))
    w[f < 4] = 0
    w[f > fmax] = 0
    return w * np.exp(-np.pi * f * ts)


def est_lnW(img, ref, tap):
    """ln |sum_k I_k R_k^*| - ln sum_k |R_k|^2 (cross-spectral wavelet amplitude, known reflectivity)."""
    I = np.fft.rfft(tap * img, NF, axis=-1)
    Rf = np.fft.rfft(tap * ref, NF, axis=-1)
    return np.log(np.abs(np.sum(I * np.conj(Rf), axis=0))) - np.log(np.sum(np.abs(Rf) ** 2, axis=0))


def est_lnA(img, tap):
    """statistical: half log of the multitaper power spectrum (reflectivity unknown, assumed white)."""
    I = np.fft.rfft(tap * img, NF, axis=-1)
    return 0.5 * np.log(np.sum(np.abs(I) ** 2, axis=0))


def f5_refl():
    a = np.fromfile(os.path.join(fld, 'well_F-5.bin'), np.float32).reshape(-1, 4).astype(float)
    zw, r = a[:, 0], a[:, 3]
    v = np.fromfile(os.path.join(fld, 'vp.bin'), np.float32).reshape(401, 671).astype(float)
    x = -3600 + (np.arange(671) - 25) * 10.0
    z = (np.arange(401) - 25) * 10.0
    vc = v[:, np.argmin(abs(x - np.median(a[:, 1])))]
    vw = np.interp(zw, z, vc)
    tw = np.r_[0, np.cumsum(2 * 10 / vw[:-1])]
    t1 = np.arange(0, tw[-1], DT)
    out = np.zeros_like(t1)
    idx = np.clip(np.rint(tw / DT).astype(int), 0, t1.size - 1)
    out[idx] = r                                                           # spikes at the interface times
    return out


RF5 = f5_refl()


def boot_refl(n, block=40):
    """block bootstrap of the F-5 reflectivity in time (40 ms blocks, as specslope_q_volve), unit variance."""
    s = []
    while len(s) < n:
        i = rng.integers(0, RF5.size - block)
        s.extend(RF5[i:i + block])
    s = np.array(s[:n])
    return s / s.std()


def mc(band, snr15_db, K, nrep=300, mode='well', refl='white'):
    fmax = band[1] + 2 if band[1] < 30 else band[1] + 5                   # source band slightly wider than the fit band
    slopes = []
    for _ in range(nrep):
        y = {}
        for T, ts, nm in ((Ta, ts_a, 'a'), (Tr, ts_a + dt_true, 'r')):
            n = int(round(T / DT))
            tap = sine_tapers(n, 3)[:, None, :] if True else None
            W = wavelet_spec(fax, fmax, ts)
            if refl == 'white':
                r = rng.standard_normal((K, n))
            else:
                r = np.array([boot_refl(n) for _ in range(K)])
            if mode == 'stat':
                r[:] = r[0]                                                # layered: same reflectivity in all columns
            sig = np.fft.irfft(np.fft.rfft(r, NF, axis=-1) * W, NF, axis=-1)[:, :n]
            # noise: white, power set so that SNR at 15 Hz = snr15 (per-frequency power ratio, r unit variance)
            w15 = np.interp(15.0, fax, W)
            sn = w15 / np.sqrt(10 ** (snr15_db / 10))
            img = sig + sn * rng.standard_normal((K, n))
            if mode == 'well':
                y[nm] = np.log(np.abs(np.sum(np.fft.rfft(tap * img[None], NF, axis=-1) *
                                             np.conj(np.fft.rfft(tap * r[None], NF, axis=-1)), axis=(0, 1)))) - \
                        np.log(np.sum(np.abs(np.fft.rfft(tap * r[None], NF, axis=-1)) ** 2, axis=(0, 1)))
            else:
                y[nm] = 0.5 * np.log(np.sum(np.abs(np.fft.rfft(tap * img[None], NF, axis=-1)) ** 2, axis=(0, 1)))
        m = (fax >= band[0]) & (fax <= band[1])
        slopes.append(wls_slope(fax[m], (y['r'] - y['a'])[m])[0])
    slopes = np.array(slopes)
    d = -slopes / np.pi
    return dict(mean_dt_res_s=float(d.mean()), sd_s=float(d.std()), bias_s=float(d.mean() - dt_true),
                ratio_sd_to_target30=float(d.std() / tgt_30), ratio_sd_to_target_any=float(d.std() / tgt_any))


M = {}
for b in BANDS:
    for s in (15, 20):
        for K in (1, 5):
            M['well %d-%d Hz SNR15 %d dB K=%d' % (b[0], b[1], s, K)] = mc(b, s, K)
    M['stat F-5 colour %d-%d Hz SNR15 20 dB K=20 layered' % b] = mc(b, 20, 20, mode='stat', refl='f5')
    M['stat white %d-%d Hz SNR15 20 dB K=20 layered' % b] = mc(b, 20, 20, mode='stat', refl='white')
R['monte_carlo'] = M
json.dump(R, open(out, 'w'), indent=1)
print('dt_true(Q120) %.4f s, target any %.5f s, target 30%% %.5f s, Ta %.3f s, Tr %.3f s, t*_a %.4f s' %
      (dt_true, tgt_any, tgt_30, Ta, Tr, ts_a))
for k, v in C.items():
    print('closed %-28s SE %.5f s (flat %.5f)  /target30 %5.1f  /target_any %4.1f  K30 %6.0f  Kany %5.1f' %
          (k, v['se_K1_s'], v['se_K1_flatSNR_s'], v['ratio_to_target30'], v['ratio_to_target_any'], v['K_needed_30'], v['K_needed_any']))
print('tie-limited rho %.2f snr_eff %.2f:' % (rho, st), {k: round(v['se_K1_s'], 5) for k, v in R['tie_limited'].items() if isinstance(v, dict)})
for k, v in M.items():
    print('MC %-48s mean %.5f sd %.5f bias %+.5f  sd/target30 %5.1f  sd/target_any %4.1f' %
          (k, v['mean_dt_res_s'], v['sd_s'], v['bias_s'], v['ratio_sd_to_target30'], v['ratio_sd_to_target_any']))
