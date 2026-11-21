"""Residual Delta t* metric (PLAN_ZH.md section 2), synthetic implementation.

For an image I(z, x) and a reference image Ref(z, x) (synthetic: elastic data migrated with v100 = the
'well' known at every column; field: well reflectivity at the well only), per column j and window w:
  1. prep() of score_stage2.m (remove 40 m Gaussian low-wavenumber trend, zero z < 140 m) on both;
  2. depth -> two-way time with the migration velocity of the image (tau = 2 cumsum dz / v), linear
     resampling to 1 ms between tau(z_top) and tau(z_bot) of the window;
  3. K sine multitapers; lag tau_j of the I-Ref cross-correlation peak (|lag| <= 10 ms, parabolic),
     cross spectrum S_IR = sum_k I_k Ref_k^* exp(+i 2 pi f tau_j), S_RR = sum_k |Ref_k|^2;
  4. wavelet amplitude of the window over a column set: |W_w(f)| = |sum_j S_IR| / sum_j S_RR;
  5. y(f) = ln|W_r| - ln|W_a|, OLS fit over the band: dt_res = -slope / pi (linear) and the slope at the
     band centre of c0 + c1 (f - fc) + c2 (f - fc)^2 (curvature as nuisance);
  6. SE: delete-one jackknife over blocks of 25 columns (100 m); single-well arm: over the K tapers.
usage: python -I metric.py <img_<model>.mat> <out.json> [--wa SPEC] [--wr SPEC] [--arms main,fieldlike]
  SPEC = anchor:<d1>:<d2> (z_a(x)+d1 .. z_a(x)+d2) or abs:<z1>:<z2>; defaults --wa abs:820:1100 --wr abs:1250:1550
  (PLAN_ZH.md section 4.1: the first smoke used anchor:10:310 / abs:820:1550, kept as smoke_D2_v1)"""
import sys, os, json
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'novelty_evidence_20261025', 'cache', 'deps'))
import scipy.io as sio
from common import sine_tapers

DT, NF = 0.001, 2048
FAX = np.fft.rfftfreq(NF, DT)


def prep(y, z):
    s = 10
    h = np.exp(-.5 * (np.arange(-30, 31) / s) ** 2)
    h /= h.sum()
    yp = np.vstack([np.repeat(y[:1], 30, 0), y, np.repeat(y[-1:], 30, 0)])
    sm = np.apply_along_axis(lambda c: np.convolve(c, h, 'valid'), 0, yp)
    out = y - sm
    out[z < 140] = 0
    return out


def to_time(col, z, tau, zlo, zhi):
    t0, t1 = np.interp(zlo, z, tau), np.interp(zhi, z, tau)
    tt = np.arange(t0, t1, DT)
    return np.interp(tt, tau, col)


def lag(a, b, maxlag=10):
    c = np.correlate(a - a.mean(), b - b.mean(), 'full')
    mid = b.size - 1
    lo, hi = mid - maxlag, mid + maxlag + 1
    seg = c[lo:hi]
    k = int(np.argmax(seg))
    d = 0.0
    if 0 < k < seg.size - 1:
        y0, y1, y2 = seg[k - 1], seg[k], seg[k + 1]
        den = y0 - 2 * y1 + y2
        d = 0.5 * (y0 - y2) / den if den != 0 else 0.0
    return (k - maxlag + d) * DT                                             # a(t) ~ b(t - lag)


def spectra(img, ref, z, tau, wins, cols, K=3, ref_noise=None, rng=None):
    """per column: S_IR (complex) and S_RR for each window; wins: list of (zlo(x), zhi(x)) arrays."""
    out = []
    for (zlo, zhi) in wins:
        SIR = np.zeros((len(cols), FAX.size), complex)
        SRR = np.zeros((len(cols), FAX.size))
        PK = []
        for n, j in enumerate(cols):
            a = to_time(img[:, j], z, tau[:, j], zlo[j], zhi[j])
            b = to_time(ref[:, j], z, tau[:, j], zlo[j], zhi[j])
            if ref_noise is not None:
                b = b + ref_noise(b, rng)
            L = lag(a, b)
            tp = sine_tapers(a.size, K)
            A = np.fft.rfft(tp * a, NF, axis=-1)
            B = np.fft.rfft(tp * b, NF, axis=-1)
            ph = np.exp(1j * 2 * np.pi * FAX * L)
            SIR[n] = np.sum(A * np.conj(B), 0) * ph
            SRR[n] = np.sum(np.abs(B) ** 2, 0)
            PK.append((A * np.conj(B) * ph, np.abs(B) ** 2))
        out.append((SIR, SRR, PK))
    return out


def fit(y, band):
    m = (FAX >= band[0]) & (FAX <= band[1])
    f = FAX[m]
    s1 = np.polyfit(f, y[m], 1)[0]
    fc = 0.5 * (band[0] + band[1])
    c = np.polyfit(f - fc, y[m], 2)
    return -s1 / np.pi, -c[1] / np.pi, c[0]


def dtres(sp, idx, band):
    (Sa, Ra, _), (Sr, Rr, _) = sp
    ya = np.log(np.abs(Sa[idx].sum(0))) - np.log(Ra[idx].sum(0))
    yr = np.log(np.abs(Sr[idx].sum(0))) - np.log(Rr[idx].sum(0))
    return fit(yr - ya, band)


def jack_blocks(sp, ncol, band, block=25):
    nb = ncol // block
    full = dtres(sp, np.arange(ncol), band)
    reps = []
    for b in range(nb):
        keep = np.r_[0:b * block, (b + 1) * block:ncol]
        reps.append(dtres(sp, keep, band))
    reps = np.array(reps)
    se = np.sqrt((nb - 1) / nb * np.sum((reps - reps.mean(0)) ** 2, 0))
    return full, se


def jack_tapers(sp, band):
    """single-well arm: delete-one over tapers (pooled over the column set)."""
    (_, _, Pa), (_, _, Pr) = sp
    K = Pa[0][0].shape[0]

    def est(keep):
        Sa = sum(p[0][keep].sum(0) for p in Pa); Ra = sum(p[1][keep].sum(0) for p in Pa)
        Sr = sum(p[0][keep].sum(0) for p in Pr); Rr = sum(p[1][keep].sum(0) for p in Pr)
        return fit((np.log(np.abs(Sr)) - np.log(Rr)) - (np.log(np.abs(Sa)) - np.log(Ra)), band)

    full = est(np.arange(K))
    reps = np.array([est(np.r_[0:k, k + 1:K]) for k in range(K)])
    se = np.sqrt((K - 1) / K * np.sum((reps - reps.mean(0)) ** 2, 0))
    return full, se


def tstar_vert(v, Q, dz):
    return np.vstack([np.zeros((1, v.shape[1])), np.cumsum(2 * dz / (v[:-1] * Q[:-1]), 0)])


def main():
    fn, out = sys.argv[1:3]
    opt = dict(zip(sys.argv[3::2], sys.argv[4::2]))
    arms = opt.get('--arms', 'main,fieldlike').split(',')
    S = sio.loadmat(fn, squeeze_me=True, struct_as_record=False)
    z, x, dz = S['z'].astype(float), S['x'].astype(float), float(S['dx'])
    vb = S['vb'].astype(float)
    tau = np.vstack([np.zeros((1, vb.shape[1])), np.cumsum(2 * dz / vb[:-1], 0)])
    za = np.atleast_1d(S['za']).astype(float)
    ms = np.atleast_1d(S['m']).astype(float)
    meths = [str(q) for q in np.atleast_1d(S['methods'])]
    ref = prep(S['ref'].astype(float), z)
    cols = np.nonzero((x >= 300) & (x <= 2100))[0]
    def win(spec):
        k, a, b = spec.split(':')
        if k == 'anchor':
            return (za + float(a), za + float(b))
        return (np.full_like(x, float(a)), np.full_like(x, float(b)))
    wa, wr = win(opt.get('--wa', 'abs:820:1100')), win(opt.get('--wr', 'abs:1250:1550'))
    # vertical Delta t* (true Q, v100), between the time centres of the windows, median over columns
    tsv = tstar_vert(S['v'].astype(float), S['Q'].astype(float), dz)
    dv = []
    for j in cols:
        ca = 0.5 * (np.interp(wa[0][j], z, tau[:, j]) + np.interp(wa[1][j], z, tau[:, j]))
        cr = 0.5 * (np.interp(wr[0][j], z, tau[:, j]) + np.interp(wr[1][j], z, tau[:, j]))
        dv.append(np.interp(np.interp(cr, tau[:, j], z), z, tsv[:, j]) - np.interp(np.interp(ca, tau[:, j], z), z, tsv[:, j]))
    R = dict(model=str(S['name']), m=ms.tolist(), dtstar_true_vert_s=float(np.median(dv)), arms={},
             wa=opt.get('--wa', 'abs:820:1100'), wr=opt.get('--wr', 'abs:1250:1550'))
    imgs = {'R0': [prep(S['R0'].astype(float), z)]}
    for q in meths:
        A = S[q]
        imgs[q] = [prep(A[:, :, r].astype(float), z) for r in range(A.shape[2])]
    BANDS = {'b8_50': (8, 50), 'b6_23': (6, 23), 'b8_40': (8, 40)}
    if 'main' in arms:
        T = {}
        for q, L in imgs.items():
            sp = [spectra(I, ref, z, tau, [wa, wr], cols) for I in L]
            T[q] = {}
            for bn, b in BANDS.items():
                res = [jack_blocks(s, len(cols), b) for s in sp]
                T[q][bn] = dict(dt_lin_s=[float(r[0][0]) for r in res], se_lin_s=[float(r[1][0]) for r in res],
                                dt_quad_s=[float(r[0][1]) for r in res], se_quad_s=[float(r[1][1]) for r in res],
                                curv=[float(r[0][2]) for r in res])
        R['arms']['main'] = T
    if 'fieldlike' in arms:
        rng = np.random.default_rng(20261121)
        T = {}
        for x0 in (600.0, 1200.0, 1800.0):
            cw = np.nonzero(abs(x - x0) <= 20)[0]
            for tie in ('exact', 'rho063'):
                def rn(b, g, tie=tie):
                    # independent 'tie mismatch' with the spectrum of the reference window itself (a white
                    # mismatch would bias S_RR where the band-limited reference is weak), correlation 0.63
                    if tie == 'exact':
                        return 0
                    e = np.fft.irfft(np.fft.rfft(g.standard_normal(b.size)) * np.abs(np.fft.rfft(b)), b.size)
                    return e / e.std() * b.std() * np.sqrt(1 / 0.63 ** 2 - 1)
                key = 'x%d_%s' % (x0, tie)
                T[key] = {}
                for q, L in imgs.items():
                    rs = []
                    for I in L:
                        g = np.random.default_rng(int(x0) + (0 if tie == 'exact' else 7))   # same mismatch for every image
                        sp = spectra(I, ref, z, tau, [wa, wr], cw, K=5, ref_noise=rn, rng=g)
                        rs.append(jack_tapers(sp, (6, 23)))
                    T[key][q] = dict(dt_lin_s=[float(r[0][0]) for r in rs], se_lin_s=[float(r[1][0]) for r in rs])
        R['arms']['fieldlike_b6_23_single_well'] = T
    json.dump(R, open(out, 'w'), indent=1)
    print('model', R['model'], 'dt*_true(vertical) %.5f s' % R['dtstar_true_vert_s'], 'm', ms)
    if 'main' in R['arms']:
        for q, T in R['arms']['main'].items():
            for bn in BANDS:
                t = T[bn]
                print('%-9s %-6s lin ' % (q, bn) + ' '.join('%+.5f(%.5f)' % (a, b) for a, b in zip(t['dt_lin_s'], t['se_lin_s'])) +
                      ' | quad ' + ' '.join('%+.5f(%.5f)' % (a, b) for a, b in zip(t['dt_quad_s'], t['se_quad_s'])))
    if 'fieldlike_b6_23_single_well' in R['arms']:
        for key, T in R['arms']['fieldlike_b6_23_single_well'].items():
            for q, t in T.items():
                print('%-14s %-9s ' % (key, q) + ' '.join('%+.5f(%.5f)' % (a, b) for a, b in zip(t['dt_lin_s'], t['se_lin_s'])))


if __name__ == '__main__':
    main()
