"""Operator-only predictions of the residual Delta t* metric (PLAN_ZH.md sections 2 and 6).

No image is read. For a window pair (anchor window a, reservoir window r) the compensation that a
method applies between the windows is the log-spectral slope of |H(f, x_r)| / |H(f, x_a)|, x = the t*
the method acts on (absolute t*/m, or anchored Delta t*/m for VA, VAfix). Converted to seconds:
applied = slope / pi. Residual = Delta t*_true - applied. In the high-SNR linear regime this is exactly
what the metric measures on a compensated image minus what it measures on R0 (PLAN section 2.4).

Vertical two-way t* profiles: field = F-5 column of the stage-4 velocity (uq_rtm_20261117/cache/fld/vp.bin),
Q = 120 below the seabed (stage-4 nominal); synthetic = macro layers of make_model.m (fine-layer dv
ignored, mean 0), true Q incl. chimney / lens / qgrad. Ray obliquity is ignored, so for the synthetic
gates only the FRACTION applied / Delta t*_true is used (PLAN section 5.4).

usage: python -I op_pred.py <uq_rtm_20261117/cache/fld> <out.json>"""
import sys, os, json
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import H, wls_slope

fld, out = sys.argv[1:3]
METH = ['PH100', 'UQW1_100', 'UQW1_25', 'VA', 'PH25', 'VAfix']
ANCH = {'UQW1_100': False, 'UQW1_25': False, 'VA': True, 'PH25': False, 'VAfix': True, 'PH100': False}


def applied(meth, band, xa, xr, inv_snr=0.0):
    f = np.arange(band[0], band[1] + 1e-9, 0.5)
    y = np.array([np.log(abs(H(meth, fi, np.array([xr]), inv_snr)[0])) - np.log(abs(H(meth, fi, np.array([xa]), inv_snr)[0]))
                  for fi in f])
    return wls_slope(f, y)[0] / np.pi


# ---------------- field: F-5 column ----------------
v = np.fromfile(os.path.join(fld, 'vp.bin'), np.float32).reshape(401, 671).astype(float)
x = -3600 + (np.arange(671) - 25) * 10.0
z = (np.arange(401) - 25) * 10.0
j = np.argmin(abs(x - (-1200.0)))
vc = v[:, j]
Q = np.where(vc < 1500, np.inf, 120.0)
ok = z >= 0
zz, vv, QQ = z[ok], vc[ok], Q[ok]
tau = np.r_[0, np.cumsum(2 * 10 / vv[:-1])]
ts120 = np.r_[0, np.cumsum(2 * 10 / (vv[:-1] * QQ[:-1]))]
win_a, win_r, z_anchor = (2200.0, 2500.0), (2600.0, 3200.0), 2300.0


def tcentre(w):
    m = (zz >= w[0]) & (zz <= w[1])
    tc = 0.5 * (tau[m][0] + tau[m][-1])
    return np.interp(tc, tau, zz), tc


za_c, ta_c = tcentre(win_a)
zr_c, tr_c = tcentre(win_r)
ts_a, ts_r = np.interp(za_c, zz, ts120), np.interp(zr_c, zz, ts120)
ts_anchor = np.interp(z_anchor, zz, ts120)
dts = {'Q120': ts_r - ts_a}
for qv in (60, 80, 100, 135):
    dts['Q%d' % qv] = (ts_r - ts_a) * 120 / qv
R = {'field': dict(column_line_coord_m=float(x[j]), win_a_m=win_a, win_r_m=win_r, anchor_m=z_anchor,
                   win_a_time_centre_depth_m=float(za_c), win_r_time_centre_depth_m=float(zr_c),
                   twt_a_s=[float(np.interp(win_a[0], zz, tau)), float(np.interp(win_a[1], zz, tau))],
                   twt_r_s=[float(np.interp(win_r[0], zz, tau)), float(np.interp(win_r[1], zz, tau))],
                   tstar_Q120_at_anchor_s=float(ts_anchor), tstar_Q120_a_s=float(ts_a), tstar_Q120_r_s=float(ts_r),
                   dtstar_true_by_Q=dts)}
mgrid = [0.5, 0.6, 0.7, 0.85, 1.0, 1.2, 1.4, 1.7, 2.0]
for band in ((6, 23), (6, 40), (6, 55)):
    for isnr in (0.0, 10 ** (-15 / 10)):
        key = 'band_%d_%d_invsnr_%s' % (band[0], band[1], '0' if isnr == 0 else '15dB')
        T = {}
        for meth in METH:
            a = []
            for m in mgrid:
                if ANCH[meth]:
                    xa, xr = max(ts_a - ts_anchor, 0) / m, (ts_r - ts_anchor) / m
                else:
                    xa, xr = ts_a / m, ts_r / m
                a.append(applied(meth, band, xa, xr, isnr))
            T[meth] = dict(applied_s=a, applied_over_dtstar_Q120=[ai / dts['Q120'] for ai in a])
        R['field'][key] = T
R['field']['m'] = mgrid

# ---------------- synthetic: macro models (make_model.m) ----------------
MODELS = {
    'D1': dict(iface=[lambda x: 380 + .05 * (x - 1200), lambda x: 760 + .02 * (x - 1200), lambda x: 1120 - .03 * (x - 1200)],
               vp=[2050, 2450, 2850, 3250], Qp=[60, 45, 75, 110], lens=[(700, 520, 250, 50, 20, -.04), (1700, 900, 300, 60, 25, -.04)],
               anchor=0),
    'D2': dict(iface=[lambda x: 450 - 120 * np.exp(-((x - 1200) / 500) ** 2), lambda x: 820 - 150 * np.exp(-((x - 1200) / 550) ** 2),
                      lambda x: 1180 - 100 * np.exp(-((x - 1200) / 600) ** 2)],
               vp=[2150, 2550, 2950, 3350], Qp=[70, 50, 60, 90], chimney=(1200, 120, 300, 800, 22, -.03), anchor=0),
    'D3': dict(iface=[lambda x: 400 + .03 * (x - 1200), lambda x: 780 + .03 * (x - 1200), lambda x: 1150 + .03 * (x - 1200)],
               vp=[2100, 2500, 2900, 3300], Qp=[60, 45, 70, 100], qgrad=(1.6, .6), anchor=0),
    'D4': dict(iface=[lambda x: 260 + 0 * x, lambda x: 330 + 0 * x, lambda x: 600 + .25 * (x - 1200), lambda x: 1200 - .15 * (x - 1200)],
               vp=[2000, 2250, 2700, 3200, 3400], Qp=[20, 60, 80, 110, 120], anchor=1),
}
# D1 fault (x0 1300, dip 70, throw 40) is ignored here (macro t* changes by < 2%).


def column(spec, xc, dz=4.0):
    zg = np.arange(0, 1600, dz)
    vcol = np.full(zg.size, 1850.0)
    qcol = np.full(zg.size, 80.0)
    tops = [lambda x: 120 + 0 * x] + spec['iface']
    for k in range(len(spec['vp'])):
        top = tops[k](xc)
        bot = tops[k + 1](xc) if k + 1 < len(tops) else np.inf
        m = (zg >= top) & (zg < bot)
        vcol[m] = spec['vp'][k]
        qcol[m] = spec['Qp'][k]
    if 'qgrad' in spec:
        g = spec['qgrad'][0] + (spec['qgrad'][1] - spec['qgrad'][0]) * min(max(xc / 2400, 0), 1)
        qcol[zg > 120] *= g
    for L in spec.get('lens', []):
        e = ((xc - L[0]) / L[2]) ** 2 + ((zg - L[1]) / L[3]) ** 2 <= 1
        vcol[e] *= 1 + L[5]
        qcol[e] = L[4]
    if 'chimney' in spec:
        C = spec['chimney']
        e = (abs(xc - C[0]) <= C[1]) & (zg >= C[2]) & (zg <= C[3])
        vcol[e] *= 1 + C[5]
        qcol[e] = C[4]
    tau = np.r_[0, np.cumsum(2 * dz / vcol[:-1])]
    ts = np.r_[0, np.cumsum(2 * dz / (vcol[:-1] * qcol[:-1]))]
    return zg, tau, ts


# synthetic window pair (PLAN_ZH.md section 4.1, after the D2 smoke): both inside W_deep
SYN_WA, SYN_WR = (820.0, 1100.0), (1250.0, 1550.0)
SYN = {}
mg5 = [0.5, 0.7, 1.0, 1.4, 2.0]
for name, spec in MODELS.items():
    xs = np.arange(300, 2101, 50.0)
    rows = []
    for xc in xs:
        zg, tau, ts = column(spec, xc)
        za = spec['iface'][spec['anchor']](xc)
        wa, wr = SYN_WA, SYN_WR
        ca = np.interp(0.5 * (np.interp(wa[0], zg, tau) + np.interp(wa[1], zg, tau)), tau, zg)
        cr = np.interp(0.5 * (np.interp(wr[0], zg, tau) + np.interp(wr[1], zg, tau)), tau, zg)
        rows.append([np.interp(ca, zg, ts), np.interp(cr, zg, ts), np.interp(za, zg, ts),
                     np.interp(wa[1], zg, tau) - np.interp(wa[0], zg, tau), np.interp(wr[1], zg, tau) - np.interp(wr[0], zg, tau)])
    rows = np.array(rows)
    ts_a, ts_r, ts_z = np.median(rows[:, 0]), np.median(rows[:, 1]), np.median(rows[:, 2])
    d = ts_r - ts_a
    S = dict(tstar_a_s=float(ts_a), tstar_r_s=float(ts_r), tstar_anchor_s=float(ts_z), dtstar_true_vert_s=float(d),
             twt_len_a_s=float(np.median(rows[:, 3])), twt_len_r_s=float(np.median(rows[:, 4])), m=mg5)
    for band in ((8, 50), (8, 40), (6, 23)):
        T = {}
        for meth in METH:
            fr = []
            for m in mg5:
                if ANCH[meth]:
                    xa, xr = max(ts_a - ts_z, 0) / m, (ts_r - ts_z) / m
                else:
                    xa, xr = ts_a / m, ts_r / m
                fr.append(applied(meth, band, xa, xr) / d)
            T[meth] = dict(applied_fraction=fr, predicted_residual_fraction=[1 - q for q in fr])
        S['band_%d_%d' % band] = T
    SYN[name] = S
R['synthetic'] = SYN
R['synthetic_windows_m'] = dict(a=SYN_WA, r=SYN_WR)
json.dump(R, open(out, 'w'), indent=1)


def show(d, ind=''):
    for k, v in d.items():
        if isinstance(v, dict):
            print(ind + k)
            show(v, ind + '  ')
        elif isinstance(v, list) and v and isinstance(v[0], float):
            print(ind + k, ' '.join('%.3g' % q for q in v))
        else:
            print(ind + k, v)


show(R)
