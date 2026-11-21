"""Development diagnostic on D2 only (PLAN_ZH.md section 5): per-window log-spectral slopes of R0 and
VA (m = 1) against the elastic reference (t* from the surface, not the difference), the vertical t* of the
window centres, the spectral ratio curves y(f), and the effect of the lag alignment.
usage: python -I diag_d2.py <img_D2.mat> <out.json>"""
import sys, os, json
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import metric as MT

fn, out = sys.argv[1:3]
S = MT.sio.loadmat(fn, squeeze_me=True, struct_as_record=False)
z, x, dz = S['z'].astype(float), S['x'].astype(float), float(S['dx'])
vb = S['vb'].astype(float)
tau = np.vstack([np.zeros((1, vb.shape[1])), np.cumsum(2 * dz / vb[:-1], 0)])
za = np.atleast_1d(S['za']).astype(float)
ref = MT.prep(S['ref'].astype(float), z)
cols = np.nonzero((x >= 300) & (x <= 2100))[0]
wa = (za + 10, za + 310)
wr = (np.full_like(x, 820.0), np.full_like(x, 1550.0))
tsv = MT.tstar_vert(S['v'].astype(float), S['Q'].astype(float), dz)
R = {}
for nm, img in (('R0', S['R0']), ('VA_m1', S['VA'][:, :, 2]), ('UQW1_25_m1', S['UQW1_25'][:, :, 2])):
    I = MT.prep(img.astype(float), z)
    sp = MT.spectra(I, ref, z, tau, [wa, wr], cols)
    d = {}
    for w, (Sx, Rx, _) in zip(('a', 'r'), sp):
        y = np.log(np.abs(Sx.sum(0))) - np.log(Rx.sum(0))
        for bn, b in (('b8_50', (8, 50)), ('b6_23', (6, 23)), ('b8_30', (8, 30))):
            d['%s_%s_tstar_from_surface_s' % (w, bn)] = float(MT.fit(y, b)[0])
        sel = (MT.FAX >= 4) & (MT.FAX <= 56)
        d['%s_y' % w] = y[sel][::2].round(4).tolist()
    yr = np.log(np.abs(sp[1][0].sum(0))) - np.log(sp[1][1].sum(0)) - (np.log(np.abs(sp[0][0].sum(0))) - np.log(sp[0][1].sum(0)))
    sel = (MT.FAX >= 4) & (MT.FAX <= 56)
    d['ratio_y'] = yr[sel][::2].round(4).tolist()
    # coherence-type check: |sum S_IR| / sqrt(sum |I|^2 sum |R|^2) not available here; record incoherent amplitude ratio instead
    R[nm] = d
R['f'] = MT.FAX[(MT.FAX >= 4) & (MT.FAX <= 56)][::2].round(3).tolist()
ca, cr = [], []
for j in cols:
    t_a = 0.5 * (np.interp(wa[0][j], z, tau[:, j]) + np.interp(wa[1][j], z, tau[:, j]))
    t_r = 0.5 * (np.interp(wr[0][j], z, tau[:, j]) + np.interp(wr[1][j], z, tau[:, j]))
    ca.append(np.interp(np.interp(t_a, tau[:, j], z), z, tsv[:, j]))
    cr.append(np.interp(np.interp(t_r, tau[:, j], z), z, tsv[:, j]))
R['tstar_vert_a_s'] = float(np.median(ca))
R['tstar_vert_r_s'] = float(np.median(cr))
json.dump(R, open(out, 'w'), indent=1)
for k in ('R0', 'VA_m1', 'UQW1_25_m1'):
    print(k, {q: round(v, 5) for q, v in R[k].items() if q.endswith('_s')})
print('vertical t*: a %.5f r %.5f' % (R['tstar_vert_a_s'], R['tstar_vert_r_s']))
f = np.array(R['f'])
for k in ('R0', 'VA_m1', 'UQW1_25_m1'):
    print(k, 'ratio y at f=', ' '.join('%g:%+.2f' % (a, b) for a, b in zip(f[::4], np.array(R[k]['ratio_y'])[::4])))
    print(k, 'anchor y      ', ' '.join('%g:%+.2f' % (a, b) for a, b in zip(f[::4], np.array(R[k]['a_y'])[::4])))
