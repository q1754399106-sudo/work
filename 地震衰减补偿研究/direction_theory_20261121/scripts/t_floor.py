"""DERIVATION_ZH.md section 2.4: smallest correlation difference a Q error can produce under a
correlation metric (constant gain and the best lag absorbed),
    loss ~= 1/2 Var_w[A(f)] + 1/2 E_w[phi_perp(f)^2],
    A = d * pi f t*  (log-amplitude tilt),  phi = d * 2 f t* ln(f/fr), phi_perp = phi minus its best
    weighted linear fit in f (the lag search removes the linear part), d = error of u (relative t*).
Weights w(f): flat in the evaluation band (the field tie uses a box band-pass).
Also checks the small-error formula against the exact correlation of the two spectra.
"""
import json, os
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
cases = [  # name, band (Hz), t* (s) of the compensated part, fr (Hz)
    ('Volve RTM F-5 tie, 4-25 Hz, two-way t* 0.033 s', (4, 25), 0.033, 100),
    ('Volve RTM, same band, fr 25 Hz', (4, 25), 0.033, 25),
    ('post-stack Volve tie, 10-62 Hz, t* below anchor 0.0153 s', (10, 62), 0.0153, 100),
    ('post-stack, hypothetical 10-62 Hz, t* 0.033 s', (10, 62), 0.033, 100),
    ('synthetic RTM deep window, 5-50 Hz, t* 0.016 s', (5, 50), 0.016, 100),
]
out = []
for name, (f1, f2), ts, fr in cases:
    f = np.linspace(f1, f2, 400); w = np.ones_like(f) / f.size
    for d in [np.log(1 / 0.77), np.log(1.3), np.log(2.0)]:
        A = d * np.pi * f * ts
        phi = d * 2 * f * ts * np.log(f / fr)
        X = np.vstack([np.ones_like(f), f]).T
        c, *_ = np.linalg.lstsq(X * np.sqrt(w)[:, None], phi * np.sqrt(w), rcond=None)
        pp = phi - X @ c
        varA = (w * A ** 2).sum() - (w * A).sum() ** 2
        loss_approx = 0.5 * varA + 0.5 * (w * pp ** 2).sum()
        # exact: correlation between spectra 1 and exp(A + i pp), flat reflectivity, best constant gain
        g = np.exp(A + 1j * pp)
        rho = (w * g).sum().real / np.sqrt((w * abs(g) ** 2).sum())
        out.append(dict(case=name, d=float(d), sd_A=float(np.sqrt(varA)), rms_phi_perp_deg=float(np.degrees(np.sqrt((w * pp ** 2).sum()))),
                        loss_approx=float(loss_approx), loss_exact=float(1 - rho)))
        print(f'{name:58s} d {d:.3f}: sd A {np.sqrt(varA):.3f}  phi_perp {np.degrees(np.sqrt((w*pp**2).sum())):5.1f} deg  '
              f'loss approx {loss_approx:.4f} exact {1-rho:.4f}')
    # bulk delay from a wrong dispersion reference: tau = t* ln(fr_used / fr_true) / pi
    print(f'   delay if the true reference is 25 Hz and {fr} Hz is used: {1e3*ts*np.log(fr/25)/np.pi:.1f} ms')
json.dump(out, open(os.path.join(ROOT, 'analysis', 't_floor.json'), 'w', encoding='utf-8'), indent=1)
