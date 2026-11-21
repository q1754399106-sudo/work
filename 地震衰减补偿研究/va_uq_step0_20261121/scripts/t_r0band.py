"""How much can ANY compensation raise a box-band well tie whose band is the uncompensated image's own
-10 dB band?  (stage-4 metric, fld_eval.m: well reflectivity and image both box-filtered to R0's -10 dB band)
Ideal case: white well reflectivity, perfect coherence, no noise; the image spectrum a0(k) is the only
imperfection. rho0 = <a0>/sqrt(<a0^2>) over the band; flattening a0 gives 1. Upper bound on the gain:
R0 * (1/rho0 - 1) if the other losses are frequency independent.
Shapes: Ricker-type amplitude x^2 exp(1-x^2) and Gaussian, both with -10 dB points at the band edges.
"""
import numpy as np
from scipy.optimize import brentq
k1, k2, R0 = 3.515625, 12.890625, 0.626
def ricker(x): return x**2*np.exp(1-x**2)
lo = brentq(lambda x: ricker(x)-10**(-0.5), 0.05, 1.0); hi = brentq(lambda x: ricker(x)-10**(-0.5), 1.0, 4.0)
print(f'Ricker -10 dB roots x = {lo:.3f}, {hi:.3f}, ratio {hi/lo:.2f} (band ratio {k2/k1:.2f})')
k = np.linspace(k1, k2, 4001)
out = {}
# Ricker shape: choose kp so that -10 dB points are as close as possible to band (ratio differs slightly)
kp = np.sqrt(k1*k2/(lo*hi))
a = ricker(k/kp)
rho = a.mean()/np.sqrt((a**2).mean()); out['ricker'] = rho
# Gaussian in k with -10 dB at both edges
kc = 0.5*(k1+k2); w = (k2-k1)/2/np.sqrt(np.log(10**0.5)); g = np.exp(-((k-kc)/w)**2)
rho_g = g.mean()/np.sqrt((g**2).mean()); out['gauss'] = rho_g
# one-sided: monotonically falling from 0 dB at k1 to -10 dB at k2 (worst attenuation-like case)
e = 10**(-0.5*(k-k1)/(k2-k1)); rho_e = e.mean()/np.sqrt((e**2).mean()); out['falling'] = rho_e
for n, r in out.items():
    print(f'{n:8s}: rho0 = {r:.3f}; max relative gain {1/r-1:+.3f}; on R0 = {R0}: +{R0*(1/r-1):.3f}')
