"""Check the solver against the 2D Green's function (i/4) H0(k r) in a homogeneous constant-Q medium."""
from common import *
from scipy.special import hankel1
c0 = np.full((NX, NZ), 2500.0); qi = np.full((NX, NZ), 1 / 30.0)
for f in (10.0, 25.0, 42.5):
    for kap in (1.0, 0.0, -1.0):
        lu = factor(laplacian(f), c0, qi, f, 1.0, kap)
        b = np.zeros(((NX + 2 * NPML) * (NZ + 2 * NPML), 1), complex); b[gidx(150, 20)] = -1 / H ** 2
        p = interior(lu.solve(b))[..., 0]
        k = 2 * np.pi * f * slowness(2500.0, 1 / 30.0, f, 1.0, kap)
        out = []
        for (ix, iz) in ((150, 100), (150, 170), (60, 120), (290, 160)):
            r = H * np.hypot(ix - 150, iz - 20); g = 0.25j * hankel1(0, k * r)
            out.append("r=%4.0f amp %.3f dph %+.3f" % (r, abs(p[ix, iz]) / abs(g), np.angle(p[ix, iz] / g)))
        print("f=%.1f kappa=%+.0f |" % (f, kap), " | ".join(out), flush=True)
