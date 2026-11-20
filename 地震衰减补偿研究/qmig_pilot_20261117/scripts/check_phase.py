"""Code checks: imaging-filter phase against the true reflectivity; sign of the dispersion correction."""
import kirch
from kirch import *
from scipy.ndimage import gaussian_filter1d
d = np.load(CACHE / "data.npz"); Dv, Dr = d["visco"], d["ref"]; T, A = tables(); M = model()
amax = 2 * float(A.max()); xk = np.linspace(0, amax, KBIN); t0 = time.time()
I0 = migrate(Dr, law("NONE", xk), 1.0, xk, T, A, phase=1.0); I90 = migrate(Dr, law("NONE", xk), 1.0, xk, T, A, phase=1j)
print("two migrations %.1f s" % (time.time() - t0))
r = np.gradient(np.log(M["c"][::2]), axis=1); r = gaussian_filter1d(r, 1.5, axis=1)
sx, sz = ZONES["all"]
for k in range(8):
    ph = k * np.pi / 4; I = np.cos(ph) * I0 + np.sin(ph) * I90
    print("phase %5.3f  corr with d ln c / dz: %+.3f" % (ph, np.corrcoef(I[sx, sz].ravel(), r[sx, sz].ravel())[0, 1]))
R = I0
for sgn in (+1, -1):
    kirch.KF = sgn * 2 * FREQ * np.log(FREQ / F0)
    for name in ("NONE", "PH", "WIEN"):
        print("kf sign %+d %-5s" % (sgn, name), score(migrate(Dv, law(name, xk), 1.0, xk, T, A, phase=1.0), R))
