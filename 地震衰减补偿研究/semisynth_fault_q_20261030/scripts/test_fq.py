"""Checks of fq.py: (1) synth with a laterally constant t* equals semi.model (experiment-8 modelling);
(2) a t* half-way between grid nodes; (3) the fault mapping gives missing section, never repeated section."""
import time
import fq
from fq import np, semi, lib

cm = semi.curtain_cm(); f = np.fft.rfftfreq(8192, fq.FINE); fref = float(cm["fref"])
rng = np.random.default_rng(1); nx = 7
R = np.zeros((len(fq.TF), nx)); idx = rng.choice(len(fq.TF), 60, replace=False); R[idx] = rng.normal(0, 0.1, (60, nx))
S = np.exp(-((f - 35) / 25) ** 2) + 0j
ts = semi.tstar_true(cm, 0.9)
t0 = time.time(); A = semi.model(R, ts, S, f, fref); t1 = time.time()
B = fq.synth(R, np.repeat(ts[:, None], nx, axis=1), S, f, fref); t2 = time.time()
print("constant t*: rel err %.2e  (semi %.1fs, fq %.1fs)" % (np.linalg.norm(A - B) / np.linalg.norm(A), t1 - t0, t2 - t1))
A0 = semi.model(R, None, S, f, fref); B0 = fq.synth(R, None, S, f, fref)
print("unattenuated: rel err %.2e" % (np.linalg.norm(A0 - B0) / np.linalg.norm(A0)))
ts2 = np.full(len(fq.TF), 0.03005)                                     # exactly between nodes 150 and 151
A2 = semi.model(R, ts2, S, f, fref); B2 = fq.synth(R, np.repeat(ts2[:, None], nx, axis=1), S, f, fref)
print("mid-node t*: rel err %.2e" % (np.linalg.norm(A2 - B2) / np.linalg.norm(A2)))
geo = fq.draw_geology("A", np.random.default_rng(fq.GEO_SEED["A"]))
tau = np.zeros((len(fq.TW), fq.NX)); tg, code, dist = fq.mapping(geo, tau)
d = np.diff(tg, axis=0); dc = np.diff(code, axis=0) != 0
print("tg steps: min %.4f ms, max jump %.1f ms, jumps %d" % (d.min() * 1e3, d.max() * 1e3, int(dc.sum())))
