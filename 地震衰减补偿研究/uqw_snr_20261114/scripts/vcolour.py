"""Variance of the log-amplitude departure r(f) of uqwp.update that comes from reflectivity colour alone (no noise, no Q):
r computed on the ideal sections (stationary wavelet, no attenuation: dx_hat = 0, snr infinite) of the 24 development sets,
windows centred 1.9-2.9 s, anchor 1.2-1.6 s. Stand-in for a statistic that would come from well-log reflectivity.
-> results/vcolour.json {f, vc}: variance over windows and datasets, smoothed over 21 bins."""
import json
from env import np, lib, TW, load, DEV, ROOT
import uqwp as U
rs = []
for did in DEV:
    I = load(did)[0]["ideal"]; z = np.zeros((len(TW), 513)); A0, w0 = U.anchor_semi(I, TW)
    _, log = U.update(I, I, TW, np.zeros(len(TW)), z, A0, np.ones_like(w0), 0.447, mode="const", ret=True)
    rs += [r for t, r, lm in log if 1.9 <= t <= 2.9]
f = np.fft.rfftfreq(U.NF, U.DT); vc = lib.smooth(np.var(np.array(rs), axis=0), 21)
(ROOT / "results" / "vcolour.json").write_text(json.dumps(dict(f=f.tolist(), vc=vc.tolist())), encoding="utf-8")
print("windows", len(rs), "sd of r at 10..70 step 10 Hz:", " ".join("%.2f" % np.sqrt(vc[int(round(q / f[1]))]) for q in range(10, 71, 10)))
