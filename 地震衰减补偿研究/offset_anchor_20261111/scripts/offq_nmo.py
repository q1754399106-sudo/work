"""Offset spectral ratio on FLATTENED gathers: reverse the moveout correction first (exploratory synthetic pilot).

Field CIP gathers (Volve ST10010) are delivered moveout-corrected, so the stretch is already in the data. Test: the
recorded synthetic gather of offq.py is NMO-corrected with the true vrms (stand-in for the migration), then the moveout is
restored with the processing velocity vrms * vscale, and offq.estimates runs as before (windows and predicted extra path
also with vscale). Interpolation: 8 x Fourier oversampling, then linear. "direct" = offq.py without the round trip.

usage: offq_nmo.py run
"""
import sys, json
import numpy as np
import offq
from offq import N, DT, T, OFF, vrms, OUT, SEEDS

UP = 8; TU = np.arange(N * UP) * DT / UP


def resample(g, tq):
    """g (n_t, n_cmp) read at times tq (n_t), zero outside the record."""
    G = np.fft.rfft(g, axis=0); Z = np.zeros((N * UP // 2 + 1, g.shape[1]), complex); Z[:G.shape[0]] = G; u = np.fft.irfft(Z, N * UP, axis=0) * UP
    out = np.stack([np.interp(tq, TU, u[:, k], left=0.0, right=0.0) for k in range(g.shape[1])], axis=1); out[tq > T[-1]] = 0.0
    return out


def round_trip(vscale):
    def post(g, i):
        if OFF[i] == 0:
            return g
        flat = resample(g, np.sqrt(T ** 2 + (OFF[i] / vrms) ** 2)); tx = np.sqrt(T ** 2 + (OFF[i] / (vrms * vscale)) ** 2)
        return resample(flat, np.interp(T, tx, T, left=0.0))
    return post


CASES = [dict(tag="noise-free white", snr=0, beta=0.0, v=1.0), dict(tag="blue", snr=30, beta=0.7, v=1.0), dict(tag="blue", snr=18, beta=0.7, v=1.0),
         dict(tag="blue, velocity +2%", snr=18, beta=0.7, v=1.02), dict(tag="blue, velocity -2%", snr=18, beta=0.7, v=0.98),
         dict(tag="blue, velocity +5%", snr=18, beta=0.7, v=1.05), dict(tag="blue, velocity -5%", snr=18, beta=0.7, v=0.95)]

if __name__ == "__main__":
    rows = []
    for case in CASES:
        for c in (0.7, 1.0, 1.4):
            lo = {"direct": [], "trip": []}
            for seed in SEEDS:
                for name, post in (("direct", None), ("trip", round_trip(case["v"]))):
                    rng = np.random.default_rng(1000 * seed + 7); r = offq.reflectivity(rng, 0.0, case["beta"])
                    co, _, n = offq.estimates(r, c, case["snr"], case["v"], rng, post); rows.append(dict(case, c=c, seed=seed, method=name, c_off=co, n=n))
                    lo[name].append(np.log(co / c) if co > 0 else np.nan)
            print("%-20s snr %2d c %.1f | direct %+.3f sd %.3f | round trip %+.3f sd %.3f (bad %d)" % (case["tag"], case["snr"], c, np.nanmean(lo["direct"]), np.nanstd(lo["direct"]),
                  np.nanmean(lo["trip"]), np.nanstd(lo["trip"]), int(np.isnan(lo["trip"]).sum())), flush=True)
    (OUT / "offq_nmo_rows.json").write_text(json.dumps(rows), encoding="utf-8")
