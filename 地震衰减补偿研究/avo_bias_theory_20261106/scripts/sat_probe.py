"""Noise-free probe, M1, events_aligned, fixed lambda: where does the remaining bias of saturated damping come from?
Per event and per band: the amplitude transfer at each offset, measured from the output by windowing around the event."""
import numpy as np
import gate1 as g
from gate1 import C, NX, n, I, DT, T, a0_ax, atten, avo, MODELS, reflectors, synth, SIN2, FS
import step2, theory as th, stable, satdamp as sd
from gate1b import EV
from diag_aligned import bp

Q = MODELS["M1"]; a0, ax = a0_ax(Q); lam = 1e-3
idx, R = reflectors("events_aligned", np.random.default_rng(0)); dp, tru = synth(idx, R, Q, False); S = dp[:n]; keep = np.ones((n, NX), bool)


def out(kind):
    O = np.zeros((n, NX))
    if kind == "split":
        G0 = atten(a0); F0 = np.linalg.solve(G0.T @ G0 + lam * I, G0.T)
        for i in range(NX):
            G = atten(ax[:, i]); O[:, i] = F0 @ np.linalg.solve(G.T @ G + 1e-3 * I, G.T @ S[:, i])
        return O
    fs = np.array([th.fc_star(a0[k], lam) for k in range(n)]); wn = stable.bound(a0, ax, lam, 6.0); fb = sd.fb_of(a0, lam, 0.05)
    for i in range(NX):
        G = atten(a0 + ax[:, i])
        if kind == "method":
            D = np.diag(np.maximum(np.exp(-2 * np.pi * fs * ax[:, i]), wn[:, i]))
        elif kind == "sat":
            L = sd.Lmat(ax[:, i], fb); D = L.T @ L
        elif kind == "satg":
            L = sd.Lmat(ax[:, i], sd.fb_floor(fb, ax[:, i], wn[:, i])); D = L.T @ L
        elif kind == "satnoin":      # stationary-frozen check: L built with a_x frozen at each event handled below
            L = sd.Lmat(ax[:, i], fb); D = L.T @ L
        O[:, i] = np.linalg.solve(G.T @ G + lam * D, G.T @ S[:, i])
    return O


print("a0 at events", np.round(a0[EV], 4), "ax far", np.round(ax[EV, -1], 4), "fb(.05)", np.round(sd.fb_of(a0, lam, 0.05)[EV], 0), "fc*", [round(th.fc_star(a0[k], lam), 1) for k in EV])
print("exp(-2 pi fb ax) floor wn far:", np.round(stable.bound(a0, ax, lam, 6.0)[EV, -1], 4), " exp(-2pi fb ax):", np.round(np.exp(-2 * np.pi * sd.fb_of(a0, lam, .05)[EV] * ax[EV, -1]), 4))
f = np.fft.rfftfreq(2048, DT); win = 0.08
for kind in ("split", "method", "sat", "satg"):
    O = out(kind); print("\n", kind)
    Ob = bp(O, (10, 60)); Tb = bp(tru, (10, 60)); A, B, ok = avo(Ob, keep); At, Bt, _ = avo(Tb, keep); r = ok
    print("  overall bias 10-60:", round(float((np.sum(B * Bt * r) / np.sum(Bt * Bt * r)) / (np.sum(A * At * r) / np.sum(At * At * r)) - 1), 3))
    for j, k0 in enumerate(EV):
        mask = np.abs(T - T[k0]) <= win; Rr = R[j]
        X = np.fft.rfft(O * mask[:, None], 2048, axis=0); Xt = np.fft.rfft(tru * mask[:, None], 2048, axis=0)
        rows = []
        for lo in (10, 20, 30, 40, 50, 60, 70, 80, 90, 100):
            k = (f >= lo) & (f < lo + 10); g_ = np.sum((X[k] * np.conj(Xt[k])).real, axis=0) / np.sum(np.abs(Xt[k]) ** 2, axis=0); rows.append(g_)
        rows = np.array(rows)   # (bands, NX): transfer at each offset; true AVO removed
        print(f"  event {j} R near {Rr[0]:+.3f} far {Rr[-1]:+.3f} sin2far {SIN2[k0, -1]:.3f} ; transfer far/near per 10 Hz from 10:", " ".join(f"{x:.2f}" for x in rows[:, -1] / rows[:, 0]))
        print(f"           transfer near:", " ".join(f"{x:.2f}" for x in rows[:, 0]))
