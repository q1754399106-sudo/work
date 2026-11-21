"""POST HOC, descriptive (written after S2): does pure NMO-stretch kinematics explain phi, without noise or lambda?
Per event and offset, m = dt_x/dt_0 (< 1). After NMO the event spectrum is RIC(f'/m)/m and the Q loss acts as exp(-pi f' a_true / m).
Variants of the noise-free readout weight w(f') (reference At = RIC(f') bp, output bandpassed bp):
  K1  operator compensates a(x) in the t0 domain (as gate A2):  w = RIC(f') bp^2 RIC(f'/m)/m exp[-pi f' a (1/m - 1)]
  K2  operator compensates a(x)/m (stretch-aware):               w = RIC(f') bp^2 RIC(f'/m)/m
  K3  K2 with a stretched reference RIC(f'/m):                    w = bp^2 RIC(f'/m)^2       (pure reparametrisation)
Q sensitivity of ln Gamma: -pi (a/m) fbar_w ; aligned: -pi a fbar(RIC^2 bp^2). phi = 1 - slope_nmo / slope_aligned over kept offsets.
Also SH: K1 read through the amplitude-domain Shuey fit y_x = co_x Gamma_x (exact Zoeppritz co, stretch mute 1.3) -> G/A,
which is what gate A2's per-event G/A does; gives the noise-free stretch bias and Q sensitivity d(G/A)/d ln u in both scenes.
usage: s2_ideal.py -> results/s2_ideal.json"""
import json, sys
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT.parent / "offset_q_avo_20261118" / "scripts" / "legacy2b"))
import gate1 as g1  # noqa: E402

MODELS = {"M24": [42, 70, 55, 95], "M25": [75, 45, 30, 65], "M26": [35, 40, 85, 50], "M27": [95, 60, 45, 35]}
EV = np.round(g1.C.event_t / g1.DT).astype(int); F = np.arange(1, 1001) * 0.1


def ric(f, fp):
    return f ** 2 * np.exp(-(f / fp) ** 2)


def bp(f, fp):
    return np.clip((f - 5) / 5, 0, 1) * np.clip((fp + 8 - f) / 8, 0, 1)


out = []
for mn, Q in MODELS.items():
    a0, ax = g1.a0_ax(np.array(Q, float))
    for fp in (40.0, 60.0):
        for j, k0 in enumerate(EV):
            a = a0[k0] + ax[k0]; s = g1.SIN2[k0]
            m = np.array([np.gradient(g1.GEO[i]["tx"], g1.DT)[k0] for i in range(g1.NX)])
            keep = np.array([g1.GEO[i]["tx"][k0] <= g1.MUTE * g1.T[k0] for i in range(g1.NX)])
            cen = lambda w: (w * F).sum(-1) / w.sum(-1)  # noqa: E731
            sens_al = a * cen(ric(F, fp) ** 2 * bp(F, fp) ** 2)
            R0, B2 = ric(F, fp), bp(F, fp) ** 2
            res = {}
            for name in ("K1", "K2", "K3"):
                w = []
                for i in range(g1.NX):
                    st = ric(F / m[i], fp) / m[i]
                    w.append(R0 * B2 * st * np.exp(-np.pi * F * a[i] * (1 / m[i] - 1)) if name == "K1" else (R0 * B2 * st if name == "K2" else B2 * st ** 2))
                res[name] = (a / m) * cen(np.array(w))
            sk, sc = s[keep], s[keep] - s[keep].mean()
            sl = lambda y: float((y[keep] * sc).sum() / (sc ** 2).sum())  # noqa: E731
            co = g1.rpp(j, np.array([g1.GEO[i]["p"][k0] * g1.C.vp[j] for i in range(g1.NX)])); X = np.c_[np.ones(keep.sum()), s[keep]]

            def ga(eps, nmo):
                y = []
                for i in range(g1.NX):
                    if nmo:
                        w = R0 * B2 * ric(F / m[i], fp) / m[i] * np.exp(-np.pi * F * a[i] * (1 / m[i] - 1)) * np.exp(-np.pi * F * eps * a[i] / m[i])
                    else:
                        w = R0 * B2 * R0 * np.exp(-np.pi * F * eps * a[i])
                    y.append(co[i] * w.sum() / (R0 ** 2 * B2).sum())
                A_, B_ = np.linalg.lstsq(X, np.array(y)[keep], rcond=None)[0]; return B_ / A_
            t_ = np.linalg.lstsq(X, co[keep], rcond=None)[0]; tg = t_[1] / t_[0]
            Sal = (ga(0.01, 0) - ga(-0.01, 0)) / 0.02; Sn = (ga(0.01, 1) - ga(-0.01, 1)) / 0.02
            out.append(dict(model=mn, fp=fp, event=j, n_keep=int(keep.sum()), m_far=float(m[keep][-1]), **{"phi_" + k: 1 - sl(v) / sl(sens_al) for k, v in res.items()},
                            SH_bias_nmo=float(ga(0, 1) - tg), SH_S_al=float(Sal), SH_S_nmo=float(Sn), phi_SH=float(1 - Sn / Sal)))
summ = {k: float(np.median([r[k] for r in out])) for k in ("phi_K1", "phi_K2", "phi_K3", "phi_SH", "SH_S_al", "SH_S_nmo", "m_far")}
(ROOT / "results" / "s2_ideal.json").write_text(json.dumps(dict(rows=out, median=summ), indent=1), encoding="utf-8")
print(json.dumps(summ, indent=1))
for r in out:
    print(r["model"], r["fp"], r["event"], r["n_keep"], round(r["m_far"], 3), round(r["phi_K1"], 2), round(r["phi_K2"], 2), round(r["phi_K3"], 2), round(r["SH_bias_nmo"], 3), round(r["SH_S_al"], 3), round(r["SH_S_nmo"], 3), round(r["phi_SH"], 2))
