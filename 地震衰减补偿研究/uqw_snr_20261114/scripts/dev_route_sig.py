"""Descriptive companion of dev_route.py cal: UQW-S on the SR and QScan routes with a smaller declared sigma (0.25, the SR
error model of DESIGN_ZH.md 2.3: SR within about +-25%; and 0.1), development data, calibrated noise.
-> results/dev_route/cal_sigma.json"""
import json
from env import np, semi, score, TW, snr_rows, CW, H, DEV, ROOT, F10, load
import uqws as S
E17 = ROOT.parent / "mu_tune_20261031" / "results" / "dev"; cm = semi.curtain_cm(); Q = np.asarray(cm["Q"], float); rec = []
for did in DEV:
    arr, meta = load(did); D = arr["D"]; r17 = json.loads((E17 / f"{did}.json").read_text(encoding="utf-8")); inv = snr_rows(D, TW, F10); w_t, fs_r = S.target(D, cm)
    for route, c in (("SR", r17["c_sr"]), ("QScan", r17["qscan"])):
        for sig in (0.1, 0.25):
            X = CW.uqw(D, cm, Q, c, sig, TW, 1.4, H.GAMMA, inv_snr=inv)
            rec.append(dict(id=did, route=route, sig=sig, UQW=score(X, arr["ideal"])["S"], UQWS=score(S.shape(X, w_t, inv, fs_r), arr["ideal"])["S"]))
(ROOT / "results" / "dev_route" / "cal_sigma.json").write_text(json.dumps(rec), encoding="utf-8")
cal = {(r["id"], r["route"]): r for r in json.loads((ROOT / "results" / "dev_route" / "cal.json").read_text(encoding="utf-8"))}
for route in ("SR", "QScan"):
    for sig in (0.1, 0.25):
        d = np.array([r["UQWS"] - cal[(r["id"], route)]["SCa"] for r in rec if r["route"] == route and r["sig"] == sig])
        u = np.mean([r["UQWS"] for r in rec if r["route"] == route and r["sig"] == sig])
        print("%-5s sigma %.2f: UQWS mean %.3f, UQWS - SCa %+.3f %d/24" % (route, sig, u, d.mean(), (d > 0).sum()))
