"""UQW-S against SC-a when Q is estimated from the data (development data, exploratory; NOTE_ZH.md section 10).

cal: the 24 development sets with the SR and QScan factors of experiment 17 (mu_tune_20261031/results/dev, archived with
     SR + SC-a and QScan + SC-a at a = 2^-13, mu 0.25). Here: UQW-S, UQ-W and VDBin at c_SR and at c_QScan.
nd:  the 8 development sets of dev_sc.py at 4x noise (res, gau): SR re-estimated on the noisy data (lib17.estimates' SR
     call, nboot 0); UQW-S, UQ-W, VDBin and SC-a (robust wavelets, as dev_sc.py) at c_SR. QScan is not run at 4x noise
     (about 13 SC solves per dataset).
UQW-S keeps the frozen sigma 0.447 (Q known only within x2) on both routes.
usage: dev_route.py cal | dev_route.py nd ID [ID ...]  -> results/dev_route/
"""
import sys, json, time
from env import np, lib, semi, score, sc11, TW, snr_rows, CW, lib12, SEL, H, DEV, ROOT, BASE, F10, noise_of, fd, load
from stage0 import A_TUNED, MU
import lib17
import uqws as S

OUT = ROOT / "results" / "dev_route"; OUT.mkdir(parents=True, exist_ok=True); E17 = BASE / "mu_tune_20261031" / "results" / "dev"
cm = semi.curtain_cm(); Q = np.asarray(cm["Q"], float)


def single(D, c, w_t, fs_r, inv):
    X = CW.uqw(D, cm, Q, c, H.SIG, TW, 1.4, H.GAMMA, inv_snr=inv)
    return dict(UQW=X, UQWS=S.shape(X, w_t, inv, fs_r), VDBin=lib12.tvw(fd.iq_t(D, cm, Q, None, c, 0, TW), w_t, SEL["B3"]["eps"], win=SEL["B3"]["win"]))


if sys.argv[1] == "cal":
    rec = []
    for did in DEV:
        arr, meta = load(did); D = arr["D"]; r17 = json.loads((E17 / f"{did}.json").read_text(encoding="utf-8"))
        inv = snr_rows(D, TW, F10); w_t, fs_r = S.target(D, cm)
        for route, c, sc_key in (("SR", r17["c_sr"], "SR_mu0.25"), ("QScan", r17["qscan"], "QT_mu0.25")):
            row = dict(id=did, route=route, c=c, c_true=meta["c_eff_cum"], SCa=r17["scores"][sc_key]["S"])
            row.update({k: score(Y, arr["ideal"])["S"] for k, Y in single(D, c, w_t, fs_r, inv).items()}); rec.append(row)
        print(did, flush=True)
    (OUT / "cal.json").write_text(json.dumps(rec), encoding="utf-8")
else:
    f = np.fft.rfftfreq(lib17.sq.NFFT, lib17.sq.DT); wins = lib17.sq.windows(); T = lib17.sq.mw_table(cm, Q, f, wins)
    sys.path.insert(0, str(BASE / "fault_q_confirm_20261030" / "scripts")); sys.path.insert(0, str(BASE / "sca_ceiling_20261030" / "scripts"))
    import libc  # noqa: E402
    for did in sys.argv[2:]:
        path = OUT / f"nd_{did}.json"
        if path.exists():
            continue
        t0 = time.monotonic(); arr, meta = load(did); rec = []
        for st, kind in (("res4", "res"), ("gau4", "gau")):
            D = arr["signal"] + noise_of(arr, kind, 4.0, 20261114000 + DEV.index(did))
            c = float(lib17.sq.estimate(D, cm, Q, table=T, wins=wins, nboot=0, seed=meta["rc_seed"])["c"])
            inv = snr_rows(D, TW, F10); w_t, fs_r = S.target(D, cm); p = sc11.prepare_sc(D, cm); wv = libc.wavelets(D, cm, True)
            row = dict(id=did, state=st, route="SR", c=c, c_true=meta["c_eff_cum"])
            row.update({k: score(Y, arr["ideal"])["S"] for k, Y in single(D, c, w_t, fs_r, inv).items()})
            X, info, _ = lib.solve(cm, Q * c, wv, D, p["sigma"], p["K"], p["live"], A_TUNED, MU); row["SCa"] = score(X, arr["ideal"])["S"]; row["SCa_conv"] = bool(info["converged"])
            rec.append(row)
        path.write_text(json.dumps(rec), encoding="utf-8"); print(did, "%.0fs" % (time.monotonic() - t0), flush=True)
