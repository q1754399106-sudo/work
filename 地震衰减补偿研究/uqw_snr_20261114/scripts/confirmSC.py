"""Confirmation of the UQW-S / SC division rule (PLAN_CONFIRM_SC_ZH.md) on the UQW-S confirmation data (second use):
semi-synthetic batch 6 (states cal, res4, gau4 with the noise and seeds of confirmS.py) and FD H13-H18 (30 and 14 dB, the
three realisations and seeds of confirmS.py).
semi: SC-a (lib.solve, a 2^-13, mu_relative 0.25, robust wavelets libc.wavelets(robust=True), noise sigma of sc11.prepare_sc)
      and UQW-S (uqws.shape after uqw, frozen) at m in M5 (c = m c_eff_cum) and at the SR factor estimated on the noisy data
      (lib17.sq.estimate, nboot 0, the dataset's rc_seed). UQW-S at M5 is checked against results/confirmS (same data).
fd:   SC (benchmark.solve with the FD dictionary, a 0.125, mu_relative 10, noise sigma = true noise std per trace) at M5;
      UQW-S read from results/confirmS.
usage: confirmSC.py smoke        one used development set (cal) and FD H1 30 dB, one m: code path and timing, no scores
       confirmSC.py run K N      worker K of N over 36 semi (dataset x state) + 36 FD jobs -> results/confirmSC/ (resumable)
"""
import sys, json, time
from env import np, lib, semi, score, sc11, TW, snr_rows, CW, H, ROOT, BASE, F10, fd
from stage0 import A_TUNED, MU
import lib17
import uqws as S
import confirmS as CS
sys.path.insert(0, str(BASE / "fault_q_confirm_20261030" / "scripts")); sys.path.insert(0, str(BASE / "sca_ceiling_20261030" / "scripts"))
import libc  # noqa: E402
from scipy.io import loadmat  # noqa: E402

M5 = (0.5, 0.7, 1.0, 1.4, 2.0); OUT = ROOT / "results" / "confirmSC"; bm = lib.bm
cm = semi.curtain_cm(); Q = np.asarray(cm["Q"], float)


def semi_unit(arr, meta, D, ms=M5, sr=True):
    f = np.fft.rfftfreq(lib17.sq.NFFT, lib17.sq.DT); wins = lib17.sq.windows(); T = lib17.sq.mw_table(cm, Q, f, wins)
    inv = snr_rows(D, TW, F10); w_t, fs_r = S.target(D, cm); p = sc11.prepare_sc(D, cm); wv = libc.wavelets(D, cm, True); rec = []
    cs = [(f"m{m:g}", m * meta["c_eff_cum"]) for m in ms]
    if sr:
        c_sr = float(lib17.sq.estimate(D, cm, Q, table=T, wins=wins, nboot=0, seed=meta["rc_seed"])["c"]); cs.append(("SR", c_sr))
    for tag, c in cs:
        U = S.shape(CW.uqw(D, cm, Q, c, H.SIG, TW, 1.4, H.GAMMA, inv_snr=inv), w_t, inv, fs_r)
        X, info, _ = lib.solve(cm, Q * c, wv, D, p["sigma"], p["K"], p["live"], A_TUNED, MU)
        rec.append(dict(point=tag, c=c, UQWS=score(U, arr["ideal"])["S"], SCa=score(X, arr["ideal"])["S"], SCa_conv=bool(info["converged"]), SCa_s=float(info["seconds"])))
    return rec


def fd_unit(D0, T, cfg_path, tw, db, seed, ms=M5):
    E = fd.noise(D0, db, seed); D = D0 + E; sig = E.std(axis=0); live = D.any(axis=0); K, _, _ = bm.predictor(D, live)
    cfg = dict(loadmat(cfg_path, simplify_cells=True)["config"]); Qf = np.asarray(cfg["Q"], float); rec = []
    for m in ms:
        c = dict(cfg); c["Q"] = Qf * m; B, W, P, scale = bm.dictionary(c); X, info = bm.solve(B, P, D, sig, scale, K, 0.125, 10.0, live)
        rec.append(dict(point=f"m{m:g}", SC=float(fd.score_fd(X, T, tw)[0]), SC_conv=bool(info["converged"]), SC_s=float(info["seconds"])))
    return rec


if __name__ == "__main__" and sys.argv[1] == "smoke":
    from env import load
    t0 = time.monotonic(); arr, meta = load("ZA_c0.75_r0"); r = semi_unit(arr, meta, arr["D"], ms=(1.0,))
    print("semi smoke: points %s, finite %s, SC-a seconds %s, %.0fs" % ([x["point"] for x in r], all(np.isfinite([x["UQWS"], x["SCa"]]).all() for x in r), [round(x["SCa_s"]) for x in r], time.monotonic() - t0))
    t0 = time.monotonic(); D0, T, cmf, tw, fp = fd.load_fd("H1"); r = fd_unit(D0, T, fd.FD / "models" / "H1.mat", tw, 30.0, 1, ms=(1.0,))
    print("FD smoke: %d record, finite %s, %.0fs" % (len(r), bool(np.isfinite(r[0]["SC"])), time.monotonic() - t0))

if __name__ == "__main__" and sys.argv[1] == "run":
    sys.path.insert(0, str(CS.D6 / "scripts")); import fq6  # noqa: E402
    k, n = int(sys.argv[2]), int(sys.argv[3]); OUT.mkdir(parents=True, exist_ok=True)
    jobs = [("semi", (d, st)) for d in CS.IDS6 for st, _, _ in CS.STATES] + [("fd", (h, lv, r)) for h in CS.MODELS3 for lv in CS.LEVELS3 for r in CS.R3]
    for j, (kind, key) in enumerate(jobs):
        name = "S6_%s_%s" % key if kind == "semi" else "FD3_%s_%s_r%d" % key; path = OUT / f"{name}.json"
        if j % n != k or path.exists():
            continue
        t0 = time.monotonic()
        if kind == "semi":
            did, st = key; arr, meta = fq6.load(did); arr = {a: v.astype(float) for a, v in arr.items()}
            _, kind_n, x = [s for s in CS.STATES if s[0] == st][0]; D = arr["signal"] + CS.noise6(arr, kind_n, x, CS.IDS6.index(did))
            out = dict(id=name, records=semi_unit(arr, meta, D))
        else:
            h, lv, r = key; D0, T, cmf, tw, fp = CS.load_fd3(h); s = CS.fd_seed(h, lv, r)
            out = dict(id=name, seed=s, records=fd_unit(D0, T, CS.FD3 / "models" / f"{h}.mat", tw, CS.LEVELS3[lv], s))
        out["seconds"] = time.monotonic() - t0; path.write_text(json.dumps(out), encoding="utf-8"); print(name, "%.0fs" % out["seconds"], flush=True)
