"""UQW-S against the literature method SC (development data, exploratory; NOTE_ZH.md section 9).

semi: SC-a (lib.solve: a = 2^-13, mu_relative 0.25 = experiment 17 / stage 0; noise sigma from sc11.prepare_sc) at 4x noise,
      res and gau (dev9/dev10 noise and seeds), variant O, 5 Q multipliers, with the robust wavelets (libc, as UQW-S's target;
      the frozen wavelets collapse at 4x noise, sca_ceiling mechanism 1). Calibrated noise: the archived stage-0 SC-a results.
fd:   SC (benchmark.solve with the FD dictionary, a 0.125, mu_relative 10: the selection of novelty_evidence_20261025 and
      uqp_selfconsistent pilot_sc.py; the Volve-tuned parameters do not converge on FD) on the noisy gathers of dev10.py (same
      seeds), noise sigma = true noise standard deviation per trace (as pilot_sc.py; favourable to SC). 5 Q multipliers.
usage: dev_sc.py semi ID [ID ...] | dev_sc.py fd K N  -> results/dev_sc/
"""
import sys, json, time
from env import np, lib, semi, score, sc11, TW, MS, DEV, ROOT, BASE, noise_of, fd
from stage0 import A_TUNED, MU
sys.path.insert(0, str(BASE / "fault_q_confirm_20261030" / "scripts")); sys.path.insert(0, str(BASE / "sca_ceiling_20261030" / "scripts"))
import libc  # noqa: E402
from scipy.io import loadmat  # noqa: E402

OUT = ROOT / "results" / "dev_sc"; OUT.mkdir(parents=True, exist_ok=True); bm = lib.bm

if sys.argv[1] == "semi":
    cm = semi.curtain_cm(); Q = np.asarray(cm["Q"], float)
    from env import load
    for did in sys.argv[2:]:
        path = OUT / f"semi_{did}.json"
        if path.exists():
            continue
        t0 = time.monotonic(); arr, meta = load(did); rec = []
        for st, kind in (("res4", "res"), ("gau4", "gau")):
            D = arr["signal"] + noise_of(arr, kind, 4.0, 20261114000 + DEV.index(did)); p = sc11.prepare_sc(D, cm); wv = libc.wavelets(D, cm, True)
            for m in MS:
                X, info, _ = lib.solve(cm, Q * m * meta["c_eff_cum"], wv, D, p["sigma"], p["K"], p["live"], A_TUNED, MU)
                rec.append(dict(state="O" + st, method="SCa", m=m, S=score(X, arr["ideal"])["S"], conv=bool(info["converged"]), s=float(info["seconds"])))
        path.write_text(json.dumps(rec), encoding="utf-8"); print(did, "%.0fs" % (time.monotonic() - t0), flush=True)
else:
    k, n = int(sys.argv[2]), int(sys.argv[3]); jobs = [(h, j, lv, db) for h in fd.MODELS for j, (lv, db) in enumerate((("30dB", 30.0), ("14dB", 14.0)))]
    for i, (h, j, lv, db) in enumerate(jobs):
        path = OUT / f"fd_{h}_{lv}.json"
        if i % n != k or path.exists():
            continue
        t0 = time.monotonic(); D0, T, cm_, tw, fp = fd.load_fd(h); E = fd.noise(D0, db, 20261113900 + 10 * fd.MODELS.index(h) + j); D = D0 + E
        sig = E.std(axis=0); live = D.any(axis=0); K, _, _ = bm.predictor(D, live); rec = []
        cfg = dict(loadmat(fd.FD / "models" / f"{h}.mat", simplify_cells=True)["config"]); Qf = np.asarray(cfg["Q"], float)
        for m in MS:
            c = dict(cfg); c["Q"] = Qf * m; B, W, P, scale = bm.dictionary(c)
            X, info = bm.solve(B, P, D, sig, scale, K, 0.125, 10.0, live)
            rec.append(dict(state=lv, method="SC", m=m, S=float(fd.score_fd(X, T, tw)[0]), conv=bool(info["converged"]), s=float(info["seconds"])))
        path.write_text(json.dumps(rec), encoding="utf-8"); print(h, lv, "%.0fs" % (time.monotonic() - t0), flush=True)
