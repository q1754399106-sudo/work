"""Gate 1b (GATE1B_ZH.md): first step correcting offset absorption and NMO stretch together.
usage: gate1b.py smoke | run M | stress   -> results/gate1b_<M>.json.gz, results/gate1b_stress.json"""
import gzip, json, sys, time
from pathlib import Path
import numpy as np
import gate1 as g
from gate1 import C, GEO, SIN2, NX, n, T, I, DT, NSYN, FS, RIC, MODELS, LAMS, LX, QF, a0_ax, atten, ricker, avo, bandpass, solve_all, build_data, MUTE, TW

OUT = Path(__file__).resolve().parents[1] / "results"; EV = np.round(C.event_t / DT).astype(int)
B_ = -np.pi * FS + 2j * FS * np.log(np.maximum(FS, 1e-300) / C.fref); B_[0] = 0; M4 = 4; Nu = M4 * NSYN; TU = np.arange(Nu) * DT / M4


def txs(i, vs=1.0):
    return np.sqrt(T ** 2 + (GEO[i]["tx"] ** 2 - T ** 2) / vs ** 2)


def nmo_cols(d, tx, keep):
    X = np.fft.rfft(d, axis=0); Xu = np.zeros((Nu // 2 + 1, d.shape[1]), complex); Xu[:X.shape[0]] = X; du = np.fft.irfft(Xu, Nu, axis=0) * M4
    pos = np.clip(tx / (DT / M4), 0, Nu - 2); i0 = pos.astype(int); fr = pos - i0
    return (du[i0, :] * (1 - fr)[:, None] + du[i0 + 1, :] * fr[:, None]) * keep[:, None]


def ricker_spec(fp):
    return np.abs(np.fft.rfft(np.roll(ricker(np.arange(NSYN) * DT - NSYN * DT / 2, fp), -NSYN // 2)))


def op_columns(a, tk, ric):
    return np.fft.irfft((ric[:, None] * np.exp(np.outer(B_, a))) * np.exp(-2j * np.pi * np.outer(FS, tk)), NSYN, axis=0)


def operators(Q, vs=1.0, fp=C.fpeak):
    """F_s, T_s per offset, F_0, W, normalised by the largest singular value of F_0 (returned)."""
    ric = ricker_spec(fp); a0, ax = a0_ax(Q); F0 = op_columns(a0, T, ric)[:n]; W = op_columns(0 * a0, T, ric)[:n]; sig = np.linalg.norm(F0, 2); Fs, Ts, keeps = [], [], []
    for i in range(NX):
        tx = txs(i, vs); keep = tx <= MUTE * np.maximum(T, DT); keeps.append(keep)
        Fs.append(nmo_cols(op_columns(a0 + ax[:, i], tx, ric), tx, keep) / sig); Ts.append(nmo_cols(op_columns(0 * a0, tx, ric), tx, keep) / sig)
    return dict(Fs=Fs, Ts=Ts, F0=F0, W=W, sig=sig, a0=a0, ax=ax, keep=np.stack(keeps, axis=1))


def ridge_apply(F, lam, Y):
    return np.linalg.solve(F.T @ F + lam * I, F.T @ Y)


def template(D, keep):
    amp = np.zeros((3, NX))
    for j, te in enumerate(C.event_t):
        k = np.abs(T - te) <= 0.07; w = ricker(T[k] - te); amp[j] = (w @ D[k]) / (w @ w)
    return np.array([np.polyfit(SIN2[EV[j]][keep[EV[j]]], amp[j][keep[EV[j]]], 1)[::-1] for j in range(3)])


def score(O, Tru, keep, At, Bt, ok, events):
    mask = (keep & TW[:, None])[:, :, None]; wave = np.sqrt(np.sum(((O - Tru) * mask) ** 2, axis=(0, 1)) / np.sum((Tru * mask) ** 2, axis=(0, 1)))
    A, B, _ = avo(bandpass(O), keep); r = ok[:, None]
    ba = np.sum(A * At * r, axis=0) / np.sum(At ** 2 * r, axis=0); bg = np.sum(B * Bt * r, axis=0) / np.sum(Bt ** 2 * r, axis=0)
    cg = np.sum(B * Bt * r, axis=0) / np.sqrt(np.sum(B ** 2 * r, axis=0) * np.sum(Bt ** 2 * r, axis=0))
    eg = np.sqrt(np.sum(((B / ba - Bt) * r) ** 2, axis=0) / np.sum((Bt * r) ** 2, axis=0))
    out = dict(wave=wave, eg=eg, mg=np.abs(bg / ba - 1), ma=np.abs(ba - 1), mb=np.abs(bg - 1), cg=cg, ba=ba, bg=bg)
    if events:
        t = []
        for c in range(O.shape[2]):
            fT = template(Tru[:, :, c], keep); fE = template(O[:, :, c], keep); t.append(np.mean(np.abs(fE[:, 1] / fE[:, 0] - fT[:, 1] / fT[:, 0])))
        out["tga"] = np.array(t)
    return out


def run_cell(label, op, Q_a, scene, cols, rows, extra=None, plain_only=False):
    a0, ax = op["a0"], op["ax"]; G0 = atten(a0); Gx = [atten(ax[:, i]) for i in range(NX)]; Gd = [atten(a0 + ax[:, i]) for i in range(NX)]; keep = cols[0][4]; nc = len(cols); ev = scene.startswith("events")
    Dat = np.stack([c[2] for c in cols], axis=2); Tru = np.stack([c[3] for c in cols], axis=2); At, Bt, ok = avo(bandpass(Tru), keep); sig = op["sig"]; meta = [(c[0], c[1]) for c in cols]

    def rec(m, lam, O):
        s = score(O, Tru, keep, At, Bt, ok, ev)
        for c in range(nc):
            rows.append(dict(label, scene=scene, snr=int(meta[c][0]), seed=int(meta[c][1]), m=m, lam=float(lam), **{k: float(v[c]) for k, v in s.items()}))

    Bx = [np.linalg.solve(G.T @ G + LX * I, G.T) for G in Gx]; Y1 = np.stack([Bx[i] @ Dat[:, i, :] for i in range(NX)], axis=1)
    if not plain_only:
        R1 = [ridge_apply(op["Fs"][i], LX, Dat[:, i, :] / sig) for i in range(NX)]; Y2 = np.stack([op["F0"] @ R1[i] for i in range(NX)], axis=1)
        D3 = np.stack([op["W"] @ ridge_apply(op["Ts"][i], LX, Dat[:, i, :] / sig) for i in range(NX)], axis=1); Y3 = np.stack([Bx[i] @ D3[:, i, :] for i in range(NX)], axis=1)
    for lam in LAMS:
        F = np.linalg.solve(G0.T @ G0 + lam * I, G0.T); O = np.zeros_like(Dat)
        for i in range(NX):
            O[:, i, :] = F @ Y1[:, i, :]
        rec("split", lam, O)
        for i in range(NX):
            O[:, i, :] = solve_all(Gd[i], np.ones(n), lam, Dat[:, i, :])
        rec("direct", lam, O)
        if plain_only:
            continue
        for i in range(NX):
            O[:, i, :] = F @ Y2[:, i, :]
        rec("split_st", lam, O)
        for i in range(NX):
            O[:, i, :] = op["W"] @ ridge_apply(op["Fs"][i], lam, Dat[:, i, :] / sig)
        rec("direct_st", lam, O)
        for i in range(NX):
            O[:, i, :] = F @ Y3[:, i, :]
        rec("dest_split", lam, O)


def smoke():
    t0 = time.monotonic(); rows = []; op = operators(MODELS["M1"]); print("operators", round(time.monotonic() - t0, 1), "s", flush=True)
    cols = build_data("M1", "dense_nmo")[:4]; run_cell(dict(M="M1", qf=1.0), op, None, "dense_nmo", cols, rows)
    print("rows", len(rows), "finite", bool(np.all(np.isfinite([r["eg"] for r in rows]))), "elapsed", round(time.monotonic() - t0, 1), "s; no scores printed", flush=True)
    # operator check against the forward model on a noiseless single reflector: relative difference only
    idx = np.array([EV[2]]); R = np.ones((1, NX)); dp, _ = g.synth(idx, R, MODELS["M1"], True); y, keep = g.to_t0(dp, True); d = np.stack([op["Fs"][i][:, idx[0]] * op["sig"] for i in range(NX)], axis=1)
    print("operator vs forward, single reflector, relative difference", float(np.linalg.norm((y - d) * keep) / np.linalg.norm(y * keep)), flush=True)


def run(Mname):
    rows = []; path = OUT / f"gate1b_{Mname}.json.gz"; Q0 = MODELS[Mname]
    for scene in ("events_nmo", "thin_nmo", "dense_nmo"):
        cols = build_data(Mname, scene)
        for qf in QF:
            t0 = time.monotonic(); run_cell(dict(M=Mname, qf=qf), operators(Q0 * qf), None, scene, cols, rows); print(Mname, scene, qf, f"{time.monotonic() - t0:.0f}s", flush=True)
            gzip.open(path, "wt", encoding="utf-8").write(json.dumps(rows))
    for scene in ("events_aligned", "dense_aligned"):
        cols = build_data(Mname, scene); run_cell(dict(M=Mname, qf=1.0), operators(Q0), None, scene, cols, rows, plain_only=True); print(Mname, scene, flush=True)
        gzip.open(path, "wt", encoding="utf-8").write(json.dumps(rows))


def stress():
    rows = []; Q0 = MODELS["M1"]
    for scene in ("events_nmo", "dense_nmo"):
        cols = [c for c in build_data("M1", scene) if c[0] == 18]
        for name, kw in (("none", {}), ("v0.98", dict(vs=0.98)), ("v1.02", dict(vs=1.02)), ("fp0.9", dict(fp=0.9 * C.fpeak)), ("fp1.1", dict(fp=1.1 * C.fpeak))):
            op = operators(Q0, **kw); c2 = cols
            if "vs" in kw:      # data NMO-corrected with the same wrong velocity: rebuild from the forward model
                c2 = []
                for snr, seed, y, tru, keep, sg in cols:
                    rng = np.random.default_rng(100 * ["events_aligned", "thin_aligned", "dense_aligned", "events_nmo", "thin_nmo", "dense_nmo"].index(scene) + seed); idx, R = g.reflectors(scene, rng); dp, tr = g.synth(idx, R, Q0, True)
                    d = dp + sg[:, None].T * rng.standard_normal(dp.shape); out = np.stack([nmo_cols(d[:, [i]], txs(i, kw["vs"]), op["keep"][:, i])[:, 0] for i in range(NX)], axis=1); c2.append((snr, seed, out, tr, op["keep"], sg))
            run_cell(dict(M="M1", qf=1.0, stress=name), op, None, scene, c2, rows); print(scene, name, flush=True)
    (OUT / "gate1b_stress.json").write_text(json.dumps(rows), encoding="utf-8")


if __name__ == "__main__":
    m = sys.argv[1]
    if m == "smoke":
        smoke()
    elif m == "stress":
        stress()
    else:
        run(sys.argv[2])
