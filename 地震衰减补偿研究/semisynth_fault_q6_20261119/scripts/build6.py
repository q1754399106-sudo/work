"""Build the sixth batch (DATA_ZH.md): python build6.py build [P R S T U V]   (as build5.py, Z only, plus signal_TV / signal_ST)
                                    python build6.py check    data-level checks (no compensation): curvature of the TV / K0 / ST
                                                              log spectral ratio 2.3-2.8 s over 1.2-1.6 s (specform.py definition),
                                                              mean over the 24 datasets, calibrated noise
                                    python build6.py seal     SEALED_SHA256.txt (every cache/data/*.npz and analysis/datasets.json)
No compensation and no Q estimation is run here."""
import sys, time, json, hashlib
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import fq6  # noqa: E402
from fq6 import fq, np, ROOT, BASE, TW, DT, lib  # noqa: E402
sys.path.insert(0, str(BASE / "semisynth_fault_q_20261030" / "scripts"))
import build as b1  # noqa: E402  (gas_shadow, json helpers, unchanged)

DATA = ROOT / "cache" / "data"; ANA = ROOT / "analysis"


def tvfilt():
    sys.path.insert(0, str(BASE / "uqw_improve_20261113" / "scripts"))
    import probe2                                   # noqa: E402  (tvfilt unchanged)
    return probe2.tvfilt


def build(keys):
    DATA.mkdir(parents=True, exist_ok=True); ANA.mkdir(exist_ok=True); tv = tvfilt()
    meta, geo_all, real = (b1.load_json(ANA / n) for n in ("datasets.json", "geology.json", "realism.json"))
    for key in keys:
        t0 = time.monotonic(); ctx = fq.prepare(key); b = fq6.bend(key)
        geo_all[key] = dict(fq.jsonable(ctx["geo"]), lines=ctx["lines"], f_knee=fq6.f_knee(key), b=b.tolist(),
                            b_curvature=fq6.curvature(b))
        real[key] = dict(inline=ctx["lines"][0], decay_slope_real=fq.decay_slope(ctx["D"]))
        print(f"{key}: prepared lines {ctx['lines']} f_knee {fq6.f_knee(key):.1f} in {time.monotonic() - t0:.0f}s", flush=True)
        for c0 in fq6.C0:
            for r in fq6.RS:
                o = fq.build(ctx, c0, r)["Z"]; m = o["meta"]; a = o["arrays"]; m["gas_shadow_30_60Hz"] = b1.gas_shadow(a)
                a["signal_TV"] = tv(a["signal"], b); a["signal_ST"] = fq6.st_filter(a["signal"], b)
                np.savez_compressed(DATA / f"{m['id']}.npz", **{k: v.astype(np.float32) for k, v in a.items()})
                meta[m["id"]] = m
                print(f"  {m['id']}: c_eff {m['c_eff_cum']:.3f}/{m['c_eff_dif']:.3f} level {m['level']:.2f} "
                      f"f_true {m['f_true']:.0f} slope {m['decay_slope']:.4f} shadow {m['gas_shadow_30_60Hz']}", flush=True)
        for n, obj in (("datasets.json", meta), ("geology.json", geo_all), ("realism.json", real)):
            b1.save_json(ANA / n, obj)
        print(f"{key}: done in {time.monotonic() - t0:.0f}s", flush=True)


def amp(D, a, c):
    k = (TW >= a) & (TW <= c); live = np.abs(D).sum(0) > 0; F = fq6.F5
    return np.sqrt(lib.smooth(np.mean(np.abs(np.fft.rfft(D[k][:, live] * np.hanning(k.sum())[:, None], fq6.NF, axis=0)) ** 2, axis=1), 5))


def check():
    out = {}
    for var, key in (("K0", "signal"), ("TV", "signal_TV"), ("ST", "signal_ST")):
        r = np.array([np.log(amp(D, 2.3, 2.8) / amp(D, 1.2, 1.6)) for D in
                      (lambda a: a[key] + a["noise"] for a in (fq6.load(i)[0] for i in fq6.IDS))])
        m = r.mean(0); FIT = fq6.FIT; p1 = np.polyfit(fq6.F5[FIT], m[FIT], 1)
        out[var] = dict(curvature=fq6.curvature(m), dtstar_fit=float(-p1[0] / np.pi),
                        curvature_per_set=[fq6.curvature(q) for q in r])
        print(f"{var}: curvature of the mean ratio {out[var]['curvature']:+.3f}, dt* fit {out[var]['dtstar_fit']:.4f} s", flush=True)
    ok = 0.45 <= out["TV"]["curvature"] <= 0.90; out["TV_in_range_0.45_0.90"] = ok
    b1.save_json(ANA / "check_curvature.json", out); print("TV curvature in [0.45, 0.90]:", ok)


def seal():
    files = sorted(DATA.glob("*.npz")) + [ANA / "datasets.json"]
    lines = [f"{hashlib.sha256(p.read_bytes()).hexdigest()}  {p.relative_to(ROOT).as_posix()}" for p in files]
    (ROOT / "SEALED_SHA256.txt").write_text("\n".join(lines) + "\n", encoding="utf-8"); print(len(lines), "files sealed")


if __name__ == "__main__":
    mode = sys.argv[1]
    if mode == "build":
        build(sys.argv[2:] or list(fq6.STRUCT))
    elif mode == "check":
        check()
    elif mode == "seal":
        seal()
