"""Build the datasets of DATA_ZH.md: python build.py [A B C]. Arrays -> cache/data/<id>.npz (float32),
per-dataset metadata -> analysis/datasets.json, geology -> analysis/geology.json, realism -> analysis/realism.json."""
import sys, json, time
import fq
from fq import np, semi

ROOT = fq.semi.ROOT.parent / "semisynth_fault_q_20261030"
DATA = ROOT / "cache" / "data"; ANA = ROOT / "analysis"


def load_json(p):
    return json.loads(p.read_text(encoding="utf-8")) if p.exists() else {}


def save_json(p, obj):
    p.write_text(json.dumps(obj, indent=1, ensure_ascii=False), encoding="utf-8")


def gas_shadow(arr):
    """Ratio of the 30-60 Hz rms of the signal below the gas lens (0.1-0.4 s under its base) to the same
    times on traces outside the lens; descriptive."""
    g = arr["gas"]; x = semi.bandpass(arr["signal"], 30, 60)
    inside = g.max(axis=0) > 0.5
    if not inside.any():
        return None
    base = int(np.max(np.where(g[:, inside].max(axis=1) > 0.5)[0]))
    rows = slice(base + 25, min(base + 100, len(fq.TW)))
    return float(np.sqrt(np.mean(x[rows][:, inside] ** 2) / np.mean(x[rows][:, ~inside] ** 2)))


def main(keys):
    DATA.mkdir(parents=True, exist_ok=True); ANA.mkdir(exist_ok=True)
    meta, geo_all, real = load_json(ANA / "datasets.json"), load_json(ANA / "geology.json"), load_json(ANA / "realism.json")
    for key in keys:
        t0 = time.monotonic(); ctx = fq.prepare(key)
        geo_all[key] = dict(fq.jsonable(ctx["geo"]), lines=ctx["lines"])
        real[key] = dict(inline=ctx["lines"][0], decay_slope_real=fq.decay_slope(ctx["D"]))
        print(f"{key}: prepared lines {ctx['lines']} in {time.monotonic() - t0:.0f}s", flush=True)
        for c0 in fq.C0:
            for r in (0, 1):
                out = fq.build(ctx, c0, r)
                for var, o in out.items():
                    m = o["meta"]; m["gas_shadow_30_60Hz"] = gas_shadow(o["arrays"])
                    np.savez_compressed(DATA / f"{m['id']}.npz", **{k: v.astype(np.float32) for k, v in o["arrays"].items()})
                    meta[m["id"]] = m
                    print(f"  {m['id']}: c_eff {m['c_eff_cum']:.3f}/{m['c_eff_dif']:.3f} level {m['level']:.2f} "
                          f"f_true {m['f_true']:.0f} slope {m['decay_slope']:.4f} shadow {m['gas_shadow_30_60Hz']}", flush=True)
        save_json(ANA / "datasets.json", meta); save_json(ANA / "geology.json", geo_all); save_json(ANA / "realism.json", real)
        print(f"{key}: done in {time.monotonic() - t0:.0f}s", flush=True)


if __name__ == "__main__":
    main(sys.argv[1:] or list(fq.STRUCT))
