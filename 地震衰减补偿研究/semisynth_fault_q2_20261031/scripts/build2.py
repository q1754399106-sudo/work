"""Build the datasets of DATA_ZH.md: python build2.py [D E F] (as build.py of semisynth_fault_q_20261030, Z only)."""
import sys, time
import fq2
from fq2 import fq, np, ROOT
import build as b1                                     # gas_shadow, json helpers (unchanged)

DATA = ROOT / "cache" / "data"; ANA = ROOT / "analysis"


def main(keys):
    DATA.mkdir(parents=True, exist_ok=True); ANA.mkdir(exist_ok=True)
    meta, geo_all, real = (b1.load_json(ANA / n) for n in ("datasets.json", "geology.json", "realism.json"))
    for key in keys:
        t0 = time.monotonic(); ctx = fq.prepare(key)
        geo_all[key] = dict(fq.jsonable(ctx["geo"]), lines=ctx["lines"])
        real[key] = dict(inline=ctx["lines"][0], decay_slope_real=fq.decay_slope(ctx["D"]))
        print(f"{key}: prepared lines {ctx['lines']} in {time.monotonic() - t0:.0f}s", flush=True)
        for c0 in fq2.C0:
            for r in (0, 1):
                o = fq.build(ctx, c0, r)["Z"]; m = o["meta"]; m["gas_shadow_30_60Hz"] = b1.gas_shadow(o["arrays"])
                np.savez_compressed(DATA / f"{m['id']}.npz", **{k: v.astype(np.float32) for k, v in o["arrays"].items()})
                meta[m["id"]] = m
                print(f"  {m['id']}: c_eff {m['c_eff_cum']:.3f}/{m['c_eff_dif']:.3f} level {m['level']:.2f} "
                      f"f_true {m['f_true']:.0f} slope {m['decay_slope']:.4f} shadow {m['gas_shadow_30_60Hz']}", flush=True)
        for n, obj in (("datasets.json", meta), ("geology.json", geo_all), ("realism.json", real)):
            b1.save_json(ANA / n, obj)
        print(f"{key}: done in {time.monotonic() - t0:.0f}s", flush=True)


if __name__ == "__main__":
    main(sys.argv[1:] or list(fq2.STRUCT))
