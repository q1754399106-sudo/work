"""Parameter selection on V1-V3 (PLAN_ZH.md "开发").

  choose.py           summarise the 16-cell grid; if the primary choice lies on an outer edge, write
                      analysis/extension.json (one extension step, at most once) and stop.
  choose.py --final   summarise grid (+ extension), write and freeze analysis/selection.json.
Primary: mean of the 30 rows; only candidates with all 30 rows converged. Ties: larger a, then smaller mu.
Secondary (report only): mean of the 18 rows with Q error in {-10, 0, +10}.
"""
import sys, json, hashlib  # noqa: E401
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(Path(__file__).parent))
from run import GRID, A_GRID, M_GRID, DEV_SCENES, DEV_ERRORS, NOISES, tag  # noqa: E402

AN = ROOT / "analysis"
AN.mkdir(exist_ok=True)
final = "--final" in sys.argv
ext_file = AN / "extension.json"
cands = list(GRID)
if final and ext_file.exists():
    cands += [tuple(x) for x in json.loads(ext_file.read_text())["added"]]


def summarise(cands):
    table = []
    for a, m in cands:
        rows = []
        for s in DEV_SCENES:
            for e in DEV_ERRORS:
                for n in NOISES:
                    f = ROOT / "cache/dev" / f"{s}_{n}_Q{e:+d}_SC_{tag(a, m)}.json"
                    rows.append(json.loads(f.read_text(encoding="utf-8")))
        assert len(rows) == 30
        err = [r["rel_error"] for r in rows]
        mid = [r["rel_error"] for r in rows if r["q_error"] in (-10, 0, 10)]
        q0 = [r["rel_error"] for r in rows if r["q_error"] == 0]
        table.append(dict(a=a, mu_relative=m, mean_error=sum(err) / 30, mean_error_q_mid=sum(mid) / len(mid),
                          mean_error_q0=sum(q0) / len(q0), n_converged=sum(r["converged"] for r in rows),
                          per_scene={s: sum(r["rel_error"] for r in rows if r["scene"] == s) / 10 for s in DEV_SCENES},
                          max_iterations=max(r["iterations"] for r in rows)))
    return table


def pick(table, key):
    ok = [t for t in table if t["n_converged"] == 30]
    if not ok:  # PLAN 失败处理
        top = max(t["n_converged"] for t in table)
        ok = [t for t in table if t["n_converged"] == top]
    return min(ok, key=lambda t: (t[key], -t["a"], t["mu_relative"]))


table = summarise(cands)
prim, sec = pick(table, "mean_error"), pick(table, "mean_error_q_mid")
print(f"{'a':>9} {'mu':>7} {'mean30':>8} {'mid18':>8} {'Q0':>8} conv  maxit")
for t in sorted(table, key=lambda t: (t["a"], t["mu_relative"])):
    mark = " <- primary" if t is prim else ""
    mark += " <- secondary" if t is sec else ""
    print(f"{t['a']:9.6g} {t['mu_relative']:7g} {t['mean_error']:8.5f} {t['mean_error_q_mid']:8.5f} "
          f"{t['mean_error_q0']:8.5f} {t['n_converged']:4d} {t['max_iterations']:6d}{mark}")

if not final:
    assert not ext_file.exists(), "extension already decided"
    a_new = [A_GRID[0] / 2] if prim["a"] == A_GRID[0] else [A_GRID[-1] * 2] if prim["a"] == A_GRID[-1] else []
    m_new = [M_GRID[0] / 10] if prim["mu_relative"] == M_GRID[0] else \
        [M_GRID[-1] * 10] if prim["mu_relative"] == M_GRID[-1] else []
    added = [(a, m) for a in a_new for m in M_GRID] + [(a, m) for a in A_GRID for m in m_new] + \
            [(a, m) for a in a_new for m in m_new]
    ext_file.write_text(json.dumps(dict(trigger=dict(a=prim["a"], mu_relative=prim["mu_relative"]),
                                        added=added), indent=2), encoding="utf-8")
    print("extension:", added if added else "none (interior choice)")
else:
    sel = dict(primary=dict(a=prim["a"], mu_relative=prim["mu_relative"], mean_error=prim["mean_error"]),
               secondary=dict(a=sec["a"], mu_relative=sec["mu_relative"], mean_error_q_mid=sec["mean_error_q_mid"]),
               original=dict(a=0.125, mu_relative=10.0),
               original_on_dev=next(t for t in table if (t["a"], t["mu_relative"]) == (0.125, 10.0)),
               candidates=table, n_candidates=len(table), frozen=True)
    p = AN / "selection.json"
    assert not p.exists(), "selection already frozen"
    p.write_text(json.dumps(sel, indent=2), encoding="utf-8")
    print("frozen:", sel["primary"], "secondary:", sel["secondary"],
          "sha256", hashlib.sha256(p.read_bytes()).hexdigest())
