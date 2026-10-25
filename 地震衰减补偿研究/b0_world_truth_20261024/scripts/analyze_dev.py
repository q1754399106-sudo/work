"""Apply the declared development rule (DEV_PLAN_ZH.md) to the 44 seen Q0 gathers."""
from pathlib import Path
import json
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
PREV = ROOT / 'cache' / 'inputs' / 'qscan_bias_diagnostic_20261024'
VARIANTS = {'B0': 'B0_rel_delta_mse', 'B0M': 'B0M_rel_delta_mse', 'B0C': 'B0C_rel_delta_mse'}


def rho(x, y):
    x, y = pd.Series(np.asarray(x)).rank(), pd.Series(np.asarray(y)).rank()
    return float(np.corrcoef(x, y)[0, 1])


def main():
    dev = pd.concat([pd.read_csv(f) for f in sorted((ROOT / 'analysis').glob('dev_w*.csv'))])
    man = pd.read_csv(PREV / 'predictions' / 'manifest.csv')
    exp = pd.read_csv(PREV / 'analysis' / 'expected_truth.csv')
    d = man.merge(dev, on='id', validate='one_to_one').merge(
        exp[['set', 'model', 'condition', 'rel_delta', 'decisive']], on=['set', 'model', 'condition'], how='left', validate='many_to_one')
    assert len(d) == 44 and d.rel_delta.notna().all()
    reproduced = bool(np.array_equal(d.B0_rel_delta_mse.to_numpy(), d.sealed_B0.to_numpy()))
    d['truth_worse'] = d.rel_delta > 0
    f1 = (d.set == 'bae') & (d.model == 'F1') & (d.condition == 'Moderate')
    f2 = (d.set == 'bae') & (d.model == 'F2') & (d.condition == 'Moderate')
    dec = d[d.decisive.astype(bool)]
    cond = d.groupby(['set', 'model', 'condition']).agg(rel_delta=('rel_delta', 'first'),
                                                        **{v: (c, 'mean') for v, c in VARIANTS.items()}).reset_index()
    out = dict(B0_reproduces_sealed=reproduced, gathers=len(d), decisive_gathers=len(dec), variants={})
    for v, col in VARIANTS.items():
        r1 = bool((d.loc[f1, col] > 0).all() and (d.loc[f2, col] <= 0).all())
        acc = float(((dec[col] > 0) == dec.truth_worse).mean())
        better, worse = dec.loc[~dec.truth_worse, col], dec.loc[dec.truth_worse, col]
        lo, hi = float(better.max()), float(worse.min())
        sep = lo < hi
        out['variants'][v] = dict(R1_F1M_F2M=r1, decisive_accuracy=acc, spearman_16_conditions=rho(cond[v], cond.rel_delta),
                                  F1M=[round(x, 4) for x in d.loc[f1, col]], F2M=[round(x, 4) for x in d.loc[f2, col]],
                                  max_better=lo, min_worse=hi, separable=bool(sep),
                                  threshold=(lo + hi) / 2 if sep else None,
                                  relative_gap=(hi - lo) / (abs(hi) + abs(lo)) if sep else None)
    passing = [v for v in VARIANTS if v != 'B0' and out['variants'][v]['R1_F1M_F2M']]
    if passing:
        best = max(passing, key=lambda v: (out['variants'][v]['decisive_accuracy'], out['variants'][v]['spearman_16_conditions']))
        out['selected'] = dict(variant=best, threshold=0.0, rule='steps 1-3 (threshold 0)')
    else:
        sep = [v for v in VARIANTS if out['variants'][v]['separable']]
        if sep:
            best = max(sep, key=lambda v: out['variants'][v]['relative_gap'])
            out['selected'] = dict(variant=best, threshold=out['variants'][best]['threshold'], rule='step 5 (calibrated threshold)')
        else:
            out['selected'] = dict(variant=None, threshold=None, rule='stop: no variant passes step 1 and none is separable')
    (ROOT / 'analysis' / 'dev_selection.json').write_text(json.dumps(out, indent=2), encoding='utf-8')
    d.to_csv(ROOT / 'analysis' / 'dev_scored.csv', index=False)
    cond.to_csv(ROOT / 'analysis' / 'dev_conditions.csv', index=False)
    print(json.dumps(out, indent=2))


if __name__ == '__main__':
    main()
