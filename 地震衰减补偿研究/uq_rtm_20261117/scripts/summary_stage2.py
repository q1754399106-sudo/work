"""Stage-2 summary tables (descriptive) from analysis/stage2_*.json -> analysis/stage2_summary.json and stdout."""
import json, glob, os, math
A = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'analysis')
L = lambda v: v if isinstance(v, list) else [v]
g = lambda T, k: (T.get(k) if T.get(k) is not None else float('nan'))
out = {}
print('== E1, W_deep / W_body: M1 M2 M3 (UQW0 vs B*) | (UQW1 vs B+) | M4 min over baselines')
for m in ['D1', 'D2', 'D3', 'D4']:
    R = json.load(open(os.path.join(A, f'stage2_{m}.json'), encoding='utf-8'))
    for ds in ['m30', 'm40', 'g30']:
        for w in ['deep', 'body']:
            T = R[ds][w]['E1']
            row = dict(M1=g(T, 'M1'), M2=g(T, 'M2'), M3=g(T, 'M3'), M1u1=g(T, 'M1_u1'), M2u1=g(T, 'M2_u1'), M3u1=g(T, 'M3_u1'),
                       M4=min(v for v in T['M4'].values() if v is not None),
                       M4u1=min(v for v in T.get('M4_u1', {'x': float('nan')}).values() if v is not None))
            out[f'{m}/{ds}/{w}'] = row
            print(f"{m} {ds} {w:4s}  UQW0 {row['M1']:+.3f} {row['M2']:+.3f} {row['M3']:+.3f} | UQW1 {row['M1u1']:+.3f} {row['M2u1']:+.3f} {row['M3u1']:+.3f} | M4 {row['M4']:+.3f} / {row['M4u1']:+.3f}")
print('\n== other Q-error forms (m30, W_deep): UQW0-B* and UQW1-B+ at listed m')
for m in ['D1', 'D2', 'D3', 'D4']:
    R = json.load(open(os.path.join(A, f'stage2_{m}.json'), encoding='utf-8'))
    for f in ['E2', 'E3', 'E4']:
        T = R['m30']['deep'][f]
        print(f"{m} {f} m={L(T['m'])}  UQW0 {' '.join(f'{x:+.3f}' for x in L(T['d_Bstar']))} | UQW1 {' '.join(f'{x:+.3f}' for x in L(T['d1_Bplus']))}")
print('\n== velocity error (m30, E1, m = 0.7 / 1 / 1.4): UQW0-B* | UQW1-B+')
for m in ['D1', 'D3']:
    for v in ['v096', 'v104']:
        R = json.load(open(os.path.join(A, f'stage2_{m}_{v}.json'), encoding='utf-8'))
        for w in ['deep', 'body']:
            T = R['m30'][w]['E1']
            print(f"{m} {v} {w:4s} UQW0 {' '.join(f'{x:+.3f}' for x in L(T['d_Bstar']))} | UQW1 {' '.join(f'{x:+.3f}' for x in L(T['d1_Bplus']))}")
json.dump(out, open(os.path.join(A, 'stage2_summary.json'), 'w'), indent=1)
