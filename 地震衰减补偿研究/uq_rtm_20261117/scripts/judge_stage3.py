"""Stage-3 verdict (PLAN_STAGE3_ZH.md section 4) from analysis/stage2_C1..C4.json -> analysis/stage3_verdict.json."""
import json, os
A = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'analysis')
units = []
for m in ['C1', 'C2', 'C3', 'C4']:
    R = json.load(open(os.path.join(A, f'stage2_{m}.json'), encoding='utf-8'))
    for ds in ['m30', 'm40']:
        T = R[ds]['deep']['E1']
        m4 = min(v for v in T['M4_u1'].values() if v is not None)
        units.append(dict(unit=f'{m}/{ds}', M1=T['M1_u1'], M2=T['M2_u1'], M3=T['M3_u1'], M4=m4))
M1 = [u['M1'] for u in units]
C1 = sum(M1) / len(M1) >= 0.05 and sum(x >= 0.03 for x in M1) >= 6
C2 = all(u['M2'] >= -0.04 for u in units)
C3 = all(u['M3'] >= -0.06 for u in units)
C4 = sum(u['M4'] > 0 for u in units) >= 7
verdict = ('confirmed' if all([C1, C2, C3, C4]) else
           'narrowed: robust against severe over-compensation only' if C1 else 'not confirmed')
out = dict(units=units, M1_mean=sum(M1) / len(M1), M1_ge_003=sum(x >= 0.03 for x in M1),
           C1=C1, C2=C2, C3=C3, C4=C4, verdict=verdict)
json.dump(out, open(os.path.join(A, 'stage3_verdict.json'), 'w'), indent=1)
for u in units:
    print(f"{u['unit']:8s} M1 {u['M1']:+.3f} M2 {u['M2']:+.3f} M3 {u['M3']:+.3f} M4 {u['M4']:+.3f}")
print(f"M1 mean {out['M1_mean']:+.3f}, units >= 0.03: {out['M1_ge_003']}/8 | C1 {C1} C2 {C2} C3 {C3} C4 {C4} -> {verdict}")
