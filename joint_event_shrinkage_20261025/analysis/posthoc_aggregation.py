"""Read-only aggregation of frozen outputs; no new candidate or selection.

Kept outside scripts/ because this is a post-hoc reporting operation and is
not part of the source set frozen before diagnostic evaluation.
"""
from pathlib import Path
import json,hashlib,sys,csv
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from common import np
def load(name):return json.loads((ROOT/'analysis'/name).read_text(encoding='utf-8'))
def write(name,value):(ROOT/'analysis'/name).write_text(json.dumps(value,indent=2,ensure_ascii=False),encoding='utf-8')
checks=load('frozen_checks.json');numerical=load('numerical_checks.json');timings=load('timings.json');summary=load('summary.json')
stats={};alltr=[]
for key,check in checks.items():
    traces=check['traces'];alltr.extend(traces)
    out=np.load(ROOT/'outputs'/f'{key}.npz')
    stats[key]=dict(traces_with_events=len(traces),median_groups=float(np.median([r['groups'] for r in traces])),fallback_traces=sum(r['fallback'] for r in traces),
        G_offdiag_relative_median=float(np.median([r['offdiag_G_relative'] for r in traces])),G_offdiag_relative_max=max(r['offdiag_G_relative'] for r in traces),
        N_offdiag_relative_median=float(np.median([r['offdiag_N_relative'] for r in traces])),N_offdiag_relative_max=max(r['offdiag_N_relative'] for r in traces),
        Joint_Standard_output_relative_difference=float(np.linalg.norm(out['Joint']-out['Standard'])/np.linalg.norm(out['Standard'])))
write('mechanism_posthoc.json',dict(scope='Post-hoc descriptive aggregation of existing frozen checks only; no selection or new algorithm.',
    relative_measure='Frobenius norm of off-diagonal matrix / Frobenius norm of full matrix per trace; medians are unweighted across nonempty traces.',cases=stats))
numerical_summary=dict(max_E_relative_error=max(t['energy_relative_error'] for t in alltr),max_V_relative_error=max(t['variance_relative_error'] for t in alltr),
    max_ridge_normal_residual=max(t['normal_residual'] for t in alltr),min_G_eigen_relative=min(t['min_eig_G_relative'] for t in alltr),min_N_eigen_relative=min(t['min_eig_N_relative'] for t in alltr),
    max_eta0_Standard_output_difference=max(c['eta0_standard_relative_output_difference'] for c in checks.values()),
    max_full_reconstruction=max(c['full_reconstruction'] for c in checks.values()),max_soft_reconstruction=max(c['soft_reconstruction'] for c in checks.values()),
    max_Standard_archived_platform_difference=max(c['standard_archived_output_difference'] for c in checks.values()),
    max_frozen_candidate_KKT=max(t['max_kkt'] for t in timings if t['method'] in ['Joint','GOnly','NOnly']),
    gradient_relative_error=numerical['max_gradient_relative_error'],orthogonal_negative_control=numerical['orthogonal_negative_control_gain_difference'])
write('numerical_summary.json',numerical_summary)
time_stats={}
for method in ['Joint','GOnly','NOnly','Standard','SC','SharedPreparation']:
    vals=[t['seconds'] for t in timings if t['method']==method]
    time_stats[method]=dict(total_seconds=sum(vals),median_seconds=float(np.median(vals)),max_seconds=max(vals),cases=len(vals))
write('timing_summary.json',dict(scope='Candidate solver excludes shared preparation; SC solver times archived, not rerun; Standard includes coherence and four reference rules.',methods=time_stats))
compact=[dict(case=t['case'],method=t['method'],seconds=t.get('seconds'),preparation_seconds=t.get('preparation_seconds'),max_kkt=t.get('max_kkt'),timing_scope=t.get('timing_scope','candidate box QP only')) for t in timings]
with (ROOT/'analysis/timings.csv').open('w',encoding='utf-8-sig',newline='') as f:
    w=csv.DictWriter(f,fieldnames=['case','method','seconds','preparation_seconds','max_kkt','timing_scope']);w.writeheader();w.writerows(compact)
write('posthoc_provenance.json',dict(code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),input_sha256={name:hashlib.sha256((ROOT/'analysis'/name).read_bytes()).hexdigest() for name in ['frozen_checks.json','numerical_checks.json','timings.json','summary.json']},frozen_selection_sha256=summary['selection_sha256']))
print(json.dumps(dict(numerical=numerical_summary,timing=time_stats,mechanism=stats),ensure_ascii=False,indent=2))
