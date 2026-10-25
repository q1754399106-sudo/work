"""Recompute submitted metrics from saved waveforms, not reported summaries."""
from pathlib import Path
import sys, json, csv, hashlib
ROOT=Path(__file__).resolve().parent
TASK=ROOT.parent/'joint_event_shrinkage_20261025'
sys.path.insert(0,str(TASK/'scripts'))
from common import load_scene
import numpy as np

def digest(p): return hashlib.sha256(p.read_bytes()).hexdigest()
selection=json.loads((TASK/'analysis/selection.json').read_text())
assert digest(TASK/'PLAN_ZH.md')==selection['plan_sha256']
assert {str(p.relative_to(TASK)).replace('\\','/'):digest(p) for p in sorted((TASK/'scripts').glob('*.py'))}==selection['code_hashes']
for rel, info in selection['development_input_hashes'].items():
    assert digest(TASK/rel)==info['sha256']

table=list(csv.DictReader((TASK/'analysis/full_waveform.csv').open(encoding='utf-8-sig')))
expected={(r['case'],r['method']):float(r['error']) for r in table}
methods=['Standard','FullCalRefit','SC','Joint','GOnly','NOnly']
values={m:[] for m in methods}
weakratios={}
maxdiff=0.
for name in ['F3','F4','V4','V5']:
    for noise in ['Moderate','Strong']:
        key=f'{name}_{noise}'
        cm,D,sigma,st,tr=load_scene(name,noise)
        target=np.asarray(next(x for x in tr['tr'] if x['name']==noise)['target'])
        outputs=np.load(TASK/'outputs'/f'{key}.npz')
        for method in methods:
            actual=float(np.sqrt(np.sum((outputs[method]-target)**2)/np.sum(target**2)))
            maxdiff=max(maxdiff,abs(actual-expected[key,method]))
            values[method].append(actual)
        if name=='V4':
            errors={}
            for method in ['Standard','Joint']:
                residuals=[];truths=[]
                for j in range(target.shape[1]):
                    for e in np.asarray(tr['weak']).ravel().astype(int)-1:
                        event=tr['events'][j,e]
                        if not np.isfinite(event): continue
                        k=int(event)-1
                        amplitude=tr['reflectivity'][k,j]
                        if amplitude==0: continue
                        ix=np.arange(max(0,k-10),min(target.shape[0],k+11))
                        z=np.pi*cm['fpeak']*(ix-k)*cm['dt']
                        wave=(1-2*z*z)*np.exp(-z*z)
                        estimate=np.dot(outputs[method][ix,j],wave)/np.dot(wave,wave)
                        residuals.append(estimate-amplitude);truths.append(amplitude)
                errors[method]=float(np.linalg.norm(residuals)/np.linalg.norm(truths))
            weakratios[key]=errors['Joint']/errors['Standard']
means={m:float(np.mean(v)) for m,v in values.items()}
ratio=np.asarray(values['Joint'])/values['Standard']
gates={
    'development_improvement_ge_2pct':selection['relative_improvement']>=.02,
    'diagnostic_improvement_ge_3pct':1-means['Joint']/means['Standard']>=.03,
    'worst_case_degradation_le_5pct':bool(np.max(ratio)<=1.05),
    'V4_weak_degradation_le_5pct':max(weakratios.values())<=1.05}
summary=json.loads((TASK/'analysis/summary.json').read_text())
assert maxdiff<1e-12
assert gates==summary['gates']
assert all(abs(means[m]-summary['means']['diagnostic'][m])<1e-12 for m in methods)
assert all(abs(weakratios[k]-summary['V4_weak_ratios'][k])<1e-12 for k in weakratios)
for filename, field in [('solver_review.json','submitted_joint_sha256'),('covariance_review.json','joint_sha256')]:
    assert json.loads((ROOT/filename).read_text())[field]==selection['code_hashes']['scripts/joint.py']
result=dict(passed=True,code_and_development_inputs_unchanged=True,
    source_sha256=selection['code_hashes'],selection_sha256=digest(TASK/'analysis/selection.json'),
    maximum_metric_difference=maxdiff,diagnostic_means=means,
    diagnostic_improvement=1-means['Joint']/means['Standard'],
    worst_case_ratio=float(np.max(ratio)),V4_weak_ratios=weakratios,gates=gates,
    accepts_candidate_for_promotion=all(gates.values()))
(ROOT/'results_review.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
print(json.dumps(result,indent=2))
