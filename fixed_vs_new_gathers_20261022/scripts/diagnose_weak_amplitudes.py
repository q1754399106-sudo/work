"""Post-hoc signed template-amplitude diagnosis; no refitting or inverse changes."""
from pathlib import Path
import json,csv,statistics,collections
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
def main():
 run=Path(json.loads((ROOT/'LATEST_RUN.json').read_text(encoding='utf-8'))['run'])
 source=run/'source/core_p1'
 extra=json.loads((source/'inputs/cached_shrinkage/trial_plan.json').read_text())['additional_models']
 configs={m['name']:m for m in extra}
 configs['OriginalReference']={'event_t':[.38,.65,.90],'fpeak':60}
 metrics=list(csv.DictReader((run/'matlab/metrics.csv').open(encoding='utf-8-sig')))
 rows=[];maxcheck=0
 for p in sorted((run/'matlab/cases').glob('*_WeakSection_*.json')):
  d=json.loads(p.read_text(encoding='utf-8'));meta=d['meta'];cfg=configs[meta['model']]
  t=np.asarray(d['t']).reshape(-1);true=[];old=[];new=[]
  for j in range(11):
   times=np.asarray(cfg['event_t'],float).copy();times[1]+=.012*np.sin(j*np.pi/10)
   z=np.pi*cfg['fpeak']*(t[:,None]-times[None,:]);W=(1-2*z*z)*np.exp(-z*z)
   for key,result in [('truth',true),('fixed',old),('new',new)]:
    arr=np.asarray(d[key],float);a=np.linalg.lstsq(W,arr[:,j],rcond=None)[0];result.append(a[1])
    if key=='truth':assert np.linalg.norm(W@a-arr[:,j])<1e-11
  ref=np.asarray(true)
  for method,values in [('FixedSplit',old),('DecoupledWeak',new)]:
   est=np.asarray(values);gain=float(ref@est/(ref@ref));err=float(np.linalg.norm(est-ref)/np.linalg.norm(ref))
   remainder=float(np.linalg.norm(est-gain*ref)/np.linalg.norm(ref))
   expected=float(next(r for r in metrics if r['CaseID']==meta['id'] and r['Method']==method)['WeakEventError'])
   maxcheck=max(maxcheck,abs(err-expected));assert abs(err-expected)<1e-9
   assert abs(err**2-((gain-1)**2+remainder**2))<1e-10
   rows.append({'CaseID':meta['id'],'Model':meta['model'],'SNR_dB':meta['snr_db'],'Seed':meta['seed'],'Method':method,
    'WeakEventError':err,'SignedCoherentGain':gain,'LateralAmplitudeDistortion':remainder,
    'PolarityMismatchCount':int(np.sum(ref*est<0)),'TraceCount':11,'TrueAmplitudeMean':float(np.mean(ref)),
    'RecoveredAmplitudeMean':float(np.mean(est))})
 out=ROOT/'analysis/weak_signed_amplitudes.csv'
 with out.open('w',encoding='utf-8-sig',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
 g=collections.defaultdict(list)
 for r in rows:g[r['Model'],r['Method']].append(r)
 summary=[]
 for (model,method),rs in g.items():
  summary.append({'Model':model,'Method':method,'Cases':len(rs),'MeanWeakError':statistics.mean(r['WeakEventError'] for r in rs),
   'MeanSignedGain':statistics.mean(r['SignedCoherentGain'] for r in rs),
   'MeanLateralDistortion':statistics.mean(r['LateralAmplitudeDistortion'] for r in rs),
   'PolarityMismatches':sum(r['PolarityMismatchCount'] for r in rs),'Traces':sum(r['TraceCount'] for r in rs)})
 result={'score_replay_max_absolute_difference':maxcheck,'interpretation':'Signed gain=1 is correct coherent amplitude; 0 is missing coherent component; negative is reversed. Distortion is the orthogonal remainder. Post-hoc diagnostic, not a causal ablation.','summary':summary}
 (ROOT/'analysis/weak_amplitude_diagnosis.json').write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf-8')
 print(json.dumps(result,ensure_ascii=False,indent=2))
if __name__=='__main__':main()
