from pathlib import Path
import csv,json,sys,hashlib
import numpy as np

ROOT=Path(__file__).resolve().parent
def writecsv(path,rows):
 with path.open('w',encoding='utf-8-sig',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def masks(target,t,live):
 local=np.sqrt(np.apply_along_axis(lambda a:np.convolve(a*a,np.ones(15)/15,mode='same'),0,target));region=(t[:,None]>=.3)&(t[:,None]<=1.6)&live[None,:]
 valid=region&(local>.001*local.max());lo,hi=np.quantile(local[valid],[.20,.45]);return region&(local>=lo)&(local<=hi)
def score(s,target,t,live,low):
 e=s-target;deep=(t>=1.05)&(t<=1.6);a=target[:,live].ravel();b=s[:,live].ravel();aa=a-a.mean();bb=b-b.mean()
 ratio=lambda x,y:float(np.linalg.norm(x)/max(np.linalg.norm(y),1e-30))
 return dict(RelativeError=ratio(e,target),LiveTraceError=ratio(e[:,live],target[:,live]),LowEnergyRegionError=ratio(e[low],target[low]),DeepLiveError=ratio(e[deep,:][:,live],target[deep,:][:,live]),CorrelationLive=float(np.dot(aa,bb)/max(np.linalg.norm(aa)*np.linalg.norm(bb),1e-30)),AmplitudeProjectionGain=float(np.dot(a,b)/max(np.dot(a,a),1e-30)))
def main():
 stage=sys.argv[1];spec=json.loads((ROOT/'SPEC.json').read_text(encoding='utf-8'));run=Path(json.loads((ROOT/f'{stage}_run.json').read_text())['run']);status=json.loads((run/'execution.json').read_text());assert status['status']=='COMPLETED' and status['exit_code']==0
 manifest=json.loads((run/'SOURCE_MANIFEST.json').read_text());assert all(hashlib.sha256((run/'source'/p).read_bytes()).hexdigest()==h for p,h in manifest.items())
 methods=['FullCalRefit','PrunedFullCalRefit','QAwareRidge','QAwareCluster','Trust05','Trust10','Trust20'] if stage=='development' else spec['methods']
 rows=[];modes=[];traces=[];checks=[]
 for case in status['cases']:
  truth=np.load(ROOT/'truth'/f'{case}.npz');t=truth['t'];live=truth['live'];check=json.loads((run/case/'checks.json').read_text());checks.append(check)
  for condition in ['Moderate','Strong']:
   target=truth[condition+'_target'];low=masks(target,t,live);result=json.loads((run/case/(condition+'.json')).read_text())
   for method in ['Observed']+methods:
    s=np.array(result[method]);assert s.shape==target.shape and np.isfinite(s).all();rows.append(dict(Stage=stage,Case=case,Condition=condition,Method=method,**score(s,target,t,live,low)))
    for j,index in enumerate(spec['indices_zero_based']):traces.append(dict(Case=case,Condition=condition,Method=method,OriginalTrace=index+1,Error=float(np.linalg.norm(s[:,j]-target[:,j])/np.linalg.norm(target[:,j]))))
   for method,di in result['diagnostics'].items():
    for j,d in enumerate(di):
     path=np.atleast_1d(d['objective_path']);assert np.all(np.diff(path)<=1e-9*max(path[0],1e-30));assert d['normal_residual']<1e-8
     if 'kappa' in d:assert d['max_protected_change_sd']<=d['kappa']+1e-7
     modes.append(dict(Case=case,Condition=condition,Method=method,OriginalTrace=spec['indices_zero_based'][j]+1,RankBefore=d['rank_before'],RankAfter=d['rank_after'],Removed=d['rank_before']-d['rank_after'],ProtectedModes=d.get('protected_modes',''),MaxProtectedChangeSD=d.get('max_protected_change_sd',''),GuardRejections=d.get('guard_rejections',''),ReferenceOutputChange=d.get('reference_change_norm',''),NormalResidual=d['normal_residual']))
 agg={m:{k:float(np.mean([r[k] for r in rows if r['Method']==m])) for k in ['LiveTraceError','LowEnergyRegionError','DeepLiveError','CorrelationLive','AmplitudeProjectionGain']} for m in methods}
 gates={};selection=None
 if stage=='development':
  names=['Trust05','Trust10','Trust20'];ks=dict(Trust05=.5,Trust10=1.,Trust20=2.);feasible=[m for m in names if agg[m]['LowEnergyRegionError']<=1.02*agg['FullCalRefit']['LowEnergyRegionError']]
  selected=min(feasible or names,key=lambda m:(agg[m]['LiveTraceError'],ks[m]))
  selection=dict(method=selected,kappa=ks[selected],feasible=bool(feasible),all_development_means=agg,rule='Among low-energy feasible candidates minimize development waveform error; tie: smaller kappa. No validation scores used.')
  selected_path=ROOT/'SELECTED.json'
  if selected_path.exists():
   prior=json.loads(selected_path.read_text(encoding='utf-8'));assert prior['method']==selected and prior['kappa']==ks[selected],'Replay would change frozen development selection'
  else:selected_path.write_text(json.dumps(selection,indent=2),encoding='utf-8')
 else:
  a=agg['QAwareOutputTrust'];b=agg['FullCalRefit'];gates['overall_2pct']=a['LiveTraceError']<=.98*b['LiveTraceError']
  gates['each_model_within_2pct']=all(np.mean([r['LiveTraceError'] for r in rows if r['Case']==c and r['Method']=='QAwareOutputTrust'])<=1.02*np.mean([r['LiveTraceError'] for r in rows if r['Case']==c and r['Method']=='FullCalRefit']) for c in status['cases'])
  gates['better_than_Qridge_1pct']=a['LiveTraceError']<=.99*agg['QAwareRidge']['LiveTraceError'];gates['better_than_NoQ_05pct']=a['LiveTraceError']<=.995*agg['OutputTrustNoQ']['LiveTraceError']
  gates['low_energy_within_2pct']=a['LowEnergyRegionError']<=1.02*b['LowEnergyRegionError'];new=[r for r in modes if r['Method']=='QAwareOutputTrust'];changed=sum(r['Removed']>0 for r in new)/len(new)
  gates['nontrivial_deletions_10pct']=changed>=.10;gates['numerical']=all(c['optimization_flags']==0 and c['max_normal']<1e-8 and c['max_guard_ratio']<=1+1e-6 and c['controls']['passed'] for c in checks);gates['all_pass']=all(gates.values())
 writecsv(ROOT/'analysis'/f'{stage}_metrics.csv',rows);writecsv(ROOT/'analysis'/f'{stage}_modes.csv',modes);writecsv(ROOT/'analysis'/f'{stage}_trace_metrics.csv',traces)
 output=dict(stage=stage,run=str(run),means=agg,gates=gates,checks=checks,selection=selection,scope='Descriptive means over one development geology or three new same-family geologies; two dependent observations per geology.')
 (ROOT/'analysis'/f'{stage}_summary.json').write_text(json.dumps(output,indent=2),encoding='utf-8');print(json.dumps(dict(stage=stage,means=agg,gates=gates,selected=selection and selection['kappa']),indent=2))
if __name__=='__main__':main()
