"""Post-hoc regional diagnosis using the original scoring definitions.

Does not change the frozen primary gates or select a new winning candidate.
"""
from common import *
import csv
from scipy.ndimage import convolve1d
sel=json.loads((ROOT/'analysis/selection.json').read_text(encoding='utf-8'))
records=[]
for name in ['F3','F4','V4','V5']:
    for ni,noise in enumerate(['Moderate','Strong']):
        cm,D,sigma,st,tr=load_scene(name,noise)
        target=next(x for x in tr['tr'] if x['name']==noise)['target']
        path=ROOT/'cache/solved'/f"{name}_{noise}_{sel['a']:g}_{sel['mu_relative']:g}.npz"
        baseline=np.load(path)['out']
        arc=archived_mat(f'round3_20261024/results/{name}_Q+0_{noise}.mat')['gathers']
        methods={'SC':baseline,'StandardArchived':arc['Standard'],'SoftArchived':arc['SoftGroupRefit'],'FullCalRefitArchived':arc['FullCalRefit']}
        env=np.sqrt(convolve1d(target**2,np.ones(21),axis=0,mode='constant')/convolve1d(np.ones((len(target),1)),np.ones(21),axis=0,mode='constant'))
        quiet=env<.02*np.max(env)
        for method,A in methods.items():
            measures={'overall':float(np.linalg.norm(A-target)/np.linalg.norm(target)),'quiet_leak':float(np.linalg.norm(A[quiet])/np.linalg.norm(target))}
            for region,M in tr.get('regions',{}).items():
                M=np.asarray(M,bool)
                if M.ndim==3:M=M[:,:,ni]
                if M.shape==target.shape and np.any(M) and np.linalg.norm(target[M])>0:
                    measures[region]=float(np.linalg.norm((A-target)[M])/np.linalg.norm(target[M]))
            if 'events' in tr:
                est=[];truth=[]
                for j in range(A.shape[1]):
                    for e in vec(tr['weak'],int)-1:
                        k=tr['events'][j,e]
                        if not np.isfinite(k):continue
                        k=int(k)-1
                        if tr['reflectivity'][k,j]==0:continue
                        ii=np.arange(max(0,k-10),min(len(A),k+11));z=np.pi*cm['fpeak']*(ii-k)*cm['dt'];w=(1-2*z*z)*np.exp(-z*z)
                        est.append(float(A[ii,j]@w/(w@w)));truth.append(float(tr['reflectivity'][k,j]))
                if truth:
                    est=np.asarray(est);truth=np.asarray(truth)
                    measures['weak_projection_error']=float(np.linalg.norm(est-truth)/np.linalg.norm(truth))
                    measures['weak_projection_gain']=float(est@truth/(truth@truth))
            records += [dict(case=f'{name}_{noise}',method=method,metric=k,value=v) for k,v in measures.items()]
with (ROOT/'analysis/regional_posthoc.csv').open('w',newline='',encoding='utf-8-sig') as f:
    w=csv.DictWriter(f,fieldnames=['case','method','metric','value']);w.writeheader();w.writerows(records)
writejson(ROOT/'analysis/regional_posthoc.json',dict(scope='Post-hoc diagnostic, not a registered gate and not a rescue of failed primary comparison.',source_metric='round3/core/score_gather.m and score_weak.m',rows=records))
print('Post-hoc diagnostic rows:',len(records))
