from pathlib import Path
import argparse,csv,hashlib,json
import numpy as np

R=Path(__file__).resolve().parent
ERRORS=[-30,-10,0,10,30]
METHODS=['fixed','unpruned','hard','soft']

def weak_metrics(a,primary):
    rr=np.asarray(primary['reflectivity']);dt=primary['config']['dt'];freq=primary['config']['fpeak']
    estimate=[];truth=[]
    for j,events in enumerate(primary['event_indices_zero_based']):
        for e in [1,3,5,7]:
            k=events[e];ii=np.arange(max(0,k-10),min(len(a),k+11));z=np.pi*freq*(ii-k)*dt
            w=(1-2*z*z)*np.exp(-z*z)
            estimate.append(np.dot(a[ii,j],w)/np.dot(w,w));truth.append(rr[k,j])
    estimate=np.array(estimate);truth=np.array(truth)
    return dict(weak_amplitude_error=float(np.linalg.norm(estimate-truth)/np.linalg.norm(truth)),
                weak_projection_gain=float(np.dot(estimate,truth)/np.dot(truth,truth)))

def analyze(final=False):
    primary=json.loads((R/'inputs/primary.json').read_text(encoding='utf-8'));rows=[];scenes=[]
    for model in ['primary','complex']:
        reference=R/'results'/model/'reference.json'
        if not reference.exists():continue
        ref=json.loads(reference.read_text(encoding='utf-8'))
        for bias in ERRORS:
            path=R/'results'/model/f'Q_{bias:+d}.json'
            if not path.exists():continue
            result=json.loads(path.read_text(encoding='utf-8'))
            for data,record in zip(ref['datasets'],result['records']):
                assert data['name']==record['name'];target=np.array(data['target']);obs=np.array(data['D'])
                live=np.linalg.norm(obs,axis=0)>0
                scene=dict(model=model,condition=data['name'],Q_percent=bias)
                scenes.append(scene)
                for method in METHODS:
                    a=np.array(record[method]);assert a.shape==target.shape and a.shape[1]==99 and np.isfinite(a).all()
                    if model=='complex':assert np.max(abs(a[:,[22,77]]))==0
                    row=dict(**scene,method=method,relative_error=float(np.linalg.norm(a-target)/np.linalg.norm(target)),
                        live_error=float(np.linalg.norm((a-target)[:,live])/np.linalg.norm(target[:,live])),
                        weak_amplitude_error=None,weak_projection_gain=None)
                    if model=='primary':row.update(weak_metrics(a,primary))
                    rows.append(row)
    if final:
        checks=json.loads((R/'checks.json').read_text(encoding='utf-8'))
        assert checks['complete'] and checks['blocks']==180 and checks['optimization_flags']==0
        assert checks['max_cached_baseline_difference']<1e-6
        assert int((R/'matlab_exit_code.txt').read_text(encoding='utf-8-sig'))==0
        assert len(rows)==80 and len(scenes)==20
        for file,h in json.loads((R/'UPSTREAM_HASHES.json').read_text()).items():
            assert hashlib.sha256((R/file).read_bytes()).hexdigest()==h,file
        for file,h in json.loads((R/'TRIAL_CODE_HASHES.json').read_text()).items():
            assert hashlib.sha256((R/file).read_bytes()).hexdigest()==h,file
        old=R.parent/'q_dense_gathers_20261023'
        for file,h in json.loads((R/'BASELINE_HASHES.json').read_text()).items():
            assert hashlib.sha256((old/file).read_bytes()).hexdigest()==h,file
    lookup={(r['model'],r['condition'],r['Q_percent'],r['method']):r for r in rows}
    comparisons=[]
    for scene in scenes:
        key=(scene['model'],scene['condition'],scene['Q_percent']);soft=lookup[key+('soft',)]
        for m in ['fixed','unpruned','hard']:
            other=lookup[key+(m,)]
            comparisons.append(dict(**scene,baseline=m,
                relative_error_change_percent=100*(soft['relative_error']/other['relative_error']-1),
                weak_error_change_percent=100*(soft['weak_amplitude_error']/other['weak_amplitude_error']-1) if scene['model']=='primary' else None))
    summary=dict(complete=final,scenes=len(scenes),metrics=rows,comparisons=comparisons,
        overall_wins={m:sum(r['relative_error_change_percent']<0 for r in comparisons if r['baseline']==m) for m in ['fixed','unpruned','hard']},
        weak_wins={m:sum(r['weak_error_change_percent']<0 for r in comparisons if r['baseline']==m and r['weak_error_change_percent'] is not None) for m in ['fixed','unpruned','hard']})
    name='summary.json' if final else 'provisional.json'
    (R/'analysis'/name).write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding='utf-8')
    if final:
        with (R/'analysis/metrics.csv').open('w',encoding='utf-8-sig',newline='') as f:
            w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
        with (R/'analysis/comparisons.csv').open('w',encoding='utf-8-sig',newline='') as f:
            w=csv.DictWriter(f,fieldnames=list(comparisons[0]));w.writeheader();w.writerows(comparisons)
    print(json.dumps({k:v for k,v in summary.items() if k not in ['metrics','comparisons']},ensure_ascii=False))
    return summary

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--final',action='store_true');analyze(p.parse_args().final)
