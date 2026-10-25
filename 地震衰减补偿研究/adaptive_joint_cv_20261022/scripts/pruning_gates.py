import math,statistics
KEYS=['RelativeError','EventAmplitudeError','WeakEventError','AVOGradientError','NoiseWindowRMS','FaultNeighborhoodError','PinchoutTailError']
def summarize(rows):
    methods=sorted({r['Method'] for r in rows});models=sorted({r['Model'] for r in rows})
    def avg(method,key,model=None):
        v=[float(r[key]) for r in rows if r['Method']==method and (model is None or r['Model']==model) and math.isfinite(float(r[key]))]
        return statistics.mean(v) if v else None
    s={'cases':len({r['CaseID'] for r in rows}),'means':{m:{k:avg(m,k) for k in KEYS} for m in methods},'per_model':{model:{m:{k:avg(m,k,model) for k in KEYS} for m in methods} for model in models}}
    baseline={r['CaseID']:r for r in rows if r['Method']=='DecoupledWeak'}
    s['case_wins']=sum(float(r['RelativeError'])<float(baseline[r['CaseID']]['RelativeError']) for r in rows if r['Method']=='PrunedFullCalRefit')
    return s
def gates(s,g,checks):
    a=s['means']['PrunedFullCalRefit'];b=s['means']['DecoupledWeak'];old=s['means']['FixedSplit'];ratio=lambda k:a[k]/b[k]
    d={'overall':ratio('RelativeError')<=g['overall_max_ratio'],'each_model':all(v['PrunedFullCalRefit']['RelativeError']<=g['each_model_max_ratio']*v['DecoupledWeak']['RelativeError'] for v in s['per_model'].values()),'weak':ratio('WeakEventError')<=g['weak_max_ratio'],'weak_vs_original':a['WeakEventError']<=g['weak_vs_original_max_ratio']*old['WeakEventError'],'avo':ratio('AVOGradientError')<=g['avo_max_ratio'],'noise':ratio('NoiseWindowRMS')<=g['noise_max_ratio'],'fault':ratio('FaultNeighborhoodError')<=g['fault_max_ratio'],'pinchout':ratio('PinchoutTailError')<=g['pinchout_max_ratio'],'solver':checks['optimization_flags']==0 and checks['max_normal']<1e-8}
    d['all_pass']=all(d.values());return d
