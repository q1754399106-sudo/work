"""Frozen SC solver; reuse verified upstream caches, or solve without reading truth."""
import os
for _key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):
    os.environ[_key]='1'
from pathlib import Path
import sys,json,hashlib,time,argparse
from concurrent.futures import ProcessPoolExecutor,as_completed
ROOT=Path(__file__).resolve().parent
RESEARCH=ROOT.parent
sys.path.insert(0,str(RESEARCH/'novelty_evidence_20261025/cache/deps'))
sys.path.insert(0,str(ROOT/'source/sc/scripts'))
import benchmark as bm
import numpy as np
from scipy.io import loadmat,savemat

def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def aslist(x):return x if isinstance(x,list) else [x]
def savejson(p,x):Path(p).write_text(json.dumps(x,indent=2,ensure_ascii=False,allow_nan=False),encoding='utf-8')

def descriptions(stage):
    ds=[]
    for n in ['F1','F2','F3','F4','V1','V2','V4','V5']+[f'U{k}' for k in range(1,5)]:
        method='QScan+SC' if n=='V5' else 'SC'
        ds.append(dict(case=n,model=n,family='Kennett 细层' if n[0] in 'FU' else '复杂地层' if n=='V5' else '稀疏反射模型',method=method))
    for k in range(1,7):
        for mode in ['SC','QSC']:
            ds.append(dict(case=f'H{k}'+('_Qminus30' if mode=='QSC' else ''),model=f'H{k}',family='二维有限差分',method='QScan+SC' if mode=='QSC' else 'SC'))
    if stage=='archived':return ds
    return [dict(case=p.stem,model=p.stem,family='四层模型',method='SC') for p in sorted((ROOT/'cache/cases').glob('*.mat')) if p.stem!='V3_qscan']+[dict(case='V3',model='V3',family='复杂地层',method='QScan+SC')]

def observation(task):
    """Only observed data/configuration. Does not open a target/truth file."""
    n=task['model'];path=None;cached=None;cachedinfo=None;extras={};sources=[]
    def read(p):
        p=Path(p);sources.append(dict(path=str(p.relative_to(RESEARCH)).replace('\\','/'),sha256=digest(p)))
        return loadmat(p,simplify_cells=True)
    inp=ROOT/'cache/inputs'
    if n.startswith('H') and n[1:].isdigit():
        level='Q070_14dB' if task['method']=='QScan+SC' else 'Q100_14dB'
        v=read(inp/f'fd2d_20261025/models/{n}.mat');o=read(RESEARCH/f'sc_neutral_20261026/cache/noisy/{n}_L1_r1.mat')
        cm=dict(v['config']);D=o['D'];sigma=bm.vec(o['sig']);noise='14 dB, replicate 1';seed=1400000000+10000000*int(n[1:])+100000+1
        folder='qscan_sc_20261026' if task['method']=='QScan+SC' else 'sc_neutral_20261026'
        lab='QSC' if task['method']=='QScan+SC' else 'SC'
        cached=RESEARCH/f'{folder}/cache/solved/{n}_{level}_r1_{lab}.npz'
        cachedinfo=cached.with_suffix('.json')
        info=json.loads(cachedinfo.read_text(encoding='utf-8'));factor=float(info.get('qscan_factor',1))
        if lab=='QSC':cm['Q']=np.asarray(cm['Q'])*.7*factor
    elif n in ['F1','F2','F3','F4','V4','V5','V1','V2','V3']:
        base=inp/'round3_20261024/dev' if n in ['F1','F2','F3','F4','V4','V5'] else inp/'coherence_soft_20261023/inputs/validation'
        v=read(base/f'{n}.mat');o=next(x for x in aslist(v['obs']) if x['name']=='Moderate');D=o['D'];sigma=bm.vec(o['sigma']);cm=dict(v['config']);noise='Moderate, 14 dB';seed=None;factor=1
        if n in ['F1','F2','F3','F4','V4']:
            cached=RESEARCH/f'novelty_evidence_20261025/cache/solved/{n}_Moderate_0.125_10.npz';cachedinfo=cached.with_suffix('.json')
        if n=='V5':
            res=read(inp/f'round3_20261024/results/{n}_Q+0_Moderate.mat');factor=float(res['scan']['factor']);cm=dict(res['scanned'])
        if n=='V3':
            res=read(ROOT/'cache/cases/V3_qscan.mat');cm=dict(res['config']);factor=float(res['factor'])
    elif n.startswith('U'):
        v=read(inp/f'b0_world_truth_20261024/inputs/{n}.mat');o=next(x for x in aslist(v['obs']) if float(x['snr_db'])==14 and int(x['replicate'])==1)
        D=o['D'];sigma=bm.vec(o['sigma']);cm=dict(v['config']);noise='SNR14, replicate 1';seed=int(o['seed']);factor=1
    else:
        v=read(ROOT/f'cache/cases/{n}.mat');D=v['D'];sigma=bm.vec(v['sigma']);cm=dict(v['config']);seed=int(v['seed']);noise='new display scene, white noise 14 dB';factor=1
    if cached is not None and cached.exists() and cachedinfo.exists():
        sources += [dict(path=str(p.relative_to(RESEARCH)).replace('\\','/'),sha256=digest(p)) for p in [cached,cachedinfo]]
    else:cached=None;cachedinfo=None
    return cm,np.asarray(D,float),sigma,cached,cachedinfo,dict(noise=noise,seed=seed,qscan_factor=factor,sources=sources)

def solve_one(task):
    dest=ROOT/'results'/f'{task["case"]}.mat';meta=dest.with_suffix('.json')
    if dest.exists() and meta.exists():return json.loads(meta.read_text(encoding='utf-8'))
    cm,D,sigma,cached,cachedinfo,context=observation(task)
    if cached:
        X=np.load(cached)['out'];info=json.loads(cachedinfo.read_text(encoding='utf-8'));reused=True
    else:
        B,W,P,scale=bm.dictionary(cm);live=D.any(axis=0);K,_,adj=bm.predictor(D,live)
        X,info=bm.solve(B,P,D,sigma,scale,K,.125,10,live);info['predictor_adjoint_error']=adj;reused=False
    assert X.shape==D.shape and np.isfinite(X).all()
    # Evaluation starts only after the solver/reuse phase is finished.
    n=task['model'];inp=ROOT/'cache/inputs';extras={}
    if n.startswith('H') and n[1:].isdigit():targetfile=inp/f'fd2d_20261025/truth/{n}_truth.mat'
    elif n.startswith('U'):targetfile=inp/f'b0_world_truth_20261024/truth/{n}_truth.mat'
    elif n in ['F1','F2','F3','F4','V4','V5']:targetfile=inp/f'round3_20261024/dev/{n}_truth.mat'
    elif n in ['V1','V2','V3']:targetfile=inp/f'coherence_soft_20261023/truth/{n}_truth.mat'
    else:targetfile=ROOT/f'cache/truth/{n}.mat'
    tr=loadmat(targetfile,simplify_cells=True)
    T=np.asarray(tr['target'] if 'target' in tr else next(x for x in aslist(tr['tr']) if x['name']=='Moderate')['target'],float)
    if n in ['V4','V5']:
        res=loadmat(inp/f'round3_20261024/results/{n}_Q+0_Moderate.mat',simplify_cells=True)
        extras={k:np.asarray(v) for k,v in res['gathers'].items() if k in ['Standard','QScan_Standard']}
        sc=RESEARCH/f'novelty_evidence_20261025/cache/solved/{n}_Moderate_0.125_10.npz'
        extras['SC']=np.load(sc)['out']
    elif n=='V3':extras['QScan_Standard']=np.asarray(loadmat(ROOT/'cache/cases/V3_qscan.mat',simplify_cells=True)['outcome']['gather'])
    live=D.any(axis=0);den=np.linalg.norm(T);error=lambda A:float(np.linalg.norm(A-T)/den)
    assert T.shape==D.shape
    record=dict(**task,**context,reused_upstream=reused,solver=info,input_error=error(D),output_error=error(X),
                output_live_error=float(np.linalg.norm((X-T)[:,live])/np.linalg.norm(T[:,live])),
                dead_traces=(np.where(~live)[0]+1).tolist(),target_sha256=digest(targetfile),
                snr_db=14,parameters=dict(a=.125,mu_relative=10,maxiter=5000,kkt_tolerance=.02),
                extra_errors={key:error(value) for key,value in extras.items()})
    dest.parent.mkdir(parents=True,exist_ok=True)
    savemat(dest,dict(D=D,target=T,compensated=X,t=np.asarray(cm['t']).reshape(-1,1),sigma=sigma,config=cm,**extras),do_compression=True)
    savejson(meta,record)
    print(f"{task['case']} {task['method']} reused={reused} error {record['input_error']:.4f} -> {record['output_error']:.4f}, KKT {info['kkt_ratio']:.4f}",flush=True)
    return record

def main():
    ap=argparse.ArgumentParser();ap.add_argument('stage',choices=['archived','new','index']);ap.add_argument('--workers',type=int,default=3);args=ap.parse_args()
    if args.stage!='index':
        tasks=descriptions(args.stage)
        with ProcessPoolExecutor(max_workers=args.workers) as pool:
            futures=[pool.submit(solve_one,x) for x in tasks]
            for f in as_completed(futures):f.result()
    records=[json.loads(p.read_text(encoding='utf-8')) for p in sorted((ROOT/'results').glob('*.json'))]
    savejson(ROOT/'CASE_INDEX.json',records)
    import csv
    with (ROOT/'metrics.csv').open('w',encoding='utf-8-sig',newline='') as f:
        w=csv.DictWriter(f,fieldnames=['case','model','family','method','noise','input_error','output_error','output_live_error','reused_upstream','converged','kkt_ratio']);w.writeheader()
        for r in records:w.writerow({**{key:r[key] for key in w.fieldnames if key not in ['converged','kkt_ratio']},'converged':r['solver']['converged'],'kkt_ratio':r['solver']['kkt_ratio']})
    print('COMPLETE CASES',len(records),flush=True)

if __name__=='__main__':main()
