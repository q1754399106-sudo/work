"""Independent implementation of Ma et al. (2023), equation 17.

Convex objective solved by FISTA, not the paper's IRLS implementation.
No author code was available. Physical-reflectivity L1 penalty is retained
under reversible column scaling. No test-set parameter selection.
"""
from common import *
import time,csv,argparse
from scipy.linalg import eigh

NOISES=['Moderate','Strong'];DEV=['F1','F2'];TEST=['F3','F4','V4','V5']
AGRID=[.25,.5,1.];MGRID=[.01,.1,1.]
MAXITER=5000;TOL=.02

def predictor(D,live):
    n,nx=D.shape
    tt,jj=np.meshgrid(np.arange(2,n-2),np.arange(2,nx-2),indexing='ij')
    tt=tt.ravel();jj=jj.ravel()
    keep=np.logical_and.reduce([live[jj+j] for j in [-2,-1,0,1,2]])
    tt,jj=tt[keep],jj[keep]
    offsets=[(i,j) for j in [-2,-1,1,2] for i in range(-2,3)]
    X=np.column_stack([D[tt+i,jj+j] for i,j in offsets])
    b=np.linalg.lstsq(X,D[tt,jj],rcond=None)[0]
    rows=np.repeat(np.arange(len(tt)),21)
    cols=np.column_stack([tt*nx+jj]+[(tt+i)*nx+jj+j for i,j in offsets]).ravel()
    vals=np.tile(np.r_[-1.,b],len(tt))
    K=sparse.csr_matrix((vals,(rows,cols)),shape=(len(tt),n*nx))
    rng=np.random.default_rng(27001);x=rng.normal(size=n*nx);y=rng.normal(size=len(tt))
    adj=abs(np.dot(K@x,y)-np.dot(x,K.T@y))/max(1,abs(np.dot(K@x,y)))
    assert adj<1e-11
    return K,b,float(adj)

def solve(B,P,D,sigma,scale,K,a,mu_rel,live,start=None):
    n,nx=D.shape;gram=B.T@B;rhs=B.T@D
    Ld=float(eigh(gram,subset_by_index=[n-1,n-1],eigvals_only=True)[0])
    def sh(x):return P.T@((K.T@(K@(P@x).ravel())).reshape(n,nx))
    rng=np.random.default_rng(27002);v=rng.normal(size=D.shape);v[:,~live]=0;v/=np.linalg.norm(v)
    for _ in range(35):
        w=sh(v);w[:,~live]=0;nv=np.linalg.norm(w)
        if nv==0:break
        v=w/nv
    Ls=float(np.sum(v*sh(v)))
    mu=mu_rel*Ld/max(Ls,1e-30)
    L=1.05*(Ld+mu*Ls)
    pen=a*np.sqrt(2*np.log(n))*sigma[None,:]*np.median(scale)/scale[:,None]
    pen=np.maximum(pen,1e-14)
    x=np.zeros_like(D) if start is None else start.copy();x[:,~live]=0
    y=x.copy();theta=1.;kkt=np.inf;backtracks=0;beg=time.monotonic()
    def h(z):return gram@z+mu*sh(z) if mu else gram@z
    def obj(z):
        r=B@z-D;s=K@(P@z).ravel()
        return float(.5*np.sum(r*r)+np.sum(pen*np.abs(z))+.5*mu*np.dot(s,s))
    for it in range(1,MAXITER+1):
        hy=h(y);grad=hy-rhs
        while True:
            z=y-grad/L;xx=np.sign(z)*np.maximum(np.abs(z)-pen/L,0);xx[:,~live]=0
            delta=xx-y
            # Exact quadratic majorization test, sparse and dense operators.
            if it==1 or it%25==0:
                hd=h(delta)
                if np.sum(delta*hd)>L*np.sum(delta*delta)*(1+1e-10):
                    L*=2;backtracks+=1;continue
            break
        newtheta=(1+np.sqrt(1+4*theta*theta))/2
        yy=xx+(theta-1)/newtheta*(xx-x)
        if np.sum((y-xx)*(xx-x))>0:yy=xx.copy();newtheta=1.
        x,y,theta=xx,yy,newtheta
        if it%25==0:
            gg=h(x)-rhs;rr=np.maximum(np.abs(gg)-pen,0);active=x!=0
            rr[active]=np.abs(gg[active]+(pen*np.sign(x))[active]);rr[:,~live]=0
            kkt=float(np.max(rr/pen))
            if kkt<=TOL:break
    out=P@x
    return out,dict(iterations=it,kkt_ratio=kkt,converged=kkt<=TOL,seconds=time.monotonic()-beg,
                    a=a,mu_relative=mu_rel,mu=mu,Ldata=Ld,Lstructure=Ls,L=L,backtracks=backtracks,objective=obj(x))

def proposed(W,cm,st,do_align=True):
    R=dense(st['R0']);supports=[vec(s,int)-1 for s in st['supports']]
    E=[vec(v) for v in st['E']];V=[vec(v) for v in st['V']];live=vec(st['live'],bool)
    n,nx=R.shape;gap=max(1,int(np.floor(.5/(cm['fpeak']*cm['dt'])+.5)))
    maxlag=int(np.ceil(.25/(cm['fpeak']*cm['dt']))) if do_align else 0
    groups=[];centers=[];contrib=[]
    for j,s in enumerate(supports):
        g=np.cumsum(np.r_[0,np.diff(s)>gap]);groups.append(g)
        cc=[];ww=[]
        for k in range(g.max()+1 if len(g) else 0):
            ii=s[g==k];w=np.abs(R[ii,j]);cc.append(float(np.dot(ii,w)/max(w.sum(),np.finfo(float).tiny)))
            ww.append(W[:,ii]@R[ii,j])
        centers.append(np.asarray(cc));contrib.append(ww)
    cs=[];fallback=[]
    for j in range(nx):
        inside=np.array([i for i in range(max(0,j-2),min(nx,j+3)) if i!=j]);ng=len(centers[j])
        c=np.zeros(ng);fb=np.ones(ng,bool)
        if live[j] and np.all(live[inside]) and len(inside)>=2:
            fb[:]=False
            for k,x in enumerate(contrib[j]):
                nxx=np.dot(x,x)
                if nxx<=0:continue
                scores=[]
                for m in inside:
                    bound=maxlag*abs(m-j);best=0.
                    for h,cy in enumerate(centers[m]):
                        if abs(cy-centers[j][k])>bound+gap:continue
                        yy=contrib[m][h];ny=np.dot(yy,yy)
                        if ny<=0:continue
                        for lag in range(-bound,bound+1):
                            xy=np.dot(x[lag:],yy[:-lag]) if lag>0 else np.dot(x[:lag],yy[-lag:]) if lag<0 else np.dot(x,yy)
                            if xy>0:best=max(best,float(xy/np.sqrt(nxx*ny)*min(xy/nxx,nxx/xy)))
                    scores.append(best)
                scores=np.asarray(scores);left=scores[inside<j];right=scores[inside>j]
                c[k]=min(1.,max(left.mean() if len(left) else 0,right.mean() if len(right) else 0))
        cs.append(c);fallback.append(fb)
    outs={}
    rules={'Soft':(1,1),'v1':(1,8),'Standard':(.5,16),'Uniform':(1.951,1.951)}
    for name,(lo,hi) in rules.items():
        rr=np.zeros_like(R)
        for j,s in enumerate(supports):
            lam=lo+(hi-lo)*(1-cs[j])**2;lam[fallback[j]]=1
            den=E[j]+lam*V[j];g=np.ones_like(den);np.divide(E[j],den,out=g,where=den>0)
            if len(s):rr[s,j]=R[s,j]*g[groups[j]]
        outs[name]=W@rr
    return outs,cs

def case(name,noise):
    cm,D,sigma,st,tr=load_scene(name,noise)
    target=np.asarray(next(t for t in tr['tr'] if t['name']==noise)['target'],float)
    B,W,P,scale=dictionary(cm);live=vec(st['live'],bool)
    den=max(np.linalg.norm(target),1e-30)
    methods,cs=proposed(W,cm,st)
    archived=archived_mat(f'round3_20261024/results/{name}_Q+0_{noise}.mat')['gathers']
    checks={'full_reconstruction':float(np.linalg.norm(W@dense(st['R0'])-dense(st['S0']))/max(np.linalg.norm(dense(st['S0'])),1e-30)),
            'soft_reconstruction':float(np.linalg.norm(methods['Soft']-dense(st['soft']))/max(np.linalg.norm(dense(st['soft'])),1e-30)),
            'standard_platform_difference':float(np.linalg.norm(methods['Standard']-archived['Standard'])/max(np.linalg.norm(archived['Standard']),1e-30))}
    assert checks['full_reconstruction']<1e-10 and checks['soft_reconstruction']<1e-10,checks
    assert checks['standard_platform_difference']<.02,checks
    methods['FullCalRefit']=dense(st['S0'])
    methods['StandardArchived']=archived['Standard']
    # QScan is reported only; it changes Q and is not a same-Q structural ablation.
    methods['QScan_StandardArchived']=archived['QScan_Standard']
    if name in TEST:
        na,_=proposed(W,cm,st,False);methods['NoAlign']=na['Standard']
    K,b,adj=predictor(D,live);checks['predictor_adjoint_error']=adj
    return cm,D,sigma,st,target,B,W,P,scale,K,live,methods,checks

def error(out,target):return float(np.linalg.norm(out-target)/max(np.linalg.norm(target),1e-30))

def run(mode):
    rows=[];allchecks={};selpath=ROOT/'analysis/selection.json'
    if mode=='test':
        sel=json.loads(selpath.read_text(encoding='utf-8'));pairs=[(sel['a'],sel['mu_relative']),(sel['no_structure_a'],0.)];names=TEST
    else:pairs=[(a,m) for a in AGRID for m in MGRID]+[(a,0.) for a in AGRID];names=DEV
    for name in names:
        for noise in NOISES:
            key=f'{name}_{noise}';beg=time.monotonic();print('CASE',mode,key,flush=True)
            cm,D,sigma,st,target,B,W,P,scale,K,live,methods,checks=case(name,noise)
            allchecks[key]=checks
            for method,out in methods.items():rows.append(dict(case=key,method=method,error=error(out,target)))
            for a,m in pairs:
                p=ROOT/'cache/solved'/f'{key}_{a:g}_{m:g}.npz';p.parent.mkdir(parents=True,exist_ok=True)
                infofile=p.with_suffix('.json')
                if p.exists() and infofile.exists():
                    saved=np.load(p);out=saved['out'];info=json.loads(infofile.read_text())
                else:
                    out,info=solve(B,P,D,sigma,scale,K,a,m,live)
                    np.savez_compressed(p,out=out);writejson(infofile,info)
                label=f'SC_a{a:g}_m{m:g}'
                er=error(out,target);rows.append(dict(case=key,method=label,error=er))
                print(key,label,'err',round(er,6),'KKT',round(info['kkt_ratio'],5),'iter',info['iterations'],'sec',round(info['seconds'],1),flush=True)
                allchecks[f'{key}/{label}']=info
                writejson(ROOT/f'analysis/{mode}_checks.json',allchecks)
                writejson(ROOT/f'analysis/{mode}_rows.json',rows)
            print('DONE',key,round(time.monotonic()-beg,1),flush=True)
    writejson(ROOT/f'analysis/{mode}_rows.json',rows)
    if mode=='dev':
        candidates=[]
        for a,m in pairs:
            label=f'SC_a{a:g}_m{m:g}';v=[r['error'] for r in rows if r['method']==label]
            ok=all(allchecks[f'{n}_{s}/{label}']['converged'] for n in names for s in NOISES)
            candidates.append(dict(a=a,mu_relative=m,mean_error=float(np.mean(v)),all_converged=ok))
        valid=[x for x in candidates if x['all_converged']]
        structured=[x for x in valid if x['mu_relative']>0];unstructured=[x for x in valid if x['mu_relative']==0]
        if not structured or not unstructured:raise RuntimeError('No fully converged candidate; do not evaluate test set.')
        best=min(structured,key=lambda x:x['mean_error']);base=min(unstructured,key=lambda x:x['mean_error'])
        writejson(selpath,dict(**best,no_structure_a=base['a'],all_candidates=candidates,frozen_before_test=True))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['smoke','dev','test']);o=p.parse_args()
    if o.mode=='smoke':
        data=case('F1','Moderate');print(data[-1],flush=True)
        _,D,sigma,st,target,B,W,P,scale,K,live,methods,checks=data
        out,info=solve(B,P,D,sigma,scale,K,.5,.1,live)
        print(info,'error',error(out,target),flush=True)
        writejson(ROOT/'analysis/smoke.json',dict(checks=checks,solver=info,error=error(out,target)))
    else:run(o.mode)
