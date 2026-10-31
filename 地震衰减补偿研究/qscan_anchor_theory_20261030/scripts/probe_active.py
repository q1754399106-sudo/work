"""Exploratory (THEORY_ZH.md 5h): kappa with weights = share of non-zero lasso coefficients per time sample.

Development seed 9940, c 0.7 / 1.0 / 1.4, anchor and scan at the truth, lambda 1/2; prints per profile the
weighted kappa and the active share in 0.3 s bins from 1.2 s. Output kept in logs/probe_active.log.
"""
import sys; sys.path.insert(0,'.')
from env15 import np, lib, semi, sc11, bm, gen
from qprofile import context
bm.MAXITER=20000
base=gen.prepare(10227,[10149,10273,10355]); TW=semi.TW
PROFS={'base':[78.5,78.5,135,135],'S':[30,200,200,200],'D':[250,250,40,40]}
for name,Qp in PROFS.items():
    ctx=dict(base) if name=='base' else None
    if name=='base':
        ctx=dict(base)
    else:
        import qprofile; qprofile.PROFILES[name]=Qp; ctx=context(base,name)
    cm=ctx['cm']; Q=np.asarray(cm['Q'],float); a=lib.absorption(cm,TW); a0=lib.absorption(cm,[1.4])[0]
    res=[]
    for c in (0.7,1.0,1.4):
        g=gen.signal(ctx,c,9940); D=g['x']+0.5*g['level']*ctx['noise'][0]; pp=sc11.prepare_sc(D,cm)
        wv=lib.wavelets(D,TW,dict(cm,Q=Q*c)); cc=dict(cm,Q=Q*c); B,_,scale=lib.dictionary(cc,wv['w_in'],wv['w_out'])
        X,info=bm.solve(B,np.eye(D.shape[0]),D,pp['sigma'],scale,pp['K'],0.5,0.0,pp['live'])
        nz=(np.abs(X[:,pp['live']])>1e-8).mean(axis=1)            # active fraction per time sample
        en=(X[:,pp['live']]**2).mean(axis=1)*0+np.abs(X[:,pp['live']]).mean(axis=1)
        k1=a0*np.sum(nz*a)/np.sum(nz*a**2)
        # smoothed active fraction per 0.2 s
        res.append((c, round(float(k1),3), np.round([nz[(TW>=t)&(TW<t+0.3)].mean() for t in np.arange(1.2,3.0,0.3)],3).tolist()))
    print(name,'uniform %.3f'%(a0*a.mean()/(a**2).mean()), res, flush=True)
