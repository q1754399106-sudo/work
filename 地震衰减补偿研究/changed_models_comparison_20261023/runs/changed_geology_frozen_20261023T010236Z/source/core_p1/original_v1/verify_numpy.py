"""Independent numerical cross-check, NOT execution of the MATLAB files.
Requires NumPy only; writes compact JSON. Run: python verify_numpy.py
"""
import json
from pathlib import Path
import numpy as np

def ray(t, x, vp, Q, bounds):
    a=np.zeros(len(t)); b=a.copy(); small=a.copy()
    for k,t0 in enumerate(t):
        tau=np.maximum(0,np.minimum(t0,bounds[1:])-bounds[:-1])
        use=tau>0
        if not use.any():
            b[k]=x/(vp[0]*Q[0]); continue
        tt=tau[use]; v=vp[use]; q=Q[use]; a[k]=sum(tt/q)
        if x==0: continue
        lo=0.; hi=(1-1e-12)/max(v)
        for _ in range(70):
            p=(lo+hi)/2
            if sum(tt*v*v*p/np.sqrt(1-(p*v)**2))<x: lo=p
            else: hi=p
        p=(lo+hi)/2
        b[k]=sum(tt/q/np.sqrt(1-(p*v)**2))-a[k]
        small[k]=.5*(x/sum(tt*v*v))**2*sum(tt*v*v/q)
    return a,b,small

def matrix(a,dt,N):
    n=len(a); f=np.arange(N//2+1)/(N*dt)
    z=-np.pi*f+2j*f*np.log(np.maximum(f,np.finfo(float).tiny)/100)
    z[0]=0
    H=np.exp(z[:,None]*a[None,:]); H[-1]=H[-1].real
    H=np.concatenate((H,H[-2:0:-1].conj()))
    out=np.fft.ifft(H*np.fft.fft(np.eye(n),n=N,axis=0),axis=0)
    assert np.max(abs(out.imag))<1e-12
    return out[:n].real

def solve(A,y,target):
    U,z,Vh=np.linalg.svd(A,full_matrices=False); uy=U.T@y
    grid=np.logspace(-7,0,41)
    candidates=(Vh.T[:,:,None]*(z[:,None]/(z[:,None]**2+grid)*uy[:,None])[None,:,:]).sum(axis=1)
    norms=np.linalg.norm(A@candidates-y[:,None],axis=0)
    k=np.argmin(abs(norms-target))
    return candidates[:,k],float(grid[k])

def main():
    # Same dt and duration as MATLAB, subset of offsets to limit verification cost.
    dt=.002; t=np.arange(601)*dt; n=len(t); N=4096; I=np.eye(n)
    vp=np.array([2100.,2250,2450,2800]); Q=np.array([40.,50,70,100])
    bounds=np.array([0.,.38,.65,.90,np.inf])
    a,_,_=ray(t,0,vp,Q,bounds); A=matrix(a,dt,N)
    zero=np.linalg.norm(matrix(np.zeros(n),dt,N)-I)/np.sqrt(n)
    homogeneous_t=np.array([.3,.7,1.]); aa,bb,_=ray(homogeneous_t,900,np.full(4,2200.),np.full(4,60.),bounds)
    hom=np.linalg.norm(aa+bb-np.sqrt(homogeneous_t**2+(900/2200)**2)/60)
    assert zero<1e-11 and hom<1e-12
    B0=np.linalg.solve(A.T@A+.001*I,A.T)
    rng=np.random.default_rng(20261017); rows=[]
    for x in (0,750,1500):
        print('NumPy check offset',x,flush=True)
        _,b,small=ray(t,x,vp,Q,bounds); X=matrix(b,dt,N); G=matrix(a+b,dt,N)
        # Chosen amplitudes for independent verification, not MATLAB Shuey AVO.
        s=np.zeros(n)
        for tj,r in zip((.38,.65,.90),(.12,-.10,.15)):
            z=np.pi*60*(t-tj); s+=r*(1-2*z*z)*np.exp(-z*z)
        clean=G@s; sigma=np.linalg.norm(clean)/np.sqrt(n)/10**(.9)
        noise=sigma*rng.standard_normal(n); d=clean+noise
        B=np.linalg.solve(X.T@X+.001*I,X.T); y=B@d; fixed=B0@y
        L=np.linalg.cholesky(B@B.T+1e-6*I)
        cov,lc=solve(np.linalg.solve(L,A),np.linalg.solve(L,y),sigma*np.sqrt(n))
        direct,ld=solve(G,d,sigma*np.sqrt(n))
        unweighted,lu=solve(A,y,sigma*np.linalg.norm(B))
        identity=np.linalg.norm(fixed-s-((B0@A-I)@s+B0@(y-A@s)))/np.linalg.norm(s)
        assert identity<1e-10
        mask=(t>=.25)&(t<=1.05)
        rows.append(dict(offset=x,composition_error=float(np.linalg.norm(G-X@A)/np.linalg.norm(G)),
            small_angle_error=float(np.linalg.norm(small[mask]-b[mask])/max(np.linalg.norm(b[mask]),1e-15)),
            identity_residual=float(identity),lambda_cov=lc,lambda_direct=ld,lambda_unweighted=lu,
            errors={name:float(np.linalg.norm(v-s)/np.linalg.norm(s)) for name,v in
                    [('observed',d),('fixed',fixed),('unweighted',unweighted),('covariance',cov),('direct',direct)]}))
    report=dict(engine='Independent NumPy implementation; MATLAB execution not implied',
                zero_loss_error=float(zero),homogeneous_ray_error=float(hom),rows=rows)
    p=Path(__file__).parent/'verification'; p.mkdir(exist_ok=True)
    (p/'numpy_checks.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    print(json.dumps(report,indent=2))

if __name__=='__main__': main()
