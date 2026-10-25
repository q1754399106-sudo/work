"""Fixed support ridge noise propagation and convex box quadratic shrinkage."""
from common import *
from scipy.linalg import svd,eigvalsh
from scipy.optimize import minimize
import reference,time

def prepare(cm,D,sigma,st):
    begin=time.perf_counter();B,W,P,scale=dictionary(cm);dictionary_seconds=time.perf_counter()-begin
    coherence_begin=time.perf_counter()
    original,cs=reference.proposed(W,cm,st)
    coherence_reference_seconds=time.perf_counter()-coherence_begin
    propagation_begin=time.perf_counter()
    R=dense(st['R0']);supports=[vec(s,int)-1 for s in st['supports']];live=vec(st['live'],bool)
    gap=max(1,int(np.floor(.5/(cm['fpeak']*cm['dt'])+.5)))
    traces=[];checks=[]
    for j,s in enumerate(supports):
        if not len(s):
            traces.append(None);continue
        groups=np.cumsum(np.r_[0,np.diff(s)>gap]);ng=int(groups.max()+1)
        Y=np.column_stack([W[:,s[groups==k]]@R[s[groups==k],j] for k in range(ng)])
        G=Y.T@Y;C=B[:,s];U,z,Vh=svd(C,full_matrices=False,check_finite=False)
        top=float(z[0]**2);bottom=float(z[-1]**2)
        ridge=max(1e-8*top,(top-1e4*bottom)/(1e4-1))
        shape=Vh.T*(z/(z*z+ridge))[None,:];cov=shape@shape.T
        source=P[:,s];source_gram=source.T@source
        # M is a PSD Schur product; grouping is the congruence Z' M Z.
        Z=np.eye(ng)[groups];N=sigma[j]**2*(Z.T@(source_gram*cov)@Z)
        N=(N+N.T)/2;G=(G+G.T)/2
        E=vec(st['E'][j]);V=vec(st['V'][j]);lam=.5+15.5*(1-cs[j])**2
        inside=np.array([m for m in range(max(0,j-2),min(R.shape[1],j+3)) if m!=j],int)
        fallback=not (live[j] and np.all(live[inside]) and len(inside)>=2)
        if fallback:lam[:]=1
        a=R[s,j]*scale[s]
        normal=np.linalg.norm(C.T@(C@a-D[:,j])+ridge*a)/max(np.linalg.norm(C.T@D[:,j]),1e-30)
        ere=float(np.linalg.norm(np.diag(G)-E)/max(np.linalg.norm(E),1e-30))
        vre=float(np.linalg.norm(np.diag(N)-V)/max(np.linalg.norm(V),1e-30))
        # Preserve the already archived diagonal exactly; only roundoff correction.
        assert ere<1e-10 and vre<1e-7, (j,ere,vre)
        # Congruence scaling makes exact diag N=V while retaining PSD and the
        # full derived correlation matrix. This is solely roundoff alignment.
        dn=np.diag(N);q=np.sqrt(np.divide(V,dn,out=np.ones_like(V),where=dn>0));N=q[:,None]*N*q[None,:]
        dg=np.diag(G);q=np.sqrt(np.divide(E,dg,out=np.ones_like(E),where=dg>0));G=q[:,None]*G*q[None,:]
        mineg=float(eigvalsh(G,check_finite=False)[0]);minen=float(eigvalsh(N,check_finite=False)[0])
        gscale=max(float(np.linalg.norm(G,2)),1e-30);nscale=max(float(np.linalg.norm(N,2)),1e-30)
        assert mineg>=-1e-10*gscale and minen>=-1e-10*nscale
        checks.append(dict(trace=j,groups=ng,support=len(s),ridge=ridge,normal_residual=float(normal),energy_relative_error=ere,variance_relative_error=vre,
            min_eig_G=mineg,min_eig_N=minen,min_eig_G_relative=mineg/gscale,min_eig_N_relative=minen/nscale,fallback=fallback,
            offdiag_G_relative=float(np.linalg.norm(G-np.diag(E))/max(np.linalg.norm(G),1e-30)),offdiag_N_relative=float(np.linalg.norm(N-np.diag(V))/max(np.linalg.norm(N),1e-30))))
        traces.append(dict(Y=Y,G=G,N=N,E=E,V=V,lam=lam))
    return dict(traces=traces,original=original,W=W,preparation_seconds=time.perf_counter()-begin,checks=checks,
        dictionary_seconds=dictionary_seconds,coherence_reference_seconds=coherence_reference_seconds,propagation_seconds=time.perf_counter()-propagation_begin)

def matrices(tr,eta,mode='Joint'):
    eg=eta if mode in ['Joint','GOnly'] else 0.;en=eta if mode in ['Joint','NOnly'] else 0.
    G=(1-eg)*np.diag(tr['E'])+eg*tr['G'];N=(1-en)*np.diag(tr['V'])+en*tr['N']
    root=np.sqrt(tr['lam']);H=G+root[:,None]*N*root[None,:];b=G@np.ones(len(root))
    return G,N,H,b

def solve_trace(tr,eta,mode='Joint'):
    G,N,H,b=matrices(tr,eta,mode);den=tr['E']+tr['lam']*tr['V']
    start=np.divide(tr['E'],den,out=np.ones_like(den),where=den>0)
    sc=max(float(np.max(np.abs(H))),float(np.max(np.abs(b))),1e-30);HH=H/sc;bb=b/sc
    if eta==0:
        g=start;iterations=0
    else:
        res=minimize(lambda x: .5*x@HH@x-bb@x,start,jac=lambda x:HH@x-bb,
            bounds=[(0,1)]*len(b),method='L-BFGS-B',options={'ftol':1e-15,'gtol':1e-11,'maxiter':3000,'maxls':50,'maxcor':30})
        g=np.clip(res.x,0,1);iterations=int(res.nit)
        # Deterministic exact coordinate polishing; meaningful projected KKT.
        for sweep in range(20000):
            grad=HH@g-bb
            pg=g-np.clip(g-grad,0,1)
            if np.max(np.abs(pg))<1e-10:break
            for k in range(len(g)):
                if HH[k,k]<=0:continue
                delta=np.clip(g[k]-grad[k]/HH[k,k],0,1)-g[k]
                if delta:g[k]+=delta;grad+=HH[:,k]*delta
        iterations+=sweep
    grad=HH@g-bb;kkt=float(np.max(np.abs(g-np.clip(g-grad,0,1))))
    assert kkt<1e-8,(mode,eta,kkt)
    obj=float((g-1)@G@(g-1)+(np.sqrt(tr['lam'])*g)@N@(np.sqrt(tr['lam'])*g))
    return g,dict(kkt_projected_normalized=kkt,iterations=iterations,objective=obj)

def solve(prepared,eta,mode='Joint'):
    begin=time.perf_counter();W=prepared['W'];out=np.zeros((W.shape[0],len(prepared['traces'])));checks=[]
    gains=[]
    for j,tr in enumerate(prepared['traces']):
        if tr is None:gains.append([]);continue
        g,info=solve_trace(tr,eta,mode);out[:,j]=tr['Y']@g;checks.append(info);gains.append(g.tolist())
    return out,dict(seconds=time.perf_counter()-begin,max_kkt=max((x['kkt_projected_normalized'] for x in checks),default=0.),
        objective=sum(x['objective'] for x in checks),max_iterations=max((x['iterations'] for x in checks),default=0),gains=gains)

def gradient_check(tr):
    G,N,H,b=matrices(tr,1.);x=np.linspace(.2,.8,len(b));v=np.cos(np.arange(len(b)));h=1e-5
    def f(z):return (z-1)@G@(z-1)+(np.sqrt(tr['lam'])*z)@N@(np.sqrt(tr['lam'])*z)
    fd=(f(x+h*v)-f(x-h*v))/(2*h);analytic=2*v@(H@x-b)
    return float(abs(fd-analytic)/max(1,abs(fd),abs(analytic)))
