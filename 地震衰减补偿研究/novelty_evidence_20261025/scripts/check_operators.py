"""Numerical derivative and adjoint checks, without any comparison scoring."""
from common import *
from benchmark import predictor
cm,D,sigma,st,_=load_scene('F1','Moderate')
B,W,P,scale=dictionary(cm);live=vec(st['live'],bool);K,b,adj=predictor(D,live)
rng=np.random.default_rng(20261025);x=rng.normal(size=D.shape)*.001;v=rng.normal(size=D.shape);v/=np.linalg.norm(v)
mu=.01
def smooth(x):
    r=B@x-D;s=K@(P@x).ravel()
    return .5*np.sum(r*r)+mu/2*np.dot(s,s)
gradient=B.T@(B@x-D)+mu*P.T@((K.T@(K@(P@x).ravel())).reshape(D.shape))
eps=1e-5;fd=(smooth(x+eps*v)-smooth(x-eps*v))/(2*eps);exact=float(np.sum(v*gradient))
relative=float(abs(fd-exact)/max(abs(exact),1e-12))
u=rng.normal(size=D.shape);m=x/scale[:,None]
physical_error=float(abs(smooth(x)-(.5*np.sum((B@(scale[:,None]*m)-D)**2)+mu/2*np.sum((K@(W@m).ravel())**2))))
checks=dict(predictor_adjoint=adj,gradient_relative_error=relative,gradient_finite_difference=float(fd),gradient_analytic=exact,physical_coordinate_objective_difference=physical_error,
            all_pass=bool(adj<1e-11 and relative<1e-6 and physical_error<1e-10))
writejson(ROOT/'analysis/operator_checks.json',checks);print(checks)
assert checks['all_pass']
