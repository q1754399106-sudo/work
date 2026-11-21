function T=va_tab(f,xg,sig,fv)
% VA-UQ expectations on the grid xg of the anchored two-way Delta t* (PLAN_ZH.md section 3).
% Forward model G = exp[(-pi f + i 2 f ln(f/f_v)) x], x = u * x_hat, u lognormal (sigma sig, 41-point
% Gauss-Hermite in the log domain, as uq_tab.m), ln f_v uniform on [fv(1), fv(2)]. The f_v expectation is
% closed form: E_fv[exp(-i 2 f x ln(f/f_v))] = exp(-i kbar x) sinc(f x Lam), kbar = 2 f ln(f/fbar),
% fbar = sqrt(fv(1) fv(2)), Lam = ln(fv(2)/fv(1)), sinc(y) = sin(y)/y. The sinc factor depends on u and
% is therefore applied inside the u quadrature (exact for the stated model). |G| does not depend on f_v.
%   T.A = E[G*]/E|G|^2,  T.C = 1/E|G|^2,  T.P = E[exp(-i k x)]  (same fields and use as uq_tab.m)
% fv scalar: f_v known (no sinc) = uq_tab(f,xg,sig,fv).
xg=xg(:)';
if isscalar(fv),fb=fv;Lam=0;else,fb=sqrt(fv(1)*fv(2));Lam=log(fv(2)/fv(1));end
k=2*f*log(f/fb);
if sig==0,e=0;lw=0;else,[e,w]=gh41();lw=log(w/sum(w));end
xs=xg.*exp(sig*e);                                            % (nq x n)
y=f*xs*Lam;sc=ones(size(y));nz=y~=0;sc(nz)=sin(y(nz))./y(nz);
l2=-2*pi*f*xs+lw;mx=max(l2,[],1);L2=mx+log(sum(exp(l2-mx),1)); % log E|G|^2
T.A=sum(exp((-pi*f-1i*k)*xs+lw-L2).*sc,1);T.C=exp(-L2);T.P=sum(exp(-1i*k*xs+lw).*sc,1);
end

function [x,w]=gh41()
% probabilists' Gauss-Hermite nodes and weights (Golub-Welsch), identical to uq_tab.m
n=41;b=sqrt(1:n-1);J=diag(b,1)+diag(b,-1);[V,D]=eig(J);x=diag(D);w=sqrt(2*pi)*V(1,:)'.^2;
[x,i]=sort(x);w=w(i);
end
