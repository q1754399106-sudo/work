function T=uq_tab(f,xg,sig,fref)
% Expectations over an uncertain two-way t* (anchor at the surface): x = x_hat exp(sig e), e ~ N(0,1),
% 41-point Gauss-Hermite quadrature in the log domain (as qunknown_comp_20261109/scripts/uqw.py with
% dx_hat = x_hat). Forward model G = exp[(-pi f + i k) x], k = 2 f ln(f/fref) (fft convention).
%   T.A = E[G*]/E|G|^2,  T.C = 1/E|G|^2,  T.P = E[exp(-i k x)]   on the grid xg (row vectors).
% The operator for noise term rho = s + E|G|^2/snr is
%   H = (A + gamma (s P C + P/snr)) / (1 + s C + 1/snr)
% sig = 0, gamma = 1, 1/snr = 0 gives the stabilised inverse Q exp(-ikx)(b+s)/(b^2+s), b = exp(-pi f x).
k=2*f*log(f/fref);xg=xg(:)';
if sig==0
 T.A=exp((pi*f-1i*k)*xg);T.C=exp(2*pi*f*xg);T.P=exp(-1i*k*xg);return;
end
[e,w]=gh41();lw=log(w/sum(w));xs=xg.*exp(sig*e);             % (41 x n)
l2=-2*pi*f*xs+lw;mx=max(l2,[],1);L2=mx+log(sum(exp(l2-mx),1)); % log E|G|^2
T.A=sum(exp((-pi*f-1i*k)*xs+lw-L2),1);T.C=exp(-L2);T.P=sum(exp(-1i*k*xs+lw),1);
end

function [x,w]=gh41()
% probabilists' Gauss-Hermite nodes and weights (Golub-Welsch), as numpy hermegauss(41)
n=41;b=sqrt(1:n-1);J=diag(b,1)+diag(b,-1);[V,D]=eig(J);x=diag(D);w=sqrt(2*pi)*V(1,:)'.^2;
[x,i]=sort(x);w=w(i);
end
