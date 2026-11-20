function R=ou_sigma(P,Qc,vmig,sloc,ou)
% Path-integral variance of the two-way t* under a local random field du (PLAN_ZH.md section 3).
%   t* = int a (1 + du) ds over the source leg and the specular-receiver leg (x_r = 2x - x_s,
%   clamped to the shot range as in two_way.m), a = 1/(Qc v), straight legs.
%   Cov(du(p),du(p')) = s(p) s(p') exp(-|dz|/Lz - |dx|/Lx).
% For L << leg length:  Var(leg) = int s^2 a^2 lc(theta) ds,  lc = 2/(|cos|/Lz + |sin|/Lx);
% cross term of the two legs: same depth fraction tau, horizontal separation d = tau |x_r - x_s|,
% Cov = int s s' a a' sqrt(lc lc') exp(-d/Lx) sqrt(ds ds')  (exact for coincident legs).
% Output (interior nodes lin, all shots): R.r = Var(t*)/t*^2 (single, nin x ns), R.L2 = log(1 + r).
% Computed on a grid decimated by ou.dec and interpolated (r is smooth).
if nargin<5,ou=P.ou;end
a=1./(Qc.*vmig);G=griddedInterpolant({P.z,P.x},a,'linear','nearest');
Gs=griddedInterpolant({P.z,P.x},sloc,'linear','nearest');
zi=P.z(P.in{1});xi=P.x(P.in{2});zd=zi(1:ou.dec:end);xd=xi(1:ou.dec:end);
if zd(end)<zi(end),zd(end+1)=zi(end);end
if xd(end)<xi(end),xd(end+1)=xi(end);end
[ZD,XD]=ndgrid(zd,xd);zp=ZD(:);xp=XD(:);np=numel(zp);
ns=numel(P.ixs);xs=P.x(P.ixs);zs=P.z(P.izs);n=ou.nseg;tau=((1:n)-.5)/n;
[ZI,XI]=ndgrid(zi,xi);R.r=zeros(numel(zi)*numel(xi),ns,'single');
for k=1:ns
 xr=min(max(2*xp-xs(k),xs(1)),xs(end));
 [ts,vs,ws,ls,dss,ss]=leg(xs(k)+0*xp);[tr,vr,wr,lr,dsr,sr]=leg(xr);
 d=abs(xr-xs(k))*tau;                                                  % np x n
 cv=sum(ss.*sr.*ws.*wr.*sqrt(ls.*lr).*exp(-d/ou.Lx),2).*sqrt(dss.*dsr);
 V=vs+vr+2*cv;t=ts+tr;r=V./max(t,eps).^2;
 F=griddedInterpolant({zd,xd},reshape(r,numel(zd),numel(xd)),'linear');R.r(:,k)=single(F(ZI(:),XI(:)));
end
R.L2=log1p(R.r);R.ou=ou;
 function [t,V,w,lc,ds,s]=leg(xe)
  L=hypot(zp-zs,xp-xe);ds=L/n;c=abs(zp-zs)./max(L,eps);sn=abs(xp-xe)./max(L,eps);
  lc=2./(c/ou.Lz+sn/ou.Lx);
  zq=zp+(zs-zp)*tau;xq=xp+(xe-xp)*tau;w=G(zq,xq);s=Gs(zq,xq);
  t=sum(w,2).*ds;V=sum(s.^2.*w.^2,2).*lc.*ds;
 end
end
