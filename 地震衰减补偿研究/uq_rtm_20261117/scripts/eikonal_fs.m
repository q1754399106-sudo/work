function T=eikonal_fs(s,h,iz,ix)
% First-arrival traveltime |grad T| = s from the node (iz,ix): first-order Godunov fast sweeping
% (Zhao 2005); four orderings, repeated until the largest change is below 1e-7 s.
[nz,nx]=size(s);T=inf(nz,nx);
[J,I]=meshgrid(max(1,ix-2):min(nx,ix+2),max(1,iz-2):min(nz,iz+2));
for k=1:numel(I),T(I(k),J(k))=s(iz,ix)*h*hypot(I(k)-iz,J(k)-ix);end
fz=T;fixed=isfinite(fz);
ord={{1:nz,1:nx},{nz:-1:1,1:nx},{1:nz,nx:-1:1},{nz:-1:1,nx:-1:1}};
for it=1:8
 dmax=0;
 for o=1:4
  zz=ord{o}{1};xx=ord{o}{2};
  for j=xx
   for i=zz
    if fixed(i,j),continue;end
    if i>1,a=T(i-1,j);else,a=inf;end
    if i<nz&&T(i+1,j)<a,a=T(i+1,j);end
    if j>1,b=T(i,j-1);else,b=inf;end
    if j<nx&&T(i,j+1)<b,b=T(i,j+1);end
    if isinf(a)&&isinf(b),continue;end
    f=s(i,j)*h;
    if abs(a-b)>=f,t=min(a,b)+f;else,t=.5*(a+b+sqrt(2*f*f-(a-b)^2));end
    if t<T(i,j),d=T(i,j)-t;if d>dmax,dmax=d;end;T(i,j)=t;end
   end
  end
 end
 if dmax<1e-7&&it>1,break;end
end
end
