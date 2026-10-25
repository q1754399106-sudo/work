function rec=fd2d_visco(M,q,opt)
% 2-D viscoacoustic finite-difference modelling, pressure-velocity staggered grid.
%   8th-order space, 2nd-order (leapfrog) time; GSLS attenuation with three mechanisms
%   (gsls_design, Crank-Nicolson memory update); C-PML on all four sides (no free surface).
%   Single precision throughout.
% M.v, M.rho, M.Q  nz-by-nx: phase velocity at opt.fref (m/s), density, Q (Inf = elastic)
% q                source injection rate, one value per time step (added to dp/dt at the source)
% opt.dx, opt.dt, opt.nt, opt.fref, opt.npml, opt.src [iz ix], opt.rec [iz ix; ...]
% rec              nt-by-nrec pressure at the receivers after each step (time n*dt)
[nz,nx]=size(M.v);dx=opt.dx;dt=opt.dt;np=opt.npml;
c=[1225/1024 -245/3072 49/5120 -5/7168];w=single([-fliplr(c) c]/dx);
Kxp=fliplr(w);Kxm=fliplr([w 0]);Kzp=Kxp.';Kzm=Kxm.';     % forward (i+1/2) and backward (i-1/2) staggering
if any(isfinite(M.Q(:)))
 [uq,~,iq]=unique(M.Q(:));G=gsls_design(uq,opt.fref);tau=reshape(G.tau(iq),nz,nx);vfac=reshape(G.vfac(iq),nz,nx);
else
 tau=zeros(nz,nx);vfac=ones(nz,nx);G=struct('ts',[]);
end
MR=M.rho.*M.v.^2.*vfac;MU=MR.*(1+numel(G.ts)*tau);
cmax=max(sqrt(MU(:)./M.rho(:)));assert(cmax*dt/dx*sqrt(2)*sum(abs(c))<.95,'CFL violated.');
bx=1./(.5*(M.rho+M.rho(:,[2:end end])));bz=1./(.5*(M.rho+M.rho([2:end end],:)));   % buoyancy at vx, vz nodes
% C-PML only on the boundary strips: columns sx (x layers), rows sz (z layers).
sx=[1:np nx-np+1:nx];sz=[1:np nz-np+1:nz];
[ax_i,bx_i,ax_h,bx_h]=cpml(nx,np,dx,dt,cmax);[az_i,bz_i,az_h,bz_h]=cpml(nz,np,dx,dt,cmax);
ax_i=single(ax_i(sx)');bx_i=single(bx_i(sx)');ax_h=single(ax_h(sx)');bx_h=single(bx_h(sx)');
az_i=single(az_i(sz));bz_i=single(bz_i(sz));az_h=single(az_h(sz));bz_h=single(bz_h(sz));
p=zeros(nz,nx,'single');vx=p;vz=p;ppx=zeros(nz,2*np,'single');pvx=ppx;ppz=zeros(2*np,nx,'single');pvz=ppz;
dtbx=single(dt*bx);dtbz=single(dt*bz);dtMU=single(dt*MU);
L=numel(G.ts);z=cell(1,L);A=zeros(1,L,'single');C=A;hMRt=single(.5*dt*MR.*tau);
for l=1:L,e=dt/(2*G.ts(l));A(l)=(1-e)/(1+e);C(l)=2*e/(1+e);z{l}=p;end   % zeta_l = r_l/(M_R*tau)
isrc=sub2ind([nz nx],opt.src(1),opt.src(2));irec=sub2ind([nz nx],opt.rec(:,1),opt.rec(:,2));
rec=zeros(opt.nt,numel(irec),'single');qs=single(dt*q(:)/dx^2);
for n=1:opt.nt
 d=conv2(p,Kxp,'same');ppx=bx_h.*ppx+ax_h.*d(:,sx);d(:,sx)=d(:,sx)+ppx;vx=vx-dtbx.*d;
 d=conv2(p,Kzp,'same');ppz=bz_h.*ppz+az_h.*d(sz,:);d(sz,:)=d(sz,:)+ppz;vz=vz-dtbz.*d;
 div=conv2(vx,Kxm,'same');pvx=bx_i.*pvx+ax_i.*div(:,sx);div(:,sx)=div(:,sx)+pvx;
 d=conv2(vz,Kzm,'same');pvz=bz_i.*pvz+az_i.*d(sz,:);d(sz,:)=d(sz,:)+pvz;div=div+d;
 if L
  s=0;
  for l=1:L,zn=A(l)*z{l}-C(l)*div;s=s+zn+z{l};z{l}=zn;end
  p=p-dtMU.*div-hMRt.*s;
 else
  p=p-dtMU.*div;
 end
 p(isrc)=p(isrc)+qs(n);rec(n,:)=p(irec);
end
rec=double(rec);
end

function [a_i,b_i,a_h,b_h]=cpml(n,np,dx,dt,cmax)
% C-PML coefficients (Komatitsch & Martin 2007, kappa = 1) at integer and half nodes.
Lp=np*dx;d0=-3*cmax*log(1e-4)/(2*Lp);amax=pi*30;
x=(1:n)';xh=x+.5;
dist=@(y) min(max(max(np+1-y,y-(n-np)),0),np)*dx;   % distance into the layer
[a_i,b_i]=coef(dist(x));[a_h,b_h]=coef(dist(xh));
 function [a,b]=coef(s)
  dd=d0*(s/Lp).^2;al=amax*(1-s/Lp);al(s==0)=0;b=exp(-(dd+al)*dt);a=zeros(size(s));k=dd>0;a(k)=dd(k)./(dd(k)+al(k)).*(b(k)-1);
 end
end
