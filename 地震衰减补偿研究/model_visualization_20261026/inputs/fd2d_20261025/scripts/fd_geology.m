function Gm=fd_geology(sp,P)
% 2-D geology on the finite-difference grid from the fine-layer generator of
% bae_q_20261024/make_independent.m (fine_log, lateral and upscale_q below are verbatim copies
% from qscan_offset_20261025/scripts/make_dev.m). Trace j sits at x = (j-1)*trace_dx; the column
% at any x is lateral() at the fractional trace index (clamped to 1..99 outside the line). A
% vertical fault at sp.fault_trace+0.5 drops everything below the cap by sp.throw on its
% right-hand side (hanging wall). Each grid cell holds the normal-incidence (Backus) average of
% the fine layers it covers: harmonic mean of the complex modulus rho*v^2*(1+i/Q), arithmetic
% mean of density.
% Gm.x, Gm.z   grid coordinates (m); z = 0 is the source/receiver datum
% Gm.v, Gm.rho, Gm.Q   nz-by-nx (PML rows/columns included, properties extended)
% Gm.config    public solver configuration (macro model, upscaled Q at trace 50)
fd=P.fd;rng(sp.seed,'twister');base=fine_log(sp);nx99=fd.ntraces;
x=(-fd.half_window-fd.npml*fd.dx:fd.dx:(nx99-1)*fd.trace_dx+fd.half_window+fd.npml*fd.dx);
z=(-fd.npml*fd.dx:fd.dx:fd.zmax+fd.npml*fd.dx)';
v=zeros(numel(z),numel(x));rho=v;Q=v;
zs=(-fd.dx/2+fd.sub/2:fd.sub:fd.dx/2)';                       % sub-samples within a cell
xf=(sp.fault_trace-.5)*fd.trace_dx;
for ix=1:numel(x)
 j=min(max(1+x(ix)/fd.trace_dx,1),nx99);[cv,crho,ch,cQ]=lateral(base,sp,j,nx99);
 if x(ix)>xf,ch(1)=ch(1)+sp.throw;end
 zi=[0 cumsum(ch)];
 zz=max(z'+zs,0);k=discretize(zz,[-Inf zi(2:end) Inf]);     % layer index of each sub-sample
 Mc=crho(k).*cv(k).^2.*(1+1i./cQ(k));Me=1./mean(1./Mc,1);re=mean(crho(k),1);
 v(:,ix)=sqrt(real(Me)./re);rho(:,ix)=re;Q(:,ix)=real(Me)./imag(Me);
end
cm=P.cm;config=cm;config.vp=sp.macro_vp;config.h=sp.macro_h;
[cv,~,ch,cQ]=lateral(base,sp,50,nx99);config.Q=upscale_q(cv,ch,cQ,sp.macro_h);config.solver_offset=0;
[cv,crho,ch,cQ]=lateral(base,sp,50,nx99);
Gm=struct('x',x,'z',z,'v',v,'rho',rho,'Q',Q,'base',base,'config',config,'trace_x',(0:nx99-1)*fd.trace_dx, ...
 'cap_v',base.cap_v,'cap_Q',base.cap_Q,'cap_rho',.31*base.cap_v^.25,'col50',struct('v',cv,'rho',crho,'h',ch,'Q',cQ));
end

% The three helpers below are verbatim copies of qscan_offset_20261025/scripts/make_dev.m.
function base=fine_log(sp)
% Fine layers below a homogeneous 0.23-0.24 km cap (keeps the 0.02-0.25 s noise window signal-free).
z=sp.macro_h(1);bounds=[0 cumsum(sp.macro_h) Inf];zmax=sum(sp.macro_h)+1800;
hs=[];vs=[];Qs=[];macro=[];e=0;
while z<zmax
 th=3+12*rand;m=find(z>=bounds(1:end-1)&z<bounds(2:end),1);
 e=.75*e+.25*randn;jump=(rand<.12)*sign(randn)*.07;
 hs(end+1)=th;vs(end+1)=sp.macro_vp(m)*(1+.05*e+jump);Qs(end+1)=sp.macro_Q(m)*exp(.25*randn);macro(end+1)=m; %#ok<AGROW>
 z=z+th;
end
gas=find(macro==2,1)+8;vs(gas)=.86*vs(gas);Qs(gas)=18;hs(gas)=20;
if sp.carbonate,c=find(macro==2,1)+2;vs(c)=4200;Qs(c)=150;hs(c)=25;end
base=struct('h',hs,'v',vs,'Q',Qs,'macro',macro,'gas',gas,'cap_v',sp.macro_vp(1),'cap_Q',sp.macro_Q(1),'cap_h',sp.macro_h(1));
end

function [v,rho,h,Q]=lateral(base,sp,j,nx)
% Smooth thickness stretch per macro unit (dip) and a gas-sand wedge thinning to 0 by trace 90.
x=(j-1)/(nx-1);stretch=[1+.06*x 1+.05*sin(2*pi*x) 1-.05*x 1];hs=base.h.*stretch(base.macro);
hs(base.gas)=base.h(base.gas)*min(max((90-j)/50,0),1);keep=hs>0;
v=[base.cap_v base.v(keep)];Q=[base.cap_Q base.Q(keep)];h=[base.cap_h hs(keep)];
v(end+1)=v(end)*1.05;Q(end+1)=Q(end);rho=.31*v.^.25;h=h(1:end);   % lower half-space
end

function q=upscale_q(v,h,Q,macro_h)
% Traveltime-weighted harmonic mean of fine Q in each public macro interval.
z=[0 cumsum(h)];tt=h./v(1:numel(h));bounds=[0 cumsum(macro_h) Inf];q=zeros(1,4);
for m=1:4
 in=z(1:end-1)>=bounds(m)&z(1:end-1)<bounds(m+1);q(m)=sum(tt(in))/sum(tt(in)./Q(in));
end
end
