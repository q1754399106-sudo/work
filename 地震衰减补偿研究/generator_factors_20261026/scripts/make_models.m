function make_models()
% Generator-factor models Z<s>_<abc> (PLAN_ZH.md "数据"): A reflectivity (0 dense fine log,
% 1 sparse blocky), B lateral (0 near-flat, 1 dips + two faults), C forward (0 exact Kennett,
% 1 the solver's own operator). fine_log, the B=0 branch of lateral and upscale_q are verbatim
% from round3_20261024/make_F34.m; kennett_trace, ray_absorption, attenuation_matrix and
% estimate_trace_snr come from the round3 snapshot in cache/runner. Solver inputs go to models/,
% truth to truth/; refuses to overwrite.
root=fileparts(fileparts(mfilename('fullpath')));runner=fullfile(root,'cache','runner');
addpath(fullfile(runner,'core'),fullfile(runner,'source'));maxNumCompThreads(4);
mdir=fullfile(root,'models');tdir=fullfile(root,'truth');
assert(~isfolder(mdir)&&~isfolder(tdir),'models exist; refusing to regenerate.');mkdir(mdir);mkdir(tdir);
cm=struct('dt',.004,'nt',451,'nfft',4096,'fref',100,'fpeak',32,'lambda_x',.001,'lambda_fixed',.001);cm.t=(0:cm.nt-1)'*cm.dt;
z=pi*cm.fpeak*(cm.t-cm.t');W=(1-2*z.^2).*exp(-z.^2);nx=99;
seeds=[2026102611 2026102621 2026102631];
for s=1:numel(seeds)
 rng(seeds(s),'twister');r=rand(1,10);
 sp=struct('seed',seeds(s),'carbonate',false, ...
  'macro_vp',[1750+60*r(1) 2300+100*r(2) 2700+250*r(3) 3300+150*r(4)], ...
  'macro_h',[210+40*r(5) 300+40*r(6) 320+40*r(7)], ...
  'macro_Q',[40+35*r(8) 40+45*r(9) 35+25*r(10) 85+40*rand]);
 rng(seeds(s)+1,'twister');dense=fine_log(sp);
 rng(seeds(s)+2,'twister');sparse=sparse_log(sp);
 cond=struct('name',{'Moderate','Strong'},'snr_db',{14,8},'seed',{sp.seed+101,sp.seed+102});
 for a=0:1
  if a,base=sparse;else,base=dense;end
  for b=0:1
   for c=0:1
    name=sprintf('Z%d_%d%d%d',s,a,b,c);config=cm;config.vp=sp.macro_vp;config.h=sp.macro_h;config.solver_offset=0;
    if c
     config.Q=sp.macro_Q;G=attenuation_matrix(ray_absorption(cm.t,0,config),cm);
    else
     [v,~,h,Q]=lateral_any(base,sp,50,nx,b);config.Q=upscale_q(v,h,Q,sp.macro_h);
    end
    D0=zeros(cm.nt,nx);T=D0;
    for j=1:nx
     [v,rho,h,Q]=lateral_any(base,sp,j,nx,b);
     if c
      T(:,j)=W*primaries(v,rho,h,cm);D0(:,j)=G*T(:,j);
     else
      [D0(:,j),T(:,j)]=kennett_trace(v,rho,h,Q,cm,cm.fpeak);
     end
    end
    obs=struct('name',{},'D',{},'sigma',{});tr=struct('name',{},'target',{});
    for k=1:numel(cond)
     rng(cond(k).seed,'twister');u=randn(cm.nt,nx);for i=2:cm.nt,u(i,:)=u(i,:)+.45*u(i-1,:);end
     u=.9*u+.1*mean(u,2);N=u*(norm(D0,'fro')/(10^(cond(k).snr_db/20)*norm(u,'fro')));D=D0+N;
     sigma=zeros(1,nx);for j=1:nx,ni=estimate_trace_snr(D(:,j),cm.t);sigma(j)=ni.sigma;end
     obs(k)=struct('name',cond(k).name,'D',D,'sigma',sigma);tr(k)=struct('name',cond(k).name,'target',T);
    end
    factors=struct('A_sparse',a,'B_structured',b,'C_solver_operator',c);clean=D0; %#ok<NASGU>
    save(fullfile(mdir,[name '.mat']),'config','obs','-v7');
    save(fullfile(tdir,[name '_truth.mat']),'tr','clean','base','sp','factors','-v7');
    fprintf('%s: %d layers, Q %s, clean rms %.3g\n',name,numel(base.v),mat2str(config.Q,4),rms(D0(:)));
   end
  end
 end
end
end

function base=fine_log(sp)
z=sp.macro_h(1);bounds=[0 cumsum(sp.macro_h) Inf];zmax=sum(sp.macro_h)+1800;hs=[];vs=[];Qs=[];macro=[];e=0;
while z<zmax
 th=3+12*rand;m=find(z>=bounds(1:end-1)&z<bounds(2:end),1);e=.75*e+.25*randn;jump=(rand<.12)*sign(randn)*.07;
 hs(end+1)=th;vs(end+1)=sp.macro_vp(m)*(1+.05*e+jump);Qs(end+1)=sp.macro_Q(m)*exp(.25*randn);macro(end+1)=m; %#ok<AGROW>
 z=z+th;
end
gas=find(macro==2,1)+8;vs(gas)=.86*vs(gas);Qs(gas)=18;hs(gas)=20;
if sp.carbonate,c=find(macro==2,1)+2;vs(c)=4200;Qs(c)=150;hs(c)=25;end
base=struct('h',hs,'v',vs,'Q',Qs,'macro',macro,'gas',gas,'cap_v',sp.macro_vp(1),'cap_Q',sp.macro_Q(1),'cap_h',sp.macro_h(1));
end

function base=sparse_log(sp)
% Thick blocks (60-160 m); velocity jumps weak (1.2-3%, p=.6) or strong (6-12%), pulled back
% toward the macro velocity once the deviation exceeds 5%. Q = macro Q. Gas as in fine_log.
z=sp.macro_h(1);bounds=[0 cumsum(sp.macro_h) Inf];zmax=sum(sp.macro_h)+1800;hs=[];vs=[];Qs=[];macro=[];d=0;
while z<zmax
 th=60+100*rand;m=find(z>=bounds(1:end-1)&z<bounds(2:end),1);
 if rand<.6,step=.012+.018*rand;else,step=.06+.06*rand;end
 if abs(d)>.05,sg=-sign(d);else,sg=1-2*(rand<.5);end
 d=d+sg*step;
 hs(end+1)=th;vs(end+1)=sp.macro_vp(m)*(1+d);Qs(end+1)=sp.macro_Q(m);macro(end+1)=m; %#ok<AGROW>
 z=z+th;
end
gas=find(macro==2,1)+1;vs(gas)=.86*vs(gas);Qs(gas)=18;hs(gas)=20;
base=struct('h',hs,'v',vs,'Q',Qs,'macro',macro,'gas',gas,'cap_v',sp.macro_vp(1),'cap_Q',sp.macro_Q(1),'cap_h',sp.macro_h(1));
end

function [v,rho,h,Q]=lateral_any(base,sp,j,nx,b) %#ok<INUSL>
% b=0: make_F34.m lateral, unchanged. b=1 adds a dipping top layer and two faults.
x=(j-1)/(nx-1);stretch=[1+.06*x 1+.05*sin(2*pi*x) 1-.05*x 1];cap=1;
if b,stretch(2)=stretch(2)*(1+.12*(j>=35));stretch(3)=stretch(3)*(1-.12*(j>=75));cap=1+.25*x+.03*sin(3*pi*x);end
hs=base.h.*stretch(base.macro);
hs(base.gas)=base.h(base.gas)*min(max((90-j)/50,0),1);keep=hs>0;
v=[base.cap_v base.v(keep)];Q=[base.cap_Q base.Q(keep)];h=[base.cap_h*cap hs(keep)];
v(end+1)=v(end)*1.05;Q(end+1)=Q(end);rho=.31*v.^.25;h=h(1:end);
end

function R=primaries(v,rho,h,cm)
% Normal-incidence primaries on the sample grid (no multiples, no transmission loss).
Z=rho.*v;r=(Z(2:end)-Z(1:end-1))./(Z(2:end)+Z(1:end-1));tau=cumsum(2*h./v(1:numel(h)));
ii=round(tau/cm.dt)+1;ok=ii<=cm.nt;R=accumarray(ii(ok)',r(ok)',[cm.nt 1]);
end

function q=upscale_q(v,h,Q,macro_h)
z=[0 cumsum(h)];tt=h./v(1:numel(h));bounds=[0 cumsum(macro_h) Inf];q=zeros(1,4);
for m=1:4,in=z(1:end-1)>=bounds(m)&z(1:end-1)<bounds(m+1);q(m)=sum(tt(in))/sum(tt(in)./Q(in));end
end
