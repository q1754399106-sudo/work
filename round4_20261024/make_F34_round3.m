function make_F34()
% New independent confirmation scenes F3/F4: bae_q_20261024/make_independent.m
% unchanged except for new seeds and macro models (registered in
% PLAN_ADDENDUM_GUARD.json before generation). Exact 1D fine-layered
% viscoacoustic (Kennett) forward, offset 0; solver inputs and truth apart.
root=fileparts(mfilename('fullpath'));addpath(fullfile(root,'core'),fullfile(root,'source'));maxNumCompThreads(4);
for n={'F3','F4'},assert(~isfile(fullfile(root,'dev',[n{1} '.mat'])),'F3/F4 exist; refusing to regenerate.');end
cm=struct('dt',.004,'nt',451,'nfft',4096,'fref',100,'fpeak',32,'lambda_x',.001,'lambda_fixed',.001);cm.t=(0:cm.nt-1)'*cm.dt;
specs=struct('name',{'F3','F4'},'seed',{2026102471,2026102481}, ...
 'macro_vp',{[1780 2350 2750 3350],[1750 2350 2900 3400]},'macro_h',{[235 310 340],[220 320 330]}, ...
 'macro_Q',{[70 45 38 90],[45 80 55 120]},'carbonate',{false,true}, ...
 'description',{'fine-layered log, gas sand wedge, Q decreasing with depth then recovering', ...
  'fine-layered log, high-velocity carbonate, gas sand wedge, Q reversal in layer 3'});
for s=1:numel(specs)
 sp=specs(s);rng(sp.seed,'twister');nx=99;base=fine_log(sp);D0=zeros(cm.nt,nx);T=D0;
 for j=1:nx,[v,rho,h,Q]=lateral(base,sp,j,nx);[D0(:,j),T(:,j)]=kennett_trace(v,rho,h,Q,cm,cm.fpeak);end
 config=cm;config.vp=sp.macro_vp;config.h=sp.macro_h;[v,~,h,Q]=lateral(base,sp,50,nx);config.Q=upscale_q(v,h,Q,sp.macro_h);config.solver_offset=0;
 cond=struct('name',{'Moderate','Strong'},'snr_db',{14,8},'seed',{sp.seed+1,sp.seed+2});
 obs=struct('name',{},'D',{},'sigma',{});tr=struct('name',{},'target',{});
 for k=1:numel(cond)
  rng(cond(k).seed,'twister');u=randn(cm.nt,nx);for i=2:cm.nt,u(i,:)=u(i,:)+.45*u(i-1,:);end
  u=.9*u+.1*mean(u,2);N=u*(norm(D0,'fro')/(10^(cond(k).snr_db/20)*norm(u,'fro')));D=D0+N;
  sigma=zeros(1,nx);for j=1:nx,ni=estimate_trace_snr(D(:,j),cm.t);sigma(j)=ni.sigma;end
  obs(k)=struct('name',cond(k).name,'D',D,'sigma',sigma);tr(k)=struct('name',cond(k).name,'target',T);
 end
 save(fullfile(root,'dev',[sp.name '.mat']),'config','obs','-v7');description=sp.description;clean=D0; %#ok<NASGU>
 save(fullfile(root,'dev',[sp.name '_truth.mat']),'tr','clean','base','description','-v7');
 fprintf('%s: %d fine layers, upscaled Q %s\n',sp.name,numel(base.v),mat2str(config.Q,4));
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

function [v,rho,h,Q]=lateral(base,sp,j,nx) %#ok<INUSL>
x=(j-1)/(nx-1);stretch=[1+.06*x 1+.05*sin(2*pi*x) 1-.05*x 1];hs=base.h.*stretch(base.macro);
hs(base.gas)=base.h(base.gas)*min(max((90-j)/50,0),1);keep=hs>0;
v=[base.cap_v base.v(keep)];Q=[base.cap_Q base.Q(keep)];h=[base.cap_h hs(keep)];
v(end+1)=v(end)*1.05;Q(end+1)=Q(end);rho=.31*v.^.25;h=h(1:end);
end

function q=upscale_q(v,h,Q,macro_h)
z=[0 cumsum(h)];tt=h./v(1:numel(h));bounds=[0 cumsum(macro_h) Inf];q=zeros(1,4);
for m=1:4,in=z(1:end-1)>=bounds(m)&z(1:end-1)<bounds(m+1);q(m)=sum(tt(in))/sum(tt(in)./Q(in));end
end
