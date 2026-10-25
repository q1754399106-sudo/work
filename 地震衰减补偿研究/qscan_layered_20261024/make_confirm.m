function make_confirm()
% New confirmation scenes F5/F6 for the depth-segmented study, same exact 1D
% viscoacoustic forward as bae_q_20261024/make_independent but new geology/seeds.
root=fileparts(mfilename('fullpath'));up=fullfile(root,'upstream');
addpath(fullfile(root,'core'),fullfile(up,'bae_core'),fullfile(up,'bae_source'));maxNumCompThreads(4);
out=fullfile(root,'inputs','confirm');truth=fullfile(root,'truth','confirm');
assert(~isfolder(out)&&~isfolder(truth),'Confirmation scenes exist; refusing to regenerate.');mkdir(out);mkdir(truth);test_kennett();
cm=struct('dt',.004,'nt',451,'nfft',4096,'fref',100,'fpeak',32,'lambda_x',.001,'lambda_fixed',.001);cm.t=(0:cm.nt-1)'*cm.dt;
specs=struct('name',{'F5','F6'},'seed',{2026102491,2026102501},'macro_vp',{[1820 2280 2650 3150],[1780 2150 2900 3400]}, ...
 'macro_h',{[260 300 340],[210 330 360]},'macro_Q',{[45 30 80 60],[38 65 50 110]},'carbonate',{true,false}, ...
 'description',{'fine-layered, low-Q layer 2 and Q reversal in layer 4, carbonate multiples, gas sand wedge','fine-layered, alternating Q with depth, gas sand wedge'});
for s=1:numel(specs)
 sp=specs(s);rng(sp.seed,'twister');nx=99;base=fine_log(sp);D0=zeros(cm.nt,nx);T=D0;
 for j=1:nx,[v,rho,h,Q]=lateral(base,sp,j,nx);[D0(:,j),T(:,j)]=kennett_trace(v,rho,h,Q,cm,cm.fpeak);end
 config=cm;config.vp=sp.macro_vp;config.h=sp.macro_h;[v,~,h,Q]=lateral(base,sp,50,nx);config.Q=upscale_q(v,h,Q,sp.macro_h);config.solver_offset=0;
 cond=struct('name',{'Moderate','Strong'},'snr_db',{14,8},'seed',{sp.seed+1,sp.seed+2});obs=struct('name',{},'D',{},'sigma',{});tr=struct('name',{},'target',{});
 for k=1:2
  rng(cond(k).seed,'twister');u=randn(cm.nt,nx);for i=2:cm.nt,u(i,:)=u(i,:)+.45*u(i-1,:);end
  u=.9*u+.1*mean(u,2);N=u*(norm(D0,'fro')/(10^(cond(k).snr_db/20)*norm(u,'fro')));D=D0+N;
  sigma=zeros(1,nx);for j=1:nx,ni=estimate_trace_snr(D(:,j),cm.t);sigma(j)=ni.sigma;end
  obs(k)=struct('name',cond(k).name,'D',D,'sigma',sigma);tr(k)=struct('name',cond(k).name,'target',T);
 end
 save(fullfile(out,[sp.name '.mat']),'config','obs','-v7');description=sp.description; %#ok<NASGU>
 save(fullfile(truth,[sp.name '_truth.mat']),'tr','base','description','-v7');
 fprintf('%s: %d fine layers, upscaled Q %s\n',sp.name,numel(base.v),mat2str(config.Q,4));
end
end
% The three helpers below are verbatim copies of bae_q_20261024/make_independent.m.
function base=fine_log(sp)
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
x=(j-1)/(nx-1);stretch=[1+.06*x 1+.05*sin(2*pi*x) 1-.05*x 1];hs=base.h.*stretch(base.macro);
hs(base.gas)=base.h(base.gas)*min(max((90-j)/50,0),1);keep=hs>0;
v=[base.cap_v base.v(keep)];Q=[base.cap_Q base.Q(keep)];h=[base.cap_h hs(keep)];
v(end+1)=v(end)*1.05;Q(end+1)=Q(end);rho=.31*v.^.25;h=h(1:end);
end
function q=upscale_q(v,h,Q,macro_h)
z=[0 cumsum(h)];tt=h./v(1:numel(h));bounds=[0 cumsum(macro_h) Inf];q=zeros(1,4);
for m=1:4
 in=z(1:end-1)>=bounds(m)&z(1:end-1)<bounds(m+1);q(m)=sum(tt(in))/sum(tt(in)./Q(in));
end
end
