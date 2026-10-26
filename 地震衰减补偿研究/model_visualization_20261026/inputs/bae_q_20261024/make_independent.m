function make_independent()
% Experiment 3: confirmation scenes from an independent forward model
% (exact 1D fine-layered viscoacoustic wave solution, kennett_trace) instead of
% the ray-absorption convolutional model used by every solver dictionary.
% Solver inputs and truth are written separately; generated once.
root=fileparts(mfilename('fullpath'));addpath(fullfile(root,'core'),fullfile(root,'source'));maxNumCompThreads(4);
out=fullfile(root,'inputs','independent');truth=fullfile(root,'truth','independent');
assert(~isfolder(out)&&~isfolder(truth),'Independent scenes already exist; refusing to regenerate.');mkdir(out);mkdir(truth);
test_kennett();
cm=struct('dt',.004,'nt',451,'nfft',4096,'fref',100,'fpeak',32,'lambda_x',.001,'lambda_fixed',.001);cm.t=(0:cm.nt-1)'*cm.dt;
specs=struct('name',{'F1','F2'},'seed',{2026102451,2026102461}, ...
 'macro_vp',{[1800 2300 2800 3300],[1850 2250 2700 3200]},'macro_h',{[230 300 350],[240 290 360]}, ...
 'macro_Q',{[40 60 85 110],[35 55 42 95]},'carbonate',{false,true}, ...
 'description',{'fine-layered log, gas sand wedge, Q trend increasing', ...
  'fine-layered log, high-velocity carbonate (strong internal multiples), gas sand wedge, Q reversal in layer 3'});
for s=1:numel(specs)
 sp=specs(s);rng(sp.seed,'twister');nx=99;base=fine_log(sp);
 D0=zeros(cm.nt,nx);T=D0;
 for j=1:nx
  [v,rho,h,Q]=lateral(base,sp,j,nx);[D0(:,j),T(:,j)]=kennett_trace(v,rho,h,Q,cm,cm.fpeak);
 end
 config=cm;config.vp=sp.macro_vp;config.h=sp.macro_h;
 [v,~,h,Q]=lateral(base,sp,50,nx);config.Q=upscale_q(v,h,Q,sp.macro_h);config.solver_offset=0;
 cond=struct('name',{'Moderate','Strong'},'snr_db',{14,8},'seed',{sp.seed+1,sp.seed+2});
 obs=struct('name',{},'D',{},'sigma',{});tr=struct('name',{},'target',{});
 for k=1:numel(cond)
  rng(cond(k).seed,'twister');u=randn(cm.nt,nx);for i=2:cm.nt,u(i,:)=u(i,:)+.45*u(i-1,:);end
  u=.9*u+.1*mean(u,2);N=u*(norm(D0,'fro')/(10^(cond(k).snr_db/20)*norm(u,'fro')));D=D0+N;
  sigma=zeros(1,nx);for j=1:nx,ni=estimate_trace_snr(D(:,j),cm.t);sigma(j)=ni.sigma;end
  obs(k)=struct('name',cond(k).name,'D',D,'sigma',sigma);tr(k)=struct('name',cond(k).name,'target',T);
 end
 save(fullfile(out,[sp.name '.mat']),'config','obs','-v7');
 description=sp.description;clean=D0; %#ok<NASGU>
 save(fullfile(truth,[sp.name '_truth.mat']),'tr','clean','base','description','-v7');
 fprintf('%s: %d fine layers, upscaled Q %s, noise-free attenuated/elastic energy %.3f\n',sp.name,numel(base.v),mat2str(config.Q,4),norm(D0,'fro')/norm(T,'fro'));
end
uq_writejson(fullfile(out,'README.json'),struct('forward','exact 1D normal-incidence viscoacoustic layered response (Kjartansson constant Q), all internal multiples and transmission losses, off-grid interfaces', ...
 'target','same stack with Q=Inf (elastic), Ricker 32 Hz zero phase','solver_Q','traveltime-weighted harmonic upscaling of fine Q into the four public macro intervals at trace 50','offset',0));
end

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
