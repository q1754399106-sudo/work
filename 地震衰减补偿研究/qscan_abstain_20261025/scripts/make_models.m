function make_models(which)
% Generate the geologies of set `which`: 'K' = the new confirmation models K1-K8 (writes models/
% and truth/, refuses to overwrite); 'G' = regenerate the archived development models G1-G8 into
% cache/regen_G (generator check against qscan_offset_20261025/models). Generator statement for
% statement as qscan_offset_20261025/scripts/make_dev.m; the helpers below are verbatim copies.
P=as_setup();
if strcmp(which,'K')
 specs=P.models;out=fullfile(P.root,'models');truth=fullfile(P.root,'truth');
 assert(~isfolder(out)&&~isfolder(truth),'Confirmation models exist; refusing to regenerate.');
 assert(isfile(fullfile(P.root,'logs','PRE_RUN_HASHES.json')),'Freeze PLAN_ZH.md first.');
else
 specs=g_specs();out=fullfile(P.root,'cache','regen_G');truth=out;
end
mkdir(out);if ~isfolder(truth),mkdir(truth);end;test_kennett();
cm=struct('dt',.004,'nt',451,'nfft',4096,'fref',100,'fpeak',32,'lambda_x',.001,'lambda_fixed',.001);cm.t=(0:cm.nt-1)'*cm.dt;
for s=1:numel(specs)
 sp=specs(s);rng(sp.seed,'twister');nx=99;base=fine_log(sp);D0=zeros(cm.nt,nx);T=D0;
 for j=1:nx,[v,rho,h,Q]=lateral(base,sp,j,nx);[D0(:,j),T(:,j)]=kennett_trace(v,rho,h,Q,cm,cm.fpeak);end
 config=cm;config.vp=sp.macro_vp;config.h=sp.macro_h;[v,~,h,Q]=lateral(base,sp,50,nx);config.Q=upscale_q(v,h,Q,sp.macro_h);config.solver_offset=0;
 clean=D0;target=T; %#ok<NASGU>
 save('-v7',fullfile(out,[sp.name '.mat']),'config','clean');save('-v7',fullfile(truth,[sp.name '_truth.mat']),'target','base');
 fprintf('%s: %d fine layers, upscaled Q %s, noise-free attenuated/elastic energy %.3f\n',sp.name,numel(base.v),mat2str(config.Q,4),norm(D0,'fro')/norm(T,'fro'));
end
end

function sp=g_specs()
% G1-G8 exactly as qscan_offset_20261025/scripts/os_setup.m.
sp=struct('name',{'G1','G2','G3','G4','G5','G6','G7','G8'}, ...
 'seed',{2026102601,2026102611,2026102621,2026102631,2026102641,2026102651,2026102661,2026102671}, ...
 'macro_vp',{[1850 2300 2750 3300],[1900 2350 2800 3350],[1780 2200 2820 3380],[1830 2380 2720 3260], ...
             [1870 2260 2790 3320],[1810 2340 2760 3290],[1760 2280 2840 3360],[1890 2310 2740 3240]}, ...
 'macro_h',{[240 300 350],[230 320 330],[250 290 360],[235 305 345],[255 315 325],[245 285 365],[225 335 335],[265 275 350]}, ...
 'macro_Q',{[50 70 95 130],[90 60 45 40],[40 80 50 110],[60 50 70 90],[30 45 60 80],[120 100 90 80],[70 40 100 60],[45 90 60 140]}, ...
 'carbonate',{false,true,false,true,false,true,false,true});
end

% The three helpers below are verbatim copies of bae_q_20261024/make_independent.m.
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
