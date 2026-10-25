function make_unseen()
% Never-used confirmation geologies U1-U4: bae_q_20261024/make_independent.m unchanged except for
% new seeds, new macro models and four Frobenius SNR levels (8/11/14/17 dB). Exact 1D
% fine-layered viscoacoustic (Kennett) forward, offset 0. Solver inputs (public config and the
% diagnosed observations r=1-3 per level) and truth (target, noise-free gather, fine log) apart.
P=bw_setup();maxNumCompThreads(4);out=fullfile(P.root,'inputs');truth=fullfile(P.root,'truth');
assert(~isfolder(out)&&~isfolder(truth),'Unseen scenes exist; refusing to regenerate.');mkdir(out);mkdir(truth);test_kennett();
cm=struct('dt',.004,'nt',451,'nfft',4096,'fref',100,'fpeak',32,'lambda_x',.001,'lambda_fixed',.001);cm.t=(0:cm.nt-1)'*cm.dt;
specs=struct('name',{'U1','U2','U3','U4'},'seed',{2026102501,2026102511,2026102521,2026102531}, ...
 'macro_vp',{[1820 2280 2780 3300],[1880 2320 2760 3280],[1760 2250 2850 3400],[1800 2400 2700 3250]}, ...
 'macro_h',{[225 310 340],[245 295 355],[230 330 320],[260 280 360]}, ...
 'macro_Q',{[45 65 90 120],[80 55 40 35],[35 70 45 100],[55 45 75 95]},'carbonate',{false,true,false,true}, ...
 'description',{'fine-layered, Q increasing with depth, gas sand wedge','fine-layered, carbonate multiples, Q decreasing with depth, gas sand wedge', ...
  'fine-layered, low shallow Q, Q reversal in layer 3, gas sand wedge','fine-layered, carbonate multiples, moderate Q, gas sand wedge'});
snr=[8 11 14 17];
for s=1:numel(specs)
 sp=specs(s);rng(sp.seed,'twister');nx=99;base=fine_log(sp);D0=zeros(cm.nt,nx);T=D0;
 for j=1:nx,[v,rho,h,Q]=lateral(base,sp,j,nx);[D0(:,j),T(:,j)]=kennett_trace(v,rho,h,Q,cm,cm.fpeak);end
 config=cm;config.vp=sp.macro_vp;config.h=sp.macro_h;[v,~,h,Q]=lateral(base,sp,50,nx);config.Q=upscale_q(v,h,Q,sp.macro_h);config.solver_offset=0;
 obs=struct('name',{},'snr_db',{},'replicate',{},'seed',{},'D',{},'sigma',{});
 for level=1:numel(snr)
  for r=1:3
   seed=1200000000+10000000*s+100000*level+r;D=D0+noise_replica(D0,snr(level),seed);
   obs(end+1)=struct('name',sprintf('SNR%02d',snr(level)),'snr_db',snr(level),'replicate',r,'seed',seed,'D',D,'sigma',estimate_sigmas(D,cm.t)); %#ok<AGROW>
  end
 end
 save(fullfile(out,[sp.name '.mat']),'config','obs','-v7');description=sp.description;clean=D0;target=T; %#ok<NASGU>
 save(fullfile(truth,[sp.name '_truth.mat']),'target','clean','base','description','-v7');
 fprintf('%s: %d fine layers, upscaled Q %s, noise-free attenuated/elastic energy %.3f\n',sp.name,numel(base.v),mat2str(config.Q,4),norm(D0,'fro')/norm(T,'fro'));
end
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
