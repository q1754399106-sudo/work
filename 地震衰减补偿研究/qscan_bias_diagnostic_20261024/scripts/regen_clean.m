function [D0,T]=regen_clean(P,set,model)
% Noise-free gather and target of a scene. bae/r3 truth files store them; the qscan_criterion
% truth stores only the fine log, so D0/T are recomputed with the same Kennett forward and
% the verbatim lateral() of make_independent/make_confirm (checked in check_inputs).
switch set
 case 'bae',t=load(fullfile(P.in,'bae_q_20261024','truth','independent',[model '_truth.mat']));D0=t.clean;T=t.tr(1).target;return
 case 'r3',t=load(fullfile(P.in,'round3_20261024','dev',[model '_truth.mat']));D0=t.clean;T=t.tr(1).target;return
end
t=load(fullfile(P.in,'qscan_criterion_20261024','truth','confirm',[model '_truth.mat']));
s=load(fullfile(P.in,'qscan_criterion_20261024','inputs','confirm',[model '.mat']));cm=s.config;cm.t=cm.t(:);
nx=size(s.obs(1).D,2);D0=zeros(cm.nt,nx);T=D0;
for j=1:nx,[v,rho,h,Q]=lateral(t.base,j,nx);[D0(:,j),T(:,j)]=kennett_trace(v,rho,h,Q,cm,cm.fpeak);end
end

function [v,rho,h,Q]=lateral(base,j,nx)
% Verbatim from bae_q_20261024/make_independent.m (the unused sp argument dropped).
x=(j-1)/(nx-1);stretch=[1+.06*x 1+.05*sin(2*pi*x) 1-.05*x 1];hs=base.h.*stretch(base.macro);
hs(base.gas)=base.h(base.gas)*min(max((90-j)/50,0),1);keep=hs>0;
v=[base.cap_v base.v(keep)];Q=[base.cap_Q base.Q(keep)];h=[base.cap_h hs(keep)];
v(end+1)=v(end)*1.05;Q(end+1)=Q(end);rho=.31*v.^.25;h=h(1:end);   % lower half-space
end
