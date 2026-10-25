function [cm,D0]=load_model(P,k)
% Public config and noise-free gather (used only to synthesise observations) of model k.
mo=P.models(k);
if strcmp(mo.family,'G')
 M=load(fullfile(P.in,'qscan_offset_20261025','models',[mo.name '.mat']));cm=M.config;D0=M.clean;
else
 in=load(fullfile(P.in,'bae_q_20261024','inputs','independent',[mo.name '.mat']),'config');cm=in.config;
 t=load(fullfile(P.in,'bae_q_20261024','truth','independent',[mo.name '_truth.mat']),'clean');D0=t.clean;
end
cm.t=cm.t(:);
end
