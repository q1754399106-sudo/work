function [cm,D0,T]=load_model(P,sc)
% Public config and noise-free gather of the scene's model; the target T only when asked for
% (scoring). As load_scene.m of combo_bias_variance_20261024 for F1/F2 and run_offset.m of
% qscan_offset_20261025 for G1-G8.
if sc.family=='F'
 in=load(fullfile(P.in,'bae_q_20261024','inputs','independent',[sc.model '.mat']),'config');cm=in.config;
 t=load(fullfile(P.in,'bae_q_20261024','truth','independent',[sc.model '_truth.mat']),'tr','clean');D0=t.clean;
 if nargout>2,T=t.tr(1).target;end
else
 in=load(fullfile(P.in,'qscan_offset_20261025','models',[sc.model '.mat']));cm=in.config;D0=in.clean;
 if nargout>2,t=load(fullfile(P.in,'qscan_offset_20261025','truth',[sc.model '_truth.mat']),'target');T=t.target;end
end
cm.t=cm.t(:);
end
