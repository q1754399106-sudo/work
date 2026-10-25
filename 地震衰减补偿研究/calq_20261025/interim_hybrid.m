function interim_hybrid()
% Interim switch after GhostCal failed the development gate (layer 4 collapsed low):
% Hybrid = GhostCal layers 1-3 + current layer 4 (no new tuning), plus a truth-side
% split of which layer error matters (diagnosis only).
root=fileparts(mfilename('fullpath'));addpath(fullfile(root,'core'),fullfile(root,'source'));maxNumCompThreads(4);
dev=jsondecode(fileread(fullfile(root,'analysis','dev_ghost.json')));out=struct();
for v={'V3','V5'}
 in=load(fullfile(root,'dev',[v{1} '.mat']));TT=load(fullfile(root,'dev',[v{1} '_truth.mat']));cur=in.config;cur.t=cur.t(:);
 qc=reshape(cur.Q,1,[]);qg=reshape(dev.(v{1}).archived_noisy_Q,1,[]);qt=reshape(dev.(v{1}).true_Q,1,[]);
 configs={'current',qc;'Hybrid',[qg(1:3) qc(4)];'true123_current4',[qt(1:3) qc(4)];'current123_true4',[qc(1:3) qt(4)];'true_all',qt};rows={};
 for k=1:2
  D=in.obs(k).D;s=reshape(in.obs(k).sigma,1,[]);T=struct('target',TT.tr(k).target,'masks',struct('low_energy',TT.regions.low_energy(:,:,k),'deep',TT.regions.deep));
  for c=1:size(configs,1)
   cm=cur;cm.Q=configs{c,2};[S0,S]=solve(D,s,cm);a=score_gather(S0,T,cm);b=score_gather(S,T,cm);
   rows{end+1}=struct('condition',in.obs(k).name,'config',configs{c,1},'Q',cm.Q,'full_overall',a.overall,'std_overall',b.overall,'std_deep',b.deep,'std_low_energy',b.low_energy); %#ok<AGROW>
   fprintf('%s %-8s %-17s Q %-26s FullCal %.4f | Standard %.4f deep %.4f lowE %.4f\n',v{1},in.obs(k).name,configs{c,1},mat2str(round(cm.Q,1)),a.overall,b.overall,b.deep,b.low_energy);
  end
 end
 out.(v{1})=[rows{:}];
end
uq_writejson(fullfile(root,'analysis','interim_hybrid.json'),out);
end

function [S0,S]=solve(D,s,cm)
d=prune_dictionary(cm,450);Ws=d.W./reshape(d.scale,1,[]);d.source_gram=Ws'*Ws;dicts=repmat({d},1,11);[nt,nx]=size(D);
R0=zeros(nt,nx);S0=R0;sup=cell(1,nx);E=sup;V=sup;
for b=1:9,ids=(b-1)*11+(1:11);ini=reflectivity_refit_prepared(D(:,ids),dicts,s(ids));so=reflectivity_soft_prepared(D(:,ids),dicts,s(ids),ini);
 R0(:,ids)=ini.reflectivity;S0(:,ids)=ini.gather;sup(ids)=ini.support;E(ids)=so.group_energy;V(ids)=so.group_noise_energy;end
cs=coherence_scores(R0,sup,d,any(D~=0,1),2,.25);st=struct('R0',R0,'supports',{sup},'E',{E},'V',{V},'groups',{cs.groups},'fallback',{cs.fallback},'W',d.W);
S=apply_gain(st,cs.c,struct('lambda_min',.5,'lambda_max',16));
end
