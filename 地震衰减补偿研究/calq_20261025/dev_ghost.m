function dev_ghost()
% Development of GhostCal on V3/V5 (inputs/PLAN_ZH.md): noise-free bias, noise
% replicates (same seeds as diag_q4), and section errors with the new Q.
root=fileparts(mfilename('fullpath'));addpath(fullfile(root,'core'),fullfile(root,'source'));maxNumCompThreads(4);
policy=jsondecode(fileread(fullfile(root,'inputs','candidate_policy.json')));base=jsondecode(fileread(fullfile(root,'analysis','diag_q4.json')));out=struct();started=tic;
for v={'V3','V5'}
 in=load(fullfile(root,'dev',[v{1} '.mat']));TT=load(fullfile(root,'dev',[v{1} '_truth.mat']));model=gen_complex_truth(TT.g);p=model.public;qtrue=model.truth.calibration_Q;
 e0=prepare_ghost(model.calibration_clean,p.offset,p);
 fprintf('%s noise-free GhostCal Q %s tau %.0f ms ln %s | %.0fs\n',v{1},mat2str(e0.layer_Q,4),1e3*e0.ghost_delay,mat2str(log(e0.layer_Q./qtrue),3),toc(started));
 [en,cn]=prepare_ghost(model.calibrationD,p.offset,p);
 fprintf('%s noisy (archived) GhostCal Q %s tau %.0f ms | current %s\n',v{1},mat2str(en.layer_Q,4),1e3*en.ghost_delay,mat2str(in.config.Q,4));
 reps=zeros(8,4);clean=model.calibration_clean;
 for r=1:8
  rng(2026102500+r,'twister');D=clean;
  for j=1:11,n=randn(p.nt,1);D(:,j)=clean(:,j)+n/sqrt(mean(n.^2))*sqrt(mean(clean(:,j).^2))*10^(-18/20);end
  e=prepare_ghost(D,p.offset,p);reps(r,:)=e.layer_Q;
 end
 lr=log(reps./qtrue);old=base.(v{1}).replicate_log_bias;
 fprintf('%s replicates GhostCal ln(Q/Qtrue) mean %s std %s | current mean %s | %.0fs\n',v{1},mat2str(mean(lr,1),3),mat2str(std(lr,0,1),3),mat2str(old(:)',3),toc(started));
 cmnew=uq_inverse_config(cn,en);cmnew.t=cmnew.t(:);cmold=in.config;cmold.t=cmold.t(:);sec=section_errors(in,TT,cmold,cmnew);
 out.(v{1})=struct('true_Q',qtrue,'noise_free_Q',e0.layer_Q,'noise_free_tau',e0.ghost_delay,'archived_noisy_Q',en.layer_Q,'archived_noisy_tau',en.ghost_delay, ...
  'replicate_Q',reps,'replicate_log_bias',mean(lr,1),'replicate_log_std',std(lr,0,1),'current_replicate_log_bias',old(:)','section',sec);
end
uq_writejson(fullfile(root,'analysis','dev_ghost.json'),out);
end

function sec=section_errors(in,TT,cmold,cmnew)
rows={};
for k=1:2
 D=in.obs(k).D;s=reshape(in.obs(k).sigma,1,[]);T=struct('target',TT.tr(k).target,'masks',struct('fault',TT.regions.fault,'pinch',TT.regions.pinch,'low_energy',TT.regions.low_energy(:,:,k),'deep',TT.regions.deep));
 for c={'current',cmold;'GhostCal',cmnew}'
  [S0,Sstd]=solve(D,s,c{2});
  for m={'FullCalRefit',S0;'Standard',Sstd}'
   x=score_gather(m{2},T,c{2});rows{end+1}=struct('condition',in.obs(k).name,'calibration',c{1},'method',m{1},'overall',x.overall,'deep',x.deep,'low_energy',x.low_energy); %#ok<AGROW>
   fprintf('  %s %-8s %-12s overall %.4f deep %.4f lowE %.4f\n',in.obs(k).name,c{1},m{1},x.overall,x.deep,x.low_energy);
  end
 end
end
sec=[rows{:}];
end

function [S0,S]=solve(D,s,cm)
d=prune_dictionary(cm,450);Ws=d.W./reshape(d.scale,1,[]);d.source_gram=Ws'*Ws;dicts=repmat({d},1,11);[nt,nx]=size(D);
R0=zeros(nt,nx);S0=R0;sup=cell(1,nx);E=sup;V=sup;
for b=1:9,ids=(b-1)*11+(1:11);ini=reflectivity_refit_prepared(D(:,ids),dicts,s(ids));so=reflectivity_soft_prepared(D(:,ids),dicts,s(ids),ini);
 R0(:,ids)=ini.reflectivity;S0(:,ids)=ini.gather;sup(ids)=ini.support;E(ids)=so.group_energy;V(ids)=so.group_noise_energy;end
cs=coherence_scores(R0,sup,d,any(D~=0,1),2,.25);st=struct('R0',R0,'supports',{sup},'E',{E},'V',{V},'groups',{cs.groups},'fallback',{cs.fallback},'W',d.W);
S=apply_gain(st,cs.c,struct('lambda_min',.5,'lambda_max',16));
end
