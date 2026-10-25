function run_coarse()
% Post-hoc follow-up (added after run_trueq): dictionary with the TRUE coarse
% four-layer Q (travel-time harmonic means of the true Q, the quantity the
% calibration tries to estimate), Ricker wavelet. Separates estimation error of
% the four layer Q values from the coarseness of a four-layer laterally constant Q.
root=fileparts(mfilename('fullpath'));addpath(fullfile(root,'core'),fullfile(root,'source'));maxNumCompThreads(4);
in=load(fullfile(root,'dev','V5.mat'));TT=load(fullfile(root,'dev','V5_truth.mat'));cm=in.config;cm.t=cm.t(:);nt=cm.nt;nx=99;
model=gen_complex_truth(TT.g);cq=model.truth.calibration_Q;fprintf('calibrated Q %s | true coarse Q %s\n',mat2str(cm.Q,4),mat2str(cq,4));
cases=struct('name',{'noisy_Moderate','noisy_Strong','clean_primary'},'D',{in.obs(1).D,in.obs(2).D,model.datasets(1).primary_att},'sigma',{in.obs(1).sigma,in.obs(2).sigma,in.obs(1).sigma},'k',{1,2,1});
c2=cm;c2.Q=cq;d=prune_dictionary(c2,450);Ws=d.W./reshape(d.scale,1,[]);d.source_gram=Ws'*Ws;dicts=repmat({d},1,11);rows={};
for c=1:3
 D=cases(c).D;D(:,~model.live)=0;s=reshape(cases(c).sigma,1,[]);R0=zeros(nt,nx);S0=R0;soft=R0;sup=cell(1,nx);E=sup;V=sup;
 for b=1:9,ids=(b-1)*11+(1:11);ini=reflectivity_refit_prepared(D(:,ids),dicts,s(ids));so=reflectivity_soft_prepared(D(:,ids),dicts,s(ids),ini);
  R0(:,ids)=ini.reflectivity;S0(:,ids)=ini.gather;soft(:,ids)=so.gather;sup(ids)=ini.support;E(ids)=so.group_energy;V(ids)=so.group_noise_energy;end
 cs=coherence_scores(R0,sup,d,any(D~=0,1),2,.25);st=struct('R0',R0,'supports',{sup},'E',{E},'V',{V},'groups',{cs.groups},'fallback',{cs.fallback},'W',d.W);
 G=struct('FullCalRefit',S0,'SoftGroupRefit',soft,'Standard',apply_gain(st,cs.c,struct('lambda_min',.5,'lambda_max',16)));
 k=cases(c).k;T=struct('target',TT.tr(k).target,'masks',struct('fault',TT.regions.fault,'pinch',TT.regions.pinch,'low_energy',TT.regions.low_energy(:,:,k),'deep',TT.regions.deep));
 for m=fieldnames(G)'
  x=score_gather(G.(m{1}),T,cm);rows{end+1}=struct('dictionary','CoarseTrueQ','data',cases(c).name,'method',m{1},'overall',x.overall,'low_energy',x.low_energy,'fault',x.fault,'pinch',x.pinch,'deep',x.deep,'leak',x.leak); %#ok<AGROW>
  fprintf('CoarseTrueQ %-15s %-15s overall %.4f lowE %.4f deep %.4f\n',cases(c).name,m{1},x.overall,x.low_energy,x.deep);
 end
end
uq_writejson(fullfile(root,'analysis','coarse_true_q_metrics.json'),[rows{:}]);
end
