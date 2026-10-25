function run_sealed_validation()
% Frozen rules on sealed V4/V5. Solver phase never loads truth/.
root=fileparts(mfilename('fullpath'));base=fileparts(root);v1=fullfile(base,'coherence_soft_20261023');
addpath(fullfile(root,'core'),fullfile(v1,'core'),fullfile(v1,'source'));maxNumCompThreads(4);
F=jsondecode(fileread(fullfile(root,'inputs','FROZEN.json')));names=fieldnames(F.rules);
cdir=fullfile(root,'cache','sealed');rdir=fullfile(root,'results');if ~isfolder(cdir),mkdir(cdir);end;if ~isfolder(rdir),mkdir(rdir);end
started=tic;flags=0;maxnormal=0;
for v={'V4','V5'}
 in=load(fullfile(root,'inputs','sealed',[v{1} '.mat']));assert(~isfield(in.obs,'target'));cm=in.config;cm.t=cm.t(:);Q0=reshape(cm.Q,1,[]);
 for e=[-30 -10 0 10 30]
  used=cm;used.Q=Q0*(1+e/100);d=prune_dictionary(used,450);Ws=d.W./reshape(d.scale,1,[]);d.source_gram=Ws'*Ws;dicts=repmat({d},1,11);
  for k=1:numel(in.obs)
   D=in.obs(k).D;sig=reshape(in.obs(k).sigma,1,[]);live=any(D~=0,1);
   R0=zeros(size(D));S0=R0;hard=R0;soft=R0;sup=cell(1,99);E=sup;V=sup;
   for b=1:9
    ids=(b-1)*11+(1:11);cf=fullfile(cdir,sprintf('%s_Q%+d_%s_%02d.mat',v{1},e,in.obs(k).name,b));
    if isfile(cf),c=load(cf);else
     c.initial=reflectivity_refit_prepared(D(:,ids),dicts,sig(ids));c.pruned=reflectivity_prune_prepared(D(:,ids),dicts,sig(ids),c.initial);
     c.softout=reflectivity_soft_prepared(D(:,ids),dicts,sig(ids),c.initial);save(cf,'-struct','c','-v7');
    end
    R0(:,ids)=c.initial.reflectivity;S0(:,ids)=c.initial.gather;hard(:,ids)=c.pruned.gather;soft(:,ids)=c.softout.gather;
    sup(ids)=c.initial.support;E(ids)=c.softout.group_energy;V(ids)=c.softout.group_noise_energy;
    flags=flags+c.initial.optimization_flags+c.pruned.optimization_flags;maxnormal=max([maxnormal,c.initial.normal_residual,c.softout.initial_normal_residual]);
   end
   cs=coherence_scores(R0,sup,d,live,F.common.half_width,F.common.dip_periods);
   s=struct('R0',R0,'supports',{sup},'E',{E},'V',{V},'groups',{cs.groups},'fallback',{cs.fallback},'cm',used,'D',D,'live',live,'W',d.W);
   s=q_spread_energy(s,d,F.rules.Robust.q_factors(:)');
   one=struct('lambda_min',1,'lambda_max',1);assert(norm(apply_gain(s,cs.c,one)-soft,'fro')<=1e-10*norm(soft,'fro'));
   gathers=struct('FullCalRefit',S0,'PrunedFullCalRefit',hard,'SoftGroupRefit',soft);
   for r=1:numel(names),o=F.rules.(names{r});o.tau=0;gathers.(names{r})=apply_gain(s,cs.c,o);end
   save(fullfile(rdir,sprintf('%s_Q%+d_%s.mat',v{1},e,in.obs(k).name)),'gathers','used','-v7');
   fprintf('%s Q%+d %s done %.1fs\n',v{1},e,in.obs(k).name,toc(started));
  end
 end
end
uq_writejson(fullfile(root,'analysis','sealed_checks.json'),struct('optimization_flags',flags,'max_normal_residual',maxnormal,'elapsed_seconds',toc(started)));
score_sealed(root);
end
