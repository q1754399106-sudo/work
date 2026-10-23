function run_validation()
% Frozen validation. Solver phase reads only inputs/validation; truth/ is
% opened afterwards for scoring. Every method in a scene shares support, Q,
% wavelet, sigma and observation; only the shrinkage gain differs.
root=fileparts(mfilename('fullpath'));addpath(fullfile(root,'core'),fullfile(root,'source'));maxNumCompThreads(4);
sel=jsondecode(fileread(fullfile(root,'analysis','SELECTED.json')));assert(sel.frozen_before_validation);
test_coherence_math();started=tic;
rule=struct('half_width',sel.half_width,'dip_periods',sel.dip_periods,'lambda_min',sel.lambda_min,'lambda_max',sel.lambda_max);
noalign=rule;noalign.dip_periods=sel.noalign_dip_periods;
uniform=rule;uniform.lambda_min=sel.uniform_lambda;uniform.lambda_max=sel.uniform_lambda;
resdir=fullfile(root,'results');cachedir=fullfile(root,'cache');flags=0;maxnormal=0;blocks=0;
for name={'V1','V2','V3'}
 v=name{1};in=load(fullfile(root,'inputs','validation',[v '.mat']));
 assert(~isfield(in,'tr')&&~isfield(in.obs,'target'),'Truth leaked into solver input.');
 cm=in.config;cm.t=cm.t(:);referenceQ=reshape(cm.Q,1,[]);
 for errorpct=[-30 -10 0 10 30]
  used=cm;used.Q=referenceQ*(1+errorpct/100);dict=prune_dictionary(used,450);
  Ws=dict.W./reshape(dict.scale,1,[]);dict.source_gram=Ws'*Ws;dicts11=repmat({dict},1,11);dicts=repmat({dict},1,99);
  for k=1:numel(in.obs)
   D=in.obs(k).D;sig=reshape(in.obs(k).sigma,1,[]);nx=size(D,2);assert(nx==99);live=any(D~=0,1);
   R0=zeros(size(D));S0=R0;hard=R0;soft=R0;supports=cell(1,nx);E=supports;V=supports;
   for block=1:9
    ids=(block-1)*11+(1:11);cf=fullfile(cachedir,sprintf('%s_Q%+d_%s_%02d.mat',v,errorpct,in.obs(k).name,block));
    if isfile(cf),c=load(cf,'initial','pruned','softout');initial=c.initial;pruned=c.pruned;softout=c.softout;
    else
     initial=reflectivity_refit_prepared(D(:,ids),dicts11,sig(ids));
     pruned=reflectivity_prune_prepared(D(:,ids),dicts11,sig(ids),initial);
     softout=reflectivity_soft_prepared(D(:,ids),dicts11,sig(ids),initial);
     save(cf,'initial','pruned','softout','-v7');
    end
    R0(:,ids)=initial.reflectivity;S0(:,ids)=initial.gather;hard(:,ids)=pruned.gather;soft(:,ids)=softout.gather;
    supports(ids)=initial.support;E(ids)=softout.group_energy;V(ids)=softout.group_noise_energy;
    flags=flags+initial.optimization_flags+pruned.optimization_flags;maxnormal=max([maxnormal,initial.normal_residual,softout.initial_normal_residual]);blocks=blocks+1;
   end
   one=rule;one.lambda_min=1;one.lambda_max=1;same=coherence_soft_gather(R0,supports,dicts,E,V,live,one);
   assert(norm(same.gather-soft,'fro')<=1e-10*norm(soft,'fro'),'lambda=1 must reproduce SoftGroupRefit.');
   coh=coherence_soft_gather(R0,supports,dicts,E,V,live,rule);
   na=coherence_soft_gather(R0,supports,dicts,E,V,live,noalign);
   un=coherence_soft_gather(R0,supports,dicts,E,V,live,uniform);
   gathers=struct('FullCalRefit',S0,'PrunedFullCalRefit',hard,'SoftGroupRefit',soft,'CoherenceSoftRefit',coh.gather,'NoAlign',na.gather,'UniformLambda',un.gather);
   diagnostics=struct('coherence',{coh.coherence},'lambda',{coh.lambda},'fallback',{coh.fallback},'gain',{coh.group_gain},'soft_gain',{same.group_gain}); %#ok<NASGU>
   save(fullfile(resdir,sprintf('%s_Q%+d_%s.mat',v,errorpct,in.obs(k).name)),'gathers','diagnostics','used','-v7');
   fprintf('%s Q%+d %s solved, elapsed %.1fs\n',v,errorpct,in.obs(k).name,toc(started));
  end
 end
end
uq_writejson(fullfile(root,'analysis','validation_checks.json'),struct('blocks',blocks,'optimization_flags',flags,'max_normal_residual',maxnormal,'elapsed_seconds',toc(started)));
score_validation(root);
end
