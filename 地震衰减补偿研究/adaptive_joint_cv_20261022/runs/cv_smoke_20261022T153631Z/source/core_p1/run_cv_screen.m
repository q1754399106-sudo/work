function run_cv_screen(outdir,maximum)
if nargin<2,maximum=60;end
root=fileparts(mfilename('fullpath'));addpath(fullfile(root,'matlab_src'));
assert(~isfolder(outdir));mkdir(outdir);mkdir(fullfile(outdir,'cases'));start=tic;
plan=jsondecode(fileread(fullfile(root,'inputs','adaptive_cv_plan.json')));uq_writejson(fullfile(outdir,'executed_plan.json'),plan);
cache=plan.cache;files=dir(fullfile(cache,'cases','*.mat'));assert(numel(files)==60);
policy=jsondecode(fileread(fullfile(root,'inputs','candidate_policy.json')));
cached=jsondecode(fileread(fullfile(root,'inputs','cached_shrinkage','trial_plan.json')));
previous='';rows=struct([]);maxnormal=0;maxq=0;maxf=0;replay=0;
for k=1:min(maximum,numel(files))
 z=load(fullfile(files(k).folder,files(k).name));meta=z.meta;
 if ~strcmp(previous,meta.model)
  a=jsondecode(fileread(fullfile(cache,['calibration_' meta.model '.json'])));
  for key={'vp','h','offset','lambda_grid'},f=key{1};a.corrected.(f)=reshape(a.corrected.(f),1,[]);end
  a.corrected.t=a.corrected.t(:);previous=meta.model;dcache=cell(1,2);
  c=config_model();c.model_name=meta.model;
  if ~strcmp(meta.model,'OriginalReference')
   list=cached.additional_models;ix=find(strcmp({list.name},meta.model),1);m=list(ix);
   for key={'vp','vs','rho','Q','event_t','fpeak'},f=key{1};c.(f)=reshape(m.(f),1,[]);end
   c.h=diff([0,c.event_t]).*c.vp(1:3)/2;
  end
 end
 id=1+strcmp(meta.domain,'Section');
 if isempty(dcache{id})
  est=p1_candidate_domain_estimate(a.diagnostics.fixed_estimate,a.estimate,a.corrected,policy,meta.domain);
  cm=uq_inverse_config(a.corrected,est);cm.candidate_profile=policy.selected_profile;cm.joint_prepare_only=true;
  if id==1,map=1:11;models=prepare_models_v2(z.offsets,cm);R=avo_projector_v2(models);
  else,map=6*ones(1,11);one=prepare_models_v2(z.offsets(1),cm);models=repmat(one,1,11);R=sparse(c.nt*11,c.nt*11);end
  masks=cv_masks(c.nt,11);tm=cell(1,numel(masks));
  for fold=1:numel(masks),tm{fold}=cv_training_models(models,masks{fold});end
  dcache{id}=struct('cm',cm,'map',map,'models',{models},'R',R,'tm',{tm});
  fprintf('Prepared independent fold operators: %s %s elapsed=%.1fs\n',meta.model,meta.domain,toc(start));
 else
  d=dcache{id};cm=d.cm;map=d.map;models=d.models;R=d.R;tm=d.tm;
 end
 maxq=max(maxq,max(abs(reshape(cm.Q,1,[])-reshape(meta.Q_used,1,[]))));maxf=max(maxf,abs(cm.fpeak-meta.estimated_frequency));
 out=adaptive_cv_prepared(z.D,models,map,cm,meta.domain,R,tm);S=out.gather;maxnormal=max(maxnormal,out.max_normal);
 scene=make_stability_scene(meta.scenario,c);assert(norm(scene.truth-z.truth,'fro')<1e-12);
 score=score_stability(S,scene,c);baseScore=score_stability(z.Snew,scene,c);replay=max(replay,abs(baseScore.RelativeError-meta.new_score.RelativeError));
 row=struct('CaseID',meta.id,'Model',meta.model,'Domain',meta.domain,'Scenario',meta.scenario,'SNR_dB',meta.snr_db,'Seed',meta.seed,'Method','AdaptiveJointCV');
 f=fieldnames(score);for i=1:numel(f),row.(f{i})=score.(f{i});end
 row.Alpha=out.alpha;row.BetaFactor=out.beta_factor;row.Beta=out.beta;row.AlphaEndpoint=out.alpha_endpoint;row.BetaEndpoint=out.beta_endpoint;
 row.NormalResidual=out.max_normal;row.FinalFlag=out.flags;row.CVLoss=min(out.cv_mean(:));
 if k==1,rows=row;else,rows(k)=row;end
 selection=rmfield(out,{'gather','prior'});uq_writejson(fullfile(outdir,'cases',[meta.id '_selection.json']),selection);
 plotdata=struct('meta',meta,'selection',selection,'t',z.t,'truth',z.truth,'observed',z.D,'fixed',z.Sfixed,'baseline',z.Snew,'optimized',S);
 uq_writejson(fullfile(outdir,'cases',[meta.id '.json']),plotdata);
 t=z.t;truth=z.truth;observed=z.D;baseline=z.Snew;fixed=z.Sfixed;Q_used=cm.Q;frequency_used=cm.fpeak;
 save(fullfile(outdir,'cases',files(k).name),'S','t','truth','observed','baseline','fixed','meta','selection','Q_used','frequency_used','-v7');
 writetable(struct2table(rows),fullfile(outdir,'metrics.csv'));
 uq_writejson(fullfile(outdir,'progress.json'),struct('cases',k,'expected',maximum,'elapsed_seconds',toc(start)));
 fprintf('%d/%d %s baseline=%.5f CV=%.5f alpha=%.3g beta=%.3g elapsed=%.1fs\n',k,maximum,meta.id,baseScore.RelativeError,score.RelativeError,out.alpha,out.beta_factor,toc(start));
end
assert(maxq<1e-10&&maxf<1e-10&&replay<1e-12);
uq_writejson(fullfile(outdir,'checks.json'),struct('complete',true,'cases',k,'expected',maximum,'max_normal',maxnormal,'same_Q_max_difference',maxq,'same_frequency_max_difference',maxf,'baseline_score_replay',replay,'elapsed_seconds',toc(start)));
end
