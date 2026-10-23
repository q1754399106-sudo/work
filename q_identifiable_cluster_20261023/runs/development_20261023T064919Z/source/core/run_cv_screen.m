function run_cv_screen(outdir,maximum,pilot,checkpoint)
if nargin<2,maximum=60;end
if nargin<3,pilot=false;end
if nargin<4,checkpoint='';end
root=fileparts(mfilename('fullpath'));addpath(fullfile(root,'matlab_src'));
assert(~isfolder(outdir));mkdir(outdir);mkdir(fullfile(outdir,'cases'));start=tic;
plan=jsondecode(fileread(fullfile(root,'inputs','adaptive_cv_plan.json')));uq_writejson(fullfile(outdir,'executed_plan.json'),plan);
cache=plan.cache;files=dir(fullfile(cache,'cases','*.mat'));assert(numel(files)==60);
if pilot
 keep=false(size(files));
 for i=1:numel(files),name=files(i).name;keep(i)=contains(name,'_12dB_20261411')&&any(contains(name,{'_SmoothAVO_','_WeakSection_','_CoherentNoiseAVO_'}));end
 files=files(keep);assert(numel(files)==9);
end
policy=jsondecode(fileread(fullfile(root,'inputs','candidate_policy.json')));
cached=jsondecode(fileread(fullfile(root,'inputs','cached_shrinkage','trial_plan.json')));
previous='';rows=struct([]);maxnormal=0;maxq=0;maxf=0;replay=0;completed=0;
if ~isempty(checkpoint)
 saved=readtable(fullfile(checkpoint,'metrics.csv'),'TextType','string');assert(mod(height(saved),2)==0);completed=height(saved)/2;
 rows=table2struct(saved);
 for ii=1:numel(rows),for field={'CaseID','Model','Domain','Scenario','Method'},key=field{1};rows(ii).(key)=char(rows(ii).(key));end,end
 maxnormal=max(saved.NormalResidual);
 for ii=1:completed
  old=load(fullfile(checkpoint,'cases',files(ii).name));assert(strcmp(rows(2*ii-1).CaseID,old.meta.id));
  maxq=max(maxq,max(abs(reshape(old.Q_used,1,[])-reshape(old.meta.Q_used,1,[]))));
  maxf=max(maxf,abs(old.frequency_used-old.meta.estimated_frequency));
  replay=max(replay,abs(norm(old.baseline-old.truth,'fro')/norm(old.truth,'fro')-old.meta.new_score.RelativeError));
 end
 copyfile(fullfile(checkpoint,'cases','*'),fullfile(outdir,'cases'));
 fprintf('Resuming %d saved cases without recomputation\n',completed);
end
for k=1:min(maximum,numel(files))
 if k<=completed,continue;end
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
 for method=1:2
  label='AdaptiveJointCV';sc=score;alpha=out.alpha;betaFactor=out.beta_factor;beta=out.beta;val=min(out.cv_mean(:));
  if method==2,label='RobustJointCV';sc=score_stability(out.robust_gather,scene,c);alpha=out.robust_alpha;betaFactor=out.robust_beta_factor;beta=out.robust_beta;val=min(out.robust_cv_mean(:));end
  row=struct('CaseID',meta.id,'Model',meta.model,'Domain',meta.domain,'Scenario',meta.scenario,'SNR_dB',meta.snr_db,'Seed',meta.seed,'Method',label);
  f=fieldnames(sc);for i=1:numel(f),row.(f{i})=sc.(f{i});end
  row.Alpha=alpha;row.BetaFactor=betaFactor;row.Beta=beta;row.AlphaEndpoint=any(alpha==[plan.alpha_grid(1),plan.alpha_grid(end)]);row.BetaEndpoint=any(betaFactor==[plan.beta_grid(1),plan.beta_grid(end)]);
  row.NormalResidual=out.max_normal;row.FinalFlag=out.flags;row.CVLoss=val;
  ix=(k-1)*2+method;if ix==1,rows=row;else,rows(ix)=row;end
 end
 selection=rmfield(out,{'gather','robust_gather','prior'});uq_writejson(fullfile(outdir,'cases',[meta.id '_selection.json']),selection);
 plotdata=struct('meta',meta,'selection',selection,'t',z.t,'truth',z.truth,'observed',z.D,'fixed',z.Sfixed,'baseline',z.Snew,'optimized',S);
 uq_writejson(fullfile(outdir,'cases',[meta.id '.json']),plotdata);
 t=z.t;truth=z.truth;observed=z.D;baseline=z.Snew;fixed=z.Sfixed;Q_used=cm.Q;frequency_used=cm.fpeak;robust=out.robust_gather;
 save(fullfile(outdir,'cases',files(k).name),'S','robust','t','truth','observed','baseline','fixed','meta','selection','Q_used','frequency_used','-v7');
 writetable(struct2table(rows),fullfile(outdir,'metrics.csv'));
 uq_writejson(fullfile(outdir,'progress.json'),struct('cases',k,'expected',maximum,'elapsed_seconds',toc(start)));
 fprintf('%d/%d %s baseline=%.5f CV=%.5f alpha=%.3g beta=%.3g elapsed=%.1fs\n',k,maximum,meta.id,baseScore.RelativeError,score.RelativeError,out.alpha,out.beta_factor,toc(start));
end
assert(maxq<1e-10&&maxf<1e-10&&replay<1e-12);
uq_writejson(fullfile(outdir,'checks.json'),struct('complete',true,'cases',k,'expected',maximum,'max_normal',maxnormal,'same_Q_max_difference',maxq,'same_frequency_max_difference',maxf,'baseline_score_replay',replay,'elapsed_seconds',toc(start)));
end
