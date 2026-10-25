function run_pruning_development(outdir)
root=fileparts(mfilename('fullpath'));addpath(fullfile(root,'matlab_src'));assert(~isfolder(outdir));mkdir(outdir);mkdir(fullfile(outdir,'cases'));start=tic;
sources=jsondecode(fileread(fullfile(root,'inputs','pruning_development_sources.json')));plan=jsondecode(fileread(fullfile(root,'inputs','adaptive_cv_plan.json')));
extra=jsondecode(fileread(fullfile(root,'inputs','refit_validation_plan.json')));cached=jsondecode(fileread(fullfile(root,'inputs','cached_shrinkage','trial_plan.json')));
spec=jsondecode(fileread(fullfile(root,'inputs','p1_experiment_spec.json')));allfiles=[dir(fullfile(sources.original_development,'cases','*.mat'));dir(fullfile(sources.previous_validation,'cases','*.mat'))];assert(numel(allfiles)==110);
previous='';rows=struct([]);count=0;flags=0;maxnormal=0;
for k=1:numel(allfiles)
 z=load(fullfile(allfiles(k).folder,allfiles(k).name));meta=z.meta;oldset=~isfield(meta,'suite');
 if ~strcmp(previous,meta.model)
  calibrationRoot=plan.cache;if ~oldset,calibrationRoot=sources.previous_validation;end
  a=jsondecode(fileread(fullfile(calibrationRoot,['calibration_' meta.model '.json'])));for key={'vp','h','offset','lambda_grid'},f=key{1};a.corrected.(f)=reshape(a.corrected.(f),1,[]);end
  a.corrected.t=a.corrected.t(:);cm=uq_inverse_config(a.corrected,a.estimate);dictAll=cell(1,11);models=prepare_models_v2(a.corrected.offset,cm);Ravo=avo_projector_v2(models);
  for j=1:11,dictAll{j}=prune_dictionary(cm,a.corrected.offset(j));end
  c=config_model();c.model_name=meta.model;
  if ~strcmp(meta.model,'OriginalReference')
   list=cached.additional_models;if ~oldset,list=extra.models;end;ix=find(strcmp({list.name},meta.model),1);m=list(ix);
   for key={'vp','vs','rho','Q','event_t','fpeak'},f=key{1};c.(f)=reshape(m.(f),1,[]);end;c.h=diff([0,c.event_t]).*c.vp(1:3)/2;
  end
  previous=meta.model;fprintf('Pruning dictionaries: %s elapsed=%.1fs\n',meta.model,toc(start));
 end
 D=z.observed;
 if oldset,initial=z.diagnosticFull;initial.gather=z.full;initial.reflectivity=z.reflectivityFull;meta.suite='normal';
 else,initial=z.diagnostic;initial.gather=z.optimized;initial.reflectivity=z.reflectivity;end
 sigma=initial.lambda/(.5*sqrt(2*log(cm.nt)));Dwork=D;
 if strcmp(meta.domain,'Offset')
  dicts=dictAll;pilot=zeros(size(D));lam=sigma;
  for j=1:11,F=models{j}.direct;[~,lam(j)]=discrepancy_apply(F,D(:,j),cm.lambda_grid,sqrt(cm.nt)*sigma(j));pilot(:,j)=F.V*(F.z./(F.z.^2+4*lam(j)).*(F.U'*D(:,j)));end
  [Dwork,~,~,~,~]=p1_noise_uncertainty(D,pilot,models,sigma,lam,'Offset',Ravo,spec.p1_hyperparameters);
 else,dicts=repmat(dictAll(6),1,11);end
 out=reflectivity_prune_prepared(Dwork,dicts,sigma,initial);Snew=out.gather;flags=flags+out.optimization_flags;maxnormal=max(maxnormal,out.normal_residual);
 scene=make_stability_scene(meta.scenario,c);assert(norm(scene.truth-z.truth,'fro')<1e-12);
 for method=1:4
  label='FixedSplit';S=z.fixed;if method==2,label='DecoupledWeak';S=z.baseline;end
  if method==3,label='FullCalRefit';S=initial.gather;end
  if method==4,label='PrunedFullCalRefit';S=Snew;end
  score=score_stability(S,scene,c);row=struct('CaseID',meta.id,'Model',meta.model,'Domain',meta.domain,'Scenario',meta.scenario,'Suite',meta.suite,'SNR_dB',meta.snr_db,'Seed',meta.seed,'Method',label);
  f=fieldnames(score);for i=1:numel(f),row.(f{i})=score.(f{i});end;row.NormalResidual=out.normal_residual;row.OptimizationFlags=out.optimization_flags;
  count=count+1;if count==1,rows=row;else,rows(count)=row;end
 end
 optimized=Snew;unpruned=initial.gather;truth=z.truth;t=z.t;baseline=z.baseline;fixed=z.fixed;observed=D;diagnostic=rmfield(out,{'gather','reflectivity'});reflectivity=out.reflectivity;
 save(fullfile(outdir,'cases',allfiles(k).name),'meta','t','truth','observed','fixed','baseline','unpruned','optimized','diagnostic','reflectivity','-v7');
 uq_writejson(fullfile(outdir,'cases',[meta.id '.json']),struct('meta',meta,'t',t,'truth',truth,'observed',D,'fixed',fixed,'baseline',baseline,'unpruned',unpruned,'optimized',optimized,'diagnostic',diagnostic));
 writetable(struct2table(rows),fullfile(outdir,'metrics.csv'));uq_writejson(fullfile(outdir,'progress.json'),struct('cases',k,'expected',110,'elapsed_seconds',toc(start)));
 if mod(k,5)==0,fprintf('Pruning %d/110 elapsed=%.1fs\n',k,toc(start));end
end
uq_writejson(fullfile(outdir,'checks.json'),struct('complete',true,'cases',110,'rows',count,'max_normal',maxnormal,'optimization_flags',flags,'elapsed_seconds',toc(start),'previous_holdouts_reclassified_as_development',true));
end
