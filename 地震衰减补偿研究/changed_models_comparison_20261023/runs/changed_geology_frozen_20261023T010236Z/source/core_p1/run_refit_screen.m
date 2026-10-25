function run_refit_screen(outdir,smoke)
if nargin<2,smoke=false;end
root=fileparts(mfilename('fullpath'));addpath(fullfile(root,'matlab_src'));
assert(~isfolder(outdir));mkdir(outdir);mkdir(fullfile(outdir,'cases'));start=tic;
plan=jsondecode(fileread(fullfile(root,'inputs','adaptive_cv_plan.json')));uq_writejson(fullfile(outdir,'executed_plan.json'),plan);
cache=plan.cache;files=dir(fullfile(cache,'cases','*.mat'));assert(numel(files)==60);
if smoke
 names={files.name};keep=contains(names,'_WeakSection_12dB_20261411')&~contains(names,'OriginalReference');files=files(keep);assert(numel(files)==2);
end
policy=jsondecode(fileread(fullfile(root,'inputs','candidate_policy.json')));spec=jsondecode(fileread(fullfile(root,'inputs','p1_experiment_spec.json')));
cached=jsondecode(fileread(fullfile(root,'inputs','cached_shrinkage','trial_plan.json')));
previous='';rows=struct([]);count=0;maxnormal=0;flags=0;replay=0;maxsameq=0;maxsamef=0;
for k=1:numel(files)
 z=load(fullfile(files(k).folder,files(k).name));meta=z.meta;
 if ~strcmp(previous,meta.model)
  a=jsondecode(fileread(fullfile(cache,['calibration_' meta.model '.json'])));
  for key={'vp','h','offset','lambda_grid'},f=key{1};a.corrected.(f)=reshape(a.corrected.(f),1,[]);end
  a.corrected.t=a.corrected.t(:);previous=meta.model;dcache=cell(1,2);
  cmfull=uq_inverse_config(a.corrected,a.estimate);fullDict=cell(1,11);fullModels=prepare_models_v2(a.corrected.offset,cmfull);Ravo=avo_projector_v2(fullModels);
  for j=1:11,fullDict{j}=refit_dictionary(cmfull,a.corrected.offset(j));end
  c=config_model();c.model_name=meta.model;
  if ~strcmp(meta.model,'OriginalReference')
   list=cached.additional_models;ix=find(strcmp({list.name},meta.model),1);m=list(ix);
   for key={'vp','vs','rho','Q','event_t','fpeak'},f=key{1};c.(f)=reshape(m.(f),1,[]);end
   c.h=diff([0,c.event_t]).*c.vp(1:3)/2;
  end
  fprintf('Prepared calibrated dictionaries: %s elapsed=%.1fs\n',meta.model,toc(start));
 end
 id=1+strcmp(meta.domain,'Section');
 if isempty(dcache{id})
  est=p1_candidate_domain_estimate(a.diagnostics.fixed_estimate,a.estimate,a.corrected,policy,meta.domain);cm=uq_inverse_config(a.corrected,est);
  if id==1,dictSame=fullDict;dictFull=fullDict;
  else,one=refit_dictionary(cm,z.offsets(1));dictSame=repmat({one},1,11);dictFull=repmat(fullDict(6),1,11);end
  dcache{id}=struct('cm',cm,'same',{dictSame},'full',{dictFull});
 else
  dc=dcache{id};cm=dc.cm;dictSame=dc.same;dictFull=dc.full;
 end
 maxsameq=max(maxsameq,max(abs(reshape(cm.Q,1,[])-reshape(meta.Q_used,1,[]))));maxsamef=max(maxsamef,abs(cm.fpeak-meta.estimated_frequency));
 D=z.D;sigma=zeros(1,11);pilot=zeros(size(D));lam=sigma;
 for j=1:11
  ni=estimate_trace_snr(D(:,j),cm.t);sigma(j)=ni.sigma;
  if id==1
   F=fullModels{j}.direct;[~,lam(j)]=discrepancy_apply(F,D(:,j),cm.lambda_grid,sqrt(cm.nt)*sigma(j));
   pilot(:,j)=F.V*(F.z./(F.z.^2+4*lam(j)).*(F.U'*D(:,j)));
  end
 end
 Dwork=D;
 if id==1,[Dwork,~,~,~,~]=p1_noise_uncertainty(D,pilot,fullModels,sigma,lam,'Offset',Ravo,spec.p1_hyperparameters);end
 outSame=reflectivity_refit_prepared(Dwork,dictSame,sigma);
 if id==1,outFull=outSame;else,outFull=reflectivity_refit_prepared(Dwork,dictFull,sigma);end
 scene=make_stability_scene(meta.scenario,c);assert(norm(scene.truth-z.truth,'fro')<1e-12);
 base=score_stability(z.Snew,scene,c);replay=max(replay,abs(base.RelativeError-meta.new_score.RelativeError));
 for method=1:2
  label='SameCalRefit';out=outSame;Q=cm.Q;freq=cm.fpeak;
  if method==2,label='FullCalRefit';out=outFull;Q=cmfull.Q;freq=cmfull.fpeak;end
  score=score_stability(out.gather,scene,c);maxnormal=max(maxnormal,out.normal_residual);flags=flags+out.optimization_flags;
  row=struct('CaseID',meta.id,'Model',meta.model,'Domain',meta.domain,'Scenario',meta.scenario,'SNR_dB',meta.snr_db,'Seed',meta.seed,'Method',label);
  f=fieldnames(score);for i=1:numel(f),row.(f{i})=score.(f{i});end
  row.NormalResidual=out.normal_residual;row.OptimizationFlags=out.optimization_flags;row.MaxKKT=max(cellfun(@(x)x.kkt_ratio,out.lasso));
  row.MeanSupport=mean(out.support_after);row.Q1=Q(1);row.Q2=Q(2);row.Q3=Q(3);row.Q4=Q(4);row.Frequency=freq;
  count=count+1;if count==1,rows=row;else,rows(count)=row;end
 end
 same=outSame.gather;full=outFull.gather;t=z.t;truth=z.truth;observed=D;baseline=z.Snew;fixed=z.Sfixed;
 diagnosticSame=rmfield(outSame,{'gather','reflectivity'});diagnosticFull=rmfield(outFull,{'gather','reflectivity'});
 plotdata=struct('meta',meta,'t',t,'truth',truth,'observed',D,'fixed',fixed,'baseline',baseline,'same',same,'full',full, ...
 'same_diagnostics',diagnosticSame,'full_diagnostics',diagnosticFull);
 uq_writejson(fullfile(outdir,'cases',[meta.id '.json']),plotdata);
 reflectivitySame=outSame.reflectivity;reflectivityFull=outFull.reflectivity;
 save(fullfile(outdir,'cases',files(k).name),'same','full','t','truth','observed','baseline','fixed','meta','diagnosticSame','diagnosticFull','reflectivitySame','reflectivityFull','-v7');
 writetable(struct2table(rows),fullfile(outdir,'metrics.csv'));
 uq_writejson(fullfile(outdir,'progress.json'),struct('cases',k,'expected',numel(files),'elapsed_seconds',toc(start),'optimization_flags',flags));
 s1=score_stability(same,scene,c);s2=score_stability(full,scene,c);
 fprintf('%d/%d %s baseline=%.5f same=%.5f full=%.5f flags=%d elapsed=%.1fs\n',k,numel(files),meta.id,base.RelativeError,s1.RelativeError,s2.RelativeError,flags,toc(start));
end
assert(maxsameq<1e-10&&maxsamef<1e-10&&replay<1e-12);
uq_writejson(fullfile(outdir,'checks.json'),struct('complete',true,'cases',numel(files),'rows',count,'max_normal',maxnormal,'same_Q_max_difference',maxsameq,'same_frequency_max_difference',maxsamef,'baseline_score_replay',replay,'optimization_flags',flags,'elapsed_seconds',toc(start)));
end
