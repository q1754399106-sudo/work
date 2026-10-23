function run_refit_validation(outdir)
root=fileparts(mfilename('fullpath'));addpath(fullfile(root,'matlab_src'));
frozen=jsondecode(fileread(fullfile(root,'inputs','frozen_refit_profile.json')));assert(strcmp(frozen.method,'FullCalRefit'));
plan=jsondecode(fileread(fullfile(root,'inputs','refit_validation_plan.json')));policy=jsondecode(fileread(fullfile(root,'inputs','candidate_policy.json')));
spec=jsondecode(fileread(fullfile(root,'inputs','p1_experiment_spec.json')));
assert(~isfolder(outdir));mkdir(outdir);mkdir(fullfile(outdir,'cases'));start=tic;
uq_writejson(fullfile(outdir,'executed_plan.json'),plan);uq_writejson(fullfile(outdir,'frozen_profile.json'),frozen);
rows=struct([]);count=0;cases=0;flags=0;maxnormal=0;linearflags=0;maxlinear=0;
for im=1:numel(plan.models)
 m=plan.models(im);c=config_model();for key={'vp','vs','rho','Q','event_t','fpeak'},f=key{1};c.(f)=reshape(m.(f),1,[]);end
 c.model_name=m.name;c.nominal_frequency=m.nominal_frequency;c.h=diff([0,c.event_t]).*c.vp(1:3)/2;public=uq_public_config(c,c.nominal_frequency);
 truthModels=cell(1,11);for j=1:11,[a,b]=ray_absorption(c.t,c.offset(j),c);truthModels{j}=struct('G',attenuation_matrix(a+b,c));end
 calScene=make_stability_scene('SmoothAVO',c);[calD,~]=uq_observation(calScene,c,truthModels,1:11,plan.calibration_seed,plan.calibration_snr_db);
 reusefile=fullfile(root,'inputs','validation_calibration_cache.json');reused=false;
 if isfile(reusefile)
  reuse=jsondecode(fileread(reusefile));jp=fullfile(reuse.path,['calibration_' m.name '.json']);mp=fullfile(reuse.path,['calibration_input_' m.name '.mat']);
  if isfile(jp)&&isfile(mp)
   cachedInput=load(mp);assert(isequal(calD,cachedInput.calD));aa=jsondecode(fileread(jp));
   baseEst=aa.estimate;corrected=aa.corrected;cal=aa.diagnostics;
   for key={'vp','h','offset','lambda_grid'},f=key{1};corrected.(f)=reshape(corrected.(f),1,[]);end
   corrected.t=corrected.t(:);reused=true;
  end
 end
 if ~reused,[baseEst,corrected,cal]=p1_prepare_candidate(calD,c.offset,public,policy);end
 cmfull=uq_inverse_config(corrected,baseEst);
 save(fullfile(outdir,['calibration_input_' m.name '.mat']),'calD','public','-v7');
 uq_writejson(fullfile(outdir,['calibration_' m.name '.json']),struct('estimate',baseEst,'diagnostics',cal,'public',public,'corrected',corrected));
 fullModels=prepare_models_v2(c.offset,cmfull);Rfull=avo_projector_v2(fullModels);fullDict=cell(1,11);for j=1:11,fullDict{j}=refit_dictionary(cmfull,c.offset(j));end
 fprintf('Validation calibration complete: %s elapsed=%.1fs\n',m.name,toc(start));
 for id=1:2
  domain='Offset';offsets=c.offset;map=1:11;if id==2,domain='Section';offsets=repmat(c.offset(6),1,11);map=6*ones(1,11);end
  est=p1_candidate_domain_estimate(cal.fixed_estimate,baseEst,corrected,policy,domain);cm=uq_inverse_config(corrected,est);cm.candidate_profile=policy.selected_profile;cm.joint_prepare_only=true;
  if id==1,models=fullModels;R=Rfull;dicts=fullDict;
  else,one=prepare_models_v2(offsets(1),cm);models=repmat(one,1,11);R=sparse(c.nt*11,c.nt*11);dicts=repmat(fullDict(6),1,11);end
  original=fixed_v1_operators(cm,offsets);
  for isc=1:numel(plan.scenarios)
   name=plan.scenarios{isc};if endsWith(name,'Section')~=(id==2),continue;end
   scene=make_stability_scene(name,c);
   suites=struct([]);q=0;
   for snr=reshape(plan.snr_db,1,[]),for seed=reshape(plan.seeds,1,[])
    q=q+1;oneSuite=struct('kind','normal','snr_db',snr,'seed',seed,'ar1_rho',0);
    if q==1,suites=oneSuite;else,suites(q)=oneSuite;end
   end,end
   q=q+1;suites(q)=struct('kind','colored','snr_db',plan.colored_stress.snr_db,'seed',plan.colored_stress.seed,'ar1_rho',plan.colored_stress.ar1_rho);
   for is=1:numel(suites)
    su=suites(is);[D,~]=p1_validation_observation(scene,c,truthModels,map,su.seed,su);
    [Sfixed,~,~]=fixed_v1_apply(D,original);
    [~,prior,~,ctx]=p1_solve_branch(D,models,map,cm,domain,R,'CombinedUncertainty');outBase=joint_select(ctx,'DecoupledWeak',.3,{});Sbase=outBase.split;
    audit=uq_audit_outputs(outBase);assert(audit.flags==0&&audit.max_normal<2e-7);linearflags=linearflags+audit.flags;maxlinear=max(maxlinear,audit.max_normal);
    Dwork=D;
    if id==1
     pilot=prior.pilot;lam=prior.lambda;sigma=prior.sigma;
     % prior.pilot may have already incorporated coherent correction. Rebuild
     % the original first pilot to match the public refit API exactly.
     for j=1:11,F=models{j}.direct;pilot(:,j)=F.V*(F.z./(F.z.^2+4*lam(j)).*(F.U'*D(:,j)));end
     [Dwork,~,~,~,~]=p1_noise_uncertainty(D,pilot,models,sigma,lam,domain,R,spec.p1_hyperparameters);
    else,sigma=prior.sigma;end
    out=reflectivity_refit_prepared(Dwork,dicts,sigma);Snew=out.gather;flags=flags+out.optimization_flags;maxnormal=max(maxnormal,out.normal_residual);
    caseid=sprintf('%s_%s_%s_%ddB_%d',m.name,name,su.kind,su.snr_db,su.seed);
    for method=1:3
     label='FixedSplit';S=Sfixed;normal=NaN;oflag=0;if method==2,label='DecoupledWeak';S=Sbase;normal=audit.max_normal;end
     if method==3,label='FullCalRefit';S=Snew;normal=out.normal_residual;oflag=out.optimization_flags;end
     score=score_stability(S,scene,c);row=struct('CaseID',caseid,'Model',m.name,'Domain',domain,'Scenario',name,'Suite',su.kind,'SNR_dB',su.snr_db,'Seed',su.seed,'Method',label);
     f=fieldnames(score);for i=1:numel(f),row.(f{i})=score.(f{i});end
     row.NormalResidual=normal;row.OptimizationFlags=oflag;count=count+1;if count==1,rows=row;else,rows(count)=row;end
    end
    meta=struct('id',caseid,'model',m.name,'domain',domain,'scenario',name,'suite',su.kind,'snr_db',su.snr_db,'seed',su.seed,'Q_baseline',cm.Q,'Q_new',cmfull.Q,'frequency_baseline',cm.fpeak,'frequency_new',cmfull.fpeak);
    t=c.t;truth=scene.truth;optimized=Snew;baseline=Sbase;fixed=Sfixed;observed=D;diagnostic=rmfield(out,{'gather','reflectivity'});reflectivity=out.reflectivity;
    save(fullfile(outdir,'cases',[caseid '.mat']),'meta','t','truth','optimized','baseline','fixed','observed','diagnostic','reflectivity','-v7');
    uq_writejson(fullfile(outdir,'cases',[caseid '.json']),struct('meta',meta,'t',t,'truth',truth,'observed',observed,'fixed',fixed,'baseline',baseline,'optimized',optimized,'diagnostic',diagnostic));
    writetable(struct2table(rows),fullfile(outdir,'metrics.csv'));cases=cases+1;uq_writejson(fullfile(outdir,'progress.json'),struct('cases',cases,'expected',50,'elapsed_seconds',toc(start)));
    aScore=score_stability(Sbase,scene,c);bScore=score_stability(Snew,scene,c);fprintf('%d/50 %s baseline=%.5f new=%.5f elapsed=%.1fs\n',cases,caseid,aScore.RelativeError,bScore.RelativeError,toc(start));
   end
  end
 end
end
assert(cases==50&&count==150);
uq_writejson(fullfile(outdir,'checks.json'),struct('complete',true,'cases',cases,'rows',count,'optimization_flags',flags,'max_normal',maxnormal,'baseline_flags',linearflags,'max_baseline_normal',maxlinear,'elapsed_seconds',toc(start)));
end
