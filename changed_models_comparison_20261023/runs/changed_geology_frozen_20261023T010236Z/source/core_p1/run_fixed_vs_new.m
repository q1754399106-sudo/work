function run_fixed_vs_new(outdir)
root=fileparts(mfilename('fullpath'));oldpath=path;oldrng=rng;cleanup=onCleanup(@()restore(oldpath,oldrng));
addpath(fullfile(root,'matlab_src'));assert(~isfolder(outdir));mkdir(outdir);mkdir(fullfile(outdir,'cases'));start=tic;
plan=jsondecode(fileread(fullfile(root,'inputs','fixed_vs_new_plan.json')));
policy=jsondecode(fileread(fullfile(root,'inputs','candidate_policy.json')));
uq_writejson(fullfile(outdir,'executed_plan.json'),plan);uq_writejson(fullfile(outdir,'policy.json'),policy);
base=config_model();base.model_name='OriginalReference';base.nominal_frequency=60;
spec=jsondecode(fileread(fullfile(root,'inputs','p1_experiment_spec.json')));
cached=jsondecode(fileread(fullfile(root,'inputs','cached_shrinkage','trial_plan.json')));
configs={base};
for k=2:numel(plan.models)
 name=plan.models{k};list=cached.additional_models;ix=find(strcmp({list.name},name),1);assert(~isempty(ix));m=list(ix);
 c=config_model();for key={'vp','vs','rho','Q','event_t','fpeak'},f=key{1};c.(f)=reshape(m.(f),1,[]);end
 c.h=diff([0,c.event_t]).*c.vp(1:3)/2;c.model_name=m.name;c.nominal_frequency=m.nominal_frequency;configs{end+1}=c;
end
rows=struct([]);count=0;observations=0;flags=0;maxnormal=0;oldnormal=0;oldreplay=0;
for im=1:numel(configs)
 c=configs{im};public=uq_public_config(c,c.nominal_frequency);truthModels=cell(1,11);
 for j=1:11,[a,b]=ray_absorption(c.t,c.offset(j),c);truthModels{j}=struct('G',attenuation_matrix(a+b,c));end
 calScene=make_stability_scene('SmoothAVO',c);
 [calD,~]=uq_observation(calScene,c,truthModels,1:11,plan.calibration_seed,plan.calibration_snr_db);
 [baseEst,corrected,cal]=p1_prepare_candidate(calD,c.offset,public,policy);
 save(fullfile(outdir,['calibration_input_' c.model_name '.mat']),'calD','public','-v7');
 uq_writejson(fullfile(outdir,['calibration_' c.model_name '.json']),struct('estimate',baseEst,'diagnostics',cal,'public',public,'corrected',corrected));
 fprintf('Calibration complete: %s elapsed=%.1fs\n',c.model_name,toc(start));
 for id=1:2
  domain='Offset';offsets=c.offset;map=1:11;if id==2,domain='Section';offsets=repmat(c.offset(6),1,11);map=6*ones(1,11);end
  est=p1_candidate_domain_estimate(cal.fixed_estimate,baseEst,corrected,policy,domain);
  cm=uq_inverse_config(corrected,est);cm.candidate_profile=policy.selected_profile;cm.joint_prepare_only=true;
  if id==1,models=prepare_models_v2(offsets,cm);R=avo_projector_v2(models);
  else,one=prepare_models_v2(offsets(1),cm);models=repmat(one,1,11);R=sparse(c.nt*11,c.nt*11);end
  original=fixed_v1_operators(cm,offsets);
  for isc=1:numel(plan.scenarios)
   name=plan.scenarios{isc};if endsWith(name,'Section')~=(id==2),continue;end
   scene=make_stability_scene(name,c);
   for snr=reshape(plan.snr_db,1,[])
    for seed=reshape(plan.seeds,1,[])
     [D,trueSigma]=uq_observation(scene,c,truthModels,map,seed,snr);Dshared=D;
     [Sfixed,normOld,replay]=fixed_v1_apply(D,original);
     [~,prior,preset,ctx]=p1_solve_branch(D,models,map,cm,domain,R,'CombinedUncertainty');
     out=joint_select(ctx,plan.new_mode,plan.beta_factor,{});Snew=out.split;
     assert(isequal(D,Dshared));a=uq_audit_outputs(out);assert(a.flags==0&&a.max_normal<2e-7);
     flags=flags+a.flags;maxnormal=max(maxnormal,a.max_normal);oldnormal=max(oldnormal,normOld);oldreplay=max(oldreplay,replay);
     caseid=sprintf('%s_%s_%ddB_%d',c.model_name,name,snr,seed);
     scoreOld=score_stability(Sfixed,scene,c);scoreNew=score_stability(Snew,scene,c);
     for method=1:2
      label='FixedSplit';score=scoreOld;nr=normOld;scale=NaN;endpoint=false;mismatch=NaN;beta=NaN;
      if method==2,label='DecoupledWeak';score=scoreNew;nr=out.split_record.normal_residual;scale=out.scale;endpoint=out.endpoint;mismatch=out.residual_mismatch;beta=out.beta;end
      row=struct('CaseID',caseid,'Model',c.model_name,'Domain',domain,'Scenario',name,'SNR_dB',snr,'Seed',seed,'Method',label);
      fields=fieldnames(score);for f=1:numel(fields),key=fields{f};row.(key)=score.(key);end
      row.NormalResidual=nr;row.FinalFlag=0;row.Scale=scale;row.ScaleEndpoint=endpoint;row.TargetMismatch=mismatch;row.Beta=beta;
      row.Q1=cm.Q(1);row.Q2=cm.Q(2);row.Q3=cm.Q(3);row.Q4=cm.Q(4);row.EstimatedFrequency=cm.fpeak;
      row.VelocityScale=mean(cm.vp./public.vp);row.AllNewSelectionMaxNormal=a.max_normal;
      count=count+1;if count==1,rows=row;else,rows(count)=row;end
     end
     meta=struct('id',caseid,'model',c.model_name,'scenario',name,'domain',domain,'snr_db',snr,'seed',seed, ...
      'Q_used',cm.Q,'velocity_used',cm.vp,'estimated_frequency',cm.fpeak,'old_score',scoreOld,'new_score',scoreNew, ...
      'normal_old',normOld,'normal_new',out.split_record.normal_residual,'alpha',out.scale,'beta',out.beta,'target_mismatch',out.residual_mismatch);
     plotdata=struct('meta',meta,'t',c.t,'offsets',offsets,'lateral',scene.lateral,'truth',scene.truth,'observed',D,'fixed',Sfixed,'new',Snew);
     uq_writejson(fullfile(outdir,'cases',[caseid '.json']),plotdata);
     truth=scene.truth;t=c.t;save(fullfile(outdir,'cases',[caseid '.mat']),'meta','t','offsets','truth','D','Sfixed','Snew','-v7');
     writetable(struct2table(rows),fullfile(outdir,'metrics.csv'));
     observations=observations+1;uq_writejson(fullfile(outdir,'progress.json'),struct('observations',observations,'rows',count,'elapsed_seconds',toc(start),'max_new_normal',maxnormal));
     fprintf('%s obs=%d/%d old=%.5f new=%.5f elapsed=%.1fs\n',caseid,observations,plan.expected_observations,scoreOld.RelativeError,scoreNew.RelativeError,toc(start));
    end
   end
  end
 end
end
assert(observations==plan.expected_observations&&count==plan.expected_rows);
uq_writejson(fullfile(outdir,'checks.json'),struct('completed',true,'observations',observations,'rows',count,'new_final_solver_flags',flags, ...
 'max_new_normal_residual',maxnormal,'max_old_normal_residual',oldnormal,'max_old_composition_replay',oldreplay, ...
 'same_observation_assertions',observations,'shared_calibration',true,'elapsed_seconds',toc(start),'matlab_version',version));
end
function restore(p,r),path(p);rng(r);end
