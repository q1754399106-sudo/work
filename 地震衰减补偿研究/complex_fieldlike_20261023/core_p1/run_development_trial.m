function run_development_trial(outdir)
root=fileparts(mfilename('fullpath'));oldpath=path;oldrng=rng;cleanup=onCleanup(@()restore(oldpath,oldrng));
addpath(fullfile(root,'matlab_src'));assert(~isfolder(outdir));mkdir(outdir);start=tic;
plan=jsondecode(fileread(fullfile(root,'inputs','development_trial.json')));
uq_writejson(fullfile(outdir,'executed_plan.json'),plan);
[base,~]=p1_models('development');configs=base;
for i=1:numel(plan.additional_models)
 m=plan.additional_models(i);c=config_model();
 for key={'vp','vs','rho','Q','event_t','fpeak'},k=key{1};c.(k)=reshape(m.(k),1,[]);end
 c.h=diff([0,c.event_t]).*c.vp(1:3)/2;c.model_name=m.name;c.nominal_frequency=m.nominal_frequency;
 configs{end+1}=c;
end
allrows=struct([]);nr=0;maxnormal=0;flags=0;
for im=1:numel(configs)
 c=configs{im};public=uq_public_config(c,c.nominal_frequency);trueModels=cell(1,11);
 for j=1:11,[a,b]=ray_absorption(c.t,c.offset(j),c);trueModels{j}=struct('G',attenuation_matrix(a+b,c));end
 calScene=make_stability_scene('SmoothAVO',c);
 for cseed=reshape(plan.calibration_seeds,1,[])
  [calD,~]=uq_observation(calScene,c,trueModels,1:11,cseed,18);
  cached_multiplier=NaN;
  for isu=1:numel(plan.suites)
   su=plan.suites(isu);provided=public;provided.vp=public.vp*su.vp_multiplier;
   if ~isequal(cached_multiplier,su.vp_multiplier)
   [est,corrected,diag]=p1_joint_calibrate(calD,c.offset,provided);
   if isfield(plan,'calibration_mode') && strcmp(plan.calibration_mode,'layered')
    [est,layered]=p1_layered_calibrate(calD,c.offset,corrected,est);diag.layered=layered;diag.layered_estimate=est;
   end
   fixed=diag.fixed_estimate;
   oldcm=uq_inverse_config(provided,fixed);oldmodels=prepare_models_v2(c.offset,oldcm);oldR=avo_projector_v2(oldmodels);
   cm=uq_inverse_config(corrected,est);
   if isequal(oldcm,rmfield_if_present(cm,{'velocity_correction','kinematic_correction_applied'}))
    models=oldmodels;R=oldR;
   else
    models=prepare_models_v2(c.offset,cm);R=avo_projector_v2(models);
   end
   cached_multiplier=su.vp_multiplier;
   end
   uq_writejson(fullfile(outdir,sprintf('calibration_%s_%s_%d.json',c.model_name,su.name,cseed)),diag);
   for isc=1:numel(su.scenarios)
    sc=su.scenarios{isc};[mm,domain,map,scene]=uq_case(sc,c,models);oldmm=oldmodels(map);
    for seed=reshape(plan.test_seeds,1,[])
     [D,tsigma]=uq_observation(scene,c,trueModels,map,seed,su.snr_db);
     for ib=1:(3+numel(plan.profiles))
      localcm=cm;mm_use=mm;R_use=R;est_use=est;
      if ib<=2
       localcm=oldcm;mm_use=oldmm;R_use=oldR;est_use=fixed;
       method='Frozen';if ib==2,method='CombinedUncertainty';end
       branch=method;
      elseif ib==3
       method='CombinedJointCalibration';branch='CombinedUncertainty';
      else
       profile=plan.profiles(ib-3);method=profile.name;branch='CombinedUncertainty';localcm.candidate_profile=profile;
      end
      [out,prior,preset]=p1_solve_branch(D,mm_use,map,localcm,domain,R_use,branch);preset.name=method;
      au=uq_audit_outputs(out);flags=flags+au.flags;maxnormal=max(maxnormal,au.max_normal);
      for ir=1:2
       route='Direct';S=out.direct;if ir==2,route='Split';S=out.split;end
       row=uq_metric_row(S,scene,c,c.model_name,su.name,seed,method,route,out,prior,preset,est_use);
       row.CalibrationSeed=cseed;row.EstimatedSigmaRatio=mean(prior.sigma./tsigma);
       row.VelocityCorrection=mean(localcm.vp./provided.vp);row.AllSelectionMaxNormalResidual=au.max_normal;row.AllSelectionSolverFlags=au.flags;
       nr=nr+1;if nr==1,allrows=row;else,allrows(nr)=row;end
      end
     end
     writetable(struct2table(allrows),fullfile(outdir,'metrics.csv'));
     fprintf('%s %s %s seed=%d rows=%d elapsed=%.1fs normal=%.3e flags=%d\n',c.model_name,su.name,sc,seed,nr,toc(start),maxnormal,flags);
    end
   end
  end
 end
end
uq_writejson(fullfile(outdir,'checks.json'),struct('completed',true,'rows',nr,'solver_flags',flags, ...
 'max_normal_residual',maxnormal,'elapsed_seconds',toc(start),'development_only',true));
end
function c=rmfield_if_present(c,names)
for i=1:numel(names),if isfield(c,names{i}),c=rmfield(c,names{i});end,end
end
function restore(p,r),path(p);rng(r);end
