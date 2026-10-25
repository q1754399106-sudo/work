function run_validation_job(jobfile,outdir)
% Frozen-candidate validation. Synthesis/scoring are the only truth consumers.
root=fileparts(mfilename('fullpath'));oldpath=path;oldrng=rng;cleanup=onCleanup(@()restore(oldpath,oldrng));
addpath(fullfile(root,'matlab_src'));assert(~isfolder(outdir));mkdir(outdir);start=tic;
job=jsondecode(fileread(jobfile));policy=jsondecode(fileread(fullfile(root,'inputs','candidate_policy.json')));
uq_writejson(fullfile(outdir,'executed_job.json'),job);uq_writejson(fullfile(outdir,'executed_policy.json'),policy);
c=config_model();m=job.model;
for key={'vp','vs','rho','Q','event_t','fpeak'},k=key{1};c.(k)=reshape(m.(k),1,[]);end
c.h=diff([0,c.event_t]).*c.vp(1:3)/2;c.model_name=m.name;c.nominal_frequency=m.nominal_frequency;
public=uq_public_config(c,c.nominal_frequency);trueModels=cell(1,11);
for j=1:11,[a,b]=ray_absorption(c.t,c.offset(j),c);trueModels{j}=struct('G',attenuation_matrix(a+b,c));end
calScene=make_stability_scene('SmoothAVO',c);count=0;observations=0;replaycount=0;flags=0;maxnormal=0;
velocity_groups=[1,1.03,.97];
for ic=1:numel(job.calibration_seeds)
 cseed=job.calibration_seeds(ic);[calD,~]=uq_observation(calScene,c,trueModels,1:11,cseed,18);
 save(fullfile(outdir,sprintf('calibration_input_%d.mat',cseed)),'calD','public','-v7');
 for vg=velocity_groups
  active=find(abs([job.suites.vp_multiplier]-vg)<1e-12);if isempty(active),continue;end
  provided=public;provided.vp=public.vp*vg;
  [est,corrected,calibration]=p1_prepare_candidate(calD,c.offset,provided,policy);
  oldest=calibration.fixed_estimate;oldcm=uq_inverse_config(provided,oldest);cm=uq_inverse_config(corrected,est);
  oldmodels=prepare_models_v2(c.offset,oldcm);oldR=avo_projector_v2(oldmodels);
  if isequal(oldcm,cm),models=oldmodels;R=oldR;else,models=prepare_models_v2(c.offset,cm);R=avo_projector_v2(models);end
  sectionest=p1_candidate_domain_estimate(oldest,est,corrected,policy,'Section');sectioncm=uq_inverse_config(corrected,sectionest);
  if isequal(sectioncm,cm),sectionmodels=models;else,sectionmodels=prepare_models_v2(c.offset,sectioncm);end
  uq_writejson(fullfile(outdir,sprintf('calibration_%d_v%.2f.json',cseed,vg)),struct('calibration',calibration,'estimate',est,'section_candidate_estimate',sectionest,'provided_vp',provided.vp,'corrected_vp',corrected.vp));
  for isu=reshape(active,1,[])
   su=job.suites(isu);
   for isc=1:numel(su.scenarios)
    scname=su.scenarios{isc};[mm,domain,map,scene]=uq_case(scname,c,models);oldmm=oldmodels(map);
    for iseed=1:numel(job.test_seeds)
     seed=job.test_seeds(iseed);[D,trueSigma]=p1_validation_observation(scene,c,trueModels,map,seed,su);
     replay=job.reuse_original && vg~=.97 && ic==1 && iseed==1;
     needold=~job.reuse_original || vg==.97 || replay;
     methods={'FrozenJointCalibration','CombinedJointCalibration','CombinedUncertaintyV2'};
     if needold,methods=[{'Frozen','CombinedUncertainty'},methods];end
     block=struct([]);replays=struct([]);nb=0;np=0;shared_combined=[];shared_prior=[];shared_preset=[];
     example=struct('t',c.t,'observed',D(:,6),'truth',scene.truth(:,6),'trace_index',6,'model',c.model_name,'suite',su.name,'scenario',scname,'calibration_seed',cseed,'seed',seed);
     for imethod=1:numel(methods)
      method=methods{imethod};localcm=cm;models_use=mm;R_use=R;est_use=est;branch='CombinedUncertainty';alias='';
      if any(strcmp(method,{'Frozen','CombinedUncertainty'}))
       localcm=oldcm;models_use=oldmm;R_use=oldR;est_use=oldest;branch=method;
      elseif strcmp(method,'FrozenJointCalibration'),branch='Frozen';
      elseif strcmp(method,'CombinedUncertaintyV2')
       if strcmp(domain,'Section'),localcm=sectioncm;models_use=sectionmodels(map);est_use=sectionest;end
       localcm.candidate_profile=policy.selected_profile;
      end
      calibration_identity=~isfield(policy.selected_profile,'section_calibration_weight') || policy.selected_profile.section_calibration_weight==1;
      identity=strcmp(method,'CombinedUncertaintyV2') && (strcmp(domain,'Offset') || ...
       (calibration_identity && policy.selected_profile.confidence_floor==1 && policy.selected_profile.section_modal_blend==1 && policy.selected_profile.section_target<0 && policy.selected_profile.section_gamma==10));
      if identity
       out=shared_combined;prior=shared_prior;preset=shared_preset;alias='CombinedJointCalibration';
      else
       [out,prior,preset]=p1_solve_branch(D,models_use,map,localcm,domain,R_use,branch);
      end
      if strcmp(method,'CombinedJointCalibration'),shared_combined=out;shared_prior=prior;shared_preset=preset;end
      preset.name=method;quality=uq_audit_outputs(out);flags=flags+quality.flags;maxnormal=max(maxnormal,quality.max_normal);
      example.(method)=out.split(:,6);
      for ir=1:2
       route='Direct';S=out.direct;rec=out.direct_record;if ir==2,route='Split';S=out.split;rec=out.split_record;end
       row=uq_metric_row(S,scene,c,c.model_name,su.name,seed,method,route,out,prior,preset,est_use);
       row.CalibrationSeed=cseed;row.EstimatedSigmaRatio=mean(prior.sigma./trueSigma);row.Group=job.group;
       row.InitialFlag=rec.initial_flag;row.Iterations=rec.iterations;row.PriorSigmaSum=sum(prior.sigma);row.PriorLambdaSum=sum(prior.lambda);row.PriorHFro=norm(prior.H,'fro');
       row.CalibrationQBoundary=est_use.q_boundary;row.CalibrationFrequencyBoundary=est_use.frequency_boundary;row.CalibrationOptimizerFlag=est_use.optimizer_exitflag;
       row.AllSelectionSolverFlags=quality.flags;row.AllSelectionMaxNormalResidual=quality.max_normal;
       row.VelocityCorrection=mean(localcm.vp./provided.vp);row.LayerQ1=localcm.Q(1);row.LayerQ2=localcm.Q(2);row.LayerQ3=localcm.Q(3);row.LayerQ4=localcm.Q(4);
       row.ComputationAliasOf=alias;row.CompensationApplied=1;row.DataEnergy=sum(D(:).^2);
       if replay && any(strcmp(method,{'Frozen','CombinedUncertainty'}))
        np=np+1;if np==1,replays=row;else,replays(np)=row;end
       else
        nb=nb+1;if nb==1,block=row;else,block(nb)=row;end
       end
      end
     end
     count=append_rows(block,fullfile(outdir,'metrics.csv'),count);
     replaycount=append_rows(replays,fullfile(outdir,'replay_metrics.csv'),replaycount);
     observations=observations+1;
     if ic==1 && iseed==1 && strcmp(su.name,'Main') && any(strcmp(scname,{'WeakSection','CoherentNoiseAVO'}))
      uq_writejson(fullfile(outdir,['example_' scname '.json']),example);
     end
     uq_writejson(fullfile(outdir,'progress.json'),struct('observations',observations,'rows',count,'replay_rows',replaycount, ...
      'elapsed_seconds',toc(start),'solver_flags',flags,'max_normal_residual',maxnormal));
     fprintf('%s %s %s cal=%d seed=%d observations=%d rows=%d elapsed=%.1fs normal=%.3e flags=%d\n', ...
      c.model_name,su.name,scname,cseed,seed,observations,count,toc(start),maxnormal,flags);
    end
   end
  end
 end
end
uq_writejson(fullfile(outdir,'checks.json'),struct('completed',true,'observations',observations,'rows',count,'replay_rows',replaycount, ...
 'solver_flags',flags,'max_normal_residual',maxnormal,'elapsed_seconds',toc(start),'matlab_version',version));
end
function count=append_rows(rows,file,count)
if isempty(rows),return;end
if count==0,writetable(struct2table(rows),file);else,writetable(struct2table(rows),file,'WriteMode','append','WriteVariableNames',false);end
count=count+numel(rows);
end
function restore(p,r),path(p);rng(r);end
