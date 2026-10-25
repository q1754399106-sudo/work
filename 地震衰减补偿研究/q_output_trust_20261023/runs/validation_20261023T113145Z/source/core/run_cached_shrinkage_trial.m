function run_cached_shrinkage_trial(outdir)
% Re-screen Section calibration shrinkage using already fitted dev calibrations.
root=fileparts(mfilename('fullpath'));oldpath=path;oldrng=rng;cleanup=onCleanup(@()restore(oldpath,oldrng));
addpath(fullfile(root,'matlab_src'));assert(~isfolder(outdir));mkdir(outdir);start=tic;
input=fullfile(root,'inputs','cached_shrinkage');plan=jsondecode(fileread(fullfile(input,'trial_plan.json')));
uq_writejson(fullfile(outdir,'executed_plan.json'),plan);[configs,~]=p1_models('development');
for i=1:numel(plan.additional_models)
 m=plan.additional_models(i);c=config_model();for key={'vp','vs','rho','Q','event_t','fpeak'},k=key{1};c.(k)=reshape(m.(k),1,[]);end
 c.h=diff([0,c.event_t]).*c.vp(1:3)/2;c.model_name=m.name;c.nominal_frequency=m.nominal_frequency;configs{end+1}=c;
end
rows=struct([]);count=0;flags=0;maxnormal=0;R=sparse(601*11,601*11);
for im=1:numel(configs)
 c=configs{im};public=uq_public_config(c,c.nominal_frequency);cseed=plan.calibration_seeds(1);
 d=jsondecode(fileread(fullfile(input,sprintf('calibration_%s.json',c.model_name))));corrected=public;corrected.vp=public.vp*d.actual_velocity_scale;
 [a,b]=ray_absorption(c.t,c.offset(6),c);trueModels=repmat({struct('G',attenuation_matrix(a+b,c))},1,11);map=6*ones(1,11);
 cached_weight=NaN;
 for ip=1:numel(plan.profiles)
  profile=plan.profiles(ip);weight=profile.section_calibration_weight;
  if ~isequal(cached_weight,weight)
   est=p1_shrink_calibration(d.fixed_estimate,d.layered_estimate,corrected,weight);
   cm=uq_inverse_config(corrected,est);one=prepare_models_v2(c.offset(6),cm);mm=repmat(one,1,11);cached_weight=weight;
  end
  cm.candidate_profile=profile;
  for isu=1:numel(plan.suites)
   su=plan.suites(isu);if su.vp_multiplier~=1,continue;end
   for isc=1:numel(su.scenarios)
    scenario=su.scenarios{isc};if ~endsWith(scenario,'Section'),continue;end
    scene=make_stability_scene(scenario,c);
    for seed=reshape(plan.test_seeds,1,[])
     [D,trueSigma]=uq_observation(scene,c,trueModels,map,seed,su.snr_db);
     [out,prior,preset]=p1_solve_branch(D,mm,map,cm,'Section',R,'CombinedUncertainty');quality=uq_audit_outputs(out);preset.name=profile.name;
     flags=flags+quality.flags;maxnormal=max(maxnormal,quality.max_normal);
     for ir=1:2
      route='Direct';S=out.direct;if ir==2,route='Split';S=out.split;end
      row=uq_metric_row(S,scene,c,c.model_name,su.name,seed,profile.name,route,out,prior,preset,est);
      row.CalibrationSeed=cseed;row.VelocityCorrection=d.actual_velocity_scale;row.EstimatedSigmaRatio=mean(prior.sigma./trueSigma);
      row.AllSelectionSolverFlags=quality.flags;row.AllSelectionMaxNormalResidual=quality.max_normal;row.CalibrationBlend=weight;
      row.ComputationAliasOf='';row.RowProvenance='new_MATLAB_solve_reused_calibration';
      count=count+1;if count==1,rows=row;else,rows(count)=row;end
     end
     writetable(struct2table(rows),fullfile(outdir,'metrics_section_only.csv'));
     fprintf('%s %s %s %s rows=%d elapsed=%.1fs normal=%.3e flags=%d\n',c.model_name,profile.name,su.name,scenario,count,toc(start),maxnormal,flags);
    end
   end
  end
 end
end
uq_writejson(fullfile(outdir,'checks.json'),struct('completed',true,'section_rows',count,'solver_flags',flags,'max_normal_residual',maxnormal,'development_only',true));
end
function restore(p,r),path(p);rng(r);end
