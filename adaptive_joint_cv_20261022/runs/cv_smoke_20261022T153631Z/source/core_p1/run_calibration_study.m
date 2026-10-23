function run_calibration_study(outdir)
root=fileparts(mfilename('fullpath'));oldpath=path;oldrng=rng;cleanup=onCleanup(@()restore(oldpath,oldrng));
addpath(fullfile(root,'matlab_src'));assert(~isfolder(outdir));mkdir(outdir);start=tic;
[configs,~]=p1_models('development');records=struct([]);count=0;
for im=1:numel(configs)
 c=configs{im};public=uq_public_config(c,c.nominal_frequency);
 scene=make_stability_scene('SmoothAVO',c);trueModels=cell(1,11);
 for j=1:11,[a,b]=ray_absorption(c.t,c.offset(j),c);trueModels{j}=struct('G',attenuation_matrix(a+b,c));end
 [D,~]=uq_observation(scene,c,trueModels,1:11,20261201,18);
 for multiplier=[1,1.03,.97]
  provided=public;provided.vp=public.vp*multiplier;
  [est,corrected,diag]=p1_joint_calibrate(D,c.offset,provided);
  record=struct('model',c.model_name,'applied_velocity_error',multiplier, ...
   'estimated_velocity_correction',est.velocity_correction, ...
   'remaining_velocity_ratio',mean(corrected.vp./public.vp), ...
   'effective_Q',est.effective_Q,'frequency',est.estimated_frequency, ...
   'accepted',diag.applied,'validation_gain',diag.validation_gain, ...
   'q_boundary',est.q_boundary,'frequency_boundary',est.frequency_boundary, ...
   'optimizer_exitflag',est.optimizer_exitflag,'elapsed_seconds',toc(start));
  count=count+1;if count==1,records=record;else,records(count)=record;end
  uq_writejson(fullfile(outdir,sprintf('%s_%.2f.json',c.model_name,multiplier)),diag);
  writetable(struct2table(records),fullfile(outdir,'calibration_metrics.csv'));
  fprintf('%s provided=%.3f correction=%.6f remaining=%.6f validation_gain=%.4f Q=%.3f f=%.3f elapsed=%.1fs\n', ...
   c.model_name,multiplier,est.velocity_correction,record.remaining_velocity_ratio,diag.validation_gain,est.effective_Q,est.estimated_frequency,toc(start));
 end
end
uq_writejson(fullfile(outdir,'summary.json'),struct('completed',true,'cases',count,'elapsed_seconds',toc(start), ...
 'development_only',true,'calibration_seed',20261201,'records',records));
end
function restore(p,r),path(p);rng(r);end
