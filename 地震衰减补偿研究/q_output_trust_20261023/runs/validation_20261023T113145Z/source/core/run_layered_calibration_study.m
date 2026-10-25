function run_layered_calibration_study(outdir)
root=fileparts(mfilename('fullpath'));oldpath=path;oldrng=rng;cleanup=onCleanup(@()restore(oldpath,oldrng));
addpath(fullfile(root,'matlab_src'));assert(~isfolder(outdir));mkdir(outdir);start=tic;
[configs,~]=p1_models('development');c=configs{1};public=uq_public_config(c,c.nominal_frequency);
qs=[30,50,90,70];fp=62;test=public;test.Q=qs;[a,b]=ray_absorption(test.t,750,test);G=attenuation_matrix(a+b,test);
times=cumsum(2*public.h./public.vp(1:3));W=zeros(public.nt,3);
for k=1:3,z=pi*fp*(public.t-times(k));W(:,k)=(1-2*z.^2).*exp(-z.^2);end
B=p1_layered_calibration_templates(public,750,qs,fp);error=norm(B(:,1:2:end)-G*W,'fro')/norm(G*W,'fro');
assert(error<1e-8);uq_writejson(fullfile(outdir,'basis_check.json'),struct('passed',true,'relative_error',error));
trueModels=cell(1,11);
for j=1:11,[a,b]=ray_absorption(c.t,c.offset(j),c);trueModels{j}=struct('G',attenuation_matrix(a+b,c));end
scene=make_stability_scene('SmoothAVO',c);[calD,~]=uq_observation(scene,c,trueModels,1:11,20261201,18);
rows=struct([]);nr=0;
for multiplier=[1,.97,1.03]
 provided=public;provided.vp=public.vp*multiplier;
 [initial,corrected,joint]=p1_joint_calibrate(calD,c.offset,provided);
 [est,diag]=p1_layered_calibrate(calD,c.offset,corrected,initial);
 cm=uq_inverse_config(corrected,est);models=prepare_models_v2(c.offset,cm);R=avo_projector_v2(models);
 uq_writejson(fullfile(outdir,sprintf('calibration_%.2f.json',multiplier)),struct('joint',joint,'layered',diag,'estimate',est));
 for name={'SmoothAVO','WeakSection'}
  [mm,domain,map,scene]=uq_case(name{1},c,models);[D,~]=uq_observation(scene,c,trueModels,map,20261211,18);
  [out,~,~]=p1_solve_branch(D,mm,map,cm,domain,R,'CombinedUncertainty');score=score_stability(out.split,scene,c);
  row=struct('provided_velocity',multiplier,'scenario',name{1},'layered_applied',diag.applied, ...
   'effective_Q',est.effective_Q,'estimated_frequency',est.estimated_frequency,'validation_gain',diag.validation_gain, ...
   'relative_error',score.RelativeError,'weak_error',score.WeakEventError,'normal_residual',out.split_record.normal_residual, ...
   'elapsed_seconds',toc(start));
  nr=nr+1;if nr==1,rows=row;else,rows(nr)=row;end
  writetable(struct2table(rows),fullfile(outdir,'metrics.csv'));
  fprintf('Layered %.3f %s accepted=%d Qeff=%.3f f=%.3f RelErr=%.6f WeakErr=%.6f elapsed=%.1fs\n', ...
   multiplier,name{1},diag.applied,est.effective_Q,est.estimated_frequency,score.RelativeError,score.WeakEventError,toc(start));
 end
end
uq_writejson(fullfile(outdir,'summary.json'),struct('completed',true,'records',nr,'development_only',true,'elapsed_seconds',toc(start)));
end
function restore(p,r),path(p);rng(r);end
