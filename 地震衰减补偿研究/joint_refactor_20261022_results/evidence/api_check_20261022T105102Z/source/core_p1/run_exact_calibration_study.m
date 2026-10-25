function run_exact_calibration_study(outdir)
root=fileparts(mfilename('fullpath'));oldpath=path;oldrng=rng;cleanup=onCleanup(@()restore(oldpath,oldrng));
addpath(fullfile(root,'matlab_src'));assert(~isfolder(outdir));mkdir(outdir);start=tic;
[configs,~]=p1_models('development');c=configs{1};p=uq_public_config(c,c.nominal_frequency);
% Verify the new calibration basis against the actual full propagation matrix.
q=65;fp=55;test=p;test.Q=q*ones(1,4);offset=750;
[aa,bb]=ray_absorption(test.t,offset,test);G=attenuation_matrix(aa+bb,test);times=cumsum(2*p.h./p.vp(1:3));W=zeros(p.nt,3);
for k=1:3,z=pi*fp*(p.t-times(k));W(:,k)=(1-2*z.^2).*exp(-z.^2);end
B=p1_exact_calibration_templates(p,offset,q,fp);basis_error=norm(B(:,1:2:end)-G*W,'fro')/norm(G*W,'fro');
assert(basis_error<1e-8,'Exact forward basis mismatch');
uq_writejson(fullfile(outdir,'basis_check.json'),struct('relative_error',basis_error,'passed',true));
truthModels=cell(1,11);
for j=1:11,[a,b]=ray_absorption(c.t,c.offset(j),c);truthModels{j}=struct('G',attenuation_matrix(a+b,c));end
calScene=make_stability_scene('SmoothAVO',c);[Dcal,~]=uq_observation(calScene,c,truthModels,1:11,20261201,18);
records=struct([]);nr=0;
for multiplier=[1,.97,1.03]
 provided=p;provided.vp=p.vp*multiplier;
 [initial,corrected,joint]=p1_joint_calibrate(Dcal,c.offset,provided);
 [refined,diag]=p1_exact_calibration_refine(Dcal,c.offset,corrected,initial);
 cm=uq_inverse_config(corrected,refined);models=prepare_models_v2(c.offset,cm);R=avo_projector_v2(models);
 [mm,domain,map,scene]=uq_case('SmoothAVO',c,models);[D,~]=uq_observation(scene,c,truthModels,map,20261211,18);
 [out,~,~]=p1_solve_branch(D,mm,map,cm,domain,R,'CombinedUncertainty');
 score=score_stability(out.split,scene,c);
 row=struct('velocity_multiplier',multiplier,'correction',initial.velocity_correction,'initial_Q',initial.effective_Q, ...
  'initial_frequency',initial.estimated_frequency,'refined_Q',refined.effective_Q,'refined_frequency',refined.estimated_frequency, ...
  'relative_error',score.RelativeError,'normal_residual',out.split_record.normal_residual,'elapsed_seconds',toc(start));
 nr=nr+1;if nr==1,records=row;else,records(nr)=row;end
 writetable(struct2table(records),fullfile(outdir,'metrics.csv'));
 uq_writejson(fullfile(outdir,sprintf('exact_%.2f.json',multiplier)),struct('initial_kinematic_fit',joint,'refinement',diag,'estimate',refined));
 fprintf('Exact calibration provided=%.3f Q=%.4f f=%.4f RelErr=%.6f elapsed=%.1fs\n',multiplier,refined.effective_Q,refined.estimated_frequency,score.RelativeError,toc(start));
end
uq_writejson(fullfile(outdir,'summary.json'),struct('completed',true,'cases',nr,'development_only',true,'elapsed_seconds',toc(start),'records',records));
end
function restore(p,r),path(p);rng(r);end
