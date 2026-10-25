function out=combined_uncertainty_v2(D,offsets,public,calibrationD,calibrationOffsets,domain)
% Candidate entry. Policy must be selected on development data and frozen.
root=fileparts(mfilename('fullpath'));oldpath=path;cleanup=onCleanup(@()path(oldpath));
addpath(fullfile(root,'matlab_src'));
allowed={'dt','nt','t','offset','vp','h','fref','nfft','lambda_x','lambda_fixed','lambda_grid','cov_floor','nominal_frequency'};
assert(isempty(setdiff(fieldnames(public),allowed)),'Nonpublic input field detected');
assert(any(strcmp(domain,{'Offset','Section'})),'Explicit acquisition domain required');
assert(size(D,1)==public.nt&&size(D,2)==numel(offsets)&&all(isfinite(D(:))));
assert(size(calibrationD,1)==public.nt&&size(calibrationD,2)==numel(calibrationOffsets));
assert(numel(calibrationOffsets)==11 && size(D,2)==11,'Current implementation requires the validated 11-trace layout');
policyfile=fullfile(root,'inputs','candidate_policy.json');
assert(isfile(policyfile),'Development selection has not been completed; no policy is silently assumed');
policy=jsondecode(fileread(policyfile));
[est,corrected,calibration]=p1_joint_calibrate(calibrationD,calibrationOffsets,public);
if isfield(policy,'calibration_mode') && strcmp(policy.calibration_mode,'layered')
 [est,layered]=p1_layered_calibrate(calibrationD,calibrationOffsets,corrected,est);calibration.layered=layered;
end
cm=uq_inverse_config(corrected,est);cm.candidate_profile=policy.selected_profile;
allmodels=prepare_models_v2(public.offset,cm);R=avo_projector_v2(allmodels);
map=zeros(1,numel(offsets));
for j=1:numel(offsets),[distance,map(j)]=min(abs(public.offset-offsets(j)));assert(distance<1e-6,'Offsets must match the public model grid');end
[out,prior,preset]=p1_solve_branch(D,allmodels(map),map,cm,domain,R,'CombinedUncertainty');
out.estimated_Q=est;out.calibration_diagnostics=calibration;out.prior=prior;out.preset=preset;
out.policy=policy;out.estimated_velocity_scale=est.velocity_correction;
out.selected_method='CombinedUncertaintyV2_Split';out.gather=out.split;
out.compensation_applied=true;
end
