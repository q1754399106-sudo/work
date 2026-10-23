function out=joint_compensate(D,offsets,public,calibrationD,calibrationOffsets,domain,options)
% Research API; configuration is explicit, never silently promoted to default.
% options.mode: DecoupledWeak, Decoupled, SharedIsotropic or Physical.
% options.beta_factor is fixed while alpha is selected by discrepancy.
root=fileparts(mfilename('fullpath'));oldpath=path;cleanup=onCleanup(@()path(oldpath));
addpath(fullfile(root,'matlab_src'));
allowed={'dt','nt','t','offset','vp','h','fref','nfft','lambda_x','lambda_fixed','lambda_grid','cov_floor','nominal_frequency'};
assert(isempty(setdiff(fieldnames(public),allowed)),'Nonpublic input field detected');
assert(any(strcmp(domain,{'Offset','Section'})));
assert(size(D,1)==public.nt&&size(D,2)==11&&numel(offsets)==11&&all(isfinite(D(:))));
assert(size(calibrationD,1)==public.nt&&size(calibrationD,2)==11&&numel(calibrationOffsets)==11&&all(isfinite(calibrationD(:))));
assert(any(strcmp(options.mode,{'DecoupledWeak','Decoupled','SharedIsotropic','Physical'})));
assert(isscalar(options.beta_factor)&&isfinite(options.beta_factor)&&options.beta_factor>=0);
policy=jsondecode(fileread(fullfile(root,'inputs','candidate_policy.json')));
[base,corrected,cal]=p1_prepare_candidate(calibrationD,calibrationOffsets,public,policy);
est=p1_candidate_domain_estimate(cal.fixed_estimate,base,corrected,policy,domain);
cm=uq_inverse_config(corrected,est);allmodels=prepare_models_v2(public.offset,cm);R=avo_projector_v2(allmodels);
map=zeros(1,11);for j=1:11,[distance,map(j)]=min(abs(public.offset-offsets(j)));assert(distance<1e-6);end
cm.candidate_profile=policy.selected_profile;cm.joint_prepare_only=true;
[~,prior,preset,ctx]=p1_solve_branch(D,allmodels(map),map,cm,domain,R,'CombinedUncertainty');
ensemble={};if any(strcmp(options.mode,{'Physical','SharedIsotropic'})),ensemble=joint_physical_ensemble(cm,offsets);end
out=joint_select(ctx,options.mode,options.beta_factor,ensemble);
out.gather=out.split;out.prior=prior;out.preset=preset;out.estimated_Q=est;
out.calibration_diagnostics=cal;out.policy=policy;out.options=options;out.compensation_applied=true;
out.selected_method=[options.mode '_Split'];out.development_only=true;
end
