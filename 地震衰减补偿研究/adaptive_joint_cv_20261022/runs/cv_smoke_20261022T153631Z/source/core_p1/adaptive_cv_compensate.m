function out=adaptive_cv_compensate(D,offsets,public,calibrationD,calibrationOffsets,domain)
% Unified research entry: data-driven alpha/beta, no model-specific options.
root=fileparts(mfilename('fullpath'));oldpath=path;cleanup=onCleanup(@()path(oldpath));addpath(fullfile(root,'matlab_src'));
allowed={'dt','nt','t','offset','vp','h','fref','nfft','lambda_x','lambda_fixed','lambda_grid','cov_floor','nominal_frequency'};
assert(isempty(setdiff(fieldnames(public),allowed)),'Nonpublic input field detected');
assert(any(strcmp(domain,{'Offset','Section'})));
assert(isequal(size(D),[public.nt,11])&&numel(offsets)==11&&all(isfinite(D(:))));
assert(isequal(size(calibrationD),[public.nt,11])&&numel(calibrationOffsets)==11&&all(isfinite(calibrationD(:))));
policy=jsondecode(fileread(fullfile(root,'inputs','candidate_policy.json')));
[base,corrected,cal]=p1_prepare_candidate(calibrationD,calibrationOffsets,public,policy);
est=p1_candidate_domain_estimate(cal.fixed_estimate,base,corrected,policy,domain);cm=uq_inverse_config(corrected,est);
modelsAll=prepare_models_v2(public.offset,cm);R=avo_projector_v2(modelsAll);map=zeros(1,11);
for j=1:11,[distance,map(j)]=min(abs(public.offset-offsets(j)));assert(distance<1e-6);end
cm.candidate_profile=policy.selected_profile;models=modelsAll(map);
out=adaptive_cv_prepared(D,models,map,cm,domain,R);
out.estimated_Q=est;out.calibration_diagnostics=cal;out.policy=policy;
end
