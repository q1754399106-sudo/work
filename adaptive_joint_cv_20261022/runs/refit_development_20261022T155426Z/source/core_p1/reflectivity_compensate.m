function out=reflectivity_compensate(D,offsets,public,calibrationD,calibrationOffsets,domain)
% One full-calibration rule across all models. Research status until validated.
root=fileparts(mfilename('fullpath'));oldpath=path;cleanup=onCleanup(@()path(oldpath));addpath(fullfile(root,'matlab_src'));
allowed={'dt','nt','t','offset','vp','h','fref','nfft','lambda_x','lambda_fixed','lambda_grid','cov_floor','nominal_frequency'};
assert(isempty(setdiff(fieldnames(public),allowed)),'Nonpublic input field detected');
assert(any(strcmp(domain,{'Offset','Section'})));
assert(isequal(size(D),[public.nt,11])&&numel(offsets)==11&&all(isfinite(D(:))));
assert(isequal(size(calibrationD),[public.nt,11])&&numel(calibrationOffsets)==11&&all(isfinite(calibrationD(:))));
policy=jsondecode(fileread(fullfile(root,'inputs','candidate_policy.json')));
[est,corrected,cal]=p1_prepare_candidate(calibrationD,calibrationOffsets,public,policy);cm=uq_inverse_config(corrected,est);
dicts=cell(1,11);sigma=zeros(1,11);
for j=1:11
 previous=find(abs(offsets(1:j-1)-offsets(j))<1e-9,1);
 if isempty(previous),dicts{j}=refit_dictionary(cm,offsets(j));else,dicts{j}=dicts{previous};end
 ni=estimate_trace_snr(D(:,j),cm.t);sigma(j)=ni.sigma;
end
Dwork=D;
if strcmp(domain,'Offset')
 models=prepare_models_v2(offsets,cm);R=avo_projector_v2(models);pilot=zeros(size(D));lam=sigma;
 for j=1:11
  F=models{j}.direct;[~,lam(j)]=discrepancy_apply(F,D(:,j),cm.lambda_grid,sqrt(cm.nt)*sigma(j));
  pilot(:,j)=F.V*(F.z./(F.z.^2+4*lam(j)).*(F.U'*D(:,j)));
 end
 spec=jsondecode(fileread(fullfile(root,'inputs','p1_experiment_spec.json')));
 [Dwork,~,~,~,~]=p1_noise_uncertainty(D,pilot,models,sigma,lam,domain,R,spec.p1_hyperparameters);
end
out=reflectivity_refit_prepared(Dwork,dicts,sigma);
out.estimated_Q=est;out.calibration_diagnostics=cal;out.selected_method='FullCalRefit';out.compensation_applied=true;
end
