function out=p1_compensate(D,offsets,public,calibrationD,calibrationOffsets,domain,variant)
% P1 public entry point. Enforces public input whitelist; zero access to true Q, fpeak, rho, vs.
% Outputs Direct and Split using common priors and uncertainty-averaged empirical weights.
if nargin<7 || isempty(variant), variant='GuardedNoise'; end
allowedVariants={'Baseline','Frozen','StrongerDamping','NoiseUncertaintyOnly','QModelErrorOnly','CombinedUncertainty','CombinedUncertainty_NoStructure','TargetMatched','CombinedFixedTarget','GuardedNoise','NoCompensation'};
assert(any(strcmp(variant,allowedVariants)),'Unsupported P1 variant');
root=fileparts(mfilename('fullpath')); oldpath=path; cleanup=onCleanup(@()path(oldpath));
addpath(fullfile(root,'matlab_src'));
allowed={'dt','nt','t','offset','vp','h','fref','nfft','lambda_x','lambda_fixed','lambda_grid','cov_floor','nominal_frequency'};
assert(isempty(setdiff(fieldnames(public),allowed)),'Nonpublic input field detected');
assert(~isfield(public,'Q') && ~isfield(public,'fpeak') && ~isfield(public,'rho') && ~isfield(public,'vs'));
assert(isnumeric(D) && all(isfinite(D(:))) && size(D,1)==public.nt);
assert(size(D,2)==numel(offsets));
assert(size(calibrationD,1)==public.nt && size(calibrationD,2)==numel(calibrationOffsets));
assert(all(isfinite(calibrationD(:))) && size(calibrationD,2)>=11,'Calibration requires 11-trace layout');
est=p1_calibrate_ensemble(calibrationD,calibrationOffsets,public);
cm=uq_inverse_config(public,est);
if strcmp(variant,'GuardedNoise') && ~cm.calibration_reliable
 [out,prior,preset]=p1_solve_branch(D,{},[],cm,domain,sparse(numel(D),numel(D)),variant);
 out.estimated_Q=est; out.prior=prior; out.preset=preset;
 out.selected_method='NoCompensation_Split'; out.requested_method=variant; out.gather=out.split;
 return;
end
allmodels=prepare_models_v2(public.offset,cm);
R_avo=avo_projector_v2(allmodels);
map=ones(1,numel(offsets));
for j=1:numel(offsets)
 [~,ix]=min(abs(public.offset-offsets(j))); map(j)=ix;
end
models=allmodels(map);
[out,prior,preset]=p1_solve_branch(D,models,map,cm,domain,R_avo,variant);
out.estimated_Q=est; out.prior=prior; out.preset=preset;
actual=variant; if isfield(out,'effective_method'), actual=out.effective_method; end
out.selected_method=[actual '_Split']; out.requested_method=variant; out.gather=out.split;
end
