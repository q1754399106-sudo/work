function out=frozen_compensate(D,offsets,public,calibrationD,calibrationOffsets,domain,variant)
% Retained public entry. No true-Q input. Same prior for Direct and Split.
if nargin<7, variant='Frozen'; end
assert(any(strcmp(variant,{'Baseline','Frozen'})),'Only retained variants allowed');
root=fileparts(mfilename('fullpath')); oldpath=path; cleanup=onCleanup(@()path(oldpath));
addpath(fullfile(root,'matlab_src'));
allowed={'dt','nt','t','offset','vp','h','fref','nfft','lambda_x','lambda_fixed','lambda_grid','cov_floor','nominal_frequency'};
assert(isempty(setdiff(fieldnames(public),allowed)),'Nonpublic input field');
assert(isnumeric(D) && all(isfinite(D(:))) && size(D,1)==public.nt);
assert(size(D,2)==numel(offsets));
assert(size(calibrationD,1)==public.nt && size(calibrationD,2)==numel(calibrationOffsets));
assert(all(isfinite(calibrationD(:))) && size(calibrationD,2)>=11,'Calibration requires original 11-trace layout');
o=uq_spec(); ix=find(strcmp({o.presets.name},variant),1); assert(~isempty(ix));
out=compensate_unknown_q(D,offsets,public,calibrationD,calibrationOffsets,domain,o.presets(ix));
end
