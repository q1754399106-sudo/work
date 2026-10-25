function out=compensate_unknown_q(D,offsets,public,calibrationD,calibrationOffsets,domain,preset)
% Production entry: no true Q accepted. Known kinematic model is required.
assert(~isfield(public,'Q') && ~isfield(public,'fpeak') && ~isfield(public,'rho'));
if nargin<7 || isempty(preset), o=uq_spec(); preset=o.presets(2); end
est=estimate_q_effective(calibrationD,calibrationOffsets,public); c=uq_inverse_config(public,est);
models=prepare_models_v2(offsets,c); prior=uq_build_prior(D,models,c,domain,preset);
out=uq_select(D,models,prior,preset); out.estimated_Q=est; out.preset=preset;
out.prior=prior; out.selected_method=['UnknownQ_' preset.name '_Split']; out.gather=out.split;
end
