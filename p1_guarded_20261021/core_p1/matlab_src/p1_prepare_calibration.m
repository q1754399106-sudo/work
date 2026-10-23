function [est,cm,models,R]=p1_prepare_calibration(D,offsets,public)
% Rebuild calibrated operators and ray projectors together, on every fit.
est=p1_calibrate_ensemble(D,offsets,public);
cm=uq_inverse_config(public,est);
models=prepare_models_v2(public.offset,cm);
R=avo_projector_v2(models);
end
