function [est,corrected,diagnostic]=prepare_sharedghost(calibrationD,offsets,public)
% Same as p1_prepare_candidate (layered mode) but with the SharedGhost layered step.
[est,corrected,diagnostic]=p1_joint_calibrate(calibrationD,offsets,public);
[est,layered]=p1_layered_calibrate_sharedghost(calibrationD,offsets,corrected,est);diagnostic.layered=layered;
end
