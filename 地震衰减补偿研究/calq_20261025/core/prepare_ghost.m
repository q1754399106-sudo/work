function [est,corrected,diagnostic]=prepare_ghost(calibrationD,offsets,public)
% Same as p1_prepare_candidate (layered mode) but with the GhostCal layered step.
[est,corrected,diagnostic]=p1_joint_calibrate(calibrationD,offsets,public);
[est,layered]=p1_layered_calibrate_ghost(calibrationD,offsets,corrected,est);diagnostic.layered=layered;
end
