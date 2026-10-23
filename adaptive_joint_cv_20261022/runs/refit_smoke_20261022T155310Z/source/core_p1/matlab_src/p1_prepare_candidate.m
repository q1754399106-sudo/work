function [est,corrected,diagnostic]=p1_prepare_candidate(calibrationD,offsets,public,policy)
% Single calibration implementation shared by the API and batched validation.
[est,corrected,diagnostic]=p1_joint_calibrate(calibrationD,offsets,public);
if isfield(policy,'calibration_mode') && strcmp(policy.calibration_mode,'layered')
 [est,layered]=p1_layered_calibrate(calibrationD,offsets,corrected,est);diagnostic.layered=layered;
end
end
