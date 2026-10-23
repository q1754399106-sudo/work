function est=p1_candidate_domain_estimate(fixed,base,corrected,policy,domain)
% Acquisition domain is public. Never receives a scenario/model label or truth.
assert(any(strcmp(domain,{'Offset','Section'})));
est=base;
if strcmp(domain,'Section') && isfield(policy.selected_profile,'section_calibration_weight')
 est=p1_shrink_calibration(fixed,base,corrected,policy.selected_profile.section_calibration_weight);
 for field={'velocity_correction','kinematic_correction_applied'}
  key=field{1};if isfield(base,key),est.(key)=base.(key);end
 end
end
end
