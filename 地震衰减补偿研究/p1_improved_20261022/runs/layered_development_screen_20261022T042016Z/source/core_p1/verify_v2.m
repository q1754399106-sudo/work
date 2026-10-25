function verify_v2(outdir)
% Production API, truth-field isolation, layered forward basis and both routes.
root=fileparts(mfilename('fullpath'));oldpath=path;oldrng=rng;cleanup=onCleanup(@()restore(oldpath,oldrng));
addpath(fullfile(root,'matlab_src'));assert(~isfolder(outdir));mkdir(outdir);
assert(isfile(fullfile(root,'inputs','candidate_policy.json')),'Select candidate before production-entry verification');
[configs,~]=p1_models('development');c=configs{1};p=uq_public_config(c,c.nominal_frequency);
poison=c;poison.Q=999*c.Q;poison.fpeak=999;poison.rho=999*c.rho;poison.vs=999*c.vs;
pp=uq_public_config(poison,c.nominal_frequency);assert(isequal(p,pp));
trueModels=cell(1,11);
for j=1:11,[aa,bb]=ray_absorption(c.t,c.offset(j),c);trueModels{j}=struct('G',attenuation_matrix(aa+bb,c));end
calScene=make_stability_scene('SmoothAVO',c);[calD,~]=uq_observation(calScene,c,trueModels,1:11,74101,18);
[D,~]=uq_observation(calScene,c,trueModels,1:11,75101,18);
a=combined_uncertainty_v2(D,c.offset,p,calD,c.offset,'Offset');
b=combined_uncertainty_v2(D,c.offset,pp,calD,c.offset,'Offset');
assert(isequal(a.direct,b.direct)&&isequal(a.split,b.split),'Production output changed after private-field poisoning');
bad=p;bad.Q=c.Q;rejected=false;
try,combined_uncertainty_v2(D,c.offset,bad,calD,c.offset,'Offset');catch ex,rejected=contains(ex.message,'Nonpublic input');end
assert(rejected,'Production API accepted hidden truth fields');
allnormal=[a.direct_record.normal_residual,a.split_record.normal_residual];
allflags=[a.direct_record.flag,a.split_record.flag];
scene=make_stability_scene('WeakSection',c);map=6*ones(1,11);[sectionD,~]=uq_observation(scene,c,trueModels,map,75101,18);
sec=combined_uncertainty_v2(sectionD,c.offset(6)*ones(1,11),p,calD,c.offset,'Section');
allnormal=[allnormal,sec.direct_record.normal_residual,sec.split_record.normal_residual];
allflags=[allflags,sec.direct_record.flag,sec.split_record.flag];
assert(all(isfinite(allnormal))&&all(allnormal<=2e-7)&&all(allflags==0));
quality=uq_audit_outputs(a);sectionquality=uq_audit_outputs(sec);
assert(quality.flags==0&&sectionquality.flags==0&&max(quality.max_normal,sectionquality.max_normal)<=2e-7);
assert(a.compensation_applied&&sec.compensation_applied,'Rejection/pass-through cannot count as restoration');
checks=struct('passed',true,'public_fields_equal_after_poison',isequal(p,pp), ...
 'production_output_exact_poison_invariance',isequal(a.direct,b.direct)&&isequal(a.split,b.split), ...
 'illegal_truth_field_rejected',rejected,'offset_and_section_direct_split_normal_residuals',allnormal, ...
 'all_scale_max_residual',max(quality.max_normal,sectionquality.max_normal), ...
 'compensation_applied',true,'estimated_velocity_correction',a.estimated_velocity_scale, ...
 'kinematic_correction_applied',a.estimated_Q.kinematic_correction_applied);
uq_writejson(fullfile(outdir,'checks.json'),checks);
uq_writejson(fullfile(outdir,'calibration_diagnostics.json'),a.calibration_diagnostics);
fprintf('PASS production-entry truth isolation, finite outputs and all-route/all-scale normal equations\n');
end
function restore(p,r),path(p);rng(r);end
