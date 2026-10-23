function run_refit_api_check(outdir)
root=fileparts(mfilename('fullpath'));addpath(fullfile(root,'matlab_src'));
assert(~isfolder(outdir));mkdir(outdir);
plan=jsondecode(fileread(fullfile(root,'inputs','adaptive_cv_plan.json')));cache=plan.cache;
z=load(fullfile(cache,'cases','OriginalReference_WeakSection_12dB_20261411.mat'));
cal=load(fullfile(cache,'calibration_input_OriginalReference.mat'));
out=reflectivity_compensate(z.D,z.offsets,cal.public,cal.calD,cal.public.offset,'Section');
development=jsondecode(fileread(fullfile(root,'inputs','frozen_refit_profile.json')));
reference=load(fullfile(development.development_run,'matlab','cases','OriginalReference_WeakSection_12dB_20261411.mat'));
replay=norm(out.gather-reference.full,'fro')/norm(reference.full,'fro');assert(replay<1e-8&&out.optimization_flags==0);
bad=cal.public;bad.Q=ones(1,4);rejected=false;
try,reflectivity_compensate(z.D,z.offsets,bad,cal.calD,cal.public.offset,'Section');catch err,rejected=contains(err.message,'Nonpublic');end
assert(rejected);
uq_writejson(fullfile(outdir,'checks.json'),struct('complete',true,'api_replay_relative_error',replay,'truth_input_rejected',rejected,'optimization_flags',out.optimization_flags,'max_normal',out.normal_residual));
end
