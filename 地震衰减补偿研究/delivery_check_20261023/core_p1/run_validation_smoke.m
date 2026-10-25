function run_validation_smoke(outdir)
% Validate the full batch harness on development data only.
root=fileparts(mfilename('fullpath'));oldpath=path;cleanup=onCleanup(@()path(oldpath));
addpath(fullfile(root,'matlab_src'));assert(~isfolder(outdir));mkdir(outdir);
spec=jsondecode(fileread(fullfile(root,'inputs','p1_experiment_spec.json')));
s=struct('name','Main','kind','normal','scenarios',{{'SmoothAVO','WeakSection'}}, ...
 'snr_db',18,'vp_multiplier',1,'h_multiplier',1,'phase_deg',0,'max_nmo_shift_s',0,'ar1_rho',0);
suites=repmat(s,1,3);suites(2).name='StressVelocityError';suites(2).kind='velocity';suites(2).vp_multiplier=1.03;suites(2).scenarios={'SmoothAVO'};
suites(3).name='StressVelocityMinus3';suites(3).kind='velocity';suites(3).vp_multiplier=.97;suites(3).scenarios={'SmoothAVO'};
job=struct('model',spec.models.development(1),'group','development', ...
 'calibration_seeds',spec.calibration_seeds(1),'test_seeds',spec.test_seeds(1),'suites',suites,'reuse_original',true);
file=fullfile(outdir,'job.json');uq_writejson(file,job);
run_validation_job(file,fullfile(outdir,'job'));
checks=jsondecode(fileread(fullfile(outdir,'job','checks.json')));
assert(checks.completed && checks.observations==4 && checks.rows==28 && checks.replay_rows==12);
assert(checks.solver_flags==0 && checks.max_normal_residual<=2e-7);
uq_writejson(fullfile(outdir,'checks.json'),struct('passed',true,'development_only',true,'observations',4,'new_rows',28,'original_replay_rows',12));
end
