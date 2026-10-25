function run_p1_research(action,outdir)
% Shared-session P1 research entry. Never touches base workspace, path, or RNG permanently.
root=fileparts(mfilename('fullpath')); pkg=fileparts(root);
assert(~isfolder(outdir),'Output must be a new directory');
oldpath=path; oldrng=rng; cleanup=onCleanup(@()restore(oldpath,oldrng));
addpath(fullfile(root,'matlab_src'));
switch action
 case 'validate'
  mkdir(outdir); validate_p1_unknown_q(fullfile(outdir,'original_units')); p1_regression_checks(outdir);
 case {'regression','replay'}, run_research(action,outdir);
 case {'smoke','p1_experiment','full'}
  spec=jsondecode(fileread(fullfile(root,'inputs','p1_experiment_spec.json')));
  if strcmp(action,'smoke')
   spec.calibration_seeds=spec.calibration_seeds(1); spec.test_seeds=spec.test_seeds(1);
   spec.models.historical=spec.models.historical(1);
   spec.models.development=[]; spec.models.heldout_new_geology=[];
   spec.stress_suites=spec.stress_suites(strcmp({spec.stress_suites.name},'StressVelocityError'));
  end
  p1_run_grid(outdir,spec,strcmp(action,'smoke'));
 otherwise, error('Unknown action: %s',action);
end
end
function restore(p,r), path(p); rng(r); end
