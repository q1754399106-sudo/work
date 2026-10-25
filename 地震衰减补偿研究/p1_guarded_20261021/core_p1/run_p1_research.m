function run_p1_research(action,outdir)
% Shared-session P1 research entry. Never touches base workspace, path, or RNG permanently.
root=fileparts(mfilename('fullpath')); pkg=fileparts(root);
assert(~isfolder(outdir),'Output must be a new directory');
oldpath=path; oldrng=rng; cleanup=onCleanup(@()restore(oldpath,oldrng));
addpath(fullfile(root,'matlab_src'));
switch action
 case 'validate'
  mkdir(outdir); validate_p1_unknown_q(fullfile(outdir,'original_units')); p1_regression_checks(outdir); p1_guard_checks(outdir);
 case {'regression','replay'}, run_research(action,outdir);
 case {'smoke','guard_trial','p1_experiment','full'}
  spec=jsondecode(fileread(fullfile(root,'inputs','p1_experiment_spec.json')));
  small=any(strcmp(action,{'smoke','guard_trial'}));
  if small
   spec.calibration_seeds=spec.calibration_seeds(1); spec.test_seeds=spec.test_seeds(1);
   spec.models.historical=spec.models.historical(1);
   spec.models.development=[]; spec.models.heldout_new_geology=[];
   spec.stress_suites=spec.stress_suites(strcmp({spec.stress_suites.name},'StressVelocityError'));
  end
  if strcmp(action,'guard_trial')
   spec.test_seeds=[20261101;20261107];
   chosen={'Frozen','NoiseUncertaintyOnly','CombinedUncertainty','GuardedNoise','NoCompensation'};
   spec.methods=spec.methods(ismember({spec.methods.name},chosen));
   positive=spec.stress_suites; negative=positive;
   positive.name='StressVelocityPlus3'; positive.vp_multiplier=1.03;
   negative.name='StressVelocityMinus3'; negative.vp_multiplier=.97;
   negative.description='-3% velocity with fixed thickness; recalibrate and rebuild operators';
   spec.stress_suites=[positive;negative];
  end
  p1_run_grid(outdir,spec,small);
 otherwise, error('Unknown action: %s',action);
end
end
function restore(p,r), path(p); rng(r); end
