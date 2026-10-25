function run_research(action,outdir)
% Existing shared session only. Run as function; never touch base workspace.
root=fileparts(mfilename('fullpath')); pkg=fileparts(root);
assert(~isfolder(outdir),'Output must be new');
oldpath=path; oldrng=rng; cleanup=onCleanup(@()restore(oldpath,oldrng));
addpath(fullfile(root,'matlab_src')); 
switch action
 case 'validate', validate_unknown_q(outdir);
 case {'regression','replay'}, run_retained_grid(action,outdir,pkg);
 otherwise, error('Unknown action');
end
end
function restore(p,r), path(p); rng(r); end
function run_retained_grid(action,outdir,pkg)
mkdir(outdir); start=tic; [configs,~]=stability_models(); o=uq_spec();
rows=struct([]); count=0; inputs=0; flags=0; maxnormal=0; maxdiff=0;
old=readtable(fullfile(pkg,'evidence','historical','metrics.csv'));
for im=2:numel(configs)
 c=configs{im}; if strcmp(action,'regression') && im~=2, continue; end
 public=uq_public_config(c,o.nominal_frequency_hz(im-1));
 est=jsondecode(fileread(fullfile(pkg,'evidence','calibration',['estimate_' c.model_name '.json'])));
 cm=uq_inverse_config(public,est); invmodels=prepare_models_v2(c.offset,cm); trueModels=prepare_models_v2(c.offset,c);
 for isu=1:2
  if strcmp(action,'regression') && isu~=1, continue; end
  suite='Main'; snr=18; scenarios=o.test_scenarios;
  if isu==2, suite='LowerSNR'; snr=12; scenarios=o.stress_scenarios; end
  for ic=1:numel(scenarios)
   name=char(scenarios{ic}); if strcmp(action,'regression') && ~strcmp(name,'SmoothAVO'), continue; end
   [models,domain,map,scene]=uq_case(name,c,invmodels);
   for seed=reshape(o.test_seeds,1,[])
    if strcmp(action,'regression') && seed~=o.test_seeds(1), continue; end
    [D,trueSigma]=uq_observation(scene,c,trueModels,map,seed,snr);
    for ip=1:2
     preset=o.presets(ip); assert(any(strcmp(preset.name,{'Baseline','Frozen'})));
     prior=uq_build_prior(D,models,cm,domain,preset); out=uq_select(D,models,prior,preset);
     au=uq_audit_outputs(out); flags=flags+au.flags; maxnormal=max(maxnormal,au.max_normal);
     for ir=1:2
      route='Direct'; S=out.direct; if ir==2, route='Split'; S=out.split; end
      row=uq_metric_row(S,scene,c,c.model_name,suite,seed,preset.name,route,out,prior,preset,est);
      row.EstimatedSigmaRatio=mean(prior.sigma./trueSigma); row.CalibrationSeed=o.calibration_seed;
      count=count+1; if count==1, rows=row; else, rows(count)=row; end
      mask=strcmp(old.Model,c.model_name)&strcmp(old.Suite,suite)&strcmp(old.Scenario,name)&old.Seed==seed&strcmp(old.Method,preset.name)&strcmp(old.Route,route);
      assert(sum(mask)==1); delta=abs(row.RelativeError-old.RelativeError(mask)); maxdiff=max(maxdiff,delta);
     end
    end
    inputs=inputs+1; writetable(struct2table(rows),fullfile(outdir,'metrics_partial.csv'));
    uq_writejson(fullfile(outdir,'run_status.json'),struct('status','RUNNING','action',action,'inputs',inputs,'rows',count,'elapsed_seconds',toc(start)));
    fprintf('%s %s %s seed %d; %d rows\n',c.model_name,suite,name,seed,count);
   end
  end
 end
end
writetable(struct2table(rows),fullfile(outdir,'metrics.csv'));
checks=struct('solver_flags',flags,'max_normal_residual',maxnormal,'max_historical_waveform_metric_delta',maxdiff);
uq_writejson(fullfile(outdir,'checks.json'),checks);
assert(flags==0 && maxnormal<2e-7 && maxdiff<1e-7,'Retained regression failed');
expected=504; if strcmp(action,'regression'), expected=4; end; assert(count==expected);
uq_writejson(fullfile(outdir,'run_status.json'),struct('status','COMPLETED','action',action,'inputs',inputs,'rows',count,'elapsed_seconds',toc(start),'scope','Historical fixed-Q and fixed-seed replay; not new independent validation'));
end
