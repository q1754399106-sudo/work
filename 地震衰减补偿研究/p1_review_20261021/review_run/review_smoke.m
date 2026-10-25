function review_smoke(outdir)
% Independent review harness; original package sources are unchanged.
here=fileparts(mfilename('fullpath')); root=fileparts(here);
addpath(fullfile(root,'core_p1')); addpath(fullfile(root,'core_p1','matlab_src'));
if nargin<1, outdir=fullfile(here,'matlab_smoke_corrected'); end
assert(~isfolder(outdir)); mkdir(outdir);
oldrng=rng; cleanup=onCleanup(@()rng(oldrng));
start=tic;
fprintf('START independent MATLAB smoke %s\n',version);
validate_p1_unknown_q(fullfile(outdir,'unit_validation'));
[configs,~]=p1_models('historical'); c=configs{1};
p=uq_public_config(c,c.nominal_frequency);
trueModels=prepare_models_v2(c.offset,c);
[~,~,~,calScene]=uq_case('SmoothAVO',c,trueModels);
[calD,~]=uq_observation(calScene,c,trueModels,1:11,20261030,18);
est=p1_calibrate_ensemble(calD,c.offset,p);
uq_writejson(fullfile(outdir,'calibration.json'),est);
cm=uq_inverse_config(p,est); invModels=prepare_models_v2(c.offset,cm);
R=avo_projector_v2(invModels);
scenarios={'SmoothAVO','FaultSection','PinchoutSection','WeakSection','CoherentNoiseAVO'};
methods={'Baseline','Frozen','StrongerDamping','NoiseUncertaintyOnly','QModelErrorOnly','CombinedUncertainty','CombinedUncertainty_NoStructure'};
rows=struct([]); count=0;
for isc=1:numel(scenarios)
 [models,domain,map,scene]=uq_case(scenarios{isc},c,invModels);
 [D,trueSigma]=uq_observation(scene,c,trueModels,map,20261101,18);
 for ib=1:numel(methods)
  [out,prior,preset]=p1_solve_branch(D,models,map,cm,domain,R,methods{ib});
  for ir=1:2
   route='Direct'; S=out.direct;
   if ir==2, route='Split'; S=out.split; end
   row=uq_metric_row(S,scene,c,c.model_name,'Main',20261101,methods{ib},route,out,prior,preset,est);
   row.CalibrationSeed=20261030; row.EstimatedSigmaRatio=mean(prior.sigma./trueSigma);
   count=count+1; if count==1, rows=row; else, rows(count)=row; end
  end
  fprintf('%s %s err=%.9f normal=%.3e/%.3e elapsed=%.1fs\n',scenarios{isc},methods{ib},row.RelativeError,out.direct_record.normal_residual,out.split_record.normal_residual,toc(start));
  writetable(struct2table(rows),fullfile(outdir,'metrics.csv'));
 end
end
uq_writejson(fullfile(outdir,'status.json'),struct('status','COMPLETED','version',version,'rows',count,'elapsed_seconds',toc(start),'scope','One historical model, one calibration, one observation seed, five scenarios, seven methods, two routes; smoke only.'));
fprintf('COMPLETED %d rows %.1fs\n',count,toc(start));
end
