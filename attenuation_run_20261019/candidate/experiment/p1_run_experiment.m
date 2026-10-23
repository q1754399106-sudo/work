function p1_run_experiment(stage,outdir,specfile)
% Fresh, paired synthetic experiment. Does not connect/start/stop MATLAB.
% Call ONLY inside the user's already shared ResearchMATLAB session.
% Truth stays in p1_synthesize and post-hoc scoring; public input is whitelisted.
assert(any(strcmp(stage,{'development','test'})),'Stage must be development or test');
assert(~isfolder(outdir),'Output directory must be new');
oldpath=path; oldrng=rng; state=onCleanup(@()restore_state(oldpath,oldrng)); %#ok<NASGU>
cfg=jsondecode(fileread(specfile));
here=fileparts(mfilename('fullpath')); cfg=resolve_paths(cfg,here);
if strcmp(stage,'test'), verify_freeze(cfg,specfile); end
addpath(cfg.paths.source_matlab,'-begin'); addpath(cfg.paths.source_core,'-begin');
addpath(cfg.paths.candidate_matlab,'-begin'); addpath(here,'-begin');
assert_same_source('config_model',cfg.paths.source_matlab);
assert_same_source('prepare_models_v2',cfg.paths.source_matlab);
assert_same_source('uq_select',cfg.paths.source_matlab);
assert_same_source('p1_compensate',cfg.paths.candidate_matlab);
assert_same_source('p1_prepare_calibration',cfg.paths.candidate_matlab);
mkdir(outdir); mkdir(fullfile(outdir,'calibration')); mkdir(fullfile(outdir,'solver_logs'));
mkdir(fullfile(outdir,'examples')); mkdir(fullfile(outdir,'figures'));
mkdir(fullfile(outdir,'worst_cases'));
started=tic; rows=struct([]); nrows=0; fail=0; observationCount=0;
worstScores=struct;
manifest=source_manifest(cfg,here);
metadata=make_metadata(cfg,stage,specfile,manifest);
writejson(fullfile(outdir,'RUN_METADATA.json'),metadata);
writejson(fullfile(outdir,'experiment_spec_used.json'),cfg);
allModels=cfg.models; selected=cfg.development_model_ids;
if strcmp(stage,'test'), selected=cfg.test_model_ids; end
methods=[{'Baseline','Frozen'},cellstr(string(cfg.candidate_branches(:)))',cellstr(string(cfg.controls(:)))'];
for im=1:numel(selected)
 id=char(selected{im}); mi=find(strcmp({allModels.id},id)); assert(numel(mi)==1);
 m=allModels(mi); c=p1_model_config(m);
 % Assertion documents that independently chosen layered configurations differ
 % from historical configs. It does not prove external/new geological validity.
 assert_new_model(c,cfg.paths.source_core);
 trueModels=p1_truth_models(c); % SYNTHESIS ONLY; never passed to solver
 suites=cfg.suites;
 if strcmp(stage,'development')
  suites=suites(ismember({suites.name},cfg.development.suite_names));
  for isu=1:numel(suites)
   suites(isu).calibration_count=cfg.development.calibration_count;
   suites(isu).observation_count=cfg.development.observation_count;
   suites(isu).scenario_names=cfg.development.scenario_names;
  end
 end
 for isu=1:numel(suites)
  suite=suites(isu);
  for ik=1:suite.calibration_count
   calSeed=m.calibration_seeds(ik);
   observationSeeds=reshape(m.observation_seeds(ik,1:suite.observation_count),1,[]);
   [public,calD,calOffsets,synthesis]=p1_synthesize('calibration',c,trueModels,suite,calSeed,[]);
   calOpts=cfg.algorithm; calOpts.bootstrap_seed=mod(calSeed+cfg.seed_rules.bootstrap_offset,2^31-1);
   cache=[]; cacheFailure='';
   try
    cache=p1_prepare_calibration(public,calD,calOffsets,c.offset,calOpts);
    assert(isfield(cache,'estimate') && isfield(cache,'models'));
    calibration_record=struct('Model',id,'Suite',suite.name,'CalibrationSeed',calSeed, ...
     'bootstrap_seed',calOpts.bootstrap_seed,'estimate',cache.estimate, ...
     'truth_was_used_for_estimation',false,'public_fields',{fieldnames(public)});
    for extra={'q_frequency_draws','bootstrap_estimates','noise_summary'}
     if isfield(cache,extra{1}), calibration_record.(extra{1})=cache.(extra{1}); end
    end
   catch exception
    cacheFailure=exception_text(exception);
    calibration_record=struct('Model',id,'Suite',suite.name,'CalibrationSeed',calSeed, ...
      'Status','FAILED','FailureMessage',cacheFailure,'truth_was_used_for_estimation',false);
   end
   calibrationName=sprintf('%s_%s_cal%d',id,suite.name,calSeed);
   writejson(fullfile(outdir,'calibration',[calibrationName '.json']),calibration_record);
   calHash=p1_sha256(uint8(getByteStreamFromArray(calD)));
   for ic=1:numel(suite.scenario_names)
    name=char(suite.scenario_names{ic});
    for iseed=1:numel(observationSeeds)
     seed=observationSeeds(iseed);
     [D,scene,trueSigma]=p1_synthesize('observation',c,trueModels,suite,seed,name);
     domain='Offset'; if endsWith(name,'Section'), domain='Section'; end
     inputHash=p1_sha256(uint8(getByteStreamFromArray({D,scene.offsets,public,calD,calOffsets})));
     observationFirstRow=nrows+1; currentSolutions=struct;
     preselected=im==1 && ik==1 && iseed==1 && any(strcmp(suite.name,cfg.example_suites));
     example=struct;
     if preselected
      example=struct('Model',id,'Suite',suite.name,'Scenario',name,'Seed',seed, ...
       'CalibrationSeed',calSeed,'observation',D,'truth',scene.truth,'time',c.t, ...
       'offsets',scene.offsets,'fixed_trace_indices',cfg.fixed_trace_indices,'methods',struct);
     end
     for ip=1:numel(methods)
      method=methods{ip}; out=[]; prior=[]; est=[]; preset=[]; errorText='';
      try
       assert(isempty(cacheFailure),'Calibration preparation failed: %s',cacheFailure);
       if any(strcmp(method,{'Baseline','Frozen'}))
        % Reuse the original point estimate but retain exact historical selection.
        opts0=uq_spec(); ix=find(strcmp({opts0.presets.name},method)); preset=opts0.presets(ix);
        inverse=cache.inverse_config; models=select_cached_models(cache,c.offset,scene.offsets);
        est=cache.estimate; prior=uq_build_prior(D,models,inverse,domain,preset);
        out=uq_select(D,models,prior,preset); out.prior=prior; out.preset=preset; out.estimated_Q=est;
       else
        opts=calOpts; opts.calibration_cache=cache;
        out=p1_compensate(D,scene.offsets,public,calD,calOffsets,domain,method,opts);
        prior=out.prior; est=out.estimated_Q; preset=out.preset;
       end
       validate_output(out,size(D));
       % Hash data available for selection, and exact shared prior across routes.
       informationHash=p1_sha256(uint8(getByteStreamFromArray({inputHash,est,calOpts})));
       priorHash=p1_sha256(uint8(getByteStreamFromArray(prior)));
      catch exception
       errorText=exception_text(exception); informationHash=inputHash; priorHash='FAILED';
      end
      for ir=1:2
       route='Direct'; if ir==2, route='Split'; end
       S=[];
       row=empty_row(id,suite,name,seed,calSeed,method,route,calHash,inputHash);
       row.PriorHash=priorHash; row.SharedInformationHash=informationHash;
       if isempty(errorText)
        try
         S=out.direct; if ir==2, S=out.split; end
         scored=uq_metric_row(S,scene,c,id,suite.name,seed,method,route,out,prior,preset,est);
         keys=fieldnames(scored); for jf=1:numel(keys), row.(keys{jf})=scored.(keys{jf}); end
         row.Status='SUCCESS'; row.FailureMessage='';
         row.EstimatedSigmaRatio=mean(prior.sigma./max(trueSigma,realmin)); % SCORING ONLY
         row.QBoundary=logical(est.q_boundary); row.FrequencyBoundary=logical(est.frequency_boundary);
         row.CalibrationEstimateHash=p1_sha256(uint8(getByteStreamFromArray(est)));
         if row.FinalFlag~=0 || row.NormalResidual>cfg.acceptance.normal_residual_max
          row.Status='SOLVER_FAILED'; row.FailureMessage='Solver flag or normal residual gate failed';
         end
         if preselected, example.methods.([method '_' route])=S; end
        catch exception
         row.Status='FAILED'; row.FailureMessage=exception_text(exception);
        end
       else
        row.Status='FAILED'; row.FailureMessage=errorText;
       end
       if ~isempty(S) && all(isfinite(S(:)))
        currentSolutions.([method '_' route])=S;
       end
       if ~strcmp(row.Status,'SUCCESS'), fail=fail+1; end
       nrows=nrows+1; if nrows==1, rows=row; else, rows(nrows)=row; end
      end
      log=struct('Model',id,'Suite',suite.name,'Scenario',name,'Seed',seed, ...
        'CalibrationSeed',calSeed,'Method',method,'FailureMessage',errorText);
      if ~isempty(out)
       wanted={'records','direct_record','split_record','scale','scale_grid','residual_grid', ...
        'target','endpoint','residual_mismatch','diagnostics'};
       for jw=1:numel(wanted), if isfield(out,wanted{jw}), log.(wanted{jw})=out.(wanted{jw}); end, end
      end
      logName=sprintf('%s_%s_%s_cal%d_seed%d_%s.json',id,suite.name,name,calSeed,seed,method);
      writejson(fullfile(outdir,'solver_logs',logName),log);
      partialFile=fullfile(outdir,'metrics_partial.csv');
      methodRows=struct2table(rows(max(1,nrows-1):nrows));
      if nrows==2
       writetable(methodRows,partialFile);
      else
       writetable(methodRows,partialFile,'WriteMode','append','WriteVariableNames',false);
      end
     end
     % Fixed post-hoc rule: retain the largest finite RelativeError per
     % method/route. This scoring-only capture never feeds parameter selection.
     % Deferring capture until all methods finish includes matched Frozen data.
     for iw=observationFirstRow:nrows
      worstRow=rows(iw); caseKey=[worstRow.Method '_' worstRow.Route];
      if ~isfinite(worstRow.RelativeError) || ~isfield(currentSolutions,caseKey), continue; end
      if isfield(worstScores,caseKey) && worstRow.RelativeError<=worstScores.(caseKey), continue; end
      worstScores.(caseKey)=worstRow.RelativeError;
      captured=struct; captured.(caseKey)=currentSolutions.(caseKey);
      for frozenKey={'Frozen_Direct','Frozen_Split'}
       if isfield(currentSolutions,frozenKey{1}), captured.(frozenKey{1})=currentSolutions.(frozenKey{1}); end
      end
      worstCase=struct('Model',id,'Suite',suite.name,'Scenario',name,'Seed',seed, ...
       'CalibrationSeed',calSeed,'observation',D,'truth',scene.truth,'time',c.t, ...
       'offsets',scene.offsets,'fixed_trace_indices',cfg.fixed_trace_indices,'methods',captured, ...
       'Method',worstRow.Method,'Route',worstRow.Route,'Status',worstRow.Status, ...
       'RelativeError',worstRow.RelativeError,'FailureMessage',worstRow.FailureMessage, ...
       'capture_rule','largest finite RelativeError per method/route, including finite solver failures; scoring only');
      save(fullfile(outdir,'worst_cases',[caseKey '.mat']),'worstCase','-v7.3');
      writejson(fullfile(outdir,'worst_cases',[caseKey '.json']),worstCase);
     end
     observationCount=observationCount+1;
     if preselected
      stem=sprintf('%s_%s_%s_cal%d_seed%d',id,suite.name,name,calSeed,seed);
      save(fullfile(outdir,'examples',[stem '.mat']),'example','-v7.3');
      writejson(fullfile(outdir,'examples',[stem '.json']),example);
     end
     writejson(fullfile(outdir,'run_status.json'),struct('status','RUNNING','stage',stage, ...
       'observations',observationCount,'rows',nrows,'failed_rows',fail,'elapsed_seconds',toc(started)));
     fprintf('%s %s %s calibration %d observation %d; %d rows (%d failed)\n', ...
       id,suite.name,name,calSeed,seed,nrows,fail);
    end
   end
  end
 end
end
writetable(struct2table(rows),fullfile(outdir,'metrics.csv'));
expected=0;
for isu=1:numel(suites)
 expected=expected+numel(selected)*numel(suites(isu).scenario_names)* ...
  suites(isu).calibration_count*suites(isu).observation_count*numel(methods)*2;
end
assert(nrows==expected,'Incomplete matrix; partial CSV retained');
metadata.status='COMPLETED'; metadata.rows=nrows; metadata.failed_rows=fail;
metadata.observations=observationCount; metadata.elapsed_seconds=toc(started);
metadata.completed_with_method_failures=fail>0;
writejson(fullfile(outdir,'RUN_METADATA.json'),metadata);
writejson(fullfile(outdir,'run_status.json'),struct('status','COMPLETED','stage',stage, ...
 'rows',nrows,'failed_rows',fail,'elapsed_seconds',toc(started)));
% Engineering acceptance and hierarchical intervals are computed separately
% by strict_evaluate.py, which must refuse a successful claim on failed rows.
end

function cfg=resolve_paths(cfg,here)
for k={'source_core','source_matlab','candidate_matlab'}
 value=cfg.paths.(k{1}); if ~is_absolute(value), value=fullfile(here,value); end
 cfg.paths.(k{1})=char(java.io.File(value).getCanonicalPath());
 assert(isfolder(cfg.paths.(k{1})),'Missing source directory: %s',cfg.paths.(k{1}));
end
end
function value=is_absolute(p), value=~isempty(regexp(p,'^([A-Za-z]:[\\/]|/|\\\\)','once')); end
function restore_state(p,r), path(p); rng(r); end
function assert_same_source(name,root)
actual=which(name); expected=fullfile(root,[name '.m']);
assert(strcmpi(char(java.io.File(actual).getCanonicalPath()),char(java.io.File(expected).getCanonicalPath())), ...
 'Unexpected path resolution for %s: %s',name,actual);
end
function assert_new_model(c,core)
history=jsondecode(fileread(fullfile(core,'inputs','stability_spec.json')));
for k=1:numel(history.models)
 h=history.models(k);
 assert(~(isequal(c.vp(:),h.vp(:)) && isequal(c.Q(:),h.Q(:)) && isequal(c.event_t(:),h.event_t(:))), ...
 'Historical layered model cannot be used as independent new geology');
end
end
function models=select_cached_models(cache,offsets,needed)
models=cell(1,numel(needed));
for j=1:numel(needed), ix=find(offsets==needed(j),1); assert(~isempty(ix)); models{j}=cache.models{ix}; end
end
function validate_output(out,sz)
for f={'direct','split'}
 value=out.(f{1}); assert(isequal(size(value),sz) && all(isfinite(value(:))),'Invalid recovered gather');
end
end
function text=exception_text(e), text=regexprep(getReport(e,'extended','hyperlinks','off'),'[\r\n]+',' | '); end
function row=empty_row(id,suite,scenario,seed,calSeed,method,route,calHash,inputHash)
row=struct('Model',id,'Suite',suite.name,'Scenario',scenario,'SNRdB',suite.snr_db, ...
 'Seed',seed,'CalibrationSeed',calSeed,'Method',method,'Route',route,'Preset',method, ...
 'RelativeError',NaN,'SignalWindowError',NaN,'NoiseWindowRMS',NaN,'EventAmplitudeError',NaN, ...
 'UnidentifiableAmplitudes',NaN,'AVOInterceptError',NaN,'AVOGradientError',NaN, ...
 'WeakEventError',NaN,'FaultNeighborhoodError',NaN,'PinchoutTailError',NaN, ...
 'EffectiveQ',NaN,'EstimatedFrequency',NaN,'Gamma',NaN,'PilotMultiplier',NaN, ...
 'TargetMultiplier',NaN,'Scale',NaN,'ScaleEndpoint',NaN,'TargetMismatch',NaN, ...
 'NormalResidual',NaN,'FinalFlag',1,'Solver','NOT_RUN','HTruthRatio',NaN, ...
 'EstimatedSigmaRatio',NaN,'QBoundary',NaN,'FrequencyBoundary',NaN, ...
 'PriorHash','','SharedInformationHash','','CalibrationDataHash',calHash, ...
 'ObservationInputHash',inputHash,'CalibrationEstimateHash','','Status','FAILED','FailureMessage','');
end
function metadata=make_metadata(cfg,stage,specfile,sources)
metadata=struct('status','RUNNING','stage',stage,'actual_matlab_execution',true, ...
 'matlab_version',version,'matlab_release',version('-release'),'all_failures_retained',true, ...
 'frozen_before_test',strcmp(stage,'test'),'acceptance_thresholds_frozen_before_test',strcmp(stage,'test'), ...
 'frozen_experiment_spec_sha256',p1_sha256(specfile),'source_manifest',sources, ...
 'truth_isolation_audit',struct('status','DESIGN_ENFORCED_PENDING_EXTERNAL_AUDIT', ...
 'allowed_uses',{{'synthesis','post_hoc_scoring','HTruthRatio'}}, ...
 'public_whitelist',{{'dt','nt','t','offset','vp','h','fref','nfft','lambda_x','lambda_fixed','lambda_grid','cov_floor','nominal_frequency'}}), ...
 'new_geological_holdout',struct('status','SYNTHETIC_DESIGN_PENDING_EXTERNAL_AUDIT', ...
 'development_model_ids',{cfg.development_model_ids},'test_model_ids',{cfg.test_model_ids}, ...
 'prior_observed_test_results_used',false,'external_geological_evidence',false), ...
 'new_geology_claim_permitted',false,'paired_routes_share_prior',true, ...
 'all_methods_share_observation_and_point_calibration_estimate',true, ...
 'candidate_information_matches_Direct_and_Split',true, ...
 'number_of_geological_independent_units',numel(cfg.test_model_ids), ...
 'stress_calibration_variance_claim_permitted',false, ...
 'hyperparameter_search_budget',0, ...
 'seed_note','Calibration seeds generate fresh calibration noise; six observation seeds are nested in each calibration. Scenes/SNRs share nested seeds for pairing.', ...
 'source_hashes',struct,'source_hash_paths',struct,'observation_seeds_by_calibration',struct);
for k=1:numel(sources)
 key=sprintf('source_%04d',k); metadata.source_hashes.(key)=sources(k).sha256;
 metadata.source_hash_paths.(key)=sources(k).path;
end
for k=1:numel(cfg.models)
 m=cfg.models(k); entry=struct;
 for j=1:numel(m.calibration_seeds), entry.(sprintf('cal_%d',m.calibration_seeds(j)))=m.observation_seeds(j,:); end
 metadata.observation_seeds_by_calibration.(m.id)=entry;
end
end
function sources=source_manifest(cfg,here)
roots={cfg.paths.source_core,cfg.paths.candidate_matlab,here}; sources=struct('path',{},'sha256',{});
for r=1:numel(roots)
 files=[dir(fullfile(roots{r},'**','*.m'));dir(fullfile(roots{r},'**','*.py'));dir(fullfile(roots{r},'inputs','*.json'))];
 for k=1:numel(files)
  p=char(java.io.File(fullfile(files(k).folder,files(k).name)).getCanonicalPath());
  if ~isempty(sources) && any(strcmp({sources.path},p)), continue; end
  sources(end+1)=struct('path',p,'sha256',p1_sha256(p)); %#ok<AGROW>
 end
end
end
function verify_freeze(cfg,specfile)
assert(isfield(cfg,'freeze') && isfield(cfg.freeze,'expected_hash') && strlength(cfg.freeze.expected_hash)==64, ...
 'Test requires explicit freeze.expected_hash created after completed development');
f=cfg.freeze; snapshot=f.config_snapshot;
currentSources=source_manifest(cfg,fileparts(mfilename('fullpath')));
assert(numel(currentSources)==numel(f.sources) && isempty(setxor({currentSources.path},{f.sources.path})), ...
 'Frozen source inventory changed (including newly added source files)');
if ~is_absolute(snapshot), snapshot=fullfile(fileparts(specfile),snapshot); end
assert(strcmp(p1_sha256(snapshot),f.config_sha256),'Frozen config snapshot changed');
original=jsondecode(fileread(snapshot)); current=rmfield(cfg,'freeze');
% Resolve both sets identically before checking equality, without serializing JSON.
original=resolve_paths(original,fileparts(mfilename('fullpath')));
assert(isequaln(current,original),'Frozen spec differs from original frozen config');
identity=[f.config_sha256 newline];
for k=1:numel(f.sources)
 s=f.sources(k); assert(strcmp(p1_sha256(s.path),s.sha256),'Frozen source changed: %s',s.path);
 identity=[identity s.path '=' s.sha256 newline]; %#ok<AGROW>
end
assert(strcmp(p1_sha256(unicode2native(identity,'UTF-8')),f.expected_hash),'Freeze identity hash mismatch');
assert(strcmp(f.development_status,'COMPLETED') && f.development_actual_matlab_execution, ...
 'A completed actual development run is required before freezing');
end
function writejson(p,value)
fid=fopen(p,'w','n','UTF-8'); assert(fid>=0); cleanup=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'%s',jsonencode(value,'PrettyPrint',true));
end
