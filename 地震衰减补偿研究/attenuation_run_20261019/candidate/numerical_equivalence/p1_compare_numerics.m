function report=p1_compare_numerics(outdir)
% Development-only numerical equivalence and computational-thread benchmark.
% No test geology is synthesized; no reference truth enters either solver.
% Call only through the existing shared ResearchMATLAB Engine connection.
assert(~isfolder(outdir),'Use a new numerical-equivalence output directory');
originalPath=path; originalRng=rng; originalThreads=maxNumCompThreads;
cleanup=onCleanup(@()restore_state(originalPath,originalRng,originalThreads)); %#ok<NASGU>
here=fileparts(mfilename('fullpath')); candidate=fileparts(here); runroot=fileparts(candidate);
active=fullfile(candidate,'matlab'); legacy=fullfile(candidate,'legacy_numerics_v1');
experiment=fullfile(candidate,'experiment'); core=fullfile(runroot,'baseline_snapshot','core');
assert(isfolder(active) && isfolder(legacy));
addpath(core,'-begin'); addpath(fullfile(core,'matlab_src'),'-begin');
addpath(experiment,'-begin'); addpath(active,'-begin');
assert_source('p1_compensate',active); mkdir(outdir);
report=struct('status','RUNNING','development_only',true,'heldout_data_used',false, ...
 'original_threads',originalThreads,'matlab_version',version, ...
 'thread_control_documentation','https://www.mathworks.com/help/matlab/ref/maxnumcompthreads.html', ...
 'scope','Numerical implementation and hardware scheduling equivalence; not scientific acceptance');
report.small_matrix=small_matrix_check;
writejson(fullfile(outdir,'numerical_equivalence.json'),report);
% All model-specific truth reads below belong exclusively to synthesis.
cfg=jsondecode(fileread(fullfile(experiment,'experiment_spec.json')));
modelId='DevGradedCompact'; m=cfg.models(strcmp({cfg.models.id},modelId));
assert(numel(m)==1 && any(strcmp(cfg.development_model_ids,modelId)));
suite=cfg.suites(strcmp({cfg.suites.name},'Main')); assert(numel(suite)==1);
c=p1_model_config(m); trueModels=p1_truth_models(c);
calSeed=310101; seed=410101; scenario='SmoothAVO'; domain='Offset';
setup=tic; maxNumCompThreads(1);
[public,calD,calOffsets]=p1_synthesize('calibration',c,trueModels,suite,calSeed,[]);
[D,scene]=p1_synthesize('observation',c,trueModels,suite,seed,scenario);
calOpts=cfg.algorithm; calOpts.bootstrap_seed=mod(calSeed+cfg.seed_rules.bootstrap_offset,2^31-1);
cache=p1_prepare_calibration(public,calD,calOffsets,c.offset,calOpts);
opts=calOpts; opts.calibration_cache=cache; clear trueModels
report.case=struct('Model',modelId,'Suite','Main','Scenario',scenario,'Seed',seed, ...
 'CalibrationSeed',calSeed,'calibration_setup_seconds',toc(setup),'calibration_setup_threads',1, ...
 'input_hash',p1_sha256(uint8(getByteStreamFromArray({D,scene.offsets,public,calD,calOffsets}))), ...
 'point_estimate_hash',p1_sha256(uint8(getByteStreamFromArray(cache.estimate))), ...
 'bootstrap_draws_hash',p1_sha256(uint8(getByteStreamFromArray(cache.q_frequency_draws))), ...
 'truth_used_for_estimation',false);
maxNumCompThreads(originalThreads);
counts=unique([originalThreads 1 2 4],'stable'); counts=counts(counts<=originalThreads);
benchmark=struct([]); reference=[]; bestThreads=originalThreads; bestSeconds=Inf;
for k=1:numel(counts)
 count=counts(k); maxNumCompThreads(count); clock=tic;
 output=p1_compensate(D,scene.offsets,public,calD,calOffsets,domain,'CombinedUncertainty',opts);
 elapsed=toc(clock); if isempty(reference), reference=output; end
 check=compare_outputs(reference,output);
 entry=struct('threads',count,'elapsed_seconds',elapsed,'comparison_to_original_threads',check, ...
  'direct_normal_residual',output.direct_record.normal_residual, ...
  'split_normal_residual',output.split_record.normal_residual,'scale',output.scale);
 if isempty(benchmark), benchmark=entry; else, benchmark(end+1)=entry; end %#ok<AGROW>
 if check.pass && elapsed<bestSeconds, bestSeconds=elapsed; bestThreads=count; end
 report.thread_benchmark=benchmark; report.recommended_threads=bestThreads;
 writejson(fullfile(outdir,'numerical_equivalence.json'),report);
 fprintf('THREAD BENCHMARK threads=%d time=%.3fs relative waveform delta=(%.3e,%.3e) scale equal=%d PASS=%d\n', ...
  count,elapsed,check.direct_relative_difference,check.split_relative_difference,check.same_selected_scale,check.pass);
end
assert(all(arrayfun(@(x)x.comparison_to_original_threads.pass,benchmark)), ...
 'Computational-thread change failed numerical-equivalence gate');
% Legacy and optimized source are compared with the same cache and hardware.
maxNumCompThreads(bestThreads); comparisons=struct([]);
methods=[cellstr(string(cfg.candidate_branches(:)))' cellstr(string(cfg.controls(:)))'];
historical=readtable(fullfile(runroot,'experiments','development_v1','matlab','metrics.csv'),'TextType','string');
for k=1:numel(methods)
 method=methods{k}; addpath(legacy,'-begin'); assert_source('p1_compensate',legacy);
 clock=tic; old=p1_compensate(D,scene.offsets,public,calD,calOffsets,domain,method,opts); oldTime=toc(clock);
 addpath(active,'-begin'); assert_source('p1_compensate',active);
 clock=tic; fresh=p1_compensate(D,scene.offsets,public,calD,calOffsets,domain,method,opts); newTime=toc(clock);
 check=compare_outputs(old,fresh);
 metricComparison=compare_historical(fresh,old,historical,scene,c,modelId,seed,calSeed,method);
 entry=struct('method',method,'threads',bestThreads,'legacy_seconds',oldTime,'optimized_seconds',newTime, ...
  'comparison',check,'historical_metrics',metricComparison,'pass',check.pass && metricComparison.pass);
 if isempty(comparisons), comparisons=entry; else, comparisons(end+1)=entry; end %#ok<AGROW>
 report.branches=comparisons; writejson(fullfile(outdir,'numerical_equivalence.json'),report);
 fprintf('NUMERICAL EQUIVALENCE %s legacy %.3fs optimized %.3fs waveform=(%.3e,%.3e) PASS=%d\n', ...
  method,oldTime,newTime,check.direct_relative_difference,check.split_relative_difference,entry.pass);
end
assert(all([comparisons.pass]),'Legacy/optimized numerical equivalence failed');
maxNumCompThreads(originalThreads); path(originalPath); rng(originalRng);
report.session_preserved=struct('path',strcmp(path,originalPath),'rng',isequal(rng,originalRng), ...
 'computational_threads',maxNumCompThreads==originalThreads);
assert(all(structfun(@(x)x,report.session_preserved)));
report.status='PASS'; report.recommended_threads=bestThreads;
report.recommendation_basis='Fastest passing CombinedUncertainty timing on one development observation; no heldout outcomes used';
writejson(fullfile(outdir,'numerical_equivalence.json'),report);
fprintf('NUMERICAL EQUIVALENCE PASS; best measured thread count %d; original %d restored\n',bestThreads,originalThreads);
end

function check=compare_outputs(a,b)
de=norm(a.direct-b.direct,'fro')/max(norm(a.direct,'fro'),realmin);
se=norm(a.split-b.split,'fro')/max(norm(a.split,'fro'),realmin);
sameScale=abs(a.scale-b.scale)<=1e-12*max(abs(a.scale),realmin);
normal=max([a.direct_record.normal_residual a.split_record.normal_residual ...
 b.direct_record.normal_residual b.split_record.normal_residual]);
flags=[a.direct_record.flag a.split_record.flag b.direct_record.flag b.split_record.flag];
check=struct('direct_relative_difference',de,'split_relative_difference',se, ...
 'same_selected_scale',sameScale,'reference_selected_index',a.selected_index, ...
 'candidate_selected_index',b.selected_index,'maximum_normal_residual',normal, ...
 'solver_flags',flags,'waveform_tolerance',1e-5, ...
 'pass',de<1e-5 && se<1e-5 && sameScale && normal<=2e-7 && all(flags==0));
end

function check=compare_historical(fresh,old,T,scene,c,id,seed,calSeed,method)
% Truth is used here only for post-hoc scoring after both solves have returned.
fields={'RelativeError','SignalWindowError','NoiseWindowRMS','EventAmplitudeError', ...
 'AVOInterceptError','AVOGradientError','WeakEventError','PinchoutTailError'};
rows=struct([]); maxDelta=0;
for j=1:2
 route='Direct'; S=fresh.direct; if j==2, route='Split'; S=fresh.split; end
 scored=uq_metric_row(S,scene,c,id,'Main',seed,method,route,fresh,fresh.prior,fresh.preset,fresh.estimated_Q);
 ix=T.Model==id & T.Suite=="Main" & T.Scenario=="SmoothAVO" & T.Seed==seed & ...
  T.CalibrationSeed==calSeed & T.Method==method & T.Route==route;
 assert(sum(ix)==1,'Missing or duplicate historical development metric');
 deltas=struct;
 for k=1:numel(fields)
  field=fields{k}; saved=T.(field)(ix); current=scored.(field);
  if isfinite(saved) && isfinite(current)
   delta=abs(current-saved); maxDelta=max(maxDelta,delta); deltas.(field)=delta;
  else
   assert(isequaln(saved,current),'Historical metric finite/nonfinite status changed');
  end
 end
 entry=struct('route',route,'absolute_metric_differences',deltas);
 if isempty(rows), rows=entry; else, rows(end+1)=entry; end %#ok<AGROW>
end
check=struct('maximum_absolute_metric_difference',maxDelta,'tolerance',1e-5, ...
 'routes',rows,'pass',maxDelta<1e-5,'post_hoc_scoring_only',true);
assert(isequaln(old.estimated_Q,fresh.estimated_Q),'Point estimates changed');
assert(isequaln(old.calibration_estimates,fresh.calibration_estimates),'Bootstrap calibration estimates changed');
end

function report=small_matrix_check
% Independent dense precision and normal equations, including rank-zero case.
stream=RandStream('mt19937ar','Seed',17420); n=13; nx=3; N=n*nx;
cases=struct([]);
for rank=[0 4]
 blocks=cell(1,nx); G=blocks; denseC=zeros(N); denseG=zeros(N);
 U=.08*randn(stream,N,rank);
 for j=1:nx
  ids=(j-1)*n+(1:n); B=randn(stream,n,n);
  blocks{j}=B*B'/n+.4*eye(n); G{j}=randn(stream,n,n)/sqrt(n)+eye(n);
  denseC(ids,ids)=blocks{j}; denseG(ids,ids)=G{j};
 end
 denseC=denseC+U*U'; C=p1_covariance(blocks,U);
 H=kron(sparse([1 -1 0;0 1 -1]),speye(n)); Y=randn(stream,n,nx);
 options=p1_options(struct('cov_floor',1e-6),struct('pcg_tolerance',1e-11,'pcg_max_iterations',400));
 lambda=[.03 .05 .07]; mu=.09; scale=.73;
 ctx=p1_context(G,C,Y,H,lambda,mu,options); assert(isfield(ctx,'normal_low_rank'));
 X=randn(stream,N,3); fast=zeros(size(X));
 for j=1:nx, ids=(j-1)*n+(1:n); fast(ids,:)=ctx.block_grams{j}*X(ids,:); end
 Z=ctx.normal_low_rank;
 if rank>0, fast=fast-Z*(C.T'\(C.T\(Z'*X))); end
 fast=fast+scale*ctx.ridge.*X+scale*mu*(ctx.HH*X);
 denseA=denseG'*(denseC\denseG)+scale*diag(ctx.ridge)+scale*mu*(H'*H);
 normalError=norm(fast-denseA*X,'fro')/norm(denseA*X,'fro');
 rhsError=norm(ctx.rhs-denseG'*(denseC\Y(:)))/norm(ctx.rhs);
 [S,rec]=p1_solve(ctx,scale,zeros(n,nx)); exact=denseA\(denseG'*(denseC\Y(:)));
 solutionError=norm(S(:)-exact)/norm(exact);
 entry=struct('covariance_rank',rank,'normal_product_relative_error',normalError, ...
  'rhs_relative_error',rhsError,'rhs_rearrangement_relative_error',ctx.rhs_rearrangement_relative_error, ...
  'solution_relative_error',solutionError,'original_precision_normal_residual',rec.normal_residual, ...
  'pass',normalError<1e-11 && rhsError<1e-11 && solutionError<1e-8 && rec.normal_residual<=2e-7);
 if isempty(cases), cases=entry; else, cases(end+1)=entry; end %#ok<AGROW>
end
assert(all([cases.pass]),'Independent dense normal-equation check failed');
report=struct('status','PASS','cases',cases);
end

function assert_source(name,root)
actual=char(java.io.File(which(name)).getCanonicalPath());
expected=char(java.io.File(fullfile(root,[name '.m'])).getCanonicalPath());
assert(strcmpi(actual,expected),'Unexpected resolution for %s: %s',name,actual);
end
function restore_state(p,r,t), maxNumCompThreads(t); path(p); rng(r); end
function writejson(p,value)
fid=fopen(p,'w','n','UTF-8'); assert(fid>=0); cleanup=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'%s',jsonencode(value,'PrettyPrint',true));
end
