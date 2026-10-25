function out=p1_compensate(D,offsets,public,calibrationD,calibrationOffsets,domain,branch,opts)
% Fixed-pilot, full-gather model-error prototype with matched Direct/Split.
% All truth belongs outside this function. Requires retained core on the path.
if nargin<8, opts=struct; end
p1_check_public(public); o=p1_options(public,opts);
branches={'NoiseUncertaintyOnly','QModelErrorOnly','CombinedUncertainty', ...
 'ScalarBudgetDamping','NoiseCovarianceOnly','NoStructureCombined'};
assert(any(strcmp(branch,branches)),'Unsupported P1 branch');
assert(any(strcmp(domain,{'Offset','Section'})));
assert(isreal(D) && size(D,1)==public.nt && size(D,2)==numel(offsets) && all(isfinite(D(:))));
if isfield(opts,'calibration_cache') && ~isempty(opts.calibration_cache)
 cache=opts.calibration_cache;
 assert(strcmp(cache.version,'P1_FIXED_PILOT_V1') && isequaln(cache.public,public));
 assert(isequaln(cache.calibrationD,calibrationD) && isequal(cache.calibrationOffsets,reshape(calibrationOffsets,1,[])));
 keys={'bootstrap_replicates','bootstrap_seed','noise_bootstrap_replicates','block_length','noise_window','noise_quantile'};
 for k=1:numel(keys), assert(isequaln(cache.options.(keys{k}),o.(keys{k})),'Calibration cache options differ'); end
else
 cache=p1_prepare_calibration(public,calibrationD,calibrationOffsets,offsets,o);
end
n=size(D,1); nx=size(D,2); N=n*nx; models=cell(1,nx); map=zeros(1,nx);
for j=1:nx
 k=find(cache.all_offsets==offsets(j),1);
 assert(~isempty(k),'Requested offset is absent from calibration cache'); map(j)=k; models{j}=cache.all_models{k};
end
% This prior is frozen once; neither covariance nor H is updated during solve.
frozenPreset=struct('name','Frozen','gamma',o.gamma,'pilot_multiplier',o.pilot_multiplier,'noise_target_multiplier',1);
prior=uq_build_prior(D,models,cache.inverse_config,domain,frozenPreset);
noise=p1_noise_summary(D,public.t,o,o.bootstrap_seed+o.noise_bootstrap_seed_offset);
R=o.bootstrap_replicates; errors=zeros(N,R);
for j=1:nx
 ids=(j-1)*n+(1:n); ref=models{j}.G*prior.pilot(:,j);
 for k=1:R, errors(ids,k)=cache.operator_draws{map(j),k}*prior.pilot(:,j)-ref; end
end
b=mean(errors,2); centered=errors-b;
ceDiagonal=sum(centered.^2,2)/(R-1);
U=sqrt(1-o.covariance_shrinkage)*centered/sqrt(R-1);
ceShrinkDiagonal=o.covariance_shrinkage*ceDiagonal;
ensembleTraceCe=sum(ceDiagonal); useQ=any(strcmp(branch,{'QModelErrorOnly','CombinedUncertainty','ScalarBudgetDamping','NoStructureCombined'}));
useNoise=any(strcmp(branch,{'NoiseUncertaintyOnly','CombinedUncertainty','ScalarBudgetDamping','NoStructureCombined'}));
variance=noise.variance_reference; if useNoise, variance=noise.variance_upper; end
rho=zeros(1,nx);
if strcmp(o.noise_covariance_model,'ar1') || strcmp(branch,'NoiseCovarianceOnly'), rho=noise.rho; end
if ~useQ, b=zeros(N,1); U=zeros(N,0); ceShrinkDiagonal=zeros(N,1); end
traceCn=n*sum(variance); traceCe=double(useQ)*ensembleTraceCe; scalarAlpha=1;
if strcmp(branch,'ScalarBudgetDamping')
 scalarAlpha=(traceCn+ensembleTraceCe)/(n*sum(noise.variance_reference));
 variance=scalarAlpha*noise.variance_reference; U=zeros(N,0); ceShrinkDiagonal=zeros(N,1);
 traceCn=n*sum(variance); traceCe=0; b=zeros(N,1);
end
blocks=cell(1,nx); splitBlocks=blocks; G=blocks; splitG=blocks; splitU=zeros(size(U));
Y=reshape(D(:)-b,n,nx); splitY=zeros(n,nx); stageFloor=zeros(1,nx);
for j=1:nx
 ids=(j-1)*n+(1:n); temporal=toeplitz(rho(j).^(0:n-1));
 blocks{j}=variance(j)*temporal+diag(ceShrinkDiagonal(ids));
 G{j}=models{j}.G; M=models{j}.M;
 % Propagate the full operator, not X*A: no nonstationary composition identity.
 splitG{j}=M*G{j}; splitY(:,j)=M*Y(:,j); splitU(ids,:)=M*U(ids,:);
 stageFloor(j)=o.covariance_floor*noise.variance_reference(j);
 splitBlocks{j}=M*blocks{j}*M'+stageFloor(j)*eye(n);
end
C=p1_covariance(blocks,U); Cs=p1_covariance(splitBlocks,splitU,stageFloor);
H=prior.H; gamma=o.gamma;
if strcmp(branch,'NoStructureCombined'), H=sparse(N,N); prior.H=H; gamma=0; end
% Fixed reference-noise units; both routes and all branches share this scale.
lambda=prior.lambda./noise.variance_reference; mu=gamma*median(lambda);
directCtx=p1_context(G,C,Y,H,lambda,mu,o); splitCtx=p1_context(splitG,Cs,splitY,H,lambda,mu,o);
out=p1_select(directCtx,splitCtx,prior.pilot);
out.prior=prior; out.estimated_Q=cache.estimate;
out.preset=struct('name',branch,'gamma',gamma,'pilot_multiplier',o.pilot_multiplier,'noise_target_multiplier',1);
out.calibration_cache=cache; out.options=o; out.selected_method=[branch '_Split']; out.gather=out.split;
out.diagnostics=struct('branch',branch,'truth_used',false,'fixed_pilot_and_covariance',true, ...
 'full_gather_covariance',true,'model_error_definition','(Gdraw-Gref)*fixed_observed_pilot', ...
 'mean_correction','D-b','model_error_mean_norm',norm(b), ...
 'model_error_ensemble_trace',ensembleTraceCe,'assembled_trace_Cn',traceCn,'assembled_trace_Ce',traceCe, ...
 'assembled_trace_Ceff',traceCn+traceCe,'covariance_low_rank_columns',size(U,2), ...
 'covariance_shrinkage',o.covariance_shrinkage,'scalar_budget_alpha',scalarAlpha, ...
 'noise_variance_reference',noise.variance_reference,'noise_variance_used',variance, ...
 'noise_upper_variance',noise.variance_upper,'noise_rho_used',rho,'noise_bootstrap_seed',noise.seed, ...
 'noise_quantile',o.noise_quantile,'noise_samples',noise.noise_samples, ...
 'q_frequency_draws',cache.q_frequency_draws,'q_bootstrap_seed',o.bootstrap_seed, ...
 'direct_numerical_covariance_jitter',C.numerical_jitter,'split_numerical_covariance_jitter',Cs.numerical_jitter, ...
 'split_floor_variance',stageFloor,'direct_target_squared',C.target_squared,'split_target_squared',Cs.target_squared, ...
 'ridge_lambda_reference_units',lambda,'mu_reference_units',mu, ...
 'split_forward','M*Gref (exact full operator propagation)', ...
 'direct_split_relative_difference',norm(out.direct-out.split,'fro')/max(norm(out.direct,'fro'),realmin), ...
 'bootstrap_is_confidence_interval',false,'s_frequency_role','frequency affects joint calibration only; G acts on bandlimited s');
out.calibration_estimates=struct('reference',cache.estimate,'paired_Q_frequency_draws',cache.q_frequency_draws, ...
 'bootstrap_estimates',{cache.bootstrap_estimates});
end
