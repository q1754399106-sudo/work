function cache=p1_prepare_calibration(public,calibrationD,calibrationOffsets,offsets,opts)
% Observed-only reusable calibration. Never accepts a synthetic truth structure.
if nargin<5, opts=struct; end
p1_check_public(public); o=p1_options(public,opts);
assert(size(calibrationD,1)==public.nt && size(calibrationD,2)==numel(calibrationOffsets));
assert(size(calibrationD,2)>=11 && all(isfinite(calibrationD(:))));
est=estimate_q_effective(calibrationD,calibrationOffsets,public);
c=uq_inverse_config(public,est); uniqueOffsets=unique([reshape(calibrationOffsets,1,[]) reshape(offsets,1,[])]);
allmodels=prepare_models_v2(uniqueOffsets,c); models=cell(1,numel(offsets));
for j=1:numel(uniqueOffsets)
 [a,b,~,p]=ray_absorption(c.t,uniqueOffsets(j),c); X=attenuation_matrix(b,c);
 M=ridge_operator(X,c.lambda_x);
 if norm(X-eye(c.nt),'fro')/sqrt(c.nt)<1e-10, M=eye(c.nt); end
 allmodels{j}.M=M; allmodels{j}.travel=(a+b)*est.effective_Q; allmodels{j}.p_time=p;
end
for j=1:numel(offsets), k=find(uniqueOffsets==offsets(j),1); models{j}=allmodels{k}; end
ns=p1_noise_summary(calibrationD,public.t,o,o.bootstrap_seed);
fitted=p1_fit_templates(calibrationD,calibrationOffsets,public,est);
stream=RandStream('mt19937ar','Seed',o.bootstrap_seed+1);
R=o.bootstrap_replicates; estimates=cell(1,R); qf=zeros(R,2); draws=cell(numel(uniqueOffsets),R);
for k=1:R
 ids=p1_block_indices(size(ns.centered_window,1),public.nt,o.block_length,stream);
 replicate=fitted+ns.centered_window(ids,:);
 ek=estimate_q_effective(replicate,calibrationOffsets,public); estimates{k}=ek;
 qf(k,:)=[ek.effective_Q ek.estimated_frequency];
 for j=1:numel(uniqueOffsets)
  draws{j,k}=attenuation_matrix(allmodels{j}.travel/ek.effective_Q,c);
 end
end
cache=struct('version','P1_FIXED_PILOT_V1','public',public,'calibrationD',calibrationD, ...
 'calibrationOffsets',reshape(calibrationOffsets,1,[]),'requested_offsets',reshape(offsets,1,[]), ...
 'options',o,'estimate',est,'inverse_config',c,'models',{models},'all_offsets',uniqueOffsets, ...
 'all_models',{allmodels},'operator_draws',{draws},'bootstrap_estimates',{estimates}, ...
 'q_frequency_draws',qf,'noise_summary',rmfield(ns,'centered_window'), ...
 'bootstrap_source','fitted nuisance templates plus joint circular blocks from the observed noise window', ...
 'truth_used',false,'confidence_interval_claim',false);
end
function fitted=p1_fit_templates(D,offsets,p,est)
o=uq_spec(); times=cumsum(2*p.h./p.vp(1:3)); mask=false(p.nt,1);
for t=times, mask=mask | abs(p.t-t)<=o.event_half_window_s; end
fitted=zeros(size(D));
for j=1:size(D,2)
 ni=estimate_trace_snr(D(:,j),p.t); sigma=max(ni.sigma,sqrt(realmin));
 B=uq_templates(p,offsets(j),est.effective_Q,est.estimated_frequency);
 X=B(mask,:)/sigma; y=D(mask,j)/sigma; weights=ones(size(y)); coef=zeros(6,1);
 for k=1:o.irls_iterations
  G=X'*(weights.*X); ridge=1e-8*max(trace(G)/6,realmin);
  coef=(G+ridge*eye(6))\(X'*(weights.*y));
  r=y-X*coef; weights=min(1,o.huber_delta./max(abs(r),realmin));
 end
 fitted(:,j)=B*coef;
end
end
