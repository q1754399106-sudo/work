function validate_unknown_q(outdir)
assert(~isfolder(outdir)); mkdir(outdir); [configs,~]=stability_models(); c=configs{2}; p=uq_public_config(c,52);
poison=c; poison.Q=[999,888,777,666]; poison.fpeak=999; poison.rho=999*c.rho; pp=uq_public_config(poison,52);
assert(isequal(p,pp),'Forbidden synthetic values leaked');
o=uq_spec(); offsets=c.offset; D=zeros(c.nt,11); coefficients=[.08;.01;-.05;.02;.07;-.01];
for j=1:11, B=uq_templates(p,offsets(j),65,55); D(:,j)=B*coefficients; end
rng(8901,'twister'); D=D+randn(size(D)).*sqrt(mean(D.^2,1))*1e-3;
est=estimate_q_effective(D,offsets,p);
qerr=abs(est.effective_Q/65-1); ferr=abs(est.estimated_frequency/55-1);
assert(qerr<.05 && ferr<.05);
cm=uq_inverse_config(p,est); assert(all(cm.Q==est.effective_Q) && cm.fpeak==est.estimated_frequency);
% Fixed prior: explicit Frozen parameters reproduce inherited selection.
small=cm; small.nt=201; small.t=(0:200)'*small.dt; small.nfft=2^nextpow2(4*small.nt);
models=prepare_models_v2([0,750,1500],small); z=zeros(201,3);
for j=1:3, x=pi*50*(small.t-.31); s=(1-2*x.^2).*exp(-x.^2); z(:,j)=models{j}.G*s+.0001*randn(201,1); end
pr=uq_build_prior(z,models,small,'Offset',o.presets(2)); a=uq_select(z,models,pr,o.presets(2)); b=select_v2(z,models,pr.H,pr.lambda,pr.sigma);
e=max(norm(a.direct-b.direct,'fro')/norm(b.direct,'fro'),norm(a.split-b.split,'fro')/norm(b.split,'fro'));
assert(e<1e-10);
checks=struct('pass',true,'public_input_poison_invariance',isequal(p,pp),'stationary_template_Q_relative_error',qerr, ...
 'stationary_template_frequency_relative_error',ferr,'frozen_selection_regression_relative_error',e, ...
 'note','Stationary-template recovery is a model-consistency test, not independent inverse validation.');
uq_writejson(fullfile(outdir,'checks.json'),checks); uq_writejson(fullfile(outdir,'estimator_unit_example.json'),est);
fprintf('PASS: public-input leakage, template recovery, frozen solver regression %.3e\n',e);
end
