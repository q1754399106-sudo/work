function validate_p1_unknown_q(outdir)
% Validate P1 truth-isolation poison invariance and numerical solver quality.
assert(~isfolder(outdir)); mkdir(outdir);
[configs,~]=p1_models('historical'); c=configs{1}; p=uq_public_config(c,c.nominal_frequency);
poison=c; poison.Q=[999,888,777,666]; poison.fpeak=999; poison.rho=999*c.rho; poison.vs=999*c.vs;
pp=uq_public_config(poison,c.nominal_frequency);
assert(isequal(p,pp),'Forbidden synthetic truth values leaked into public config');
offsets=c.offset; D=zeros(c.nt,11); coefficients=[.08;.01;-.05;.02;.07;-.01];
for j=1:11, B=uq_templates(p,offsets(j),65,55); D(:,j)=B*coefficients; end
rng(8901,'twister'); D=D+randn(size(D)).*sqrt(mean(D.^2,1))*1e-3;
est=p1_calibrate_ensemble(D,offsets,p);
qerr=abs(est.effective_Q/65-1); ferr=abs(est.estimated_frequency/55-1);
assert(qerr<.05 && ferr<.05);
out_f=p1_compensate(D,offsets,p,D,offsets,'Offset','Frozen');
out_c=p1_compensate(D,offsets,p,D,offsets,'Offset','CombinedUncertainty');
assert(out_f.direct_record.normal_residual<=2e-7 && out_f.split_record.normal_residual<=2e-7);
assert(out_c.direct_record.normal_residual<=2e-7 && out_c.split_record.normal_residual<=2e-7);
checks=struct('pass',true,'public_input_poison_invariance',isequal(p,pp), ...
 'stationary_template_Q_relative_error',qerr,'stationary_template_frequency_relative_error',ferr, ...
 'combined_direct_normal_residual',out_c.direct_record.normal_residual, ...
 'combined_split_normal_residual',out_c.split_record.normal_residual);
uq_writejson(fullfile(outdir,'checks.json'),checks);
fprintf('PASS: P1 public-input poison invariance and solver residual checks\n');
end
