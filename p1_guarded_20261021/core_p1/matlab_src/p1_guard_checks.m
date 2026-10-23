function p1_guard_checks(outdir)
good=struct('fit_samples_per_trace',200,'selected_traces',[1,4,7,11], ...
 'minimum_loss',400,'effective_Q',40,'estimated_frequency',55,'q_boundary',false, ...
 'frequency_boundary',false,'optimizer_exitflag',1);
a=p1_calibration_reliability(good); assert(a.reliable);
bad=good; bad.q_boundary=true; bad.optimizer_exitflag=0; bad.minimum_loss=2000;
b=p1_calibration_reliability(bad); assert(~b.reliable && numel(b.reasons)==3);
cm=config_model(); cm.calibration_reliable=false;
rng(1081,'twister'); D=randn(cm.nt,11)*.001;
[out,~,~]=p1_solve_branch(D,{},[],cm,'Offset',sparse(numel(D),numel(D)),'GuardedNoise');
assert(isequal(D,out.direct) && isequal(D,out.split));
assert(~out.compensation_applied && out.fallback_used && strcmp(out.effective_method,'NoCompensation'));
assert(strcmp(out.pipeline_status,'REJECTED_CALIBRATION_PASSTHROUGH'));
uq_writejson(fullfile(outdir,'guard_checks.json'),struct('pass',true,'good_calibration_accepted',a.reliable, ...
 'bad_calibration_rejected',~b.reliable,'bad_reasons',{b.reasons}, ...
 'rejected_output_max_difference',max(abs(out.direct(:)-D(:))),'rejected_compensation_applied',out.compensation_applied));
fprintf('PASS: calibration rejection and explicit non-compensating pass-through\n');
end
