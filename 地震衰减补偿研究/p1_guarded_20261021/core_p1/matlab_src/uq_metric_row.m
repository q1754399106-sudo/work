function row=uq_metric_row(S,scene,c,model,suite,seed,method,route,out,prior,preset,est)
row=struct('Model',model,'Suite',suite,'Scenario',scene.name,'Seed',seed,'Method',method,'Route',route,'Preset',preset.name);
m=score_stability(S,scene,c); fields=fieldnames(m); for i=1:numel(fields), row.(fields{i})=m.(fields{i}); end
if strcmp(route,'Direct'), r=out.direct_record; else, r=out.split_record; end
row.EffectiveQ=est.effective_Q; row.EstimatedFrequency=est.estimated_frequency;
row.Gamma=preset.gamma; row.PilotMultiplier=preset.pilot_multiplier; row.TargetMultiplier=preset.noise_target_multiplier;
row.Scale=out.scale; row.ScaleEndpoint=out.endpoint; row.TargetMismatch=out.residual_mismatch;
row.NormalResidual=r.normal_residual; row.FinalFlag=r.flag; row.Solver=r.solver;
row.HTruthRatio=norm(prior.H*scene.truth(:))/norm(scene.truth(:));
row.CalibrationReliable=est.reliability.reliable;
row.CalibrationFitDeviance=est.reliability.fit_deviance;
row.FallbackUsed=false; row.CompensationApplied=true; row.EffectiveMethod=method;
row.PipelineStatus='UNGUARDED_COMPENSATION';
if isfield(out,'fallback_used'), row.FallbackUsed=out.fallback_used; end
if isfield(out,'compensation_applied'), row.CompensationApplied=out.compensation_applied; end
if isfield(out,'effective_method'), row.EffectiveMethod=out.effective_method; end
if isfield(out,'pipeline_status'), row.PipelineStatus=out.pipeline_status; end
end
