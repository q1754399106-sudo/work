function guard=p1_calibration_reliability(est)
% Uses calibration diagnostics only; no truth or reconstruction score.
root=fileparts(fileparts(mfilename('fullpath')));
policy=jsondecode(fileread(fullfile(root,'inputs','guard_policy.json')));
n=est.fit_samples_per_trace*numel(est.selected_traces);
deviance=2*est.minimum_loss/max(n,1);
reasons={};
if ~isfinite(est.effective_Q) || est.effective_Q<=0, reasons{end+1}='invalid_Q'; end
if ~isfinite(est.estimated_frequency) || est.estimated_frequency<=0, reasons{end+1}='invalid_frequency'; end
if est.q_boundary, reasons{end+1}='Q_search_boundary'; end
if est.frequency_boundary, reasons{end+1}='frequency_search_boundary'; end
if est.optimizer_exitflag<=0, reasons{end+1}='optimizer_not_converged'; end
if ~isfinite(deviance) || deviance>policy.fit_deviance_max, reasons{end+1}='calibration_fit_mismatch'; end
guard=struct('reliable',isempty(reasons),'reasons',{reasons},'fit_deviance',deviance, ...
 'fit_deviance_max',policy.fit_deviance_max,'interpretation',policy.interpretation);
end
