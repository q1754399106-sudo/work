function est=p1_calibrate_ensemble(D,offsets,p)
% Joint (Q, frequency) point fit and heuristic profile-loss uncertainty.
% This is not a bootstrap or a calibrated Bayesian posterior.
% Never reads true Q, true wavelet peak, or reflection coefficients.
assert(~isfield(p,'Q') && ~isfield(p,'fpeak') && ~isfield(p,'rho'),'Public inputs only');
est=estimate_q_effective(D,offsets,p);
% Extract joint (Q, f) uncertainty covariance from profile_losses surface
qgrid=est.q_grid(:); fgrid=est.frequency_grid(:)';
L=est.profile_losses; minL=min(L(:));
temp=max(est.near_optimal_threshold-minL,1.0);
W=exp(-max(0,L-minL)/temp); W=W/max(sum(W(:)),realmin);
[FF,QQ]=meshgrid(log(fgrid),log(qgrid));
mq=sum(W(:).*QQ(:)); mf=sum(W(:).*FF(:));
vq=sum(W(:).*(QQ(:)-mq).^2); vf=sum(W(:).*(FF(:)-mf).^2);
cov_qf=sum(W(:).*(QQ(:)-mq).*(FF(:)-mf));
est.log_q_std=sqrt(max(vq,0.01));
est.log_f_std=sqrt(max(vf,0.005));
cov_qf=max(-sqrt(vq*vf),min(sqrt(vq*vf),cov_qf));
est.log_parameter_covariance=[est.log_q_std^2,cov_qf;cov_qf,est.log_f_std^2];
est.log_qf_corr=cov_qf/(est.log_q_std*est.log_f_std);
est.uncertainty_interpretation='Heuristic profile-loss covariance with variance floors; not a confidence interval';
est.reliability=p1_calibration_reliability(est);
end
