function est=p1_calibrate_ensemble(D,offsets,p)
% Joint (Q, frequency) calibration with profile-likelihood posterior covariance.
% Never reads true Q, true wavelet peak, or reflection coefficients.
assert(~isfield(p,'Q') && ~isfield(p,'fpeak') && ~isfield(p,'rho'),'Public inputs only');
est=estimate_q_effective(D,offsets,p);
% Extract joint (Q, f) uncertainty covariance and interval span from profile_losses surface
qgrid=est.q_grid(:); fgrid=est.frequency_grid(:)';
L=est.profile_losses; minL=min(L(:));
temp=max([est.near_optimal_threshold-est.minimum_loss, 0.15*minL, 1.0]);
W=exp(-max(0,L-minL)/temp); W=W/max(sum(W(:)),realmin);
[FF,QQ]=meshgrid(log(fgrid),log(qgrid));
mq=sum(W(:).*QQ(:)); mf=sum(W(:).*FF(:));
vq=sum(W(:).*(QQ(:)-mq).^2); vf=sum(W(:).*(FF(:)-mf).^2);
cov_qf=sum(W(:).*(QQ(:)-mq).*(FF(:)-mf));
est.log_q_std=sqrt(max(vq,1e-4));
est.log_f_std=sqrt(max(vf,1e-4));
est.log_qf_corr=cov_qf/max(sqrt(vq*vf),realmin);
q_span=0.5*log(max(est.q_sensitivity_range)/max(min(est.q_sensitivity_range),1.0));
f_span=0.5*log(max(est.frequency_sensitivity_range)/max(min(est.frequency_sensitivity_range),1.0));
est.q_uncertainty=sqrt(est.log_q_std^2+0.25*q_span^2);
est.f_uncertainty=sqrt(est.log_f_std^2+0.25*f_span^2);
end
