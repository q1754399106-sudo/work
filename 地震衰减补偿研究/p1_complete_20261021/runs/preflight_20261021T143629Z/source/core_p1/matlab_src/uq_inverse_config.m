function c=uq_inverse_config(p,est)
assert(~isfield(p,'Q') && ~isfield(p,'fpeak'));
c=p; c.Q=repmat(est.effective_Q,1,4); c.fpeak=est.estimated_frequency;
c.event_t=cumsum(2*p.h./p.vp(1:3));
if isfield(est,'profile_losses') && ~isfield(est,'log_q_std')
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
end
if isfield(est,'log_q_std'), c.log_q_std=est.log_q_std; else, c.log_q_std=0.12; end
if isfield(est,'log_f_std'), c.log_f_std=est.log_f_std; else, c.log_f_std=0.12; end
if isfield(est,'log_qf_corr'), c.log_qf_corr=est.log_qf_corr; else, c.log_qf_corr=-0.75; end
if isfield(est,'q_sensitivity_range')
 q_span=0.5*log(max(est.q_sensitivity_range)/max(min(est.q_sensitivity_range),1.0));
else
 q_span=0.05;
end
if isfield(est,'frequency_sensitivity_range')
 f_span=0.5*log(max(est.frequency_sensitivity_range)/max(min(est.frequency_sensitivity_range),1.0));
else
 f_span=0.05;
end
if isfield(est,'q_uncertainty'), c.q_uncertainty=est.q_uncertainty; else, c.q_uncertainty=sqrt(c.log_q_std^2+0.25*q_span^2); end
if isfield(est,'f_uncertainty'), c.f_uncertainty=est.f_uncertainty; else, c.f_uncertainty=sqrt(c.log_f_std^2+0.25*f_span^2); end
end
