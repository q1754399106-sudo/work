function est=estimate_q_effective(D,offsets,p)
assert(~isfield(p,'Q') && ~isfield(p,'fpeak') && ~isfield(p,'rho'),'Public inputs only');
o=uq_spec(); traces=reshape(o.calibration_trace_indices,1,[]); sigma=zeros(1,numel(traces));
times=cumsum(2*p.h./p.vp(1:3)); mask=false(p.nt,1);
for t=times, mask=mask | abs(p.t-t)<=o.event_half_window_s; end
for j=1:numel(traces), info=estimate_trace_snr(D(:,traces(j)),p.t); sigma(j)=info.sigma; end
qgrid=reshape(o.q_grid,1,[]); fgrid=reshape(o.frequency_factors,1,[])*p.nominal_frequency; losses=zeros(numel(qgrid),numel(fgrid));
for iq=1:numel(qgrid), for jf=1:numel(fgrid), losses(iq,jf)=objective(qgrid(iq),fgrid(jf)); end, end
[best,idx]=min(losses(:)); [iq,jf]=ind2sub(size(losses),idx); q=qgrid(iq); fp=fgrid(jf);
opts=optimset('Display','off','MaxFunEvals',100,'MaxIter',70,'TolX',1e-4,'TolFun',1e-3);
[v,value,exitflag,output]=fminsearch(@bounded,log([q,fp]),opts);
if isfinite(value) && value<best, q=exp(v(1)); fp=exp(v(2)); best=value; end
threshold=best+o.profile_loss_increment_per_sample*sum(mask)*numel(traces);
[qi,fi]=find(losses<=threshold); qrange=[q,q]; frange=[fp,fp];
if ~isempty(qi), qrange=[min(q,min(qgrid(qi))),max(q,max(qgrid(qi)))]; frange=[min(fp,min(fgrid(fi))),max(fp,max(fgrid(fi)))]; end
est=struct('effective_Q',q,'estimated_frequency',fp,'nominal_frequency',p.nominal_frequency,'q_grid',qgrid,'frequency_grid',fgrid, ...
 'profile_losses',losses,'minimum_loss',best,'near_optimal_threshold',threshold,'q_sensitivity_range',qrange,'frequency_sensitivity_range',frange, ...
 'sensitivity_is_confidence_interval',false,'q_boundary',q<1.01*qgrid(1) || q>.99*qgrid(end), ...
 'frequency_boundary',fp<1.01*fgrid(1) || fp>.99*fgrid(end),'optimizer_exitflag',exitflag,'optimizer_evaluations',output.funcCount, ...
 'selected_traces',traces,'estimated_sigma',sigma,'fit_samples_per_trace',sum(mask),'public_input_fields',{fieldnames(p)});
 function val=bounded(v)
  qq=exp(v(1)); ff=exp(v(2));
  if qq<qgrid(1) || qq>qgrid(end) || ff<fgrid(1) || ff>fgrid(end), val=1e12+sum(v.^2); else, val=objective(qq,ff); end
 end
 function loss=objective(qq,ff)
  loss=o.source_prior_weight*(log(ff/p.nominal_frequency)/o.source_prior_log_scale)^2;
  for jj=1:numel(traces)
   j=traces(jj); B=uq_templates(p,offsets(j),qq,ff); X=B(mask,:)/max(sigma(jj),realmin); y=D(mask,j)/max(sigma(jj),realmin);
   weights=ones(size(y)); coef=zeros(6,1);
   for it=1:o.irls_iterations
    G=X'*(weights.*X); ridge=1e-8*max(trace(G)/6,realmin); coef=(G+ridge*eye(6))\(X'*(weights.*y));
    r=y-X*coef; weights=min(1,o.huber_delta./max(abs(r),realmin));
   end
   r=abs(y-X*coef); small=r<=o.huber_delta;
   loss=loss+sum(.5*r(small).^2)+sum(o.huber_delta*(r(~small)-.5*o.huber_delta));
  end
 end
end
