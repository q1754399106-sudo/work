function [est,diagnostic]=p1_exact_calibration_refine(D,offsets,public,initial)
% Q/f refinement with the exact nonstationary forward model, fixed kinematics.
assert(~any(isfield(public,{'Q','fpeak','rho','vs','event_t','truth'})),'Public inputs only');
o=uq_spec();train=reshape(o.calibration_trace_indices,1,[]);validation=setdiff(1:numel(offsets),train);
sigma=zeros(1,numel(offsets));
for j=1:numel(offsets),n=estimate_trace_snr(D(:,j),public.t);sigma(j)=n.sigma;end
times=cumsum(2*public.h./public.vp(1:3));mask=false(public.nt,1);
for t=times,mask=mask|abs(public.t-t)<=o.event_half_window_s;end
qgrid=reshape(o.q_grid,1,[]);fgrid=reshape(o.frequency_factors,1,[])*public.nominal_frequency;
losses=zeros(numel(qgrid),numel(fgrid));
for iq=1:numel(qgrid),for jf=1:numel(fgrid),losses(iq,jf)=objective(log([qgrid(iq),fgrid(jf)]),train,true);end,end
[minimum,index]=min(losses(:));[iq,jf]=ind2sub(size(losses),index);best=[qgrid(iq),fgrid(jf)];
initial_loss=objective(log([initial.effective_Q,initial.estimated_frequency]),train,true);
if initial_loss<minimum,best=[initial.effective_Q,initial.estimated_frequency];minimum=initial_loss;end
optim=optimset('Display','off','MaxFunEvals',200,'MaxIter',140,'TolX',1e-5,'TolFun',1e-3);
[x,value,exitflag,output]=fminsearch(@(x)objective(x,train,true),log(best),optim);
if isfinite(value)&&value<minimum,best=exp(x);minimum=value;end
est=initial;est.effective_Q=best(1);est.estimated_frequency=best(2);est.minimum_loss=minimum;
est.q_grid=qgrid;est.frequency_grid=fgrid;est.profile_losses=losses;
est.near_optimal_threshold=minimum+o.profile_loss_increment_per_sample*sum(mask)*numel(train);
[qi,fi]=find(losses<=est.near_optimal_threshold);
est.q_sensitivity_range=[best(1),best(1)];est.frequency_sensitivity_range=[best(2),best(2)];
if ~isempty(qi),est.q_sensitivity_range=[min(best(1),min(qgrid(qi))),max(best(1),max(qgrid(qi)))];est.frequency_sensitivity_range=[min(best(2),min(fgrid(fi))),max(best(2),max(fgrid(fi)))];end
temperature=max([est.near_optimal_threshold-minimum,.15*minimum,1]);
weight=exp(-max(0,losses-minimum)/temperature);weight=weight/sum(weight(:));
[FF,QQ]=meshgrid(log(fgrid),log(qgrid));mq=sum(weight(:).*QQ(:));mf=sum(weight(:).*FF(:));
vq=sum(weight(:).*(QQ(:)-mq).^2);vf=sum(weight(:).*(FF(:)-mf).^2);covariance=sum(weight(:).*(QQ(:)-mq).*(FF(:)-mf));
est.log_q_std=sqrt(max(vq,1e-4));est.log_f_std=sqrt(max(vf,1e-4));est.log_qf_corr=covariance/max(sqrt(vq*vf),realmin);
qs=.5*log(max(est.q_sensitivity_range)/max(min(est.q_sensitivity_range),1));fs=.5*log(max(est.frequency_sensitivity_range)/max(min(est.frequency_sensitivity_range),1));
est.q_uncertainty=sqrt(est.log_q_std^2+.25*qs^2);est.f_uncertainty=sqrt(est.log_f_std^2+.25*fs^2);
est.q_boundary=best(1)<1.01*qgrid(1)||best(1)>.99*qgrid(end);est.frequency_boundary=best(2)<1.01*fgrid(1)||best(2)>.99*fgrid(end);
est.optimizer_exitflag=exitflag;est.optimizer_evaluations=output.funcCount;
est.calibration_forward_model='nonstationary_G_times_Ricker';
diagnostic=struct('initial_Q',initial.effective_Q,'initial_frequency',initial.estimated_frequency, ...
 'initial_exact_training_loss',initial_loss,'refined_exact_training_loss',minimum, ...
 'initial_exact_validation_loss',objective(log([initial.effective_Q,initial.estimated_frequency]),validation,false), ...
 'refined_exact_validation_loss',objective(log(best),validation,false),'optimizer_exitflag',exitflag, ...
 'training_traces',train,'validation_traces',validation,'basis_tail_threshold',1e-10);
 function loss=objective(x,traces,with_prior)
  v=exp(x);qq=v(1);ff=v(2);
  if qq<qgrid(1)||qq>qgrid(end)||ff<fgrid(1)||ff>fgrid(end),loss=1e12+sum(x.^2);return;end
  loss=0;if with_prior,loss=o.source_prior_weight*(log(ff/public.nominal_frequency)/o.source_prior_log_scale)^2;end
  for j=traces
   B=p1_exact_calibration_templates(public,offsets(j),qq,ff);X=B(mask,:)/max(sigma(j),realmin);y=D(mask,j)/max(sigma(j),realmin);weight=ones(size(y));
   for iteration=1:o.irls_iterations
    G=X'*(weight.*X);ridge=1e-8*max(trace(G)/6,realmin);coef=(G+ridge*eye(6))\(X'*(weight.*y));residual=y-X*coef;weight=min(1,o.huber_delta./max(abs(residual),realmin));
   end
   residual=abs(y-X*coef);small=residual<=o.huber_delta;
   loss=loss+sum(.5*residual(small).^2)+sum(o.huber_delta*(residual(~small)-.5*o.huber_delta));
  end
 end
end
