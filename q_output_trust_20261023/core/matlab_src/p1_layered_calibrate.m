function [est,diagnostic]=p1_layered_calibrate(D,offsets,public,initial)
% Four bounded layer Qs plus wavelet frequency, with held-out trace validation.
% No access to true layer Q/frequency, clean amplitudes or reconstruction scores.
assert(~any(isfield(public,{'Q','fpeak','rho','vs','event_t','truth'})));
o=uq_spec();train=reshape(o.calibration_trace_indices,1,[]);validation=setdiff(1:numel(offsets),train);
times=cumsum(2*public.h./public.vp(1:3));mask=false(public.nt,1);
for t=times,mask=mask|abs(public.t-t)<=o.event_half_window_s;end
sigma=zeros(1,numel(offsets));
for j=1:numel(offsets),ni=estimate_trace_snr(D(:,j),public.t);sigma(j)=ni.sigma;end
qbounds=[min(o.q_grid),max(o.q_grid)];fbounds=public.nominal_frequency*[min(o.frequency_factors),max(o.frequency_factors)];
q0=initial.effective_Q;f0=initial.estimated_frequency;
starts=[q0*ones(1,4),f0;q0*[.7,1,1.3,1.6],public.nominal_frequency];
starts(:,1:4)=min(qbounds(2)*.99,max(qbounds(1)*1.01,starts(:,1:4)));
starts(:,5)=min(fbounds(2)*.99,max(fbounds(1)*1.01,starts(:,5)));
best=[q0*ones(1,4),f0];best_loss=objective(log(best),train,true);
fit=struct([]);optim=optimset('Display','off','MaxFunEvals',850,'MaxIter',600,'TolX',2e-4,'TolFun',1e-3);
for k=1:size(starts,1)
 [x,loss,flag,output]=fminsearch(@(x)objective(x,train,true),log(starts(k,:)),optim);
 rec=struct('parameters',exp(x),'loss',loss,'exitflag',flag,'evaluations',output.funcCount);
 if k==1,fit=rec;else,fit(k)=rec;end
 if isfinite(loss)&&loss<best_loss,best=exp(x);best_loss=loss;end
end
point=[q0*ones(1,4),f0];point_validation=objective(log(point),validation,false);
layer_validation=objective(log(best),validation,false);gain=1-layer_validation/max(point_validation,realmin);
accepted=gain>=.02;
est=initial;
if accepted
 est.layer_Q=best(1:4);est.estimated_frequency=best(5);
 dt=diff([0,times]);est.effective_Q=sum(dt)/sum(dt./best(1:3));
 est.minimum_loss=best_loss;est.optimizer_exitflag=fit(1).exitflag;
 [~,index]=min([fit.loss]);est.optimizer_exitflag=fit(index).exitflag;
 est.optimizer_evaluations=sum([fit.evaluations]);
 est.q_boundary=any(best(1:4)<qbounds(1)*1.01 | best(1:4)>qbounds(2)*.99);
 est.frequency_boundary=best(5)<fbounds(1)*1.01||best(5)>fbounds(2)*.99;
 est.calibration_forward_model='nonstationary_layered_Q';
 % Do not invent layerwise posterior coverage from the original effective-Q profile.
 est.uncertainty_scope='Conservative effective-Q/f initial profile; layerwise interval coverage unverified';
end
est.layered_calibration_applied=accepted;
diagnostic=struct('training_traces',train,'validation_traces',validation,'point_Q',q0,'point_frequency',f0, ...
 'layer_Q',best(1:4),'frequency',best(5),'point_validation_loss',point_validation, ...
 'layered_validation_loss',layer_validation,'validation_gain',gain,'required_gain',.02, ...
 'applied',accepted,'optimization_runs',fit,'adjacent_logQ_penalty_weight',.5, ...
 'scope','Bounded four-layer Q point estimation from known-geometry calibration; not a calibrated layerwise uncertainty interval');
 function loss=objective(x,traces,with_prior)
  value=exp(x);qs=value(1:4);fp=value(5);
  if any(qs<qbounds(1)|qs>qbounds(2))||fp<fbounds(1)||fp>fbounds(2)
   loss=1e12+sum(x.^2);return;
  end
  loss=0;
  if with_prior
   loss=o.source_prior_weight*(log(fp/public.nominal_frequency)/o.source_prior_log_scale)^2+.5*sum(diff(log(qs)).^2);
  end
  for j=traces
   B=p1_layered_calibration_templates(public,offsets(j),qs,fp);X=B(mask,:)/max(sigma(j),realmin);y=D(mask,j)/max(sigma(j),realmin);weight=ones(size(y));
   for iteration=1:o.irls_iterations
    G=X'*(weight.*X);ridge=1e-8*max(trace(G)/6,realmin);coef=(G+ridge*eye(6))\(X'*(weight.*y));residual=y-X*coef;weight=min(1,o.huber_delta./max(abs(residual),realmin));
   end
   residual=abs(y-X*coef);small=residual<=o.huber_delta;
   loss=loss+sum(.5*residual(small).^2)+sum(o.huber_delta*(residual(~small)-.5*o.huber_delta));
  end
 end
end
