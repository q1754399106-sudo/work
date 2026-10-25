function [est,diagnostic]=p1_layered_calibrate_sharedghost(D,offsets,public,initial)
% SharedGhost (inputs/PLAN_ADDENDUM_SHAREDGHOST_ZH.md): p1_layered_calibrate with every
% template column replaced by B + g*B(t - tau), one (tau, g) shared by all events and
% traces, chosen on a grid by the training loss. Same starts, bounds, prior and gate.
assert(~any(isfield(public,{'Q','fpeak','rho','vs','event_t','truth'})));
taus=(4:2:30)*1e-3;gs=-.40:.05:.40;
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
optim=optimset('Display','off','MaxFunEvals',850,'MaxIter',600,'TolX',2e-4,'TolFun',1e-3);
point=[q0*ones(1,4),f0];gh=pick(point);history=gh;
[best,best_loss,fit]=optimise(gh);gh2=pick(best);
if any(gh2~=gh),gh=gh2;history(end+1,:)=gh;[best,best_loss,fit]=optimise(gh);end
point_validation=objective(log(point),validation,false,gh);layer_validation=objective(log(best),validation,false,gh);
gain=1-layer_validation/max(point_validation,realmin);accepted=gain>=.02;
est=initial;
if accepted
 est.layer_Q=best(1:4);est.estimated_frequency=best(5);dt=diff([0,times]);est.effective_Q=sum(dt)/sum(dt./best(1:3));
 est.minimum_loss=best_loss;[~,index]=min([fit.loss]);est.optimizer_exitflag=fit(index).exitflag;est.optimizer_evaluations=sum([fit.evaluations]);
 est.q_boundary=any(best(1:4)<qbounds(1)*1.01|best(1:4)>qbounds(2)*.99);
 est.frequency_boundary=best(5)<fbounds(1)*1.01||best(5)>fbounds(2)*.99;
 est.calibration_forward_model='nonstationary_layered_Q_shared_ghost';est.ghost_delay=gh(1);est.ghost_amplitude=gh(2);
end
est.layered_calibration_applied=accepted;est.ghost_delay=gh(1);est.ghost_amplitude=gh(2);
diagnostic=struct('layer_Q',best(1:4),'frequency',best(5),'ghost',gh,'ghost_history',history, ...
 'point_validation_loss',point_validation,'layered_validation_loss',layer_validation,'validation_gain',gain,'applied',accepted,'optimization_runs',fit);

 function gh=pick(params)
  L=zeros(numel(taus),numel(gs));
  for a=1:numel(taus),for b=1:numel(gs),L(a,b)=objective(log(params),train,true,[taus(a),gs(b)]);end,end
  [~,i]=min(L(:));[a,b]=ind2sub(size(L),i);gh=[taus(a),gs(b)];
 end
 function [best,best_loss,fit]=optimise(gh)
  best=point;best_loss=objective(log(best),train,true,gh);fit=struct([]);
  for k=1:size(starts,1)
   [x,loss,flag,output]=fminsearch(@(x)objective(x,train,true,gh),log(starts(k,:)),optim);
   rec=struct('parameters',exp(x),'loss',loss,'exitflag',flag,'evaluations',output.funcCount,'ghost',gh);
   if k==1,fit=rec;else,fit(k)=rec;end
   if isfinite(loss)&&loss<best_loss,best=exp(x);best_loss=loss;end
  end
 end
 function loss=objective(x,traces,with_prior,gh)
  value=exp(x);qs=value(1:4);fp=value(5);
  if any(qs<qbounds(1)|qs>qbounds(2))||fp<fbounds(1)||fp>fbounds(2),loss=1e12+sum(x.^2);return;end
  loss=0;
  if with_prior,loss=o.source_prior_weight*(log(fp/public.nominal_frequency)/o.source_prior_log_scale)^2+.5*sum(diff(log(qs)).^2);end
  lag=round(gh(1)/public.dt);
  for j=traces
   B=p1_layered_calibration_templates(public,offsets(j),qs,fp);B=B+gh(2)*[zeros(lag,6);B(1:end-lag,:)];
   X=B(mask,:)/max(sigma(j),realmin);y=D(mask,j)/max(sigma(j),realmin);weight=ones(size(y));
   for iteration=1:o.irls_iterations
    G=X'*(weight.*X);ridge=1e-8*max(trace(G)/6,realmin);coef=(G+ridge*eye(6))\(X'*(weight.*y));residual=y-X*coef;weight=min(1,o.huber_delta./max(abs(residual),realmin));
   end
   residual=abs(y-X*coef);small=residual<=o.huber_delta;
   loss=loss+sum(.5*residual(small).^2)+sum(o.huber_delta*(residual(~small)-.5*o.huber_delta));
  end
 end
end
