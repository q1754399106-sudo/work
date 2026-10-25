function [est,corrected,diagnostic]=p1_joint_calibrate(D,offsets,public)
% Joint effective Q, frequency and small global velocity correction.
% Calibration observations/public kinematics only. No synthetic truth inputs.
% Independent calibration traces gate the extra degree of freedom.
forbidden={'Q','fpeak','rho','vs','event_t','truth','coeff'};
assert(~any(isfield(public,forbidden)),'Public inputs only');
assert(size(D,1)==public.nt && size(D,2)==numel(offsets) && all(isfinite(D(:))));
fixed=p1_calibrate_ensemble(D,offsets,public); opts=uq_spec();
train=reshape(opts.calibration_trace_indices,1,[]);
validation=setdiff(1:numel(offsets),train);
assert(numel(validation)>=3,'Independent calibration validation traces required');
sigma=zeros(1,numel(offsets));
for j=1:numel(offsets), ni=estimate_trace_snr(D(:,j),public.t); sigma(j)=ni.sigma; end
times=cumsum(2*public.h./public.vp(1:3)); mask=false(public.nt,1);
for t=times, mask=mask | abs(public.t-t)<=opts.event_half_window_s; end
qbounds=[min(opts.q_grid),max(opts.q_grid)];
fbounds=public.nominal_frequency*[min(opts.frequency_factors),max(opts.frequency_factors)];
vbounds=[.94,1.06];
seed_q=unique([fixed.effective_Q,32,80]);
seed_f=unique([fixed.estimated_frequency,.85*public.nominal_frequency,1.15*public.nominal_frequency]);
seed_v=.94:.01:1.06;
starts=zeros(numel(seed_q)*numel(seed_f)*numel(seed_v)+1,4); n=0;
for qq=seed_q, for ff=seed_f, for vv=seed_v
 if qq<qbounds(1)||qq>qbounds(2)||ff<fbounds(1)||ff>fbounds(2), continue; end
 n=n+1; starts(n,:)=[qq,ff,vv,objective(log([qq,ff,vv]),train,true)];
end, end, end
n=n+1; starts(n,:)=[fixed.effective_Q,fixed.estimated_frequency,1,objective(log([fixed.effective_Q,fixed.estimated_frequency,1]),train,true)];
starts=starts(1:n,:); [~,order]=sort(starts(:,4));
best=starts(order(1),:); optim=optimset('Display','off','MaxFunEvals',350,'MaxIter',220,'TolX',1e-5,'TolFun',1e-3);
fits=struct([]);
for k=1:min(3,numel(order))
 [x,loss,flag,output]=fminsearch(@(x)objective(x,train,true),log(starts(order(k),1:3)),optim);
 rec=struct('parameters',exp(x),'training_loss',loss,'exitflag',flag,'evaluations',output.funcCount);
 if k==1,fits=rec;else,fits(k)=rec;end
 if isfinite(loss)&&loss<best(4),best=[exp(x),loss];end
end
fixed_validation=objective(log([fixed.effective_Q,fixed.estimated_frequency,1]),validation,false);
joint_validation=objective(log(best(1:3)),validation,false);
gain=1-joint_validation/max(fixed_validation,realmin);
at_boundary=best(3)<=vbounds(1)+.0005 || best(3)>=vbounds(2)-.0005;
accepted=gain>=.05 && abs(best(3)-1)>=.002 && ~at_boundary;
corrected=public; est=fixed; actual_scale=1;
if accepted
 corrected.vp=public.vp*best(3);
 est=p1_calibrate_ensemble(D,offsets,corrected);
 refit_validation=objective(log([est.effective_Q,est.estimated_frequency,best(3)]),validation,false);
 if refit_validation>fixed_validation*.95
  accepted=false; corrected=public; est=fixed;
 else
  actual_scale=best(3);
 end
else
 refit_validation=fixed_validation;
end
est.velocity_correction=actual_scale;
est.kinematic_correction_applied=accepted;
diagnostic=struct('method','calibration_trace_cross_validation_global_velocity', ...
 'training_traces',train,'validation_traces',validation,'velocity_bounds',vbounds, ...
 'validation_gain_required',.05,'minimum_velocity_change',.002, ...
 'fixed_training_loss',fixed.minimum_loss,'fixed_validation_loss',fixed_validation, ...
 'joint_training_loss',best(4),'joint_validation_loss',joint_validation, ...
 'refit_validation_loss',refit_validation,'validation_gain',gain, ...
 'proposed_Q',best(1),'proposed_frequency',best(2),'proposed_velocity_scale',best(3), ...
 'velocity_boundary',at_boundary,'applied',accepted,'actual_velocity_scale',actual_scale, ...
 'fixed_estimate',fixed,'optimization_runs',fits, ...
 'scope','Global kinematic correction identifiable from known calibration reflectors; not layerwise velocity inversion or uncertainty coverage validation');

 function loss=objective(x,traces,with_prior)
  pp=exp(x);qq=pp(1);ff=pp(2);vv=pp(3);
  if qq<qbounds(1)||qq>qbounds(2)||ff<fbounds(1)||ff>fbounds(2)||vv<vbounds(1)||vv>vbounds(2)
   loss=1e12+sum(x.^2);return;
  end
  p=public;p.vp=public.vp*vv;loss=0;
  if with_prior,loss=opts.source_prior_weight*(log(ff/public.nominal_frequency)/opts.source_prior_log_scale)^2;end
  for jj=traces
   B=uq_templates(p,offsets(jj),qq,ff);X=B(mask,:)/max(sigma(jj),realmin);y=D(mask,jj)/max(sigma(jj),realmin);
   weights=ones(size(y));
   for it=1:opts.irls_iterations
    G=X'*(weights.*X);ridge=1e-8*max(trace(G)/6,realmin);coef=(G+ridge*eye(6))\(X'*(weights.*y));
    residual=y-X*coef;weights=min(1,opts.huber_delta./max(abs(residual),realmin));
   end
   residual=abs(y-X*coef);small=residual<=opts.huber_delta;
   loss=loss+sum(.5*residual(small).^2)+sum(opts.huber_delta*(residual(~small)-.5*opts.huber_delta));
  end
 end
end
