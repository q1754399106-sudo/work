function est=p1_shrink_calibration(fixed,layered,public,weight)
% Geometric shrinkage is a point-estimate rule, not a posterior mixture.
assert(isscalar(weight)&&weight>=0&&weight<=1);
if weight==0,est=fixed;return;end
if weight==1,est=layered;return;end
est=fixed;
if isfield(layered,'layer_Q'),q=reshape(layered.layer_Q,1,[]);else,q=layered.effective_Q*ones(1,4);end
est.layer_Q=exp((1-weight)*log(fixed.effective_Q)+weight*log(q));
est.estimated_frequency=exp((1-weight)*log(fixed.estimated_frequency)+weight*log(layered.estimated_frequency));
travel=2*public.h./public.vp(1:3);est.effective_Q=sum(travel)/sum(travel./est.layer_Q(1:3));
est.calibration_shrinkage_weight=weight;
est.calibration_forward_model='geometric_shrinkage_effective_and_layered_Q';
est.uncertainty_scope='Inherited effective-Q/f regularization profile; not a validated shrinkage or layerwise posterior interval';
o=uq_spec();est.q_boundary=any(est.layer_Q<min(o.q_grid)*1.01|est.layer_Q>max(o.q_grid)*.99);
est.frequency_boundary=est.estimated_frequency<public.nominal_frequency*min(o.frequency_factors)*1.01||est.estimated_frequency>public.nominal_frequency*max(o.frequency_factors)*.99;
end
