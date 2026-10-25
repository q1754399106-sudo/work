function o=p1_options(public,opts)
% Frozen defaults for the first, deliberately finite P1 prototype.
if nargin<2 || isempty(opts), opts=struct; end
o=struct('bootstrap_replicates',8,'bootstrap_seed',2026101901, ...
 'noise_bootstrap_replicates',64,'block_length',8,'noise_quantile',.8, ...
 'covariance_shrinkage',.1,'noise_window',[.02 .25], ...
 'noise_covariance_model','white','covariance_floor',public.cov_floor, ...
 'pcg_tolerance',1e-8,'pcg_max_iterations',1000,'preconditioner','block','gamma',10, ...
 'pilot_multiplier',4,'scale_grid',[.01 .03162277660168379 .1 .3162277660168379 1 3.16227766016838 10], ...
 'residual_tolerance',.005,'refinement_max_steps',8,'rho_limit',.9, ...
 'noise_bootstrap_seed_offset',9109);
keys=fieldnames(o);
for k=1:numel(keys), key=keys{k}; if isfield(opts,key), o.(key)=opts.(key); end, end
assert(o.bootstrap_replicates>=8 && o.bootstrap_replicates==floor(o.bootstrap_replicates));
assert(o.noise_bootstrap_replicates>=8 && o.noise_bootstrap_replicates==floor(o.noise_bootstrap_replicates));
assert(o.block_length>=1 && o.block_length==floor(o.block_length));
assert(o.noise_quantile>=.5 && o.noise_quantile<1);
assert(o.covariance_shrinkage>=0 && o.covariance_shrinkage<=1);
assert(o.covariance_floor>0 && o.pcg_tolerance>0 && o.pcg_max_iterations>=1);
assert(o.gamma>=0 && o.pilot_multiplier>0 && all(o.scale_grid>0));
assert(all(diff(o.scale_grid)>0) && numel(o.scale_grid)>=2);
assert(o.rho_limit>=0 && o.rho_limit<1);
assert(any(strcmp(o.noise_covariance_model,{'white','ar1'})));
assert(any(strcmp(o.preconditioner,{'block','diagonal'})));
end
