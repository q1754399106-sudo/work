function [S,R,gains]=apply_gain(s,c,opts)
% gain = E/(E + lambda(c) V + kappa VQ); hybrid deletion sets gain=0 when c<tau.
% Fallback groups keep lambda=1, no deletion, no VQ (existing soft rule).
% s: scene (R0, supports, E, V, groups, fallback, optional VQ, W); c: per-trace cells.
nx=numel(s.supports);R=zeros(size(s.R0));gains=cell(1,nx);
if ~isfield(opts,'tau'),opts.tau=0;end
if ~isfield(opts,'kappa'),opts.kappa=0;end
for j=1:nx
 sup=s.supports{j}(:);if isempty(sup),gains{j}=[];continue;end
 e=s.E{j}(:);v=s.V{j}(:);fb=s.fallback{j}(:);cc=c{j}(:);
 lam=opts.lambda_min+(opts.lambda_max-opts.lambda_min)*(1-cc).^2;lam(fb)=1;
 extra=zeros(size(e));
 if opts.kappa>0
  extra=opts.kappa*s.VQ{j}(:);
  if isfield(opts,'vq_coherence')&&opts.vq_coherence,extra=extra.*(1-cc).^2;end
  extra(fb)=0;
 end
 den=e+lam.*v+extra;g=ones(size(e));pos=den>0;g(pos)=e(pos)./den(pos);
 g(~fb&cc<opts.tau)=0;
 a=s.R0(sup,j).*g(s.groups{j});R(sup,j)=a;gains{j}=g;
end
S=s.W*R;
end
