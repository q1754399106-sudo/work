function s=q_spread_energy(s,d,factors)
% VQ_g = mean over alternative Q (used Q x factor) of ||W r_g(alt) - W r_g||^2,
% refitting the same fixed support with the same ridge rule as
% reflectivity_refit_prepared. Data-only: no truth or true Q.
if nargin<3,factors=[.8 1.25];end
nx=numel(s.supports);s.VQ=cell(1,nx);for j=1:nx,s.VQ{j}=zeros(numel(s.E{j}),1);end
for f=factors
 cm=s.cm;cm.Q=cm.Q*f;da=prune_dictionary(cm,450);
 for j=1:nx
  sup=s.supports{j}(:);if isempty(sup)||~s.live(j),continue;end
  C=da.B(:,sup);[U,Z,V]=svd(C,'econ');z=diag(Z);top=max(z.^2);bottom=min(z.^2);
  reg=max(1e-8*top,(top-1e4*bottom)/(1e4-1));a=V*((z./(z.^2+reg)).*(U'*s.D(:,j)));
  r=a./da.scale(sup);g=s.groups{j};
  for k=1:max(g)
   ii=g==k;dif=d.W(:,sup(ii))*(r(ii)-s.R0(sup(ii),j));s.VQ{j}(k)=s.VQ{j}(k)+(dif'*dif)/numel(factors);
  end
 end
end
end
