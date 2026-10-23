function out=reflectivity_prune_prepared(D,dicts,sigma,initial)
% Global group-evidence pruning over all sampled times. No scene labels/picks.
[n,nx]=size(D);S=zeros(n,nx);R=S;supports=cell(1,nx);removed=zeros(1,nx);normals=removed;
for j=1:nx
 d=dicts{j};support=initial.support{j};support=support(:);gap=max(1,round(.5/(d.frequency*d.dt)));
 [a,rss,nr]=fit(d.B,D(:,j),support);
 while ~isempty(support)
  group=cumsum([1;diff(support)>gap]);ids=unique(group);margin=Inf(size(ids));fits=cell(size(ids));
  for k=1:numel(ids)
   omit=group==ids(k);reduced=support(~omit);[ar,loss,normal]=fit(d.B,D(:,j),reduced);
   parameter_cost=(sum(omit)+2)*log(n)*sigma(j)^2;
   margin(k)=loss-rss-parameter_cost;fits{k}=struct('support',reduced,'a',ar,'rss',loss,'normal',normal);
  end
  [worst,k]=min(margin);if worst>=0,break;end
  best=fits{k};removed(j)=removed(j)+numel(support)-numel(best.support);support=best.support;a=best.a;rss=best.rss;nr=best.normal;
 end
 coeff=zeros(n,1);coeff(support)=a;R(:,j)=coeff./d.scale;S(:,j)=d.W*R(:,j);supports{j}=support;normals(j)=nr;
end
assert(all(isfinite(S(:)))&&max(normals)<1e-8);
out=initial;out.gather=S;out.reflectivity=R;out.support=supports;out.removed_coefficients=removed;
out.support_before_pruning=cellfun(@numel,initial.support);out.support_after=cellfun(@numel,supports);
if isfield(out,'ridge'),out=rmfield(out,'ridge');end
out.normal_residual=max(normals);out.pruning_rule='Cluster gap 0.5/f; remove group if RSS gain < sigma^2*(group coefficient count+2)*log(nt)';
out.selected_method='PrunedFullCalRefit';
end
function [a,rss,normal]=fit(B,y,support)
if isempty(support),a=[];rss=sum(y.^2);normal=0;return;end
C=B(:,support);[U,Z,V]=svd(C,'econ');z=diag(Z);top=max(z.^2);bottom=min(z.^2);
reg=max(1e-8*top,(top-1e4*bottom)/(1e4-1));a=V*((z./(z.^2+reg)).*(U'*y));
e=C*a-y;rss=sum(e.^2);normal=norm(C'*e+reg*a)/max(norm(C'*y),realmin);
end
