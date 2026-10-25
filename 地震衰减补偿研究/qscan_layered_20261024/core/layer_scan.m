function out=layer_scan(D,sigma,used,offset,spec)
% Depth-segmented Q correction. L0: common factor c0 by the original QScan
% objective on the registered grid. Then for each beta in spec.betas: two
% top-down coordinate-descent passes over per-layer factors c_l = c0*g,
% g in exp(linspace(-r,r,m)), minimising
%   J(c) + beta*nlive*sum_l log(c_l/c0)^2,
% J = sum over live traces of (0.5||D-Bx||^2 + lambda||x||_1)/sigma^2 (same lambda as FullCalRefit).
% Layers whose Q exceeds 1e5 (water) are not scanned. Final output: FullCalRefit with the chosen Q.
[n,nx]=size(D);sigma=reshape(sigma,1,[]);live=any(D~=0,1);lambda=.5*sqrt(2*log(n))*sigma;nl=nnz(live);
cache=containers.Map('KeyType','char','ValueType','double');
 function J=obj(c)
  key=sprintf('%.6f,',log(c));if isKey(cache,key),J=cache(key);return;end
  u=used;u.Q=used.Q.*c;d=prune_dictionary(u,offset);X=refit_lasso(d.Gram,d.B'*D(:,live),lambda(live),d.L);
  r=D(:,live)-d.B*X;J=sum((.5*sum(r.^2,1)+lambda(live).*sum(abs(X),1))./sigma(live).^2);cache(key)=J;
 end
grid=unique([exp(linspace(log(spec.grid_min),log(spec.grid_max),spec.grid_points)) 1]);J0=arrayfun(@(g)obj(g*ones(size(used.Q))),grid);
[~,k]=min(J0);c0=grid(k);layers=find(used.Q<=1e5);g=exp(linspace(-spec.layer_range,spec.layer_range,spec.layer_points));
variants=struct('name',{},'factors',{},'objective',{});variants(1)=struct('name','L0','factors',c0*ones(size(used.Q)),'objective',J0(k));
for b=reshape(spec.betas,1,[])
 c=c0*ones(size(used.Q));
 for pass=1:spec.passes
  for l=layers
   best=Inf;bv=c(l);
   for q=1:numel(g)
    t=c;t(l)=c0*g(q);v=obj(t)+b*nl*sum(log(t(layers)/c0).^2);
    if v<best,best=v;bv=t(l);end
   end
   c(l)=bv;
  end
 end
 if b==0,name='L1';else,name=sprintf('L2_b%g',b);end
 variants(end+1)=struct('name',name,'factors',c,'objective',obj(c)); %#ok<AGROW>
end
G=struct();flags=0;normal=0;
for v=1:numel(variants)
 u=used;u.Q=used.Q.*variants(v).factors;d=prune_dictionary(u,offset);S=zeros(n,nx);
 for first=1:11:nx
  ids=first:min(first+10,nx);o=reflectivity_refit_prepared(D(:,ids),repmat({d},1,numel(ids)),sigma(ids));
  S(:,ids)=o.gather;flags=flags+o.optimization_flags;normal=max(normal,o.normal_residual);
 end
 G.(strrep(variants(v).name,'.','p'))=single(S);
end
d=prune_dictionary(used,offset);S=zeros(n,nx);
for first=1:11:nx,ids=first:min(first+10,nx);o=reflectivity_refit_prepared(D(:,ids),repmat({d},1,numel(ids)),sigma(ids));S(:,ids)=o.gather;end
G.FullCalRefit=single(S);
out=struct('c0',c0,'variants',variants,'gathers',G,'evaluations',cache.Count,'optimization_flags',flags,'max_normal_residual',normal,'used_Q',used.Q);
end
