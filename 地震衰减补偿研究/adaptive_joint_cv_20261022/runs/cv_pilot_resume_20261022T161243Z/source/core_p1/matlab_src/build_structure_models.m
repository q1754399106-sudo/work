function models=build_structure_models(offsets,c)
% Independent columns: exact total operator versus whitened two-stage XA.
% K whitens the first-stage ridge output with a finite covariance floor.
models=cell(1,numel(offsets)); cache=containers.Map('KeyType','double','ValueType','any');
for j=1:numel(offsets)
 x=offsets(j);
 if isKey(cache,x), models{j}=cache(x); continue; end
 [a,b]=ray_absorption(c.t,x,c); A=attenuation_matrix(a,c);
 X=attenuation_matrix(b,c); G=attenuation_matrix(a+b,c);
 M=ridge_operator(X,c.lambda_x);
 if norm(X-eye(c.nt),'fro')/sqrt(c.nt)<1e-10, M=eye(c.nt); end
 K=chol(M*M'+c.cov_floor*eye(c.nt),'lower')\M;
 F=K*(X*A); fd=discrepancy_factor(G); fs=discrepancy_factor(F);
 m=struct('G',G,'K',K,'F',F,'direct',fd,'split',fs, ...
  'composition_error',norm(G-X*A,'fro')/norm(G,'fro'));
 models{j}=m; cache(x)=m;
end
end
