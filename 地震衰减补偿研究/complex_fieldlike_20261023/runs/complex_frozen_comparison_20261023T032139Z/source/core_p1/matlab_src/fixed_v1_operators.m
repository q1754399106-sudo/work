function op=fixed_v1_operators(cm,offsets)
% Original run_rebuild FixedSplit formulas, including ridge at zero offset.
assert(cm.lambda_x==1e-3 && cm.lambda_fixed==1e-3);
[a0,~]=ray_absorption(cm.t,0,cm);G0=attenuation_matrix(a0,cm);
B0=ridge_operator(G0,cm.lambda_fixed);
op=struct('B0',B0,'G0',G0,'BX',{cell(1,numel(offsets))},'GX',{cell(1,numel(offsets))},'lambda_x',cm.lambda_x,'lambda_fixed',cm.lambda_fixed);
cache=containers.Map('KeyType','double','ValueType','any');
for j=1:numel(offsets)
 x=offsets(j);
 if isKey(cache,x),z=cache(x);else
  [~,ax]=ray_absorption(cm.t,x,cm);X=attenuation_matrix(ax,cm);
  z=struct('X',X,'BX',ridge_operator(X,cm.lambda_x));cache(x)=z;
 end
 op.BX{j}=z.BX;op.GX{j}=z.X;
end
end
