function train=cv_training_models(models,mask)
% Mask before the fixed K transformation: K*M*D cannot access held-out data.
train=models;nx=numel(models);
for j=1:nx
 if j>5&&isequal(mask(:,j),mask(:,j-5))&&isequal(models{j}.G,models{j-5}.G)
  train{j}=train{j-5};continue;
 end
 m=models{j};m.G=models{j}.G.*mask(:,j);m.F=m.K*m.G;
 m.direct=discrepancy_factor(m.G);m.split=discrepancy_factor(m.F);
 m.gram_direct=m.G'*m.G;m.gram_split=m.F'*m.F;train{j}=m;
end
end
