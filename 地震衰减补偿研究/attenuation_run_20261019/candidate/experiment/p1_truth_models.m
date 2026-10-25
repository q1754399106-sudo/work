function models=p1_truth_models(c)
% SYNTHESIS ONLY: build exact G without unused inverse factorization/SVD.
models=cell(1,numel(c.offset));
for j=1:numel(c.offset)
 [a,b]=ray_absorption(c.t,c.offset(j),c);
 models{j}=struct('G',attenuation_matrix(a+b,c));
end
end
