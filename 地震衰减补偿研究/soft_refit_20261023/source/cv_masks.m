function masks=cv_masks(n,nx)
% Two deterministic validation folds, 8-sample blocks, independent of data.
% A trace-dependent shift limits simultaneous gaps across neighboring traces.
block=floor((0:n-1)'/8);masks=cell(1,2);folds=[0,2];
for k=1:2
 mask=true(n,nx);
 for j=1:nx,mask(:,j)=mod(block+2*(j-1),5)~=folds(k);end
 masks{k}=mask;
end
end
