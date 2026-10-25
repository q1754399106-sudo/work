function y=p1_forward(G,x,transposeFlag)
% Block-diagonal full nonstationary forward / adjoint.
n=size(G{1},1); nx=numel(G); x=reshape(x,n,nx); y=zeros(n,nx);
for j=1:nx
 if transposeFlag, y(:,j)=G{j}'*x(:,j); else, y(:,j)=G{j}*x(:,j); end
end
y=y(:);
end
