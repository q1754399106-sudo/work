function [S,normal,replay]=fixed_v1_apply(D,op)
% Exact two-stage arithmetic of the earliest archived FixedSplit.
[n,nx]=size(D);S=zeros(n,nx);normal=0;replay=0;
N0=op.G0'*op.G0+op.lambda_fixed*eye(n);
for j=1:nx
 y=op.BX{j}*D(:,j);S(:,j)=op.B0*y;
 X=op.GX{j};r1=(X'*X+op.lambda_x*eye(n))*y-X'*D(:,j);
 r2=N0*S(:,j)-op.G0'*y;
 normal=max([normal,norm(r1)/max(norm(X'*D(:,j)),eps),norm(r2)/max(norm(op.G0'*y),eps)]);
 if j==6
  test=(op.B0*op.BX{j})*D(:,j);replay=norm(test-S(:,j))/max(norm(S(:,j)),eps);
 end
end
assert(normal<1e-10 && replay<1e-10);
end
