function C=p1_covariance(blocks,U,addedVariance)
% C=blockdiag(blocks)+U*U'; Woodbury precision without a gather-sized matrix.
nx=numel(blocks); n=size(blocks{1},1); L=cell(1,nx); V=zeros(n*nx,size(U,2));
if nargin<3, addedVariance=zeros(1,nx); end
jitter=zeros(1,nx);
for j=1:nx
 B=(blocks{j}+blocks{j}')/2; [Lj,p]=chol(B,'lower');
 if p~=0
  jitter(j)=max(eps(max(abs(diag(B))))*n,realmin);
  for k=1:8
   [Lj,p]=chol(B+jitter(j)*eye(n),'lower');
   if p==0, break; end
   jitter(j)=10*jitter(j);
  end
  assert(p==0,'Positive covariance factorization failed');
 end
 L{j}=Lj; ids=(j-1)*n+(1:n); V(ids,:)=Lj\U(ids,:);
end
T=chol(eye(size(V,2))+V'*V,'lower');
C=struct('L',{L},'V',V,'T',T,'n',n,'nx',nx,'rank',size(U,2), ...
 'added_variance',reshape(addedVariance,1,[])+jitter,'numerical_jitter',jitter);
% Exact expected squared whitened residual for physical covariance before floor.
% tr(P*(Cphysical + floor I))=N; subtract each artificial diagonal addition.
correction=0;
for j=1:nx
 if C.added_variance(j)==0, continue; end
 ids=(j-1)*n+(1:n); B=L{j}\eye(n);
 A=L{j}'\(V(ids,:)/T');
 traceP=sum(B(:).^2)-sum(A(:).^2);
 correction=correction+C.added_variance(j)*max(traceP,0);
end
C.target_squared=max(n*nx-correction,realmin);
end
