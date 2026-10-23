function [R,info]=avo_projector_v2(models)
% R=I-P; per-time AVO subspace span[1,p(t,offset)^2]. Column-major layout.
% No clean signal, event picks, event amplitudes or true Q are inputs.
n=size(models{1}.G,1); nx=numel(models); ray=zeros(n,nx);
for j=1:nx, ray(:,j)=models{j}.p_time; end
row=zeros(n*nx*nx,1); col=row; value=row; count=0; annih=0; symmetry=0; idempotence=0;
for t=1:n
 q=ray(t,:)'.^2; q=q/max(max(q),realmin); B=[ones(nx,1),q];
 P=B*pinv(B); Rt=eye(nx)-P; ids=t+(0:nx-1)*n;
 [rr,cc]=ndgrid(ids,ids); take=count+(1:nx*nx);
 row(take)=rr(:); col(take)=cc(:); value(take)=Rt(:); count=count+nx*nx;
 annih=max(annih,norm(Rt*B,'fro')/norm(B,'fro'));
 symmetry=max(symmetry,norm(Rt-Rt','fro')); idempotence=max(idempotence,norm(Rt*Rt-Rt,'fro'));
end
R=sparse(row,col,value,n*nx,n*nx);
info=struct('annihilation_error',annih,'symmetry_error',symmetry,'idempotence_error',idempotence);
end
