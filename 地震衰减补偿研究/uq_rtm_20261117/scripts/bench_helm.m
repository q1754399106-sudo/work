function bench_helm()
maxNumCompThreads(4);
for g=[4 5]
 nz=round(1600/g)+40;nx=round(2400/g)+40;N=nz*nx;
 e=ones(N,1);h=g;f=50;v=2500;k2=(2*pi*f/v)^2;
 Dx=spdiags([-e/12 4*e/3 -5*e/2 4*e/3 -e/12]/h^2,[-2*nz -nz 0 nz 2*nz],N,N);
 Dz=spdiags([-e/12 4*e/3 -5*e/2 4*e/3 -e/12]/h^2,[-2 -1 0 1 2],N,N);
 A=Dx+Dz+spdiags(k2*(1+.02i)*e,0,N,N);
 tic;[L,U,P,Q,R]=lu(A);t1=toc;b=zeros(N,100);b(sub2ind([nz nx],22*ones(1,100),round(linspace(30,nx-30,100))))=1;
 tic;x=Q*(U\(L\(P*(R\b))));t2=toc;
 fprintf('dx=%g: N=%d, lu %.1fs, 100 solves %.1fs, nnz(L+U)=%.1fM\n',g,N,t1,t2,(nnz(L)+nnz(U))/1e6);
end
end
