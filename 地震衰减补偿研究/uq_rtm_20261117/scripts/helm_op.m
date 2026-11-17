function A=helm_op(v,f,P)
% Frequency-domain acoustic operator A = Lx + Lz + (w/v)^2 (constant density) on the full grid,
% Fourier convention of fft (X(w) = sum x(t) exp(-i w t)); A U = -F is the wave equation with source F.
% Interior: 4th-order second differences. Absorbing strips: 2nd-order stretched differences,
% s = 1 - i sigma/w, sigma = sigma0 (d/L)^2. Dirichlet outside the grid.
[nz,nx]=size(v);h=P.dx;w=2*pi*f;np=P.npml;L=np*h;s0=3*max(v(:))*log(1e4)/(2*L);
A=lap1(nz,np,h,w,s0,L);                                     % along z (rows, stride 1)
Bx=lap1(nx,np,h,w,s0,L);                                    % along x (columns, stride nz)
A=kron(speye(nx),A)+kron(Bx,speye(nz))+spdiags((w./v(:)).^2,0,nz*nx,nz*nx);
end

function D=lap1(n,np,h,w,s0,L)
d=@(y) min(max(max(np+1-y,y-(n-np)),0),np)*h;               % distance into the strip
s=@(y) 1-1i*s0*(d(y)/L).^2/w;
i=(1:n)';si=s(i);sp=s(i+.5);sm=s(i-.5);
pml=d(i)>0|d(i-1)>0|d(i+1)>0|d(i-2)>0|d(i+2)>0;          % 4th order only where all 5 nodes are unstretched
% 2nd-order stretched rows
c_m=1./(si.*sm)/h^2;c_p=1./(si.*sp)/h^2;c_0=-(c_m+c_p);
% 4th-order rows
e=ones(n,1);q=[-e/12 4*e/3 -5*e/2 4*e/3 -e/12]/h^2;
B=zeros(n,5);B(pml,2)=c_m(pml);B(pml,3)=c_0(pml);B(pml,4)=c_p(pml);B(~pml,:)=q(~pml,:);
% spdiags takes column k of the diagonal matrix from rows shifted by the offset: build by triplets
r=[];c=[];vals=[];off=-2:2;
for k=1:5
 j=i+off(k);ok=j>=1&j<=n&B(:,k)~=0;r=[r;i(ok)];c=[c;j(ok)];vals=[vals;B(ok,k)]; %#ok<AGROW>
end
D=sparse(r,c,vals,n,n);
end
