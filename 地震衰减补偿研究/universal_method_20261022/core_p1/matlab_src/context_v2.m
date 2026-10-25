function ctx=context_v2(D,models,H,lambda,gamma)
[n,nx]=size(D); rhsD=zeros(n,nx); rhsS=rhsD; uD=rhsD; uS=rhsD; Yw=rhsD;
for j=1:nx
 m=models{j}; Yw(:,j)=m.K*D(:,j); rhsD(:,j)=m.G'*D(:,j); rhsS(:,j)=m.F'*Yw(:,j);
 uD(:,j)=m.direct.U'*D(:,j); uS(:,j)=m.split.U'*Yw(:,j);
end
HH=H'*H; blocks=cell(1,nx);
for j=1:nx, ids=(j-1)*n+(1:n); blocks{j}=full(HH(ids,ids)); end
ctx=struct('D',D,'Yw',Yw,'models',{models},'H',H,'HH',HH,'HHblocks',{blocks}, ...
 'lambda',lambda,'mu_base',gamma*median(lambda),'rhsD',rhsD,'rhsS',rhsS,'uD',uD,'uS',uS,'n',n,'nx',nx);
end
