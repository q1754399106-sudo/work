function ctx=p1_context(G,C,Y,H,lambda,mu,o)
n=size(Y,1); nx=size(Y,2); gramDiagonal=zeros(n,nx); grams=cell(1,nx); HHblocks=grams;
for j=1:nx
 B=C.L{j}\G{j}; gramDiagonal(:,j)=sum(B.^2,1)'; grams{j}=B'*B;
end
HH=H'*H;
for j=1:nx, ids=(j-1)*n+(1:n); HHblocks{j}=full(HH(ids,ids)); end
ctx=struct('G',{G},'C',C,'y',Y(:),'H',H,'HH',HH,'HHdiagonal',full(diag(HH)), ...
 'ridge',repelem(lambda(:),n),'mu',mu,'n',n,'nx',nx,'options',o, ...
 'gram_diagonal',gramDiagonal(:),'block_grams',{grams},'HHblocks',{HHblocks}, ...
 'rhs',p1_forward(G,p1_precision(C,Y(:)),true));
end
