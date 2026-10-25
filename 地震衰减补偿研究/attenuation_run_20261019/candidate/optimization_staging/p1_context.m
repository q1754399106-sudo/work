function ctx=p1_context(G,C,Y,H,lambda,mu,o)
n=size(Y,1); nx=size(Y,2); gramDiagonal=zeros(n,nx); grams=cell(1,nx); HHblocks=grams;
whiteG=cell(1,nx); whiteY=zeros(n*nx,1); Z=zeros(n*nx,C.rank); rhsBase=zeros(n*nx,1);
for j=1:nx
 ids=(j-1)*n+(1:n); B=C.L{j}\G{j}; whiteG{j}=B;
 whiteY(ids)=C.L{j}\Y(:,j); rhsBase(ids)=B'*whiteY(ids);
 Z(ids,:)=B'*C.V(ids,:);
 gramDiagonal(:,j)=sum(B.^2,1)'; grams{j}=B'*B;
end
HH=H'*H;
for j=1:nx, ids=(j-1)*n+(1:n); HHblocks{j}=full(HH(ids,ids)); end
rhsFast=rhsBase;
if C.rank>0, rhsFast=rhsFast-Z*(C.T'\(C.T\(C.V'*whiteY))); end
% Retain the original precision-based RHS exactly as the scientific definition.
% Only the repeated matrix-vector normal product is algebraically rearranged.
rhsOriginal=p1_forward(G,p1_precision(C,Y(:)),true);
ctx=struct('G',{G},'C',C,'y',Y(:),'H',H,'HH',HH,'HHdiagonal',full(diag(HH)), ...
 'ridge',repelem(lambda(:),n),'mu',mu,'n',n,'nx',nx,'options',o, ...
 'gram_diagonal',gramDiagonal(:),'block_grams',{grams},'HHblocks',{HHblocks}, ...
 'white_G',{whiteG},'white_y',whiteY,'normal_low_rank',Z,'rhs_fast',rhsFast, ...
 'rhs_rearrangement_relative_error',norm(rhsFast-rhsOriginal)/max(norm(rhsOriginal),realmin), ...
 'rhs',rhsOriginal,'numerics_version','GRAM_LOW_RANK_MG_CACHE_V2');
end
