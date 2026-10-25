function y=p1_precision(C,x)
% Apply Ceff^{-1}; x may contain several right-hand sides.
assert(size(x,1)==C.n*C.nx); z=zeros(size(x));
for j=1:C.nx, ids=(j-1)*C.n+(1:C.n); z(ids,:)=C.L{j}\x(ids,:); end
if C.rank>0, z=z-C.V*(C.T'\(C.T\(C.V'*z))); end
y=zeros(size(z));
for j=1:C.nx, ids=(j-1)*C.n+(1:C.n); y(ids,:)=C.L{j}'\z(ids,:); end
end
