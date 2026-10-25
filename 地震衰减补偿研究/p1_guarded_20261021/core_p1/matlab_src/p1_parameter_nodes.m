function [nodes,weights]=p1_parameter_nodes(cm)
% Positive five-point cubature matching mean/covariance in log(Q,f).
% Nodes average empirical modal penalties, not a full data-error covariance.
C=zeros(2); if isfield(cm,'log_parameter_covariance'), C=cm.log_parameter_covariance; end
assert(isequal(size(C),[2,2]) && all(isfinite(C(:))),'Invalid log parameter covariance');
C=(C+C')/2; [V,E]=eig(C); ev=diag(E);
assert(min(ev)>=-1e-12,'Log parameter covariance must be positive semidefinite');
L=V*diag(sqrt(max(ev,0))); center=log([cm.Q(1);cm.fpeak]);
delta=sqrt(2.5)*L;
nodes=exp(center+[zeros(2,1),delta,-delta]);
weights=ones(1,5)/5;
end
