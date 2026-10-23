function F = discrepancy_factor(G)
% Cache the SVD for repeated observations; lambda multiplies ||s||^2.
[F.U,S,F.V]=svd(G,'econ'); F.z=diag(S); F.m=size(G,1);
end
