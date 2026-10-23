function B = ridge_operator(G,lambda)
% Solve the linear system; never explicitly invert G'*G.
B=(G'*G+lambda*eye(size(G,2)))\G';
end
