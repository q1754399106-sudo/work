function N=ar1_noise(phi,sig,nt,seed)
% Stationary AR(1) noise in time with per-trace marginal standard deviation sig.
rng(seed,'twister');nx=numel(sig);burn=60;u=filter(1,[1 -phi],randn(nt+burn,nx));
N=u(burn+1:end,:)*sqrt(1-phi^2).*reshape(sig,1,[]);
end
