function N=noise_replica(D0,snr_db,seed)
% Noise process of bae_q_20261024/make_independent.m, statement for statement:
% AR(1) 0.45 in time, 10% trace-common part, Frobenius SNR against the noise-free gather.
[nt,nx]=size(D0);rng(seed,'twister');u=randn(nt,nx);for i=2:nt,u(i,:)=u(i,:)+.45*u(i-1,:);end
u=.9*u+.1*mean(u,2);N=u*(norm(D0,'fro')/(10^(snr_db/20)*norm(u,'fro')));
end
