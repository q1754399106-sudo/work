function [Gs,factors]=q_sample_operators(cm,offset,s,M,rho,seed)
% Attenuation operators for Q_s=Qhat.*exp(s*(rho*xi0+sqrt(1-rho^2)*xi_k)).
% Antithetic pairs; fixed seed, so every scene uses the same normalised draws.
rng(seed,'twister');half=randn(M/2,1+numel(cm.Q));xi=[half;-half];
logf=s*(rho*xi(:,1)+sqrt(1-rho^2)*xi(:,2:end));factors=exp(logf);Gs=cell(1,M);t=cm.t(:);
for k=1:M
 c=cm;c.Q=reshape(cm.Q,1,[]).*factors(k,:);[a,b]=ray_absorption(t,offset,c);Gs{k}=attenuation_matrix(a+b,c);
end
end
