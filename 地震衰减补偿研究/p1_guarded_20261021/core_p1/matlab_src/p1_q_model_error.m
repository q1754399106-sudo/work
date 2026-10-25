function [r2_D,r2_S,eta_eff,b_eD,b_eS]=p1_q_model_error(pilot,models,cm,domain,use_q,use_noise,p1)
% Average empirical modal penalties over calibrated log(Q,f) uncertainty.
% No claim of a full forward-model-error covariance or Bayesian inversion.
[n,nx]=size(pilot); fpeak=max(cm.fpeak,10); q_eff=cm.Q(1);
low_q_weight=exp(-max(0,q_eff-p1.low_q_reference)/p1.low_q_decay_scale);
freqs=(0:floor(n/2))'/(n*cm.dt);
r2_D=cell(1,nx); r2_S=cell(1,nx);
b_eD=zeros(n,nx); b_eS=zeros(n,nx);
[nodes,weights]=p1_parameter_nodes(cm);
for j=1:nx
 for ir=1:2
  if ir==1, F=models{j}.direct; else, F=models{j}.split; end
  z=F.z; Vf=abs(fft(F.V,n,1)); Vf=Vf(1:numel(freqs),:).^2;
  f_m=(freqs'*Vf)./max(sum(Vf,1),realmin); f_m=f_m(:);
  u_f=f_m/fpeak;
  s_prior=(u_f.^2).*exp(1-min(u_f.^2,20));
  if use_q
   r2=zeros(size(z));
   for k=1:numel(weights)
    r2=r2+weights(k)*p1_modal_penalty(z,f_m,nodes(1,k),nodes(2,k),p1);
   end
  elseif use_noise
   hf_pen=0.6*max(0,u_f-p1.hf_cutoff_factor_noise).^2;
   r2=(1+hf_pen)./(1+0.06*s_prior);
  else
   r2=ones(size(z));
  end
  if ir==1, r2_D{j}=r2; else, r2_S{j}=r2; end
 end
end
if use_q
 if strcmp(domain,'Offset')
  eta_eff=p1.eta_offset_base+p1.eta_offset_lowq_add*low_q_weight;
 else
  eta_eff=p1.eta_section_base+p1.eta_section_lowq_add*low_q_weight;
 end
elseif use_noise
 if strcmp(domain,'Offset'), eta_eff=p1.eta_noise_only_offset; else, eta_eff=p1.eta_noise_only_section; end
else
 eta_eff=1.0;
end
end
