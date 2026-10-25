function [r2_D,r2_S,eta_eff,b_eD,b_eS]=p1_q_model_error(pilot,models,cm,domain,use_q,use_noise,p1)
% Compute singular-basis Q-model-error regularization weights and calibrated frequency prior R(s)
% driven directly by calibrated (Q, fpeak) posterior uncertainty (cm.q_uncertainty, cm.f_uncertainty, cm.log_qf_corr).
[n,nx]=size(pilot); fpeak=max(cm.fpeak,10); q_eff=cm.Q(1);
if isfield(cm,'effective_Q'),q_eff=cm.effective_Q;end
low_q_weight=exp(-max(0,q_eff-p1.low_q_reference)/p1.low_q_decay_scale);
if isfield(cm,'q_uncertainty'), u_q=cm.q_uncertainty; else, u_q=0.12; end
if isfield(cm,'f_uncertainty'), u_f_cal=cm.f_uncertainty; else, u_f_cal=0.14; end
if isfield(cm,'log_qf_corr'), rho_qf=cm.log_qf_corr; else, rho_qf=-0.75; end
uq_norm=min(1.35,max(0.75,u_q/0.12));
uf_norm=min(1.35,max(0.75,u_f_cal/0.14));
rho_abs=min(0.95,abs(rho_qf));
freqs=(0:floor(n/2))'/(n*cm.dt);
r2_D=cell(1,nx); r2_S=cell(1,nx);
b_eD=zeros(n,nx); b_eS=zeros(n,nx);
for j=1:nx
 for ir=1:2
  if ir==1, F=models{j}.direct; else, F=models{j}.split; end
  z=F.z; Vf=abs(fft(F.V,n,1)); Vf=Vf(1:numel(freqs),:).^2;
  f_m=(freqs'*Vf)./max(sum(Vf,1),realmin); f_m=f_m(:);
  u_f=f_m/fpeak;
  s_prior=(u_f.^2).*exp(1-min(u_f.^2,20));
  if use_q
   j_q=log(1./max(z,1e-4));
   j_f=max(0,u_f-p1.hf_cutoff_factor_q);
   q_sens=(p1.q_sensitivity_scale*low_q_weight*uq_norm*j_q).^2;
   hf_pen=(2.4*low_q_weight*uq_norm+0.38*uf_norm)*(j_f.^2) + 0.05*rho_abs*low_q_weight*uq_norm*uf_norm*(j_q.*j_f);
   pass_boost=(p1.passband_boost_lowq*low_q_weight+p1.passband_boost_base)/(uf_norm^1.2);
   r2=(1+hf_pen+q_sens)./(1+pass_boost*s_prior);
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
 uq_eta=0.75+0.25*uq_norm;
 if strcmp(domain,'Offset')
  eta_eff=p1.eta_offset_base+p1.eta_offset_lowq_add*low_q_weight*uq_eta;
 else
  eta_eff=p1.eta_section_base+p1.eta_section_lowq_add*low_q_weight*uq_eta+0.025*max(0,uf_norm-1.0);
 end
elseif use_noise
 if strcmp(domain,'Offset'), eta_eff=p1.eta_noise_only_offset; else, eta_eff=p1.eta_noise_only_section; end
else
 eta_eff=1.0;
end
end
