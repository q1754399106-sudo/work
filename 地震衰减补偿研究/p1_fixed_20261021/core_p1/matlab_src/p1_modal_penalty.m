function r2=p1_modal_penalty(z,f_m,q,fp,p1)
low_q_weight=exp(-max(0,q-p1.low_q_reference)/p1.low_q_decay_scale);
u=f_m/max(fp,10); s_prior=u.^2.*exp(1-min(u.^2,20));
q_sens=(p1.q_sensitivity_scale*low_q_weight*log(1./max(z,1e-4))).^2;
hf_pen=(2.5*low_q_weight+0.4)*max(0,u-p1.hf_cutoff_factor_q).^2;
boost=p1.passband_boost_lowq*low_q_weight+p1.passband_boost_base;
r2=(1+hf_pen+q_sens)./(1+boost*s_prior);
end
