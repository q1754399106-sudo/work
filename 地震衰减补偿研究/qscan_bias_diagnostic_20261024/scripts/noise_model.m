function [phi,sig,tend]=noise_model(D,cm)
% Truth-free AR(1) noise model for the bootstrap. Window starts as estimate_trace_snr (0.02 s)
% but ends 1.5 dominant periods before the first public macro-interface time, so the first
% reflection's wavelet does not leak into it. sig: per-trace Gaussian-calibrated MAD.
tend=min(.25,2*cm.h(1)/cm.vp(1)-1.5/cm.fpeak);m=cm.t>=.02&cm.t<=tend;Z=D(m,:);
Z=Z-median(Z,1);sig=median(abs(Z),1)/.6744897501960817;Z=Z-mean(Z,1);
phi=sum(Z(2:end,:).*Z(1:end-1,:),'all')/sum(Z.^2,'all');
end
