function phi=noise_ar1(D,t,window)
% Pooled lag-1 autocorrelation of the demeaned noise window (the window estimate_trace_snr uses).
m=t>=window(1)&t<=window(2);Z=D(m,:)-mean(D(m,:),1);
phi=sum(Z(2:end,:).*Z(1:end-1,:),'all')/sum(Z.^2,'all');
end
