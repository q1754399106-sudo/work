function [fp, f, P] = peak_freq(d, L, win)
% Dominant frequency of the data in a time window: mean amplitude spectrum over traces (Hann taper, 1 Hz padding),
% 15-point (15 Hz) smoothing, maximum between 5 and 120 Hz (Cheng et al. 2025, CJG, field example: spectral analysis of a shallow window).
i=L.t>=win(1)-1e-9 & L.t<=win(2)+1e-9; x=d(i,:); n=sum(i);
h=0.5-0.5*cos(2*pi*(0:n-1)'/(n-1)); nf=round(1/L.dt);
F=abs(fft(x.*h,nf)); F=mean(F(1:floor(nf/2),:),2); f=(0:floor(nf/2)-1)'/(nf*L.dt);
P=movmean(F,15); m=f>=5&f<=120; ff=f(m); [~,j]=max(P(m)); fp=ff(j);
end
