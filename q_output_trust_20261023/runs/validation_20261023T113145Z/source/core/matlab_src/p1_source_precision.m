function P=p1_source_precision(n,dt,frequency,floor_power)
% Reflecting-boundary stationary source covariance in an orthonormal DCT basis.
% Uses calibrated Ricker power only, not event positions or amplitudes.
assert(n>1 && dt>0 && frequency>0 && floor_power>0);
persistent cache_key cache_P
key=sprintf('%d_%.17g_%.17g_%.17g',n,dt,frequency,floor_power);
if isequal(key,cache_key),P=cache_P;return;end
k=0:n-1;U=sqrt(2/n)*cos(pi*((0:n-1)'+.5)*k/n);U(:,1)=1/sqrt(n);
f=k/(2*n*dt);u=f/frequency;power=u.^4.*exp(2*(1-u.^2));
precision=(1+floor_power)./(power+floor_power);
P=(U.*precision)*U';P=(P+P')/2;
cache_key=key;cache_P=P;
end
