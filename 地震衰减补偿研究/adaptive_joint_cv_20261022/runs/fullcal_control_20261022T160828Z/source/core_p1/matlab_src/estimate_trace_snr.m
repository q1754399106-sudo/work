function info = estimate_trace_snr(d,t,opts)
% Estimate observed (attenuated-signal) SNR from data ONLY. No clean trace/Q/SNR.
% Default: raw Gaussian-calibrated MAD in a prescribed noise-only time window.
% method='haar' uses all nonoverlapping first-level Haar details (blind).
% Window selection and white stationary Gaussian noise are explicit assumptions.
if nargin<3, opts=struct; end
if ~isfield(opts,'method'), opts.method='window'; end
if ~isfield(opts,'noise_window'), opts.noise_window=[0.02 0.25]; end
d=d(:); t=t(:);
assert(numel(d)==numel(t) && numel(d)>=16 && all(isfinite(d)) && all(isfinite(t)), ...
    'Input must be finite data/time vectors of equal length, at least 16 points.');
normalMAD=0.6744897501960817;
switch lower(string(opts.method))
 case 'window'
  w=opts.noise_window;
  assert(numel(w)==2 && w(2)>w(1),'noise_window must be [start end].');
  ii=t>=w(1) & t<=w(2); z=d(ii);
  assert(numel(z)>=32,'At least 32 samples are required in the prescribed noise window.');
  sigma=median(abs(z-median(z)))/normalMAD;
  noiseRegion=z; method='WindowMAD';
 case 'haar'
  k=floor(numel(d)/2); z=(d(2:2:2*k)-d(1:2:2*k))/sqrt(2);
  sigma=median(abs(z-median(z)))/normalMAD;
  noiseRegion=[]; method='BlindHaarMAD';
 otherwise
  error('Unknown noise estimator: %s',string(opts.method));
end
% Demean to avoid interpreting arbitrary recording DC as seismic signal energy.
totalPower=sum((d-mean(d)).^2)/(numel(d)-1);
noisePower=sigma^2; signalPower=max(totalPower-noisePower,0);
if noisePower==0
 if signalPower>0, snrdb=Inf; else, snrdb=NaN; end
elseif signalPower==0, snrdb=-Inf;
else, snrdb=10*log10(signalPower/noisePower);
end
rho1=NaN; whiteOK=NaN;
if ~isempty(noiseRegion)
 e=noiseRegion-mean(noiseRegion);
 rho1=(e(1:end-1)'*e(2:end))/max(e'*e,realmin);
 % Diagnostic flag, not a formal whiteness test or a covariance estimator.
 whiteOK=abs(rho1)<=max(0.25,2/sqrt(numel(e)));
end
info=struct('sigma',sigma,'snr_db',snrdb,'noise_power',noisePower, ...
 'estimated_signal_power',signalPower,'observed_power',totalPower, ...
 'method',method,'samples_used',numel(z),'lag1',rho1, ...
 'white_noise_assumption_ok',whiteOK,'signal_resolved',signalPower>0, ...
 'noise_window',opts.noise_window);
end
