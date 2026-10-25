function sigma=estimate_sigmas(D,t)
% Per-trace noise estimate used when the archived observations were written.
sigma=zeros(1,size(D,2));for j=1:size(D,2),ni=estimate_trace_snr(D(:,j),t);sigma(j)=ni.sigma;end
end
