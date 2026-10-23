function p=uq_public_config(c,nominal_frequency)
% Whitelist: no true Q, wavelet peak, density, Vs or reflection amplitudes.
keys={'dt','nt','t','offset','vp','h','fref','nfft','lambda_x','lambda_fixed','lambda_grid','cov_floor'};
p=struct; for i=1:numel(keys), k=keys{i}; p.(k)=c.(k); end
p.nominal_frequency=nominal_frequency;
end
