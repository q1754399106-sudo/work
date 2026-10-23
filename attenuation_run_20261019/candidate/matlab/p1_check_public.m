function p1_check_public(p)
allowed={'dt','nt','t','offset','vp','h','fref','nfft','lambda_x','lambda_fixed','lambda_grid','cov_floor','nominal_frequency'};
assert(isempty(setdiff(fieldnames(p),allowed)),'P1 requires the strict public-input whitelist');
assert(all(isfield(p,allowed)),'Missing required public field');
assert(numel(p.t)==p.nt && all(isfinite(p.t)) && p.dt>0 && p.nt>=16);
end
