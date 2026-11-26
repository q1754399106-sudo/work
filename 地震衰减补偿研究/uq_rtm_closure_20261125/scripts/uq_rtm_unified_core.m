function core = uq_rtm_unified_core()
% UQ_RTM_UNIFIED_CORE  Unified Recommended UQ-RTM Workflow: DSA-UQW1b-Anc
%
% Conflates four verified physical & statistical modules into one zero-extra-LU pipeline:
%   1) B = 300 shot-block polarity-flip blind noise prior (Theorem 3, eliminates 13-28x
%      underestimation under laterally correlated noise, band-summed +-5 Hz + spatial box);
%   2) Blind structure-tensor dip-aware active-array specular gating + Nyquist anti-aliasing
%      angle taper + dip-consistent two-way t* (DSA, eliminates ~36.5-67% out-of-array arcs);
%   3) Parametric Smooth Shallow-Anchor Spectral Equalization (PSAE) + relative dt* = max(t* - t*_anc, 0)
%      (eliminates (2*pi*f)^2 + unknown-wavelet blue tilt without injecting thin-bed notches);
%   4) Gamma=1 Gauss-Hermite Bayesian expected-risk operator (UQW1b, sigma=0.447) with
%      in-band velocity-consistent causal reference frequency f_ref = 25 Hz (default).
%
% Usage:
%   core = uq_rtm_unified_core();
%   cfg  = core.default_cfg('field_obc'); % or 'synth_10m'

core.default_cfg       = @default_cfg;
core.blind_dip         = @blind_dip;
core.build_dsa_full    = @build_dsa_full;
core.build_dsa_in      = @build_dsa_in;
core.make_block_signs  = @make_block_signs;
core.band_inv_snr      = @band_inv_snr;
core.fit_psae          = @fit_psae;
core.fit_env_anc       = @fit_env_anc;
core.tab_idx           = @tab_idx;
core.prep_image        = @prep_image;
core.steep_smile_frac  = @steep_smile_frac;
core.res_band          = @res_band;
core.bp_z              = @bp_z;
core.well_tie_full     = @well_tie_full;
core.load_well         = @load_well;
core.cos_taper         = @cos_taper;
core.smooth2d          = @smooth2d;
core.highpass_kz       = @highpass_kz;
core.score_e1          = @score_e1;
core.cc                = @cc;
end

function cfg = default_cfg(preset)
if nargin < 1, preset = 'field_obc'; end
cfg = struct();
cfg.B            = 300;
cfg.sig          = 0.447;
cfg.gamma        = 1.0;
cfg.fref         = 25.0;
cfg.dip_cap_deg  = 30.0;
cfg.snr_band_hz  = 5.0;
cfg.snr_cap      = 1e3;
cfg.w_floor      = 1e-4;
cfg.flip_seed    = 20261118;
switch preset
    case 'field_obc'
        cfg.anchor_mode   = 'psae';
        cfg.bg_sigma_m    = 40.0;
        cfg.z_mute_m      = 150.0;
        cfg.dip_sigma_m   = 60.0;
        cfg.arr_taper_m   = [-150.0, 250.0];
        cfg.theta_max_deg = 50.0;
        cfg.theta_tap_deg = 10.0;
        cfg.snr_box_m     = 400.0;
        cfg.anc_z_m       = [1000.0, 1600.0];
        cfg.anc_x_m       = [-2500.0, 1500.0];
        cfg.res_z_m       = [2600.0, 3200.0];
        cfg.res_x_m       = [-2500.0, 1500.0];
        cfg.ovb_z_m       = [1600.0, 2400.0];
    case 'synth_10m'
        cfg.anchor_mode   = 'none';
        cfg.bg_sigma_m    = 50.0;
        cfg.z_mute_m      = 160.0;
        cfg.dip_sigma_m   = 40.0;
        cfg.arr_taper_m   = [-80.0, 120.0];
        cfg.theta_max_deg = 75.0;
        cfg.theta_tap_deg = 10.0;
        cfg.snr_box_m     = 160.0;
        cfg.anc_z_m       = [190.0, 320.0];
        cfg.anc_x_m       = [240.0, 1360.0];
        cfg.res_z_m       = [560.0, 950.0];
        cfg.res_x_m       = [180.0, 1420.0];
    otherwise
        error('Unknown preset: %s', preset);
end
end

function [phi_eff, c_st, phi_raw] = blind_dip(R0_2d, dx, sig_m, cap_deg)
if nargin < 4, cap_deg = 30.0; end
sig_px = max(sig_m / dx, 1.0);
[gx, gz] = gradient(double(R0_2d), dx);
Jxx = smooth2d(gx.^2, sig_px);
Jzz = smooth2d(gz.^2, sig_px);
Jxz = smooth2d(gx.*gz, sig_px);
phi_raw = max(min(-0.5 * atan2(2*Jxz, Jzz - Jxx), deg2rad(cap_deg)), -deg2rad(cap_deg));
c_st    = sqrt((Jzz - Jxx).^2 + 4*Jxz.^2) ./ max(Jzz + Jxx, 1e-30);
phi_eff = c_st .* phi_raw;
end

function [w_dsa, TS_dsa] = build_dsa_full(xv, zv, z_src_vec, z_rec_val, xs, xr_lo, xr_hi, phi_v, F_full, zi_v, xi_v, nz, nx, cfg)
% Build DSA weight mask & dip-consistent two-way t* from full-grid eikonal tensor F_full (nz x nx x ns)
nin = numel(xv); ns = numel(xs);
[xs_sort, ord] = sort(xs(:));
w_dsa  = zeros(nin, ns, 'single');
TS_dsa = zeros(nin, ns, 'single');
dz_r   = max(zv - z_rec_val, 10.0);
for k = 1:ns
    dz_s = max(zv - z_src_vec(k), 10.0);
    th_s = atan2(xv - xs(k), dz_s);
    th_r = max(min(th_s + 2*phi_v, deg2rad(78)), -deg2rad(78));
    xmir = xv + dz_r .* tan(th_r);
    d_in = min(xmir - xr_lo(k), xr_hi(k) - xmir);
    t_arr = cos_taper(d_in, cfg.arr_taper_m(1), cfg.arr_taper_m(2));
    t_gs  = cos_taper(cfg.theta_max_deg - abs(rad2deg(th_s)), 0.0, cfg.theta_tap_deg);
    t_gr  = cos_taper(cfg.theta_max_deg - abs(rad2deg(th_r)), 0.0, cfg.theta_tap_deg);
    w_dsa(:, k) = single(max(t_arr .* t_gs .* t_gr, cfg.w_floor));

    q  = interp1(xs_sort, 1:ns, min(max(xmir, xs_sort(1)), xs_sort(end)));
    k0 = min(floor(q), ns - 1); ww = single(q - k0);
    a  = sub2ind([nz nx ns], zi_v, xi_v, ord(k0));
    b  = sub2ind([nz nx ns], zi_v, xi_v, ord(k0 + 1));
    s0 = sub2ind([nz nx ns], zi_v, xi_v, k * ones(nin, 1));
    TS_dsa(:, k) = F_full(s0) + (1 - ww).*F_full(a) + ww.*F_full(b);
end
end

function [w_dsa, TS_dsa] = build_dsa_in(xv, zv, z_src_val, xs, xr_lo, xr_hi, phi_v, Fts_in, cfg)
% Build DSA weight mask & dip-consistent two-way t* from interior eikonal matrix Fts_in (nin x ns)
nin = numel(xv); ns = numel(xs);
[xs_sort, ord] = sort(xs(:));
w_dsa  = zeros(nin, ns);
TS_dsa = zeros(nin, ns);
dz     = max(zv - z_src_val, 10.0);
row_idx = (1:nin)';
for k = 1:ns
    th_s = atan2(xv - xs(k), dz);
    th_r = max(min(th_s + 2*phi_v, deg2rad(78)), -deg2rad(78));
    xmir = xv + dz .* tan(th_r);
    d_in = min(xmir - xr_lo(k), xr_hi(k) - xmir);
    t_arr = cos_taper(d_in, cfg.arr_taper_m(1), cfg.arr_taper_m(2));
    t_gs  = cos_taper(cfg.theta_max_deg - abs(rad2deg(th_s)), 0.0, cfg.theta_tap_deg);
    t_gr  = cos_taper(cfg.theta_max_deg - abs(rad2deg(th_r)), 0.0, cfg.theta_tap_deg);
    w_dsa(:, k) = max(t_arr .* t_gs .* t_gr, cfg.w_floor);

    q  = interp1(xs_sort, 1:ns, min(max(xmir, xs_sort(1)), xs_sort(end)));
    k0 = min(floor(q), ns - 1); ww = q - k0;
    i0 = sub2ind([nin ns], row_idx, ord(k0));
    i1 = sub2ind([nin ns], row_idx, ord(k0 + 1));
    TS_dsa(:, k) = Fts_in(:, k) + (1 - ww).*Fts_in(i0) + ww.*Fts_in(i1);
end
end

function [EPS_tr, shot_sgn] = make_block_signs(shot_of, rec_ix, ns, B, flip_seed)
% Theorem 3 block-polarity signs: B=1 trace-wise, B=300 shot-wise (exact est_block / fld_run_h rule)
ntrace = numel(shot_of);
if B == 1
    rng(flip_seed, 'twister');
    EPS_tr = single(2*(rand(ntrace, 1) > 0.5) - 1);
    shot_sgn = [];
elseif B >= 300
    rng(flip_seed + B, 'twister');
    shot_sgn = single(2*(rand(ns, 1) > 0.5) - 1);
    EPS_tr = shot_sgn(shot_of(:));
else
    nper = accumarray(shot_of(:), 1, [ns 1]);
    nb = ceil(max(nper) / B);
    rng(flip_seed + B, 'twister');
    s = single(2*(rand(nb, ns) > 0.5) - 1);
    blk = zeros(ntrace, 1);
    for k = 1:ns
        t = find(shot_of == k);
        [~, o] = sort(rec_ix(t));
        r = zeros(numel(t), 1); r(o) = 1:numel(t);
        blk(t) = ceil(r / B);
    end
    EPS_tr = s(sub2ind([nb ns], blk, shot_of(:)));
    shot_sgn = [];
end
end

function iv = band_inv_snr(I_mat, If_mat, freqs, nzi, nxi, dx, band_hz, box_m, cap)
% Exact band-summed (+-band_hz) & spatial box-smoothed (box_m) blind inverse SNR estimator
if nargin < 9, cap = 1e3; end
[nin, nf] = size(I_mat);
bw = max(1, round(box_m / dx));
iv = zeros(nin, nf, 'single');
for j = 1:nf
    bnd = abs(freqs - freqs(j)) <= band_hz + 1e-9;
    s_tot = sum(double(I_mat(:, bnd)), 2);
    s_flp = sum(double(If_mat(:, bnd)), 2);
    p  = box_pow(s_tot, nzi, nxi, bw);
    pf = box_pow(s_flp, nzi, nxi, bw);
    iv(:, j) = single(min(pf ./ max(p - pf, realmin('single')), cap));
end
end

function y = box_pow(v, nzi, nxi, bw)
y = reshape(v, nzi, nxi).^2;
y = movmean(movmean(y, bw, 1), bw, 2);
y = y(:);
end

function [W_eq, beta_anc, eff_tstar_ms] = fit_psae(freqs, S_anc, f_fit_range, f_mid)
% Parametric Smooth Shallow-Slope Anchoring (PSAE): notch-free single-parameter equalization
fm = (freqs >= f_fit_range(1)) & (freqs <= f_fit_range(2));
p = polyfit(freqs(fm), log(max(S_anc(fm), 1e-30)), 1);
beta_anc = p(1);
eff_tstar_ms = -beta_anc / pi * 1e3;
W_eq = exp(-beta_anc * (freqs - f_mid));
end

function W_env = fit_env_anc(freqs, S_anc, W_taper, smooth_hz, f_mid)
% Broad-window smooth shallow envelope normalization (ablation comparison vs 1-param PSAE)
df = median(diff(freqs));
nw = max(3, 2*floor(round(smooth_hz / df) / 2) + 1);
S_sm = movmean(double(S_anc), nw);
floor_val = 0.15 * max(S_sm);
W_env = double(W_taper(:)') ./ max(S_sm, floor_val);
[~, im] = min(abs(freqs - f_mid));
W_env = W_env / max(W_env(im), 1e-30);
end

function [idx, ww] = tab_idx(ts_m, dxg, ntab_m2)
u = ts_m / dxg;
i0 = min(max(floor(u), 0), ntab_m2);
ww = u - i0;
idx = int32(i0) + 1;
end

function y = prep_image(v, nzi, nxi, z_in, dx, bg_sigma_m, z_mute_m)
% Exact fld_eval / fld_eval_h vertical Gaussian background removal + shallow water mute
y = reshape(double(v), nzi, nxi);
s = bg_sigma_m / dx;
m = ceil(3 * s);
h = exp(-0.5 * ((-m:m) / s).^2); h = h / sum(h);
yp = [repmat(y(1,:), m, 1); y; repmat(y(end,:), m, 1)];
y = y - conv2(h', 1, yp, 'valid');
y(z_in < z_mute_m, :) = 0;
end

function fr = steep_smile_frac(A_mid, dx, deg_cut, k_lo_cyckm, k_hi_cyckm)
if nargin < 3, deg_cut = 25.0; end
if nargin < 4, k_lo_cyckm = 3.0; end
if nargin < 5, k_hi_cyckm = 25.0; end
[nz, nx] = size(A_mid);
win = hann(nz) * hann(nx)';
F2 = abs(fftshift(fft2(A_mid .* win))).^2;
kz = ((-floor(nz/2)):(ceil(nz/2)-1))' / (nz * dx);
kx = ((-floor(nx/2)):(ceil(nx/2)-1))  / (nx * dx);
[KZ, KX] = ndgrid(kz, kx);
k_mag = sqrt(KX.^2 + KZ.^2);
in_band = (k_mag >= k_lo_cyckm / 1000) & (k_mag <= k_hi_cyckm / 1000);
steep = in_band & (abs(KX) > tand(deg_cut) * max(abs(KZ), 1e-9));
fr = sum(F2(steep), 'all') / max(sum(F2(in_band), 'all'), 1e-30);
end

function [bw10, cen] = res_band(A_res, dx, nfft)
if nargin < 3, nfft = 256; end
kz = (0:nfft/2) / (nfft * dx) * 1000;
s = mean(abs(fft(A_res .* hann(size(A_res, 1)), nfft, 1)), 2);
s = s(1:nfft/2+1)';
[mx, i] = max(s); t = mx * 10^(-10/20);
l = find(s(1:i) < t, 1, 'last'); if isempty(l), l = 1; end
h = i - 1 + find(s(i:end) < t, 1); if isempty(h), h = numel(s); end
bw10 = [kz(l), kz(h)];
cen = sum(kz .* (s.^2)) / max(sum(s.^2), 1e-30);
end

function y = bp_z(r, dz, lo, hi)
n = 2^nextpow2(4 * numel(r));
R = fft(r, n);
k = abs([0:n/2, -n/2+1:-1]') / (n * dz) * 1000;
y = real(ifft(R .* (k >= lo & k <= hi)));
y = y(1:numel(r));
end

function w = load_well(fn, z, x)
fid = fopen(fn); a = fread(fid, Inf, 'single'); fclose(fid);
a = reshape(a, 4, [])';
w.zw = double(a(:, 1)); w.sw = double(a(:, 2)); w.rw = double(a(:, 4));
w.dzw = median(diff(w.zw));
[~, w.iz] = min(abs(z(:) - w.zw'), [], 1);
w.ix = arrayfun(@(xv) find(abs(x - xv) == min(abs(x - xv)), 1), w.sw);
end

function [best_cc, best_lag_m, zero_cc, tr_b] = well_tie_full(A, w, rb, dx_grid, lo, hi, half_avg_m, max_lag_m)
% Exact fld_eval / fld_eval_h well tie: +-half_avg_m lateral avg, +-max_lag_m depth shift, plus L=0!
if nargin < 7, half_avg_m = 20.0; end
if nargin < 8, max_lag_m  = 30.0; end
na = round(half_avg_m / dx_grid);
nw = numel(w.zw); tr = zeros(nw, 1);
for i = 1:nw
    tr(i) = mean(A(w.iz(i), max(w.ix(i) - na, 1):min(w.ix(i) + na, size(A, 2))));
end
tr_b = bp_z(tr, w.dzw, lo, hi);
zero_cc = cc(tr_b, rb);
ns_ = round(max_lag_m / w.dzw);
best_cc = -Inf; best_lag_m = 0;
for L = -ns_:ns_
    c = cc(circshift(tr_b, L), rb);
    if c > best_cc
        best_cc = c;
        best_lag_m = L * w.dzw;
    end
end
end

function t = cos_taper(x, x0, x1)
u = min(max((x - x0) / max(x1 - x0, 1e-12), 0), 1);
t = 0.5 - 0.5 * cos(pi * u);
end

function B = smooth2d(A, sig)
r = max(1, ceil(3 * sig));
g = exp(-0.5 * ((-r:r) / sig).^2); g = g / sum(g);
Ap = [repmat(A(1,:), r, 1); A; repmat(A(end,:), r, 1)];
B1 = conv2(g', 1, Ap, 'valid');
B1p = [repmat(B1(:,1), 1, r), B1, repmat(B1(:,end), 1, r)];
B = conv2(1, g, B1p, 'valid');
end

function yh = highpass_kz(y, dx, kz_cut)
n = 2^nextpow2(2 * size(y, 1));
Y = fft(y, n, 1);
kz = abs([0:n/2, -n/2+1:-1]') / (n * dx) * 1000;
g = 0.5 - 0.5 * cos(pi * min(max((kz - (kz_cut - 2)) / 4, 0), 1));
yh = real(ifft(Y .* g, n, 1));
yh = yh(1:size(y, 1), :);
end

function s = score_e1(a, r, iz, ix)
rr = r(iz, ix);
c0 = cc(a(iz, ix), rr);
cl = c0;
nz = size(a, 1);
for L = [-2 -1 1 2]
    cl = max(cl, cc(a(min(max(iz + L, 1), nz), ix), rr));
end
s = 0.5 * (c0 + cl);
end

function c = cc(a, b)
a = a(:) - mean(a(:));
b = b(:) - mean(b(:));
c = (a' * b) / max(sqrt((a' * a) * (b' * b)), 1e-30);
end
