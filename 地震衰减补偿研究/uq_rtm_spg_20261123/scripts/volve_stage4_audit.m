function volve_stage4_audit()
% Part B: Strictly read-only quantitative audit of Volve Stage 4 (uq_rtm_20261117/cache/fld).
% Never touches fld45. Archives results/volve_stage4_audit.json.
maxNumCompThreads(2);
here = fileparts(mfilename('fullpath'));
root = fileparts(here);
out_dir = fullfile(root, 'results');
if ~isfolder(out_dir), mkdir(out_dir); end
addpath('D:\ChatGPT\研究\地震衰减补偿研究\uq_rtm_20261117\scripts');

P = fld_setup();
S1 = load(fullfile(P.cache, 'state1.mat'));
T  = load(fullfile(P.cache, 'ts.mat'));
Ia = load(fullfile(P.cache, 'images_amp.mat'));
z = P.z(P.in{1}); x = P.x(P.in{2});
[Zg, Xg] = ndgrid(z, x);

% 1. Out-of-array specular mirror fraction in reservoir window [2600-3200 m, -2500..1500 m]
rx = P.x(P.rec_ix);
xr_lo = arrayfun(@(s) min(rx(P.shot_of==s)), 1:P.nshot);
xr_hi = arrayfun(@(s) max(rx(P.shot_of==s)), 1:P.nshot);
wres = (Zg >= 2600 & Zg <= 3200) & (Xg >= -2500 & Xg <= 1500);
xs_res = Xg(wres);
out_cnt = 0; tot_cnt = 0;
for s = 1:P.nshot
    xmir = 2 * xs_res - P.src_x(s);
    out_cnt = out_cnt + sum(xmir < xr_lo(s) | xmir > xr_hi(s));
    tot_cnt = tot_cnt + numel(xmir);
end
out_frac = out_cnt / tot_cnt;

% 2. Depth-resolved R0 spectral ratio & empirical attenuation dt*_emp vs model dt*
wx = (x >= -2500 & x <= 1500);
ts_med = reshape(median(double(T.TS), 2), P.nzi, P.nxi);
zwins = [500 1000; 1000 1600; 1600 2200; 2200 2600; 2600 3200];
fb = [4 8; 8 12; 12 16; 16 20; 20 25];
depth_spec = struct();
for iw = 1:size(zwins, 1)
    wz = (z >= zwins(iw,1) & z <= zwins(iw,2));
    vals = zeros(1, size(fb,1));
    for ib = 1:size(fb,1)
        idx = P.f >= fb(ib,1) & P.f <= fb(ib,2);
        Ak = reshape(double(sum(S1.I(:, idx), 2)), P.nzi, P.nxi);
        vals(ib) = rms(Ak(wz, wx), 'all');
    end
    key = sprintf('z_%d_%d', zwins(iw,1), zwins(iw,2));
    depth_spec.(key) = struct('rms_bands', vals, 'ratio_20_25_to_8_12', vals(5)/vals(2), ...
        'med_tstar_ms', median(ts_med(wz, wx), 'all') * 1e3);
end

wz_s = (z >= 1000 & z <= 1600); wz_d = (z >= 2600 & z <= 3200);
spec_s = zeros(numel(P.f), 1); spec_d = zeros(numel(P.f), 1);
for jf = 1:numel(P.f)
    Aj = reshape(double(S1.I(:, jf)), P.nzi, P.nxi);
    spec_s(jf) = rms(Aj(wz_s, wx), 'all');
    spec_d(jf) = rms(Aj(wz_d, wx), 'all');
end
fmask = P.f >= 6 & P.f <= 24;
pfit = polyfit(P.f(fmask), log(spec_d(fmask) ./ spec_s(fmask)), 1);
dtstar_emp_ms = (-pfit(1) / pi) * 1e3;
dtstar_mod_ms = (median(ts_med(wz_d, wx), 'all') - median(ts_med(wz_s, wx), 'all')) * 1e3;
Q_eff_emp = 120.0 / (dtstar_emp_ms / dtstar_mod_ms);

% 3. Well F-5 tie: 1D post-stack UQ vs Per-shot UQW1_amp & oracle alpha scan
fid = fopen(fullfile(P.cache, 'well_F-5.bin')); a = fread(fid, Inf, 'single'); fclose(fid);
a = reshape(a, 4, [])'; zw = a(:,1); sw = a(:,2); rw = double(a(:,4)); nw = numel(zw);
[~, iz] = min(abs(z(:) - zw'), [], 1);
ix = arrayfun(@(xv) find(abs(x-xv)==min(abs(x-xv)), 1), sw);
ker = exp(-0.5*((-12:12)/4).^2); ker = ker/sum(ker);
hp = @(A) A - conv2(ker', 1, [repmat(A(1,:),12,1); A; repmat(A(end,:),12,1)], 'valid');
wtr = @(A) arrayfun(@(i) mean(A(iz(i), max(ix(i)-2,1):min(ix(i)+2,numel(x)))), 1:nw)';
rb = bp_z(rw, P.dx, 3.52, 12.89);

ts_med_vec = median(double(T.TS), 2); s70 = gain_s(70);
m_idx = [1 3 5 7 9]; m_vals = Ia.meta.m(m_idx);
cc_1d_uq = zeros(1, numel(m_idx)); cc_ps_uq = zeros(1, numel(m_idx));
cc_1d_w  = zeros(1, numel(m_idx)); cc_ps_w  = zeros(1, numel(m_idx));
for k = 1:numel(m_idx)
    im = m_idx(k); m = m_vals(k);
    I1w = zeros(size(ts_med_vec)); I1u = I1w;
    xg = linspace(0, 1.01*max(ts_med_vec)/m, 601);
    for jf = 1:numel(P.f)
        f = P.f(jf);
        T0 = uq_tab(f, xg, 0, f); T1 = uq_tab(f, xg, P.sig, f);
        A0 = interp1(xg, real(T0.A), ts_med_vec/m); C0 = interp1(xg, real(T0.C), ts_med_vec/m);
        A1 = interp1(xg, real(T1.A), ts_med_vec/m); C1 = interp1(xg, real(T1.C), ts_med_vec/m);
        I1w = I1w + ((A0 + s70*C0)./(1 + s70*C0)) .* double(S1.I(:, jf));
        I1u = I1u + ((A1 + s70*C1)./(1 + s70*C1)) .* double(S1.I(:, jf));
    end
    cc_1d_w(k)  = best_cc(bp_z(wtr(hp(reshape(I1w, P.nzi, P.nxi))), P.dx, 3.52, 12.89), rb);
    cc_ps_w(k)  = best_cc(bp_z(wtr(hp(reshape(double(Ia.W1(:,im)), P.nzi, P.nxi))), P.dx, 3.52, 12.89), rb);
    cc_1d_uq(k) = best_cc(bp_z(wtr(hp(reshape(I1u, P.nzi, P.nxi))), P.dx, 3.52, 12.89), rb);
    cc_ps_uq(k) = best_cc(bp_z(wtr(hp(reshape(double(Ia.UQW1(:,im)), P.nzi, P.nxi))), P.dx, 3.52, 12.89), rb);
end
cc_r0 = best_cc(bp_z(wtr(hp(reshape(double(sum(S1.I, 2)), P.nzi, P.nxi))), P.dx, 3.52, 12.89), rb);

alphas = [-10, -5, 0, 5, 10, 23];
cc_alpha = zeros(size(alphas));
for k = 1:numel(alphas)
    Ial = zeros(size(ts_med_vec));
    for jf = 1:numel(P.f)
        Ial = Ial + exp(pi * P.f(jf) * (alphas(k)*1e-3)) * double(S1.I(:, jf));
    end
    cc_alpha(k) = best_cc(bp_z(wtr(hp(reshape(Ial, P.nzi, P.nxi))), P.dx, 3.52, 12.89), rb);
end

out = struct( ...
    'out_of_array_reservoir_frac', out_frac, ...
    'dtstar_emp_ms_1000_to_3200', dtstar_emp_ms, ...
    'dtstar_mod_ms_Q120', dtstar_mod_ms, ...
    'Q_eff_emp', Q_eff_emp, ...
    'depth_spec', depth_spec, ...
    'cc_R0', cc_r0, ...
    'm', m_vals, 'cc_1D_UQ', cc_1d_uq, 'cc_PerShot_UQW1', cc_ps_uq, ...
    'cc_1D_W1', cc_1d_w, 'cc_PerShot_W1', cc_ps_w, ...
    'max_abs_diff_1D_vs_PerShot_UQ', max(abs(cc_1d_uq - cc_ps_uq)), ...
    'oracle_alpha_ms', alphas, 'oracle_alpha_cc', cc_alpha);

fn = fullfile(out_dir, 'volve_stage4_audit.json');
fid = fopen(fn, 'w'); fwrite(fid, jsonencode(out, 'PrettyPrint', true), 'char'); fclose(fid);
fprintf('Volve Stage 4 audit saved: out_of_array=%.1f%%, dt*_emp=%.2f ms (Q_eff=%.1f), max|1D-PS UQ|=%.4f\n', ...
    100*out_frac, dtstar_emp_ms, Q_eff_emp, out.max_abs_diff_1D_vs_PerShot_UQ);
end

function y = bp_z(r, dx, lo, hi)
n = 2^nextpow2(4 * numel(r)); R = fft(r, n);
k = abs([0:n/2, -n/2+1:-1]') / (n * dx) * 1000;
y = real(ifft(R .* (k >= lo & k <= hi))); y = y(1:numel(r));
end

function best = best_cc(tr, rb)
best = -Inf;
for L = -3:3
    a = circshift(tr, L); a = a - mean(a); b = rb - mean(rb);
    c = (a' * b) / max(sqrt((a' * a) * (b' * b)), 1e-30);
    if c > best, best = c; end
end
end
