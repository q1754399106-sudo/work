function diag_anomalies_wp1()
% WP1 Quantitative Root-Cause Diagnostics (D1-D4)
% Strictly read-only on D:\tmp_uqrtm_unified_20261124\ and uq_rtm_unified_20261124\results\.
% Outputs:
%   results/wp1_diagnostics.json
%   figures/fig_wp1_diagnostics.png

maxNumCompThreads(8);
t0 = tic;

here = fileparts(mfilename('fullpath'));
root = fileparts(here);
res_dir = fullfile(root, 'results');
fig_dir = fullfile(root, 'figures');
if ~isfolder(res_dir), mkdir(res_dir); end
if ~isfolder(fig_dir), mkdir(fig_dir); end

wt_root = fileparts(root);
addpath(fullfile(wt_root, 'uq_rtm_20261117', 'scripts'));
core = uq_rtm_unified_core();
cfg_f = core.default_cfg('field_obc');
cfg_s = core.default_cfg('synth_10m');

tmp_old  = 'D:\tmp_uqrtm_unified_20261124';
well_dir = 'D:\ChatGPT\研究\.claude\worktrees\study-research-docs-plan-c4d4c6\地震衰减补偿研究\uq_rtm_20261117\cache\fld';
fld45_cache = 'D:\ChatGPT\研究\.claude\worktrees\study-research-docs-plan-c4d4c6\地震衰减补偿研究\uq_rtm_20261117\cache\fld45';
fld25_cache = 'D:\ChatGPT\研究\地震衰减补偿研究\uq_rtm_20261117\cache\fld';

fprintf('=== WP1 Diagnostics (D1-D4) ===\n');

% Load Volve 4-45 Hz and 4-25 Hz snapshots (read-only)
S45 = load(fullfile(tmp_old, 'volve45_final_imgs.mat'));
S25 = load(fullfile(tmp_old, 'volve25_final_imgs.mat'));
z45 = S45.z_in(:); x45 = S45.x_in(:)'; dx45 = median(diff(z45)); nzi45 = numel(z45); nxi45 = numel(x45);
z25 = S25.z_in(:); x25 = S25.x_in(:)'; dx25 = median(diff(z25)); nzi25 = numel(z25); nxi25 = numel(x25);

w5_45  = core.load_well(fullfile(well_dir, 'well_F-5.bin'),  z45, x45);
w14_45 = core.load_well(fullfile(well_dir, 'well_F-14.bin'), z45, x45);
w5_25  = core.load_well(fullfile(fld25_cache, 'well_F-5.bin'),  z25, x25);
w14_25 = core.load_well(fullfile(fld25_cache, 'well_F-14.bin'), z25, x25);

prep45 = @(v) core.prep_image(v, nzi45, nxi45, z45, dx45, cfg_f.bg_sigma_m, cfg_f.z_mute_m);
prep25 = @(v) core.prep_image(v, nzi25, nxi25, z25, dx25, cfg_f.bg_sigma_m, cfg_f.z_mute_m);

diag = struct();

%% ========================================================================
%  D1: Volve 4-45 Hz Blind Well F-14 Negative Correlation Root Cause
%  ========================================================================
fprintf('\n--- D1: Volve 4-45 Hz Blind Well F-14 Diagnostics ---\n');
kfix = S45.res.kfix;         % [4.17, 20.86] cyc/km (8-40 Hz)
kr0  = S45.res.kr0;          % [3.75, 20.00] cyc/km
k25  = [3.515625, 12.890625]; % Exact 4-25 Hz R0-band
khf  = [12.890625, 20.86333]; % High-k part (25-40 Hz)

lags_d1 = -60:2.5:60;
d1_methods = {'RAW_R0', 'DSA_R0', 'DSA_R0_Eq', 'RAW_UQW1b_Unanc_100', 'DSA_UQW1b_PAnc_25', 'RAW_UQW1b_PAnc_25'};
d1 = struct('lags_m', lags_d1, 'kfix', kfix, 'kr0', kr0, 'k25', k25, 'khf', khf, ...
            'F14_zw_range', [min(w14_45.zw), max(w14_45.zw)], ...
            'F14_sw_range', [min(w14_45.sw), max(w14_45.sw)], ...
            'F5_zw_range',  [min(w5_45.zw),  max(w5_45.zw)], ...
            'F5_sw_range',  [min(w5_45.sw),  max(w5_45.sw)], ...
            'F14_n_samples', numel(w14_45.zw), 'F14_dzw', w14_45.dzw, ...
            'methods', struct());

rb14_fix = core.bp_z(w14_45.rw, w14_45.dzw, kfix(1), kfix(2));
rb14_r0  = core.bp_z(w14_45.rw, w14_45.dzw, kr0(1),  kr0(2));
rb14_k25 = core.bp_z(w14_45.rw, w14_45.dzw, k25(1),  k25(2));
rb14_khf = core.bp_z(w14_45.rw, w14_45.dzw, khf(1),  khf(2));

% Wavenumber spectra of F-14 well reflectivity and seismic traces
nfft_w = 512;
kw_axis = (0:nfft_w/2) / (nfft_w * w14_45.dzw) * 1000; % cyc/km
Rw14_spec = abs(fft(w14_45.rw .* hann(numel(w14_45.rw)), nfft_w));
d1.kw_axis = kw_axis;
d1.F14_well_spec = Rw14_spec(1:nfft_w/2+1)';

for i = 1:numel(d1_methods)
    mk = d1_methods{i};
    M = S45.Img_s.(mk);
    cidx = min(3, size(M, 2)); % m=1.0
    A = prep45(M(:, cidx));
    tr14 = extract_well_trace(A, w14_45, dx45, 20.0);

    tr_fix = core.bp_z(tr14, w14_45.dzw, kfix(1), kfix(2));
    tr_r0  = core.bp_z(tr14, w14_45.dzw, kr0(1),  kr0(2));
    tr_k25 = core.bp_z(tr14, w14_45.dzw, k25(1),  k25(2));
    tr_khf = core.bp_z(tr14, w14_45.dzw, khf(1),  khf(2));

    cc_fix_curve = lag_curve(tr_fix, rb14_fix, w14_45.dzw, lags_d1, core);
    cc_k25_curve = lag_curve(tr_k25, rb14_k25, w14_45.dzw, lags_d1, core);

    [b_fix, l_fix, l0_fix] = core.well_tie_full(A, w14_45, rb14_fix, dx45, kfix(1), kfix(2), 20.0, 30.0);
    [b_r0,  l_r0,  l0_r0 ] = core.well_tie_full(A, w14_45, rb14_r0,  dx45, kr0(1),  kr0(2),  20.0, 30.0);
    [b_k25, l_k25, l0_k25] = core.well_tie_full(A, w14_45, rb14_k25, dx45, k25(1),  k25(2),  20.0, 30.0);
    [b_khf, l_khf, l0_khf] = core.well_tie_full(A, w14_45, rb14_khf, dx45, khf(1),  khf(2),  20.0, 30.0);
    [b_fix60, idx60] = max(cc_fix_curve); l_fix60 = lags_d1(idx60);

    Tspec = abs(fft(tr14 .* hann(numel(tr14)), nfft_w));
    d1.methods.(mk) = struct( ...
        'F14_kfix_30m', b_fix, 'F14_kfix_lag30', l_fix, 'F14_kfix_L0', l0_fix, ...
        'F14_kfix_60m', b_fix60, 'F14_kfix_lag60', l_fix60, ...
        'F14_kr0_30m', b_r0, 'F14_kr0_lag30', l_r0, ...
        'F14_k25_30m', b_k25, 'F14_k25_lag30', l_k25, 'F14_k25_L0', l0_k25, ...
        'F14_khf_30m', b_khf, 'F14_khf_lag30', l_khf, ...
        'cc_kfix_curve', cc_fix_curve, 'cc_k25_curve', cc_k25_curve, ...
        'trace_spec', Tspec(1:nfft_w/2+1)');
    fprintf('  %-22s | kfix(+-30m)=%+.4f(@%+4.0fm,L0=%+.3f) | kfix(+-60m)=%+.4f(@%+4.0fm) | k25=%+.4f(@%+4.0fm) | khf=%+.4f\n', ...
        mk, b_fix, l_fix, l0_fix, b_fix60, l_fix60, b_k25, l_k25, b_khf);
end

% Compare with 4-25 Hz F-14 on exact same k25 band
A25_r0  = prep25(S25.Img_s.RAW_R0);
A25_dsa = prep25(S25.Img_s.DSA_UQW1b_PAnc_25(:, 3));
rb14_25 = core.bp_z(w14_25.rw, w14_25.dzw, k25(1), k25(2));
[d1.volve25_F14_RAW_R0, d1.volve25_F14_RAW_R0_lag] = core.well_tie_full(A25_r0, w14_25, rb14_25, dx25, k25(1), k25(2), 20.0, 30.0);
[d1.volve25_F14_DSA_UQW1b, d1.volve25_F14_DSA_UQW1b_lag] = core.well_tie_full(A25_dsa, w14_25, rb14_25, dx25, k25(1), k25(2), 20.0, 30.0);
fprintf('  [Cross-check 4-25 Hz F-14 in k25] RAW_R0=%+.4f(@%+4.0fm) | DSA_UQW1b_PAnc_25=%+.4f(@%+4.0fm)\n', ...
    d1.volve25_F14_RAW_R0, d1.volve25_F14_RAW_R0_lag, d1.volve25_F14_DSA_UQW1b, d1.volve25_F14_DSA_UQW1b_lag);
diag.D1_F14 = d1;

%% ========================================================================
%  D2: Volve 4-45 Hz RAW_UQW1b_PAnc_25 (0.5791) > DSA_UQW1b_PAnc_25 (0.5597)
%  ========================================================================
fprintf('\n--- D2: Volve 4-45 Hz RAW_UQW1b_PAnc_25 vs DSA_UQW1b_PAnc_25 ---\n');
beta_raw = S45.res.beta_anc_raw;
beta_dsa = S45.res.beta_anc_dsa;
vr = S45.res.v_r;
rb5_fix = core.bp_z(w5_45.rw, w5_45.dzw, kfix(1), kfix(2));

A_raw_p = prep45(S45.Img_s.RAW_UQW1b_PAnc_25(:, 3));
A_dsa_p = prep45(S45.Img_s.DSA_UQW1b_PAnc_25(:, 3));
tr_raw5 = extract_well_trace(A_raw_p, w5_45, dx45, 20.0);
tr_dsa5 = extract_well_trace(A_dsa_p, w5_45, dx45, 20.0);

% Swap PSAE slopes via kz-domain equivalent filter exp(-(beta_target - beta_orig)*(kz*vr/2 - 25))
tr_dsa_with_raw_slope = apply_kz_slope(tr_dsa5, w5_45.dzw, vr, -(beta_raw - beta_dsa));
tr_raw_with_dsa_slope = apply_kz_slope(tr_raw5, w5_45.dzw, vr, -(beta_dsa - beta_raw));

[cc_raw_orig, lag_raw_orig] = max_lag_cc(core.bp_z(tr_raw5, w5_45.dzw, kfix(1), kfix(2)), rb5_fix, w5_45.dzw, 30.0, core);
[cc_dsa_orig, lag_dsa_orig] = max_lag_cc(core.bp_z(tr_dsa5, w5_45.dzw, kfix(1), kfix(2)), rb5_fix, w5_45.dzw, 30.0, core);
[cc_dsa_swap, ~] = max_lag_cc(core.bp_z(tr_dsa_with_raw_slope, w5_45.dzw, kfix(1), kfix(2)), rb5_fix, w5_45.dzw, 30.0, core);
[cc_raw_swap, ~] = max_lag_cc(core.bp_z(tr_raw_with_dsa_slope, w5_45.dzw, kfix(1), kfix(2)), rb5_fix, w5_45.dzw, 30.0, core);

% Sub-band F-5 correlation breakdown: low [4.17-10], mid [10-15], high [15-20.86] cyc/km
sub_bands = {[4.17, 10.0], [10.0, 15.0], [15.0, 20.86]};
sub_names = {'low_8_19Hz', 'mid_19_29Hz', 'high_29_40Hz'};
sub_raw = zeros(1, 3); sub_dsa = zeros(1, 3);
for ib = 1:3
    b_ = sub_bands{ib};
    rb_s = core.bp_z(w5_45.rw, w5_45.dzw, b_(1), b_(2));
    [sub_raw(ib), ~] = max_lag_cc(core.bp_z(tr_raw5, w5_45.dzw, b_(1), b_(2)), rb_s, w5_45.dzw, 30.0, core);
    [sub_dsa(ib), ~] = max_lag_cc(core.bp_z(tr_dsa5, w5_45.dzw, b_(1), b_(2)), rb_s, w5_45.dzw, 30.0, core);
end

% Check F-5 shot fold & incidence angles suppressed by w_dsa from fld45_dsa_geom.mat
G45 = load(fullfile(tmp_old, 'fld45_dsa_geom.mat'), 'w_dsa', 'phi_eff', 'c_st');
M45 = jsondecode(fileread(fullfile(fld45_cache, 'meta.json')));
fid = fopen(fullfile(fld45_cache, 'src_ix.bin')); src_ix = double(fread(fid, Inf, '*int32')); fclose(fid);
fid = fopen(fullfile(fld45_cache, 'src_iz.bin')); src_iz = double(fread(fid, Inf, '*int32')); fclose(fid);
x_full45 = M45.xmin + ((1:M45.nx) - M45.npml - 1) * M45.dx;
z_full45 = ((1:M45.nz) - M45.npml - 1) * M45.dx;
xs45 = x_full45(src_ix); zs45 = z_full45(src_iz);

% Locate F-5 reservoir cells [2600-3200 m]
f5_x_med = median(w5_45.sw);
[~, ix_f5] = min(abs(x45 - f5_x_med));
iz_res45 = find(z45 >= 2600 & z45 <= 3200);
lin_f5 = sub2ind([nzi45, nxi45], iz_res45, repmat(ix_f5, numel(iz_res45), 1));
w_f5 = double(G45.w_dsa(lin_f5, :)); % nz_res x ns
fold_f5_med = median(sum(w_f5, 2));
fold_f5_active = median(sum(w_f5 > 0.5, 2));
% Angle distribution of kept vs cut shots at z=2900m, x=f5_x_med
[~, iz_2900] = min(abs(z45 - 2900));
lin_2900 = sub2ind([nzi45, nxi45], iz_2900, ix_f5);
w_2900 = double(G45.w_dsa(lin_2900, :));
th_s_2900 = abs(rad2deg(atan2(f5_x_med - xs45(:)', max(2900 - zs45(:)', 10))));
kept_idx = w_2900 >= 0.5; cut_idx = w_2900 < 0.5;

% Also check F-14 fold at z=2900m
f14_x_med = median(w14_45.sw);
[~, ix_f14] = min(abs(x45 - f14_x_med));
lin_f14_2900 = sub2ind([nzi45, nxi45], iz_2900, ix_f14);
w_f14_2900 = double(G45.w_dsa(lin_f14_2900, :));

d2 = struct( ...
    'beta_anc_raw', beta_raw, 'beta_anc_dsa', beta_dsa, 'delta_beta', beta_raw - beta_dsa, ...
    'RAW_Res_Cen', S45.res.methods.RAW_UQW1b_PAnc_25.Res_Cen(3), ...
    'DSA_Res_Cen', S45.res.methods.DSA_UQW1b_PAnc_25.Res_Cen(3), ...
    'F5_RAW_orig', cc_raw_orig, 'F5_DSA_orig', cc_dsa_orig, 'diff_orig', cc_raw_orig - cc_dsa_orig, ...
    'F5_DSA_with_RAW_slope', cc_dsa_swap, 'F5_RAW_with_DSA_slope', cc_raw_swap, ...
    'slope_swap_shift_on_DSA', cc_dsa_swap - cc_dsa_orig, ...
    'slope_swap_shift_on_RAW', cc_raw_swap - cc_raw_orig, ...
    'subband_names', {sub_names}, 'subband_RAW', sub_raw, 'subband_DSA', sub_dsa, ...
    'F5_x_m', f5_x_med, 'F5_fold_sum_med', fold_f5_med, 'F5_fold_gt05_med', fold_f5_active, ...
    'F5_n_shots_kept_2900m', sum(kept_idx), 'F5_n_shots_cut_2900m', sum(cut_idx), ...
    'F5_kept_angle_med_deg', median(th_s_2900(kept_idx)), 'F5_kept_angle_max_deg', max(th_s_2900(kept_idx)), ...
    'F5_cut_angle_med_deg', median(th_s_2900(cut_idx)), 'F5_cut_angle_min_deg', min(th_s_2900(cut_idx)), ...
    'F14_x_m', f14_x_med, 'F14_fold_sum_2900m', sum(w_f14_2900), 'F14_n_shots_kept_2900m', sum(w_f14_2900 >= 0.5));
clear G45
fprintf('  beta_anc: RAW=%.5f vs DSA=%.5f (diff=%+.5f /Hz -> RAW has stronger HF cut!)\n', beta_raw, beta_dsa, beta_raw-beta_dsa);
fprintf('  F5 m=1: RAW=%.4f, DSA=%.4f (diff=%+.4f) | Swap slope: DSA+RAW_beta=%.4f (%+.4f), RAW+DSA_beta=%.4f (%+.4f)\n', ...
    cc_raw_orig, cc_dsa_orig, cc_raw_orig-cc_dsa_orig, cc_dsa_swap, cc_dsa_swap-cc_dsa_orig, cc_raw_swap, cc_raw_swap-cc_raw_orig);
fprintf('  Sub-bands [8-19 / 19-29 / 29-40 Hz]: RAW=[%.3f, %.3f, %.3f] vs DSA=[%.3f, %.3f, %.3f]\n', ...
    sub_raw(1), sub_raw(2), sub_raw(3), sub_dsa(1), sub_dsa(2), sub_dsa(3));
fprintf('  F-5 fold at 2900m: kept=%d/%d shots (med angle %.1f deg, max %.1f deg), cut=%d shots (med angle %.1f deg, min %.1f deg)\n', ...
    d2.F5_n_shots_kept_2900m, numel(xs45), d2.F5_kept_angle_med_deg, d2.F5_kept_angle_max_deg, d2.F5_n_shots_cut_2900m, d2.F5_cut_angle_med_deg, d2.F5_cut_angle_min_deg);
diag.D2_RAW_vs_DSA = d2;

%% ========================================================================
%  D3: F-5 Zero-Lag L0 & Optimal Depth Lag Trajectory lag*(m) (25 Hz & 45 Hz)
%  ========================================================================
fprintf('\n--- D3: F-5 Depth Lag Curves corr(lag) & lag*(m) ---\n');
lags_d3 = -40:2.5:40;
d3 = struct('lags_m', lags_d3, 'volve25', struct(), 'volve45', struct());

% 4-25 Hz
m25_list = {'RAW_R0', 'DSA_R0_Eq', 'RAW_UQW1_Unanc_25', 'DSA_UQW1b_PAnc_25', 'DSA_UQW1b_PAnc_Amp', 'DSA_UQW1b_PAnc_100'};
rb5_25 = core.bp_z(w5_25.rw, w5_25.dzw, k25(1), k25(2));
d3.volve25.m = S25.res.m;
for i = 1:numel(m25_list)
    mk = m25_list{i}; M = S25.Img_s.(mk); nc = size(M, 2);
    curves = zeros(nc, numel(lags_d3)); best_lag = zeros(1, nc); best_cc = zeros(1, nc); l0_cc = zeros(1, nc);
    for c = 1:nc
        A = prep25(M(:, c));
        tr = core.bp_z(extract_well_trace(A, w5_25, dx25, 20.0), w5_25.dzw, k25(1), k25(2));
        curves(c, :) = lag_curve(tr, rb5_25, w5_25.dzw, lags_d3, core);
        [best_cc(c), best_lag(c), l0_cc(c)] = core.well_tie_full(A, w5_25, rb5_25, dx25, k25(1), k25(2), 20.0, 40.0);
    end
    d3.volve25.(mk) = struct('corr_lag', curves, 'best_lag_m', best_lag, 'best_cc', best_cc, 'L0_cc', l0_cc);
    fprintf('  [25Hz] %-22s | lag*(m)=[%s] | L0=[%s]\n', mk, num2str(best_lag, '%+4.0f '), num2str(l0_cc, '%.3f '));
end

% 4-45 Hz
m45_list = {'RAW_R0', 'DSA_R0_Eq', 'RAW_UQW1b_Unanc_100', 'DSA_UQW1b_PAnc_25', 'DSA_UQW1b_PAnc_Amp', 'DSA_UQW1b_PAnc_100'};
d3.volve45.m = S45.res.m;
for i = 1:numel(m45_list)
    mk = m45_list{i}; M = S45.Img_s.(mk); nc = size(M, 2);
    curves = zeros(nc, numel(lags_d3)); best_lag = zeros(1, nc); best_cc = zeros(1, nc); l0_cc = zeros(1, nc);
    for c = 1:nc
        A = prep45(M(:, c));
        tr = core.bp_z(extract_well_trace(A, w5_45, dx45, 20.0), w5_45.dzw, kfix(1), kfix(2));
        curves(c, :) = lag_curve(tr, rb5_fix, w5_45.dzw, lags_d3, core);
        [best_cc(c), best_lag(c), l0_cc(c)] = core.well_tie_full(A, w5_45, rb5_fix, dx45, kfix(1), kfix(2), 20.0, 40.0);
    end
    d3.volve45.(mk) = struct('corr_lag', curves, 'best_lag_m', best_lag, 'best_cc', best_cc, 'L0_cc', l0_cc);
    fprintf('  [45Hz] %-22s | lag*(m)=[%s] | L0=[%s]\n', mk, num2str(best_lag, '%+4.0f '), num2str(l0_cc, '%.3f '));
end
diag.D3_Lag_Curves = d3;

%% ========================================================================
%  D4: Synthetic Twin U4b Diagnostic (Why B=300 < B=1 / IQ70 / DSA_R0_Eq on well_cc)
%  ========================================================================
fprintf('\n--- D4: Synthetic Twin U4b Mechanism Diagnostic ---\n');
diag.D4_U4b = run_d4_u4b_diag(core, cfg_s);

%% Save JSON
diag.elapsed_s = toc(t0);
fn_json = fullfile(res_dir, 'wp1_diagnostics.json');
fid = fopen(fn_json, 'w'); fwrite(fid, jsonencode(diag, 'PrettyPrint', true), 'char'); fclose(fid);
fprintf('\nSaved %s (total %.1fs)\n', fn_json, diag.elapsed_s);

%% Plot 2x2 Diagnostic Figure
make_wp1_figure(diag, fullfile(fig_dir, 'fig_wp1_diagnostics.png'));
fprintf('Saved %s\n', fullfile(fig_dir, 'fig_wp1_diagnostics.png'));
end

function tr = extract_well_trace(A, w, dx_grid, half_avg_m)
na = round(half_avg_m / dx_grid);
nw = numel(w.zw); tr = zeros(nw, 1);
for i = 1:nw
    tr(i) = mean(A(w.iz(i), max(w.ix(i) - na, 1):min(w.ix(i) + na, size(A, 2))));
end
end

function cvec = lag_curve(tr_b, rb, dzw, lags_m, core)
cvec = zeros(1, numel(lags_m));
for i = 1:numel(lags_m)
    L = round(lags_m(i) / dzw);
    cvec(i) = core.cc(circshift(tr_b, L), rb);
end
end

function [best_c, best_lag] = max_lag_cc(tr_b, rb, dzw, max_lag_m, core)
ns_ = round(max_lag_m / dzw);
best_c = -Inf; best_lag = 0;
for L = -ns_:ns_
    c = core.cc(circshift(tr_b, L), rb);
    if c > best_c, best_c = c; best_lag = L * dzw; end
end
end

function y = apply_kz_slope(tr, dzw, vr, d_beta)
n = 2^nextpow2(4 * numel(tr));
R = fft(tr, n);
kz = abs([0:n/2, -n/2+1:-1]') / (n * dzw); % cyc/m
f_eq = kz * vr / 2.0;
H = exp(d_beta * (f_eq - 25.0));
H(f_eq > 45.0) = 1.0;
y = real(ifft(R .* H));
y = y(1:numel(tr));
end

function d4 = run_d4_u4b_diag(core, cfg)
seed = 20261304; rng(seed, 'twister');
P.dx=10; P.npml=18; P.nzi=100; P.nxi=160; P.nz=136; P.nx=196;
P.z=((1:P.nz)'-19)*10; P.x=((1:P.nx)-19)*10; P.in={19:118, 19:178}; P.izs=21;
P.ixs=19+(8:3:152); P.ixr=19+(0:1:159); ns=numel(P.ixs); nr=numel(P.ixr); N=P.nz*P.nx;
[Zi,Xi]=ndgrid(P.in{1},P.in{2}); lin=sub2ind([P.nz P.nx],Zi(:),Xi(:)); nin=numel(lin);
isrc=sub2ind([P.nz P.nx],P.izs*ones(1,ns),P.ixs); irec=sub2ind([P.nz P.nx],P.izs*ones(1,nr),P.ixr);
[Zfull,Xfull]=ndgrid(P.z,P.x); z_in=P.z(P.in{1}); x_in=P.x(P.in{2}); [Zin,Xin]=ndgrid(z_in,x_in);
xs=P.x(P.ixs); xr=P.x(P.ixr);

rec_mask = false(nr, ns); xr_lo = zeros(1,ns); xr_hi = zeros(1,ns);
for s=1:ns, rec_mask(:,s)=(xr>=xs(s)-240)&(xr<=xs(s)+860); xr_lo(s)=min(xr(rec_mask(:,s))); xr_hi(s)=max(xr(rec_mask(:,s))); end
v_bg=1800+1.45*max(Zfull,0); z_dip=Zfull+0.11*(Xfull-800); dv=zeros(size(Zfull));
ze=[175,225,285,345,415,485,555,590,630,665,705,745,785,825,860,900,945];
for il=1:numel(ze)-1, dv(z_dip>=ze(il)&z_dip<ze(il+1))=max(min(0.055*randn(),0.09),-0.09); end
v_true=v_bg.*(1+dv); v_mig=core.smooth2d(v_true,3.0);
Q_true=55*ones(P.nz,P.nx); Q_true(Zfull<80)=1000; Q_true(Zfull>=520)=35;

s_mig=1./v_mig; eps_q=0.5; F.T=zeros(P.nz,P.nx,ns,'single'); F.ts=F.T;
for k=1:ns
    F.T(:,:,k)=eikonal_fs(s_mig,10,P.izs,P.ixs(k));
    F.ts(:,:,k)=(eikonal_fs(s_mig.*(1+eps_q./Q_true),10,P.izs,P.ixs(k))-eikonal_fs(s_mig.*(1-eps_q./Q_true),10,P.izs,P.ixs(k)))/(2*eps_q);
end
TS_flat=zeros(nin,ns); Fts_in=zeros(nin,ns);
for k=1:ns, [~,tt]=two_way(F,k,P); TS_flat(:,k)=double(tt(lin)); tk=F.ts(:,:,k); Fts_in(:,k)=double(tk(lin)); end; clear F

freqs=4:1:37; nf=numel(freqs); fref=25; fpeak_true=20; sig=cfg.sig; s70=gain_s(70); ms=[0.5 0.7 1.0 1.4 2.0]; nm=numel(ms);
corr_noise_sig=5.0; snr_db=2.0; rc=ceil(3*corr_noise_sig); g_nc=exp(-0.5*((-rc:rc)'/corr_noise_sig).^2); g_nc=g_nc/sqrt(sum(g_nc.^2));
rng(seed+1,'twister'); EPS_B1=2*(rand(nr,ns)>.5)-1; rng(seed+300,'twister'); sgn_B300=2*(rand(1,ns)>.5)-1;
W_mig=ones(1,nf); a1=freqs<6; W_mig(a1)=0.5*(1-cos(pi*(freqs(a1)-4)/2)); a2=freqs>34; W_mig(a2)=0.5*(1+cos(pi*(freqs(a2)-34)/3));

W20=(20/fpeak_true)^2*exp(-(20/fpeak_true)^2); b20=zeros(N,ns); b20(sub2ind([N ns],isrc,1:ns))=-W20/100;
[Lm,Um,Pm,Qm,Rm]=lu(helm_op(v_mig,20,P)); Ub=Qm*(Um\(Lm\(Pm*(Rm\b20)))); Dbg20=Ub(irec,:); clear Ub Lm Um Pm Qm Rm
vc20=v_true./(1-(1./(pi*Q_true))*log(20/fref)-1i./(2*Q_true));
[Lv,Uv,Pv,Qv,Rv]=lu(helm_op(vc20,20,P)); Uv20=Qv*(Uv\(Lv\(Pv*(Rv\b20)))); Dcal=(Uv20(irec,:)-Dbg20).*rec_mask; clear Uv20 Lv Uv Pv Qv Rv
sig_field=sqrt(mean(abs(Dcal(rec_mask)).^2))*10^(-snr_db/20);

Cs_sig_all=zeros(nin,ns,nf,'single'); Cs_nu_all=Cs_sig_all; Cs_b1_all=Cs_sig_all; Cs_el_all=Cs_sig_all; I_RAW_R0=zeros(nin,1);
for jf=1:nf
    f=freqs(jf); w=2*pi*f; W_true=(f/fpeak_true)^2*exp(-(f/fpeak_true)^2);
    bt=zeros(N,ns); bt(sub2ind([N ns],isrc,1:ns))=-W_true/100;
    bm=zeros(N,ns); bm(sub2ind([N ns],isrc,1:ns))=-W_mig(jf)/100;
    [Lm,Um,Pm,Qm,Rm]=lu(helm_op(v_mig,f,P));
    Ub=Qm*(Um\(Lm\(Pm*(Rm\bt)))); D_bg=Ub(irec,:); clear Ub
    Um_=Qm*(Um\(Lm\(Pm*(Rm\bm)))); Us_c=conj(Um_(lin,:))*w^2; clear Um_
    [Le,Ue,Pe,Qe,Re]=lu(helm_op(v_true,f,P)); Uel=Qe*(Ue\(Le\(Pe*(Re\bt)))); D_el=(Uel(irec,:)-D_bg).*rec_mask; clear Uel Le Ue Pe Qe Re
    vc=v_true./(1-(1./(pi*Q_true))*log(max(f,1e-6)/fref)-1i./(2*Q_true));
    [Lv,Uv,Pv,Qv,Rv]=lu(helm_op(vc,f,P)); Uv_=Qv*(Uv\(Lv\(Pv*(Rv\bt)))); D_vi=(Uv_(irec,:)-D_bg).*rec_mask; clear Uv_ Lv Uv Pv Qv Rv
    N_f=conv2(g_nc,1,sig_field*(1/sqrt(2))*(randn(nr,ns)+1i*randn(nr,ns)),'same').*rec_mask;
    Vel=extract_bp(D_el,Lm,Um,Pm,Qm,Rm,irec,lin,N,ns,10);
    Vs=extract_bp(D_vi,Lm,Um,Pm,Qm,Rm,irec,lin,N,ns,10); Vn=extract_bp(N_f,Lm,Um,Pm,Qm,Rm,irec,lin,N,ns,10);
    Vb=extract_bp((D_vi+N_f).*EPS_B1,Lm,Um,Pm,Qm,Rm,irec,lin,N,ns,10); clear Lm Um Pm Qm Rm
    Cs_el_all(:,:,jf)=single(Us_c.*Vel);
    Cs_sig_all(:,:,jf)=single(Us_c.*Vs); Cs_nu_all(:,:,jf)=single(Us_c.*Vn); Cs_b1_all(:,:,jf)=single(Us_c.*Vb);
    I_RAW_R0=I_RAW_R0+real(sum(Us_c.*(Vs+Vn),2));
end

R0_2d=core.prep_image(I_RAW_R0,P.nzi,P.nxi,z_in,10,cfg.bg_sigma_m,cfg.z_mute_m);
[phi_eff,~]=core.blind_dip(R0_2d,10,cfg.dip_sigma_m,cfg.dip_cap_deg);
[w_dsa,TS_dsa]=core.build_dsa_in(Xin(:),Zin(:),20,xs,xr_lo,xr_hi,phi_eff(:),Fts_in,cfg);

I_dsa_s=zeros(nin,nf,'single'); I_dsa_n=I_dsa_s; I_dsa_t=I_dsa_s; I_dsa_f1=I_dsa_s; I_dsa_f3=I_dsa_s;
for jf=1:nf
    cs=double(Cs_sig_all(:,:,jf)); cn=double(Cs_nu_all(:,:,jf)); ct=cs+cn; cb=double(Cs_b1_all(:,:,jf));
    I_dsa_s(:,jf)=single(real(sum(w_dsa.*cs,2))); I_dsa_n(:,jf)=single(real(sum(w_dsa.*cn,2)));
    I_dsa_t(:,jf)=single(real(sum(w_dsa.*ct,2))); I_dsa_f1(:,jf)=single(real(sum(w_dsa.*cb,2)));
    I_dsa_f3(:,jf)=single(real(sum((w_dsa.*ct).*sgn_B300,2)));
end
iv_dsa_b1=core.band_inv_snr(I_dsa_t,I_dsa_f1,freqs,P.nzi,P.nxi,10,cfg.snr_band_hz,cfg.snr_box_m,cfg.snr_cap);
iv_dsa_b3=core.band_inv_snr(I_dsa_t,I_dsa_f3,freqs,P.nzi,P.nxi,10,cfg.snr_band_hz,cfg.snr_box_m,cfg.snr_cap);
iv_tr_dsa=true_inv_snr_loc(I_dsa_s,I_dsa_n,freqs,P.nzi,P.nxi,10,cfg.snr_band_hz,cfg.snr_box_m,cfg.snr_cap);

anc_mask=(Zin(:)>=190&Zin(:)<=320)&(Xin(:)>=240&Xin(:)<=1360);
res_mask=(Zin(:)>=540&Zin(:)<=940)&(Xin(:)>=240&Xin(:)<=1360);
S_anc_dsa=arrayfun(@(j) rms(double(I_dsa_t(anc_mask,j))),1:nf);
[W_eq_dsa,beta_d,~]=core.fit_psae(freqs,S_anc_dsa,[6 32],20);
tanc_d=median(TS_dsa(anc_mask,:),'all'); TS_dr=max(TS_dsa-tanc_d,0);

% Evaluate at m=1.0: DSA_R0_Eq, DSA_IQ70_PAnc_25, DSA_UQW1_PAnc_25, DSA_UQW1b_PAnc_25 vs Elastic reference (with W_eq_dsa!)
I_el_eq = zeros(nin,1); I_r0_eq = zeros(nin,1); I_iq = zeros(nin,1); I_u1 = zeros(nin,1); I_ub = zeros(nin,1);
spec_el = zeros(1,nf); spec_r0 = zeros(1,nf); spec_iq = zeros(1,nf); spec_u1 = zeros(1,nf); spec_ub = zeros(1,nf);
xmax=1.05*max(TS_dsa(:))/0.5; xg=linspace(0,xmax,1201); dxg=xg(2);
for jf=1:nf
    f=freqs(jf); wd=W_eq_dsa(jf);
    Ct=double(Cs_sig_all(:,:,jf))+double(Cs_nu_all(:,:,jf));
    Ce=double(Cs_el_all(:,:,jf));
    T0=uq_tab(f,xg,0,fref); T1=uq_tab(f,xg,sig,fref);
    ivd1=double(iv_dsa_b1(:,jf)); ivd3=double(iv_dsa_b3(:,jf));
    u=TS_dr/dxg; i0=min(floor(u),numel(xg)-2); ww=u-i0; i0=i0+1; it=@(tb) tb(i0).*(1-ww)+tb(i0+1).*ww;
    A0=it(T0.A); C0=it(T0.C); P0=it(T0.P); A1=it(T1.A); C1=it(T1.C); P1=it(T1.P);
    fe = wd * real(sum(w_dsa.*Ce,2));
    fr = wd * real(sum(w_dsa.*Ct,2));
    fiq= wd * real(sum((w_dsa.*((A0+s70*P0.*C0)./(1+s70*C0))).*Ct,2));
    fu1= wd * real(sum((w_dsa.*((A1+s70*P1.*C1+P1.*ivd1)./(1+s70*C1+ivd1))).*Ct,2));
    fub= wd * real(sum((w_dsa.*((A1+s70*P1.*C1+P1.*ivd3)./(1+s70*C1+ivd3))).*Ct,2));
    I_el_eq=I_el_eq+fe; I_r0_eq=I_r0_eq+fr; I_iq=I_iq+fiq; I_u1=I_u1+fu1; I_ub=I_ub+fub;
    spec_el(jf)=rms(fe(res_mask)); spec_r0(jf)=rms(fr(res_mask));
    spec_iq(jf)=rms(fiq(res_mask)); spec_u1(jf)=rms(fu1(res_mask)); spec_ub(jf)=rms(fub(res_mask));
end

% Compute BOTH well_cc AND amplitude recovery (deep/shallow ratio vs elastic truth) AND L2 waveform error vs elastic truth!
well_xs=[400 600 800 1000 1200]; iz_res=find(z_in>=540&z_in<=940); ix_res=find(x_in>=240&x_in<=1360);
iz_anc=find(z_in>=190&z_in<=320);
A_el = core.prep_image(I_el_eq,P.nzi,P.nxi,z_in,10,cfg.bg_sigma_m,cfg.z_mute_m);
ref_ratio = rms(A_el(iz_res,ix_res),'all') / max(rms(A_el(iz_anc,ix_res),'all'), 1e-30);

v_in=v_true(P.in{1},P.in{2}); rw_map=zeros(size(v_in)); rw_map(2:end,:)=diff(log(v_in),1,1);
names = {'DSA_R0_Eq', 'DSA_IQ70_PAnc_25', 'DSA_UQW1_PAnc_25', 'DSA_UQW1b_PAnc_25'};
vecs  = {I_r0_eq, I_iq, I_u1, I_ub};
d4_m = struct();
for k = 1:numel(names)
    A = core.prep_image(vecs{k},P.nzi,P.nxi,z_in,10,cfg.bg_sigma_m,cfg.z_mute_m);
    wcc = score_w5_loc(core, A, rw_map, x_in, iz_res, well_xs);
    e1_el = core.score_e1(A, A_el, iz_res, ix_res);
    rat = rms(A(iz_res,ix_res),'all') / max(rms(A(iz_anc,ix_res),'all'), 1e-30);
    eamp = abs(log10(rat / ref_ratio));
    rel_l2 = norm(A(iz_res,ix_res) - A_el(iz_res,ix_res), 'fro') / max(norm(A_el(iz_res,ix_res), 'fro'), 1e-30);
    d4_m.(names{k}) = struct('well_cc_m1', wcc, 'E1_vs_elastic_m1', e1_el, 'E_amp_deep_shallow_m1', eamp, 'rel_L2_vs_elastic_m1', rel_l2);
    fprintf('  [U4b m=1] %-18s | well_cc=%.4f | E1_vs_el=%.4f | E_amp=%.4f | relL2=%.4f\n', names{k}, wcc, e1_el, eamp, rel_l2);
end

iv_med_tr = median(double(iv_tr_dsa(res_mask, :)), 1);
iv_med_b1 = median(double(iv_dsa_b1(res_mask, :)), 1);
iv_med_b3 = median(double(iv_dsa_b3(res_mask, :)), 1);
fprintf('  [U4b iv med @ 15/20/25/30 Hz] true=[%.3f,%.3f,%.3f,%.3f] | B1=[%.3f,%.3f,%.3f,%.3f] | B300=[%.3f,%.3f,%.3f,%.3f]\n', ...
    iv_med_tr(12), iv_med_tr(17), iv_med_tr(22), iv_med_tr(27), ...
    iv_med_b1(12), iv_med_b1(17), iv_med_b1(22), iv_med_b1(27), ...
    iv_med_b3(12), iv_med_b3(17), iv_med_b3(22), iv_med_b3(27));

d4 = struct('freqs', freqs, 'beta_anc_dsa', beta_d, 'S_anc_dsa', S_anc_dsa, ...
    'iv_med_true', iv_med_tr, 'iv_med_b1', iv_med_b1, 'iv_med_b300', iv_med_b3, ...
    'spec_el', spec_el, 'spec_r0', spec_r0, 'spec_iq', spec_iq, 'spec_u1', spec_u1, 'spec_ub', spec_ub, ...
    'metrics_m1', d4_m);
end

function c_mean = score_w5_loc(core, A, rw_map, x_in, iz_res, well_xs)
cs=zeros(1,numel(well_xs));
for iw=1:numel(well_xs)
    [~,jx]=min(abs(x_in-well_xs(iw))); tr=mean(A(iz_res,max(jx-1,1):min(jx+1,numel(x_in))),2);
    tb=core.bp_z(tr,10,4,18); rb=core.bp_z(rw_map(iz_res,jx),10,4,18);
    best=-Inf; for L=-3:3, best=max(best,core.cc(circshift(tb,L),rb)); end; cs(iw)=best;
end
c_mean=mean(cs);
end

function iv = true_inv_snr_loc(Is_mat, In_mat, freqs, nzi, nxi, dx, band_hz, box_m, cap)
[nin, nf] = size(Is_mat); bw = max(1, round(box_m / dx)); iv = zeros(nin, nf, 'single');
for j = 1:nf
    bnd = abs(freqs - freqs(j)) <= band_hz + 1e-9;
    ps = movmean(movmean(reshape(sum(double(Is_mat(:,bnd)),2),nzi,nxi).^2,bw,1),bw,2);
    pn = movmean(movmean(reshape(sum(double(In_mat(:,bnd)),2),nzi,nxi).^2,bw,1),bw,2);
    iv(:, j) = single(min(pn(:) ./ max(ps(:), realmin('single')), cap));
end
end

function V = extract_bp(Dmat, L, U, Pp, Qq, R, irec, lin, N, ns, dx)
b = zeros(N, ns); b(irec, :) = Dmat;
Vf = -(R \ (Pp' * (L' \ (U' \ (Qq' * b))))) / dx^2; V = Vf(lin, :);
end

function make_wp1_figure(diag, out_png)
fig = figure('Visible', 'off', 'Position', [50 50 1400 950], 'Color', 'w');

% (a) D1: F-14 lag curves in kfix vs k25
subplot(2, 2, 1); hold on; grid on;
lags = diag.D1_F14.lags_m;
plot(lags, diag.D1_F14.methods.RAW_R0.cc_kfix_curve, '--k', 'LineWidth', 1.4, 'DisplayName', 'RAW\_R0 (8-40 Hz)');
plot(lags, diag.D1_F14.methods.RAW_UQW1b_Unanc_100.cc_kfix_curve, ':m', 'LineWidth', 1.5, 'DisplayName', 'Stage4H Unanc 100Hz (8-40 Hz)');
plot(lags, diag.D1_F14.methods.DSA_UQW1b_PAnc_25.cc_kfix_curve, '-r', 'LineWidth', 1.8, 'DisplayName', 'DSA\_UQW1b\_PAnc\_25 (8-40 Hz)');
plot(lags, diag.D1_F14.methods.DSA_UQW1b_PAnc_25.cc_k25_curve, '-b', 'LineWidth', 1.8, 'DisplayName', 'DSA\_UQW1b\_PAnc\_25 (k25: 3.5-12.9 cyc/km)');
xline(-30, ':k', 'HandleVisibility', 'off'); xline(30, ':k', 'HandleVisibility', 'off'); yline(0, '-k', 'HandleVisibility', 'off');
xlabel('Depth Lag (m)'); ylabel('Correlation r(lag)');
title('(a) D1: Volve 4-45 Hz Blind Well F-14 Lag Scan (kfix vs k25)');
legend('Location', 'best', 'FontSize', 8);

% (b) D2: RAW vs DSA sub-band & slope-swap breakdown
subplot(2, 2, 2); hold on; grid on;
bar_data = [diag.D2_RAW_vs_DSA.F5_RAW_orig, diag.D2_RAW_vs_DSA.F5_DSA_orig; ...
            diag.D2_RAW_vs_DSA.F5_RAW_with_DSA_slope, diag.D2_RAW_vs_DSA.F5_DSA_with_RAW_slope; ...
            diag.D2_RAW_vs_DSA.subband_RAW(1), diag.D2_RAW_vs_DSA.subband_DSA(1); ...
            diag.D2_RAW_vs_DSA.subband_RAW(2), diag.D2_RAW_vs_DSA.subband_DSA(2); ...
            diag.D2_RAW_vs_DSA.subband_RAW(3), diag.D2_RAW_vs_DSA.subband_DSA(3)];
b = bar(bar_data);
b(1).DisplayName = 'RAW\_UQW1b\_PAnc\_25'; b(2).DisplayName = 'DSA\_UQW1b\_PAnc\_25';
set(gca, 'XTick', 1:5, 'XTickLabel', {'Orig (8-40Hz)', 'Swapped \beta_{anc}', '8-19 Hz', '19-29 Hz', '29-40 Hz'});
ylabel('F-5 Correlation'); title('(b) D2: RAW vs DSA Slope-Swap & Sub-Band Breakdown');
legend('Location', 'best', 'FontSize', 8);

% (c) D3: F-5 lag*(m) & L0(m) trajectories on 4-45 Hz
subplot(2, 2, 3); hold on; grid on;
m45 = diag.D3_Lag_Curves.volve45.m;
plot(m45, diag.D3_Lag_Curves.volve45.RAW_R0.best_lag_m * ones(size(m45)), '--k', 'LineWidth', 1.5, 'DisplayName', 'RAW\_R0 lag*');
plot(m45, diag.D3_Lag_Curves.volve45.RAW_UQW1b_Unanc_100.best_lag_m, '-om', 'LineWidth', 1.8, 'DisplayName', 'Stage4H Unanc 100Hz lag*');
plot(m45, diag.D3_Lag_Curves.volve45.DSA_UQW1b_PAnc_100.best_lag_m, '-^c', 'LineWidth', 1.6, 'DisplayName', 'DSA\_PAnc 100Hz lag*');
plot(m45, diag.D3_Lag_Curves.volve45.DSA_UQW1b_PAnc_25.best_lag_m, '-sr', 'LineWidth', 2.0, 'DisplayName', 'DSA\_PAnc 25Hz lag*');
plot(m45, diag.D3_Lag_Curves.volve45.DSA_UQW1b_PAnc_Amp.best_lag_m, ':db', 'LineWidth', 1.6, 'DisplayName', 'DSA\_PAnc Amp lag*');
xlabel('Q Perturbation Multiplier m'); ylabel('Optimal Depth Lag lag*(m) [m]');
title('(c) D3: Volve 4-45 Hz F-5 Optimal Lag Trajectory vs m');
ylim([-35 15]); legend('Location', 'best', 'FontSize', 8);

% (d) D4: U4b per-frequency iv & amplitude spectrum vs elastic truth
subplot(2, 2, 4); hold on; grid on;
f4 = diag.D4_U4b.freqs;
semilogy(f4, diag.D4_U4b.spec_el / max(diag.D4_U4b.spec_el), '-k', 'LineWidth', 2.0, 'DisplayName', 'Elastic Truth (Eq)');
semilogy(f4, diag.D4_U4b.spec_r0 / max(diag.D4_U4b.spec_el), '--', 'LineWidth', 1.5, 'DisplayName', 'DSA\_R0\_Eq');
semilogy(f4, diag.D4_U4b.spec_iq / max(diag.D4_U4b.spec_el), ':g', 'LineWidth', 1.6, 'DisplayName', 'DSA\_IQ70\_PAnc\_25');
semilogy(f4, diag.D4_U4b.spec_u1 / max(diag.D4_U4b.spec_el), '-.b', 'LineWidth', 1.5, 'DisplayName', 'DSA\_UQW1\_PAnc\_25 (B=1)');
semilogy(f4, diag.D4_U4b.spec_ub / max(diag.D4_U4b.spec_el), '-r', 'LineWidth', 1.8, 'DisplayName', 'DSA\_UQW1b\_PAnc\_25 (B=300)');
xlabel('Frequency (Hz)'); ylabel('Normalized Reservoir RMS Spectrum');
title('(d) D4: U4b Reservoir Spectrum vs Elastic Truth (m=1)');
xlim([4 37]); ylim([1e-2 3]); legend('Location', 'best', 'FontSize', 8);

exportgraphics(fig, out_png, 'Resolution', 250);
close(fig);
end
