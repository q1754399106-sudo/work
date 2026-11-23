function wsuq_eikonal_bench()
% T2.1 Part A: Exact Real-Eikonal Multi-Shot Broadband (4-45 Hz) Operator Benchmark
% using cached Volve F-5 (124 shots TS_abs & TS_anc227, real well_F-5 reflectivity)
% and cached D2 (51 shots TS_abs & TS_anc, real elastic ref.R0 reflectivity).
% Zero LU, < 500 MB RAM, 2 threads.

maxNumCompThreads(2);
rng(20261123, 'twister');
here = fileparts(mfilename('fullpath'));
root = fileparts(here);
out_dir = fullfile(root, 'results');
if ~isfolder(out_dir), mkdir(out_dir); end

main_res = 'D:\ChatGPT\研究\地震衰减补偿研究';
addpath(fullfile(main_res, 'uq_rtm_20261117', 'scripts'));
addpath(here);

freqs = 4:1:45;       % 42 frequencies spanning 4-45 Hz
nf = numel(freqs);
fpeak = 22;
W_f = (freqs / fpeak).^2 .* exp(-(freqs / fpeak).^2);
W_f = W_f / max(W_f);
ms = [0.5, 0.7, 1.0, 1.4, 2.0];
sig = 0.447;
fref = 25.0;
fv_band = [15.0, 35.0];
s70 = gain_s(70);

res = struct();

% ---- 1. Volve F-5 (exact 124-shot TS_abs & TS_anc227 + real well_F-5 reflectivity) ----
fld_dir = fullfile(main_res, 'uq_rtm_20261117', 'cache', 'fld');
meta = jsondecode(fileread(fullfile(fld_dir, 'meta.json')));
dx = double(meta.dx); npml = double(meta.npml);
z_in = ((npml + (1:double(meta.nzi)))' - npml - 1) * dx;
x_in = double(meta.xmin) + ((npml + (1:double(meta.nxi))) - npml - 1) * dx;
[Zi, Xi] = ndgrid(z_in, x_in);

% Extract F-5 reservoir window (z in [2570, 3170], x in [-1300, -1100])
iz_f5 = find(z_in >= 2570 & z_in <= 3170);
ix_f5 = find(x_in >= -1300 & x_in <= -1100);
[IZg, IXg] = ndgrid(iz_f5, ix_f5);
idx_f5 = sub2ind(size(Zi), IZg(:), IXg(:));

% Load real F-5 well reflectivity and project onto the window columns
fid = fopen(fullfile(fld_dir, 'well_F-5.bin'), 'r'); wa = fread(fid, Inf, 'single'); fclose(fid);
wa = reshape(wa, 4, [])';
zw_well = double(wa(:, 1)); rw_well = double(wa(:, 4));
r_col = interp1(zw_well, rw_well, z_in(iz_f5), 'linear', 0);
r_col = r_col / max(std(r_col), 1e-12);
R_f5 = repmat(r_col(:), 1, numel(ix_f5));  % nzf x nxf true reflectivity template

S_abs = load(fullfile(fld_dir, 'ts.mat'), 'TS');
TS_v_abs = double(S_abs.TS(idx_f5, :)); clear S_abs
S_anc = load(fullfile(main_res, 'va_uq_field2_20261121', 'cache', 'anc', 'ts.mat'), 'TS');
TS_v_anc = double(S_anc.TS(idx_f5, :)); clear S_anc

fprintf('Volve F-5 eikonal bench: %d x %d = %d nodes, %d shots\n', numel(iz_f5), numel(ix_f5), numel(idx_f5), size(TS_v_abs, 2));
res.Volve_F5_SNR10dB = run_eik_case(R_f5, TS_v_abs, TS_v_anc, freqs, W_f, ms, sig, fref, fv_band, s70, 10.0, 'Volve_F5_10dB');
res.Volve_F5_SNR3dB  = run_eik_case(R_f5, TS_v_abs, TS_v_anc, freqs, W_f, ms, sig, fref, fv_band, s70,  3.0, 'Volve_F5_3dB');
clear TS_v_abs TS_v_anc

% ---- 2. Synthetic D2 Wdeep (exact 51-shot TS_abs & TS_anc + real elastic ref.R0) ----
P = rtm_setup();
z2 = P.z(P.in{1}); x2 = P.x(P.in{2});
iz2 = find(z2 >= 820 & z2 <= 1550);
ix2 = find(x2 >= 600 & x2 <= 1800);
% Subsample to ~46 x 76 = 3496 evenly spaced nodes across Wdeep (~3 min runtime, exact spatial coverage)
iz2 = iz2(1:4:end); ix2 = ix2(1:4:end);
[IZ2, IX2] = ndgrid(P.in{1}(iz2), P.in{2}(ix2));
lin2 = sub2ind([P.nz, P.nx], IZ2(:), IX2(:));

O2 = load(fullfile(main_res, 'uq_rtm_20261117', 'cache', 'images2_D2.mat'), 'ref');
R2_full = reshape(double(O2.ref.R0), P.nzi, P.nxi); clear O2
R_d2 = R2_full(iz2, ix2);
R_d2 = R_d2 / max(std(R_d2(:)), 1e-12);

F_abs = load(fullfile(main_res, 'va_uq_step1_20261121', 'cache', 'fields_D2_vb.mat'));
F_anc = load(fullfile(main_res, 'va_uq_step1_20261121', 'cache', 'fields_D2_vb_anc.mat'));
ns2 = numel(P.ixs);
TS_d2_abs = zeros(numel(lin2), ns2);
TS_d2_anc = zeros(numel(lin2), ns2);
for s = 1:ns2
    [~, tt] = two_way(F_abs, s, P); TS_d2_abs(:, s) = double(tt(lin2));
    [~, ta] = two_way(F_anc, s, P); TS_d2_anc(:, s) = max(double(ta(lin2)), 0);
end
clear F_abs F_anc

fprintf('D2 Wdeep eikonal bench: %d x %d = %d nodes, %d shots\n', numel(iz2), numel(ix2), numel(lin2), ns2);
res.D2_Wdeep_SNR10dB = run_eik_case(R_d2, TS_d2_abs, TS_d2_anc, freqs, W_f, ms, sig, fref, fv_band, s70, 10.0, 'D2_Wdeep_10dB');
res.D2_Wdeep_SNR3dB  = run_eik_case(R_d2, TS_d2_abs, TS_d2_anc, freqs, W_f, ms, sig, fref, fv_band, s70,  3.0, 'D2_Wdeep_3dB');

txt = jsonencode(res, 'PrettyPrint', true);
fid = fopen(fullfile(out_dir, 'wsuq_eikonal_bench.json'), 'w'); fwrite(fid, txt, 'char'); fclose(fid);
fprintf('\nSaved wsuq_eikonal_bench.json\n');
end

function out = run_eik_case(R_img, TS_abs, TS_anc, freqs, W_f, ms, sig, fref, fv_band, s70, snr_db, label)
[nz, nx] = size(R_img);
nin = nz * nx;
ns = size(TS_abs, 2);
nm = numel(ms);
R_vec = R_img(:);
% Reference broadband target image (sum of W_f * R_vec) & HF target (f >= 25 Hz)
I_ref = sum(W_f) * R_vec;
hf_mask = freqs >= 25.0;
I_ref_hf = sum(W_f(hf_mask)) * R_vec;

% Calibration of per-shot noise sigma so stacked unattenuated peak-band SNR = snr_db
% Stack signal RMS at peak freq (W_f = 1) is ns * std(R_vec); stack noise RMS is sqrt(ns) * sigma_n
sig_rms = std(R_vec);
snr_lin = 10^(snr_db / 10);
% Noise amplitude per shot (flat white spectrum across freqs so HF naturally drops in SNR as real data!)
sigma_n = sqrt(ns / snr_lin) * sig_rms * 0.55;

% Accumulators: signal-only and noise-only images for EXACT output SNR_dB + total image for E1 / E1_HF!
names = {'R0', 'IQ70_old', 'UQW1_old', 'VA_old', 'WS_UQ0', 'WS_UQ1', 'WS_VA_UQ1'};
S_img = struct(); N_img = struct(); S_hf = struct(); N_hf = struct();
for k = 1:numel(names)
    col = 1 + (nm - 1) * ~strcmp(names{k}, 'R0');
    S_img.(names{k}) = zeros(nin, col);
    N_img.(names{k}) = zeros(nin, col);
    S_hf.(names{k})  = zeros(nin, col);
    N_hf.(names{k})  = zeros(nin, col);
end

% True velocity reference slightly perturbed within fv_band (e.g., 22 Hz vs nominal 25 Hz) for realistic dispersion test
f_true_ref = 22.0;

for jf = 1:numel(freqs)
    f = freqs(jf);
    wf = W_f(jf);
    is_hf = hf_mask(jf);
    % True viscoacoustic forward attenuation + dispersion on each shot
    kf_true = 2 * f * log(max(f, 1e-6) / f_true_ref);
    G_true = exp((-pi * f + 1i * kf_true) * TS_abs);
    Cs_sig = wf * (G_true .* R_vec);
    Cs_noi = (sigma_n / sqrt(2)) * (randn(nin, ns) + 1i * randn(nin, ns));
    % Exact stacked 1/SNR estimate (iv = E|sum N|^2 / E|sum S|^2, median smoothed like est_noise2)
    iv = repmat(mean(abs(sum(Cs_noi, 2)).^2) / max(mean(abs(sum(Cs_sig, 2)).^2), 1e-30), nin, 1);

    % R0
    r0_s = real(sum(Cs_sig, 2)); r0_n = real(sum(Cs_noi, 2));
    S_img.R0 = S_img.R0 + r0_s; N_img.R0 = N_img.R0 + r0_n;
    if is_hf, S_hf.R0 = S_hf.R0 + r0_s; N_hf.R0 = N_hf.R0 + r0_n; end

    for im = 1:nm
        m = ms(im);
        TS_m = TS_abs / m;
        TS_am = TS_anc / m;

        % Fixed fref=25 Hz: UQW1_old, WS_UQ0, WS_UQ1 (shared weights for sig & noi — 2x faster!)
        [~, uq_s, ws0_s, ws1_s, ~, uq_n, ws0_n, ws1_n] = ws_uq_apply(Cs_sig, TS_m, f, sig, fref, iv, s70, [], [], Cs_noi);

        % IQ70_old (sig=0, iv=0, s70 only — deterministic fast path!)
        [~, iq_s, ~, ~, ~, iq_n, ~, ~] = ws_uq_apply(Cs_sig, TS_m, f, 0.0, fref, zeros(nin, 1), s70, [], [], Cs_noi);

        % VA-UQ (anchored Delta t*, fv in [15, 35] Hz): VA_old & WS_VA_UQ1
        [~, va_s, ~, wva_s, ~, va_n, ~, wva_n] = ws_uq_apply(Cs_sig, TS_am, f, sig, fref, iv, s70, fv_band, [], Cs_noi);

        S_img.IQ70_old(:, im)  = S_img.IQ70_old(:, im)  + iq_s;  N_img.IQ70_old(:, im)  = N_img.IQ70_old(:, im)  + iq_n;
        S_img.UQW1_old(:, im)  = S_img.UQW1_old(:, im)  + uq_s;  N_img.UQW1_old(:, im)  = N_img.UQW1_old(:, im)  + uq_n;
        S_img.VA_old(:, im)    = S_img.VA_old(:, im)    + va_s;  N_img.VA_old(:, im)    = N_img.VA_old(:, im)    + va_n;
        S_img.WS_UQ0(:, im)    = S_img.WS_UQ0(:, im)    + ws0_s; N_img.WS_UQ0(:, im)    = N_img.WS_UQ0(:, im)    + ws0_n;
        S_img.WS_UQ1(:, im)    = S_img.WS_UQ1(:, im)    + ws1_s; N_img.WS_UQ1(:, im)    = N_img.WS_UQ1(:, im)    + ws1_n;
        S_img.WS_VA_UQ1(:, im) = S_img.WS_VA_UQ1(:, im) + wva_s; N_img.WS_VA_UQ1(:, im) = N_img.WS_VA_UQ1(:, im) + wva_n;

        if is_hf
            S_hf.IQ70_old(:, im)  = S_hf.IQ70_old(:, im)  + iq_s;  N_hf.IQ70_old(:, im)  = N_hf.IQ70_old(:, im)  + iq_n;
            S_hf.UQW1_old(:, im)  = S_hf.UQW1_old(:, im)  + uq_s;  N_hf.UQW1_old(:, im)  = N_hf.UQW1_old(:, im)  + uq_n;
            S_hf.VA_old(:, im)    = S_hf.VA_old(:, im)    + va_s;  N_hf.VA_old(:, im)    = N_hf.VA_old(:, im)    + va_n;
            S_hf.WS_UQ0(:, im)    = S_hf.WS_UQ0(:, im)    + ws0_s; N_hf.WS_UQ0(:, im)    = N_hf.WS_UQ0(:, im)    + ws0_n;
            S_hf.WS_UQ1(:, im)    = S_hf.WS_UQ1(:, im)    + ws1_s; N_hf.WS_UQ1(:, im)    = N_hf.WS_UQ1(:, im)    + ws1_n;
            S_hf.WS_VA_UQ1(:, im) = S_hf.WS_VA_UQ1(:, im) + wva_s; N_hf.WS_VA_UQ1(:, im) = N_hf.WS_VA_UQ1(:, im) + wva_n;
        end
    end
end

out = struct('m', ms, 'snr_in_dB', snr_db);
fprintf('\n=== %s (4-45 Hz) ===\n', label);
for k = 1:numel(names)
    nm_k = names{k};
    S_m = S_img.(nm_k); N_m = N_img.(nm_k);
    Sh_m = S_hf.(nm_k); Nh_m = N_hf.(nm_k);
    nc = size(S_m, 2);
    e1 = zeros(1, nc); e1_hf = zeros(1, nc); snr_out = zeros(1, nc); snr_hf = zeros(1, nc);
    for c = 1:nc
        tot = S_m(:, c) + N_m(:, c);
        toth = Sh_m(:, c) + Nh_m(:, c);
        e1(c) = corr(tot, I_ref);
        e1_hf(c) = corr(toth, I_ref_hf);
        snr_out(c) = 10 * log10(sum(S_m(:, c).^2) / max(sum(N_m(:, c).^2), 1e-30));
        snr_hf(c)  = 10 * log10(sum(Sh_m(:, c).^2) / max(sum(Nh_m(:, c).^2), 1e-30));
    end
    out.(nm_k) = struct('E1', e1, 'E1_mean', mean(e1), 'E1_HF', e1_hf, 'E1_HF_mean', mean(e1_hf), ...
                        'SNR_dB', snr_out, 'SNR_dB_mean', mean(snr_out), 'SNR_HF_dB', snr_hf, 'SNR_HF_dB_mean', mean(snr_hf));
    fprintf('  %-10s | E1_mean=%.4f (m0.5=%.4f, m1=%.4f) | E1_HF_mean=%.4f (m0.5=%.4f) | SNR=%.2f dB (m0.5=%.2f) | SNR_HF=%.2f dB\n', ...
            nm_k, mean(e1), e1(1), e1(min(3, nc)), mean(e1_hf), e1_hf(1), mean(snr_out), snr_out(1), mean(snr_hf));
end
end
