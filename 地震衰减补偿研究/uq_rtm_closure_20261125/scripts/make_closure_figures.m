function make_closure_figures()
% WP5: Generate 4 publication-grade 300-DPI honest comparison figures
%   fig1_synth_u6_u9_full_honest.png
%   fig2_volve25_psae_vs_now2.png
%   fig3_volve45_true_inv3_and_dual_wells.png
%   fig4_volve45_2d_sections_true_closure.png

here = fileparts(mfilename('fullpath'));
root = fileparts(here);
res_dir = fullfile(root, 'results');
fig_dir = fullfile(root, 'figures');
if ~isfolder(fig_dir), mkdir(fig_dir); end
core = uq_rtm_unified_core();
cfg  = core.default_cfg('field_obc');

%% ========================================================================
%  Fig 1: Synthetic U6-U9 Dual-Metric (E1 + E_amp) & SA vs DSA Tradeoff
%  ========================================================================
models = {'U6', 'U7', 'U8', 'U9'};
titles = {'(a) U6: Gas Chimney + Gentle Dip (l30)', ...
          '(b) U7: 22^\circ Steep Reverse Fault (l30)', ...
          '(c) U8: Non-Scalar Shifted Chimney (+180m, 1.4x)', ...
          '(d) U9: Non-Scalar Opposite Q (1.5x Shallow / 0.6x Deep)'};
m_show = {'RAW_R0', 'SA_UQW1b_25', 'DSA_PH_25', 'DSA_IQ70_25', 'DSA_W1b_25', 'DSA_UQW1_25', 'DSA_UQW1b_25'};
m_lbl  = {'RAW\_R0', 'SA\_UQW1b', 'DSA\_PH', 'DSA\_IQ70', 'DSA\_W1b', 'DSA\_UQW1', 'DSA\_UQW1b'};

fig1 = figure('Visible', 'off', 'Position', [40 40 1500 950], 'Color', 'w');
for i = 1:4
    J = jsondecode(fileread(fullfile(res_dir, ['synth_' models{i} '.json'])));
    e1_vals = zeros(1, numel(m_show)); eamp_vals = zeros(1, numel(m_show));
    for k = 1:numel(m_show)
        e1_vals(k)   = J.W_deep.(m_show{k}).E1_mean;
        eamp_vals(k) = J.W_deep.(m_show{k}).E_amp_mean;
    end
    subplot(2, 2, i);
    yyaxis left;
    b1 = bar((1:numel(m_show)) - 0.18, e1_vals, 0.34, 'FaceColor', [0.22 0.48 0.78]);
    ylabel('Mean E1 (higher is better)'); ylim([0.55 0.84]);
    yyaxis right;
    b2 = bar((1:numel(m_show)) + 0.18, eamp_vals, 0.34, 'FaceColor', [0.85 0.38 0.22]);
    ylabel('Mean E_{amp} (lower is better)'); ylim([0 0.42]);
    grid on;
    set(gca, 'XTick', 1:numel(m_show), 'XTickLabel', m_lbl, 'XTickLabelRotation', 20, 'FontSize', 9);
    title(sprintf('%s | RAW SNR=%+.1f dB, B300 err=%+.2f oct', ...
        titles{i}, J.W_deep.RAW_R0.SNR_img_dB_mean, J.noise_acc.DSA_B300_log2err), 'FontSize', 10);
    if i == 1, legend([b1, b2], {'E1 (Waveform)', 'E_{amp} (Amplitude Err)'}, 'Location', 'northwest', 'FontSize', 8); end
end
exportgraphics(fig1, fullfile(fig_dir, 'fig1_synth_u6_u9_full_honest.png'), 'Resolution', 300);
close(fig1);

%% ========================================================================
%  Fig 2: Volve 4-25 Hz PSAE Attribution Control (Routes A-D)
%  ========================================================================
J25 = jsondecode(fileread(fullfile(res_dir, 'volve25_psae_ctrl.json')));
f25 = J25.f(:)';
fig2 = figure('Visible', 'off', 'Position', [40 40 1550 480], 'Color', 'w');

subplot(1, 3, 1); hold on; grid on;
semilogy(f25, J25.S_anc_dsa / max(J25.S_anc_dsa), '-b', 'LineWidth', 1.8, 'DisplayName', sprintf('With (2\\pif)^2 (\\beta=%+.4f /Hz)', J25.beta_anc_dsa));
semilogy(f25, J25.S_anc_nw2 / max(J25.S_anc_nw2), '-r', 'LineWidth', 1.8, 'DisplayName', sprintf('No (2\\pif)^2 [NoW2] (\\beta=%+.4f /Hz)', J25.beta_anc_nw2));
xlabel('Frequency (Hz)'); ylabel('Normalized Shallow [1000-1600m] RMS');
title('(a) Shallow Anchor Spectrum: (2\pif)^2 vs NoW2'); xlim([4 25]); ylim([1e-2 2]);
legend('Location', 'southwest', 'FontSize', 8);

subplot(1, 3, 2); hold on; grid on;
skeys = {'RAW_R0', 'DSA_UQW1b_Unanc_25', 'DSA_UQW1b_NoW2_Rel_25', 'DSA_UQW1b_PAnc_25', 'DSA_UQW1b_Decon_25'};
slbls = {'RAW\_R0', 'Route A: Unanc', 'Route B: NoW2\_Rel (Zero PSAE)', 'Route C: PSAE (1-param)', 'Route D: 6Hz Decon'};
cols  = {'--k', ':m', '-.c', '-r', '-g'};
for k = 1:numel(skeys)
    sv = J25.spec_res.(skeys{k});
    semilogy(f25, sv / max(sv), cols{k}, 'LineWidth', 1.7, 'DisplayName', sprintf('%s (sl=%+.4f)', slbls{k}, J25.res_slope.(skeys{k})));
end
xlabel('Frequency (Hz)'); ylabel('Normalized Reservoir [2600-3200m] RMS');
title('(b) Reservoir Spectrum Across Routes A-D'); xlim([4 25]); ylim([5e-3 2]);
legend('Location', 'southwest', 'FontSize', 7.5);

subplot(1, 3, 3); hold on; grid on;
rkeys = {'RAW_R0', 'DSA_R0_NoW2', 'DSA_UQW1b_Unanc_25', 'DSA_UQW1b_NoW2_Rel_25', 'DSA_UQW1b_NoW2_Abs_25', 'DSA_UQW1b_PAnc_25', 'DSA_UQW1b_Decon_25'};
rlbls = {'RAW\_R0', 'R0\_NoW2', 'A:Unanc', 'B2:NoW2\_Rel', 'B1:NoW2\_Abs', 'C:PSAE', 'D:Decon'};
vals = zeros(numel(rkeys), 3);
for k = 1:numel(rkeys)
    vals(k, 1) = J25.methods.(rkeys{k}).F5_R0band_mean;
    vals(k, 2) = J25.methods.(rkeys{k}).F5_Broad_mean;
    vals(k, 3) = J25.methods.(rkeys{k}).F14_R0band_mean;
end
bb = bar(vals);
bb(1).DisplayName = 'F-5 R0-band'; bb(2).DisplayName = 'F-5 Broad'; bb(3).DisplayName = 'Blind F-14';
set(gca, 'XTick', 1:numel(rkeys), 'XTickLabel', rlbls, 'XTickLabelRotation', 22, 'FontSize', 8.5);
ylabel('Well-Tie Correlation'); title('(c) Dual-Well Correlation Across Routes A-D');
legend('Location', 'northwest', 'FontSize', 8); ylim([0.15 0.72]);

exportgraphics(fig2, fullfile(fig_dir, 'fig2_volve25_psae_vs_now2.png'), 'Resolution', 300);
close(fig2);

%% ========================================================================
%  Fig 3: Volve 4-45 Hz True INV3_DSA Active Gain, F-5(m), & Dual-Well F-14
%  ========================================================================
J45 = jsondecode(fileread(fullfile(res_dir, 'volve45_closure.json')));
Jwp1 = jsondecode(fileread(fullfile(res_dir, 'wp1_diagnostics.json')));
f45 = J45.f(:)'; m45 = J45.m(:)';

fig3 = figure('Visible', 'off', 'Position', [40 40 1500 950], 'Color', 'w');

subplot(2, 2, 1); hold on; grid on;
plot(f45, 10*log10(max(J45.snr_med_b300_dsa_full, 1e-4)), '-r', 'LineWidth', 2.0, 'DisplayName', 'True Unaliased INV3\_DSA (df=0.25Hz)');
plot(f45, 10*log10(max(J45.snr_med_b300_raw_full, 1e-4)), '--b', 'LineWidth', 1.6, 'DisplayName', 'True Unaliased INV3\_RAW\_FULL');
plot(f45, 10*log10(max(J45.snr_med_b300_s4h, 1e-4)), ':k', 'LineWidth', 1.6, 'DisplayName', 'Stage 4H Aliased INV3\_S4H (df=0.50Hz)');
yline(0, '-k', 'HandleVisibility', 'off');
xlabel('Frequency (Hz)'); ylabel('Reservoir Median SNR (dB)');
title('(a) Volve 4-45 Hz Unaliased vs Aliased B=300 SNR'); xlim([6 44]); ylim([-32 22]);
legend('Location', 'southwest', 'FontSize', 8);

subplot(2, 2, 2); hold on; grid on;
s_base = J45.spec_res.DSA_R0_Eq(:)';
plot(f45, J45.spec_res.DSA_UQW1b_TrueDSAiv_PAnc_25(:)' ./ max(s_base, 1e-30), '-r', 'LineWidth', 2.0, ...
    'DisplayName', sprintf('TrueDSAiv PAnc 25Hz (G_{HF}=%.2fx)', J45.G_HF_28_40.DSA_UQW1b_TrueDSAiv_PAnc_25));
plot(f45, J45.spec_res.DSA_UQW1b_S4Hiv_PAnc_25(:)' ./ max(s_base, 1e-30), '--k', 'LineWidth', 1.6, ...
    'DisplayName', sprintf('Old Dormant S4Hiv (G_{HF}=%.2fx)', J45.G_HF_28_40.DSA_UQW1b_S4Hiv_PAnc_25));
plot(f45, J45.spec_res.DSA_W1b_TrueDSAiv_PAnc_25(:)' ./ max(s_base, 1e-30), '-.b', 'LineWidth', 1.5, ...
    'DisplayName', sprintf('W1b TrueDSAiv (G_{HF}=%.2fx)', J45.G_HF_28_40.DSA_W1b_TrueDSAiv_PAnc_25));
plot(f45, J45.spec_res.DSA_IQ70_PAnc_25(:)' ./ max(s_base, 1e-30), ':g', 'LineWidth', 1.6, ...
    'DisplayName', sprintf('IQ70 Det (G_{HF}=%.2fx)', J45.G_HF_28_40.DSA_IQ70_PAnc_25));
yline(1.0, '-k', 'HandleVisibility', 'off'); xline(28, ':k', 'HandleVisibility', 'off'); xline(40, ':k', 'HandleVisibility', 'off');
xlabel('Frequency (Hz)'); ylabel('Spectral Gain Ratio vs DSA\_R0\_Eq');
title('(b) Active Broadband Q-Compensation Gain (28-40 Hz)'); xlim([6 42]);
legend('Location', 'northwest', 'FontSize', 8);

subplot(2, 2, 3); hold on; grid on;
plot(m45, J45.methods.RAW_R0.F5_fix * ones(size(m45)), '--k', 'LineWidth', 1.5, 'DisplayName', sprintf('RAW\_R0 (%.4f)', J45.methods.RAW_R0.F5_fix_mean));
plot(m45, J45.methods.DSA_R0_Eq.F5_fix * ones(size(m45)), ':c', 'LineWidth', 1.6, 'DisplayName', sprintf('DSA\_R0\_Eq (%.4f)', J45.methods.DSA_R0_Eq.F5_fix_mean));
plot(m45, J45.methods.RAW_UQW1b_S4H_Unanc_100.F5_fix, '-^m', 'LineWidth', 1.5, 'DisplayName', sprintf('Stage4H 100Hz (%.4f)', J45.methods.RAW_UQW1b_S4H_Unanc_100.F5_fix_mean));
plot(m45, J45.methods.DSA_IQ70_PAnc_25.F5_fix, '-dg', 'LineWidth', 1.5, 'DisplayName', sprintf('DSA\_IQ70 (%.4f)', J45.methods.DSA_IQ70_PAnc_25.F5_fix_mean));
plot(m45, J45.methods.DSA_W1b_TrueDSAiv_PAnc_25.F5_fix, '-.ob', 'LineWidth', 1.6, 'DisplayName', sprintf('DSA\_W1b TrueDSAiv (%.4f)', J45.methods.DSA_W1b_TrueDSAiv_PAnc_25.F5_fix_mean));
plot(m45, J45.methods.DSA_UQW1b_TrueDSAiv_PAnc_25.F5_fix, '-sr', 'LineWidth', 2.2, 'DisplayName', sprintf('DSA\_UQW1b TrueDSAiv (%.4f)', J45.methods.DSA_UQW1b_TrueDSAiv_PAnc_25.F5_fix_mean));
plot(m45, J45.methods.DSA_UQW1b_S4Hiv_PAnc_25.F5_fix, ':ok', 'LineWidth', 1.4, 'DisplayName', sprintf('Old S4Hiv Dormant (%.4f)', J45.methods.DSA_UQW1b_S4Hiv_PAnc_25.F5_fix_mean));
xlabel('Q Perturbation Multiplier m'); ylabel('F-5 Broadband Tie (8-40 Hz)');
title('(c) Volve 4-45 Hz F-5 Well-Tie vs Q Perturbation m');
legend('Location', 'southwest', 'FontSize', 7.5);

subplot(2, 2, 4); hold on; grid on;
lags14 = Jwp1.D1_F14.lags_m;
plot(lags14, Jwp1.D1_F14.methods.RAW_R0.cc_kfix_curve, '--k', 'LineWidth', 1.4, 'DisplayName', 'RAW\_R0 (8-40Hz, peak -60m)');
plot(lags14, Jwp1.D1_F14.methods.DSA_UQW1b_PAnc_25.cc_kfix_curve, '-r', 'LineWidth', 1.8, 'DisplayName', 'DSA\_UQW1b (8-40Hz, peak -60m)');
plot(lags14, Jwp1.D1_F14.methods.DSA_UQW1b_PAnc_25.cc_k25_curve, '-b', 'LineWidth', 1.8, 'DisplayName', 'DSA\_UQW1b (k25: 3.5-12.9 cyc/km)');
xline(-30, ':k', 'HandleVisibility', 'off'); xline(30, ':k', 'HandleVisibility', 'off'); yline(0, '-k', 'HandleVisibility', 'off');
xlabel('Depth Lag (m)'); ylabel('Blind Well F-14 Correlation');
title('(d) Blind Well F-14: -60m True Tie vs Truncated \pm30m Window');
legend('Location', 'northeast', 'FontSize', 8);

exportgraphics(fig3, fullfile(fig_dir, 'fig3_volve45_true_inv3_and_dual_wells.png'), 'Resolution', 300);
close(fig3);

%% ========================================================================
%  Fig 4: Volve 4-45 Hz 2D Migration Sections (True Active Confluence)
%  ========================================================================
S45 = load('D:\tmp_uqrtm_closure_20261125\volve45_closure_imgs.mat');
z_in = S45.z_in(:); x_in = S45.x_in(:)'; dx = median(diff(z_in));
nzi = numel(z_in); nxi = numel(x_in);
well_dir = 'D:\ChatGPT\研究\.claude\worktrees\study-research-docs-plan-c4d4c6\地震衰减补偿研究\uq_rtm_20261117\cache\fld';
w5  = core.load_well(fullfile(well_dir, 'well_F-5.bin'),  z_in, x_in);
w14 = core.load_well(fullfile(well_dir, 'well_F-14.bin'), z_in, x_in);

panels = {'RAW_R0', 'RAW_UQW1b_S4H_Unanc_100', 'DSA_R0_Eq', 'DSA_UQW1b_TrueDSAiv_PAnc_25'};
ptitles = { ...
    sprintf('(a) RAW\\_R0 | F5=%.4f | Smile=%.1f%%', J45.methods.RAW_R0.F5_fix_mean, J45.methods.RAW_R0.Smile_Frac*100), ...
    sprintf('(b) Stage 4H Unanc 100Hz | F5=%.4f | Smile=%.1f%%', J45.methods.RAW_UQW1b_S4H_Unanc_100.F5_fix_mean, J45.methods.RAW_UQW1b_S4H_Unanc_100.Smile_Frac(3)*100), ...
    sprintf('(c) DSA\\_R0\\_Eq (Zero-Q Base) | F5=%.4f | Smile=%.1f%%', J45.methods.DSA_R0_Eq.F5_fix_mean, J45.methods.DSA_R0_Eq.Smile_Frac*100), ...
    sprintf('(d) DSA\\_UQW1b\\_TrueDSAiv\\_PAnc\\_25 | F5=%.4f | G_{HF}=%.2fx | Smile=%.1f%%', ...
        J45.methods.DSA_UQW1b_TrueDSAiv_PAnc_25.F5_fix_mean, J45.G_HF_28_40.DSA_UQW1b_TrueDSAiv_PAnc_25, J45.methods.DSA_UQW1b_TrueDSAiv_PAnc_25.Smile_Frac(3)*100)};

fig4 = figure('Visible', 'off', 'Position', [40 40 1600 1000], 'Color', 'w');
iz_p = find(z_in >= 1200 & z_in <= 3400);
ix_p = find(x_in >= -2800 & x_in <= 1800);
for ip = 1:4
    M = S45.Img_s.(panels{ip});
    A = core.prep_image(M(:, min(3, size(M,2))), nzi, nxi, z_in, dx, cfg.bg_sigma_m, cfg.z_mute_m);
    Ap = A(iz_p, ix_p);
    cl = prctile(abs(Ap(:)), 98.5);
    subplot(2, 2, ip);
    imagesc(x_in(ix_p), z_in(iz_p), Ap, [-cl cl]);
    colormap(gca, redblue_cmap(256)); hold on;
    plot(w5.sw, w5.zw, '-k', 'LineWidth', 1.6); text(median(w5.sw)+80, min(w5.zw)-60, 'F-5', 'Color', 'k', 'FontWeight', 'bold', 'FontSize', 9);
    plot(w14.sw, w14.zw, '--g', 'LineWidth', 1.6); text(median(w14.sw)+80, min(w14.zw)-60, 'F-14', 'Color', [0 0.5 0], 'FontWeight', 'bold', 'FontSize', 9);
    xlabel('Inline X (m)'); ylabel('Depth Z (m)'); title(ptitles{ip}, 'FontSize', 10);
end
exportgraphics(fig4, fullfile(fig_dir, 'fig4_volve45_2d_sections_true_closure.png'), 'Resolution', 300);
close(fig4);
fprintf('Generated all 4 publication-grade 300-DPI figures in %s\n', fig_dir);
end

function cmap = redblue_cmap(n)
m = floor(n/2);
r = [linspace(0, 1, m)'; ones(n-m, 1)];
g = [linspace(0, 1, m)'; linspace(1, 0, n-m)'];
b = [ones(m, 1); linspace(1, 0, n-m)'];
cmap = [r g b];
end
