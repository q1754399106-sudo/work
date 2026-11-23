function make_unified_figures()
% Publication-grade 300-DPI composite figures for uq_rtm_unified_20261124 (Geophysics / GJI style).
% Reads JSON results from results/ and 2D image snapshots from D:\tmp_uqrtm_unified_20261124\.
% Exports 4 PNG figures into figures/.

here = fileparts(mfilename('fullpath'));
root = fileparts(here);
res_dir = fullfile(root, 'results');
fig_dir = fullfile(root, 'figures');
if ~isfolder(fig_dir), mkdir(fig_dir); end
scratch_dir = 'D:\tmp_uqrtm_unified_20261124';

wt_root = fileparts(root);
addpath(fullfile(wt_root, 'uq_rtm_20261117', 'scripts'));
core = uq_rtm_unified_core();
cfg  = core.default_cfg('field_obc');

set(0, 'DefaultAxesFontName', 'Microsoft YaHei', ...
       'DefaultTextFontName', 'Microsoft YaHei', ...
       'DefaultAxesFontSize', 9.5, ...
       'DefaultLineLineWidth', 1.6);

% Fig 1: Tier 1 Synthetic Multi-Physics Confluence (U1-U5)
if isfile(fullfile(res_dir, 'synth_U1.json')) && isfile(fullfile(res_dir, 'synth_U4.json'))
    plot_fig1_synth(res_dir, fullfile(fig_dir, 'fig1_synth_theorem3_confluence.png'));
    fprintf('Exported fig1_synth_theorem3_confluence.png\n');
end

% Fig 2: Tier 2 Volve 4-25 Hz Factorial Confluence
if isfile(fullfile(res_dir, 'volve25_unified.json'))
    plot_fig2_volve25(res_dir, fullfile(fig_dir, 'fig2_volve25_factorial_closure.png'));
    fprintf('Exported fig2_volve25_factorial_closure.png\n');
end

% Fig 3 & Fig 4: Tier 3 Volve 4-45 Hz Broadband Confluence & 2D Sections
if isfile(fullfile(res_dir, 'volve45_unified.json')) && isfile(fullfile(scratch_dir, 'volve45_final_imgs.mat'))
    plot_fig3_volve45(res_dir, scratch_dir, core, cfg, fullfile(fig_dir, 'fig3_volve45_broadband_confluence.png'));
    fprintf('Exported fig3_volve45_broadband_confluence.png\n');
    plot_fig4_volve45_2d(scratch_dir, core, cfg, fullfile(fig_dir, 'fig4_volve45_2d_sections.png'));
    fprintf('Exported fig4_volve45_2d_sections.png\n');
end
end

function plot_fig1_synth(res_dir, out_png)
U1 = jsondecode(fileread(fullfile(res_dir, 'synth_U1.json')));
U2 = jsondecode(fileread(fullfile(res_dir, 'synth_U2.json')));
U3 = jsondecode(fileread(fullfile(res_dir, 'synth_U3.json')));
U4 = jsondecode(fileread(fullfile(res_dir, 'synth_U4.json')));
Us = {U1, U2, U3};

fig = figure('Visible', 'off', 'Position', [60, 60, 1280, 860], 'Color', 'w');
tl = tiledlayout(2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

% (a) Theorem 3 Noise Acc across U1-U4
nexttile(tl);
models_n = {'U1', 'U2', 'U3', 'U4'};
all4 = {U1, U2, U3, U4};
M_acc = zeros(4, 4);
for i = 1:4
    na = all4{i}.noise_acc;
    M_acc(i, :) = [na.RAW_B1_log2err, na.DSA_B1_log2err, na.RAW_B300_log2err, na.DSA_B300_log2err];
end
b = bar(M_acc);
b(1).FaceColor = [0.82 0.25 0.25]; b(2).FaceColor = [0.95 0.58 0.35];
b(3).FaceColor = [0.30 0.60 0.85]; b(4).FaceColor = [0.10 0.48 0.28];
hold on; yline(1, '--k', 'LineWidth', 1); yline(-1, '--k', 'LineWidth', 1); yline(0, ':k');
set(gca, 'XTickLabel', models_n); ylabel('log_2(iv_{est} / iv_{true})');
title('(a) 定理 3 盲估噪精度（l30 横向相关噪，目标 \pm1 档）');
legend({'RAW B=1', 'DSA B=1', 'RAW B=300', 'DSA B=300 (推荐)'}, 'Location', 'southwest', 'FontSize', 8);
grid on;

% (b) Full factorial E1 progression on U1-U3
nexttile(tl);
keys = {'RAW_UQW1b_100', 'RAW_UQW1_25', 'RAW_UQW1b_25', 'DSA_UQW1_25', 'DSA_UQW1b_25'};
lbls = {'RAW UQW1b\n(100Hz)', 'RAW UQW1\n(B=1,25Hz)', 'RAW UQW1b\n(B=300,25Hz)', 'DSA UQW1\n(B=1,25Hz)', 'DSA-UQW1b\n(合流主推)'};
E_mat = zeros(3, numel(keys));
for i = 1:3
    for j = 1:numel(keys)
        E_mat(i, j) = Us{i}.W_deep.(keys{j}).E1_mean;
    end
end
b2 = bar(E_mat');
b2(1).FaceColor = [0.22 0.45 0.72]; b2(2).FaceColor = [0.85 0.45 0.15]; b2(3).FaceColor = [0.16 0.62 0.36];
set(gca, 'XTick', 1:numel(keys), 'XTickLabel', strrep(lbls, '\n', ' ')); xtickangle(18);
ylabel('深部窗 5-m 平均 E1'); title('(b) U1–U3 四因子逐级合流 E1 增益（v_{band}=25 Hz）');
legend({'U1 气烟囱+l30', 'U2 22°逆断层+l30', 'U3 单边拖缆+l30'}, 'Location', 'southeast', 'FontSize', 8);
grid on;

% (c) m-robustness curve (U1-U3 average)
nexttile(tl);
ms = U1.m(:)';
avg_curve = @(k) mean(cell2mat(cellfun(@(u) u.W_deep.(k).E1(:)', Us, 'UniformOutput', false)'), 1);
plot(ms, avg_curve('DSA_UQW1b_25'), '-o', 'Color', [0.05 0.50 0.25], 'LineWidth', 2.2, 'MarkerFaceColor', [0.05 0.50 0.25]); hold on;
plot(ms, avg_curve('DSA_W1b_25'),   '--s', 'Color', [0.18 0.45 0.78]);
plot(ms, avg_curve('DSA_IQ70_25'),  ':^',  'Color', [0.80 0.20 0.20]);
plot(ms, avg_curve('RAW_UQW1_25'),  '-.d', 'Color', [0.60 0.40 0.15]);
r0_avg = mean([U1.W_deep.RAW_R0.E1, U2.W_deep.RAW_R0.E1, U3.W_deep.RAW_R0.E1]);
yline(r0_avg, '--k', 'RAW R0', 'LineWidth', 1.2);
xlabel('Q 模型乘子 m = Q_{used}/Q_{true}'); ylabel('U1–U3 平均 E1');
title('(c) 同门控过补偿稳健性（m=0.5 贝叶斯防崩盘）');
legend({'DSA-UQW1b (推荐)', 'DSA-W1b', 'DSA-IQ70', 'RAW-UQW1 (B=1)'}, 'Location', 'southeast', 'FontSize', 8);
grid on;

% (d) U4 Volve Twin 5-well tie
nexttile(tl);
M4 = U4.methods;
plot(ms, M4.DSA_UQW1b_PAnc_25.tie, '-o', 'Color', [0.05 0.50 0.25], 'LineWidth', 2.2, 'MarkerFaceColor', [0.05 0.50 0.25]); hold on;
plot(ms, M4.DSA_W1b_PAnc_25.tie,   '--s', 'Color', [0.18 0.45 0.78]);
plot(ms, M4.DSA_IQ70_PAnc_25.tie,  ':^',  'Color', [0.80 0.20 0.20]);
plot(ms, M4.RAW_UQW1_25.tie,       '-.x', 'Color', [0.70 0.25 0.65]);
yline(M4.RAW_R0.mean, '--k', 'RAW R0', 'LineWidth', 1.3);
yline(M4.DSA_R0_Eq.mean, ':b', 'DSA R0 Eq', 'LineWidth', 1.2);
xlabel('Q 模型乘子 m'); ylabel('5 井平均固定带相关 [4–18 cyc/km]');
title('(d) U4 Volve 孪生体（平盒 W=1 + l30 + 单边拖缆）');
legend({'DSA-UQW1b-PAnc (25Hz)', 'DSA-W1b-PAnc', 'DSA-IQ70-PAnc', 'RAW-UQW1 (未锚定)'}, 'Location', 'southeast', 'FontSize', 8);
grid on;

exportgraphics(fig, out_png, 'Resolution', 300);
close(fig);
end

function plot_fig2_volve25(res_dir, out_png)
D = jsondecode(fileread(fullfile(res_dir, 'volve25_unified.json')));
M = D.methods; ms = D.m(:)'; f = D.f(:)';

fig = figure('Visible', 'off', 'Position', [60, 60, 1280, 860], 'Color', 'w');
tl = tiledlayout(2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

% (a) F-5 R0-band & Broadband curves
nexttile(tl);
plot(ms, M.DSA_UQW1b_PAnc_25.F5_R0band,  '-o', 'Color', [0.05 0.50 0.25], 'LineWidth', 2.2, 'MarkerFaceColor', [0.05 0.50 0.25]); hold on;
plot(ms, M.DSA_UQW1b_PAnc_Amp.F5_R0band, '-s', 'Color', [0.12 0.45 0.78], 'LineWidth', 1.8);
plot(ms, M.DSA_W1b_PAnc_25.F5_R0band,    '--^', 'Color', [0.85 0.50 0.10]);
plot(ms, M.RAW_UQW1_Unanc_25.F5_R0band,  '-.x', 'Color', [0.80 0.20 0.20]);
yline(M.RAW_R0.F5_R0band, '--k', 'RAW R0 (0.626)', 'LineWidth', 1.3);
xlabel('Q 模型乘子 m'); ylabel('F-5 井旁相关（R0 带）');
title('(a) Volve 4–25 Hz 实测 F-5 井震标定（合流 vs 消融）');
legend({'DSA-UQW1b-PAnc (25Hz)', 'DSA-UQW1b-PAnc (Amp)', 'DSA-W1b-PAnc (25Hz)', 'RAW-UQW1 (未锚定25Hz)'}, 'Location', 'southeast', 'FontSize', 8);
grid on;

% (b) Shallow vs Reservoir normalized log-spectra
nexttile(tl);
fm = (f >= 5.5 & f <= 23.0);
norm_s = @(v) v(fm) / max(median(v(f>=12 & f<=17)), 1e-30);
semilogy(f(fm), norm_s(D.S_anc_dsa), '-k', 'LineWidth', 1.8); hold on;
semilogy(f(fm), norm_s(D.spec_res.RAW_UQW1b_Unanc_25), '--r', 'LineWidth', 1.6);
semilogy(f(fm), norm_s(D.spec_res.DSA_R0_Eq), ':b', 'LineWidth', 1.8);
semilogy(f(fm), norm_s(D.spec_res.DSA_UQW1b_PAnc_25), '-', 'Color', [0.05 0.50 0.25], 'LineWidth', 2.2);
xlabel('频率 f (Hz)'); ylabel('归一化 RMS 振幅谱');
title('(b) 浅层蓝倾与储层谱拉平闭环（PSAE）');
legend({'浅层锚点 DSA R0', sprintf('未锚定 RAW UQW1b (sl=%+.3f)', D.res_slope.RAW_UQW1b_Unanc_25), ...
        sprintf('均衡未补 DSA R0 Eq (sl=%+.3f)', D.res_slope.DSA_R0_Eq), ...
        sprintf('DSA-UQW1b-PAnc (sl=%+.4f)', D.res_slope.DSA_UQW1b_PAnc_25)}, 'Location', 'northwest', 'FontSize', 7.8);
grid on;

% (c) Reservoir SNR curves B=1 vs RAW B=300 vs DSA B=300
nexttile(tl);
plot(f, 10*log10(max(D.snr_med_b1, 1e-3)), '--r', 'LineWidth', 1.6); hold on;
plot(f, 10*log10(max(D.snr_med_b300_raw, 1e-3)), '-b', 'LineWidth', 1.7);
plot(f, 10*log10(max(D.snr_med_b300_dsa, 1e-3)), '-', 'Color', [0.05 0.50 0.25], 'LineWidth', 2.1);
yline(0, ':k');
xlabel('频率 f (Hz)'); ylabel('储层中位成像域 SNR (dB)');
title('(c) 实测储层盲估噪对比（B=1 vs RAW B=300 vs DSA B=300）');
legend({'RAW B=1 (Stage 4)', 'RAW B=300', 'DSA B=300 (同域)'}, 'Location', 'northeast', 'FontSize', 8);
grid on;

% (d) Smile_Frac suppression
nexttile(tl);
sk = {'RAW_R0', 'RAW_UQW1b_Unanc_25', 'RAW_UQW1b_PAnc_25', 'DSA_R0', 'DSA_UQW1b_PAnc_25', 'DSA_UQW1b_PAnc_Amp'};
slbl = {'RAW R0', 'RAW Unanc', 'RAW PAnc', 'DSA R0', 'DSA PAnc 25', 'DSA PAnc Amp'};
vals = zeros(1, numel(sk));
for i = 1:numel(sk)
    v = M.(sk{i}).Smile_Frac;
    vals(i) = v(min(3, numel(v))) * 100;
end
bh = bar(vals); bh.FaceColor = 'flat';
bh.CData = [0.75 0.25 0.25; 0.85 0.30 0.30; 0.88 0.50 0.30; 0.20 0.55 0.80; 0.05 0.52 0.28; 0.12 0.62 0.45];
set(gca, 'XTick', 1:numel(sk), 'XTickLabel', slbl); xtickangle(18);
ylabel('上覆 [1600–2400 m] 陡倾画弧能量占比 (%)');
title('(d) DSA 活排列防假频门控对上覆假频画弧的 ~10\times 压制');
grid on;

exportgraphics(fig, out_png, 'Resolution', 300);
close(fig);
end

function plot_fig3_volve45(res_dir, scratch_dir, core, cfg, out_png)
D = jsondecode(fileread(fullfile(res_dir, 'volve45_unified.json')));
S = load(fullfile(scratch_dir, 'volve45_final_imgs.mat'), 'Img_s', 'z_in', 'x_in');
M = D.methods; ms = D.m(:)'; f = D.f(:)';

fig = figure('Visible', 'off', 'Position', [60, 60, 1320, 880], 'Color', 'w');
tl = tiledlayout(2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

% (a) F-5 Broadband Fixed-Band [8-40 Hz] tie vs m
nexttile(tl);
plot(ms, M.DSA_UQW1b_PAnc_25.F5_fix,   '-o', 'Color', [0.05 0.50 0.25], 'LineWidth', 2.3, 'MarkerFaceColor', [0.05 0.50 0.25]); hold on;
plot(ms, M.DSA_UQW1b_PAnc_100.F5_fix,  '-d', 'Color', [0.55 0.20 0.72], 'LineWidth', 1.8);
plot(ms, M.DSA_UQW1b_PAnc_Amp.F5_fix,  '-s', 'Color', [0.12 0.48 0.80], 'LineWidth', 1.8);
plot(ms, M.RAW_UQW1b_Unanc_100.F5_fix, '--k', 'LineWidth', 1.6);
plot(ms, M.DSA_W1b_PAnc_25.F5_fix,     ':^', 'Color', [0.85 0.48 0.10], 'LineWidth', 1.5);
plot(ms, M.DSA_IQ70_PAnc_25.F5_fix,    ':v', 'Color', [0.82 0.18 0.18], 'LineWidth', 1.5);
yline(M.RAW_R0.F5_fix, '--', 'RAW R0 (0.474)', 'Color', [0.35 0.35 0.35], 'LineWidth', 1.4);
xlabel('Q 模型乘子 m'); ylabel('F-5 固定宽带相关 tie_{fix} (8–40 Hz)');
title('(a) Volve 4–45 Hz 宽频 F-5 井震标定（合流 vs Stage 4H）');
legend({'DSA-UQW1b-PAnc (25Hz, 主推)', 'DSA-UQW1b-PAnc (100Hz)', 'DSA-UQW1b-PAnc (Amp)', ...
        'Stage 4H RAW-UQW1b (100Hz)', 'DSA-W1b-PAnc (25Hz)', 'DSA-IQ70-PAnc (25Hz)'}, 'Location', 'southeast', 'FontSize', 7.5);
grid on;

% (b) Broadband 4-45 Hz spectral flattening
nexttile(tl);
fm = (f >= 6.0 & f <= 41.0);
nrm = @(v) v(fm) / max(median(v(f>=15 & f<=22)), 1e-30);
semilogy(f(fm), nrm(D.S_anc_dsa), '-k', 'LineWidth', 1.7); hold on;
semilogy(f(fm), nrm(D.spec_res.RAW_UQW1b_Unanc_100), '--r', 'LineWidth', 1.6);
semilogy(f(fm), nrm(D.spec_res.DSA_R0_Eq), ':b', 'LineWidth', 1.7);
semilogy(f(fm), nrm(D.spec_res.DSA_UQW1b_PAnc_25), '-', 'Color', [0.05 0.50 0.25], 'LineWidth', 2.2);
xlabel('频率 f (Hz)'); ylabel('归一化 RMS 振幅谱');
title('(b) 4–45 Hz 宽频 (2\pif)^2 蓝倾消除与储层频谱均衡');
legend({'浅层锚点 DSA R0', sprintf('Stage 4H RAW UQW1b (sl=%+.3f)', D.res_slope.RAW_UQW1b_Unanc_100), ...
        sprintf('DSA R0 Eq (sl=%+.3f)', D.res_slope.DSA_R0_Eq), ...
        sprintf('DSA-UQW1b-PAnc 25Hz (sl=%+.4f)', D.res_slope.DSA_UQW1b_PAnc_25)}, 'Location', 'northwest', 'FontSize', 7.6);
grid on;

% (c) Overburden Smile_Frac bar comparison on 4-45 Hz
nexttile(tl);
sk = {'RAW_R0', 'RAW_UQW1b_Unanc_100', 'RAW_UQW1b_PAnc_25', 'DSA_R0', 'DSA_UQW1b_PAnc_25', 'DSA_UQW1b_PAnc_100'};
slbl = {'RAW R0', 'Stage4H UQW1b', 'RAW PAnc 25', 'DSA R0', 'DSA PAnc 25', 'DSA PAnc 100'};
vals = zeros(1, numel(sk));
for i = 1:numel(sk)
    v = M.(sk{i}).Smile_Frac;
    vals(i) = v(min(3, numel(v))) * 100;
end
bh = bar(vals); bh.FaceColor = 'flat';
bh.CData = [0.72 0.22 0.22; 0.88 0.18 0.18; 0.85 0.52 0.28; 0.18 0.52 0.80; 0.05 0.52 0.28; 0.52 0.22 0.72];
set(gca, 'XTick', 1:numel(sk), 'XTickLabel', slbl); xtickangle(18);
ylabel('上覆 [1600–2400 m] 陡倾假频画弧占比 (%)');
title('(c) 4–45 Hz 宽频上覆假频画弧压制（Stage 4H 83% \rightarrow DSA）');
grid on;

% (d) F-5 Well-log waveform tie overlay
nexttile(tl);
nzi = numel(S.z_in); nxi = numel(S.x_in); dx = median(diff(S.z_in));
well_dir = 'D:\ChatGPT\研究\.claude\worktrees\study-research-docs-plan-c4d4c6\地震衰减补偿研究\uq_rtm_20261117\cache\fld';
w5 = core.load_well(fullfile(well_dir, 'well_F-5.bin'), S.z_in, S.x_in);
kfix = D.kfix(:)';
rb = core.bp_z(w5.rw, w5.dzw, kfix(1), kfix(2)); rb = rb / max(rms(rb), 1e-30);
A_r0  = core.prep_image(S.Img_s.RAW_R0, nzi, nxi, S.z_in, dx, cfg.bg_sigma_m, cfg.z_mute_m);
A_s4h = core.prep_image(S.Img_s.RAW_UQW1b_Unanc_100(:,3), nzi, nxi, S.z_in, dx, cfg.bg_sigma_m, cfg.z_mute_m);
A_h25 = core.prep_image(S.Img_s.DSA_UQW1b_PAnc_25(:,3), nzi, nxi, S.z_in, dx, cfg.bg_sigma_m, cfg.z_mute_m);
[c0, l0, ~, t0] = core.well_tie_full(A_r0,  w5, rb, dx, kfix(1), kfix(2), 20, 30);
[c1, l1, ~, t1] = core.well_tie_full(A_s4h, w5, rb, dx, kfix(1), kfix(2), 20, 30);
[c2, l2, ~, t2] = core.well_tie_full(A_h25, w5, rb, dx, kfix(1), kfix(2), 20, 30);
t0 = circshift(t0, round(l0/w5.dzw)); t0 = t0 / max(rms(t0), 1e-30);
t1 = circshift(t1, round(l1/w5.dzw)); t1 = t1 / max(rms(t1), 1e-30);
t2 = circshift(t2, round(l2/w5.dzw)); t2 = t2 / max(rms(t2), 1e-30);
plot(w5.zw, rb, 'k', 'LineWidth', 1.8); hold on;
plot(w5.zw, t0 + 3.2, '--', 'Color', [0.45 0.45 0.45], 'LineWidth', 1.2);
plot(w5.zw, t1 + 6.4, '-',  'Color', [0.82 0.25 0.25], 'LineWidth', 1.3);
plot(w5.zw, t2 + 9.6, '-',  'Color', [0.05 0.50 0.25], 'LineWidth', 1.6);
yticks([0 3.2 6.4 9.6]);
yticklabels({'F-5 测井反射率', sprintf('RAW R0 (%.3f)', c0), sprintf('Stage4H UQW1b (%.3f)', c1), sprintf('DSA-UQW1b-PAnc (%.3f)', c2)});
xlabel('深度 z (m)'); title('(d) F-5 储层段 [2600–3200 m] 8–40 Hz 井震波形叠合');
xlim([2600 3200]); grid on;

exportgraphics(fig, out_png, 'Resolution', 300);
close(fig);
end

function plot_fig4_volve45_2d(scratch_dir, core, cfg, out_png)
S = load(fullfile(scratch_dir, 'volve45_final_imgs.mat'), 'Img_s', 'res', 'z_in', 'x_in');
z = S.z_in; x = S.x_in; nzi = numel(z); nxi = numel(x); dx = median(diff(z));
well_dir = 'D:\ChatGPT\研究\.claude\worktrees\study-research-docs-plan-c4d4c6\地震衰减补偿研究\uq_rtm_20261117\cache\fld';
w5 = core.load_well(fullfile(well_dir, 'well_F-5.bin'), z, x);

iz = (z >= 800 & z <= 3400); ix = (x >= -2600 & x <= 1600);
zv = z(iz); xv = x(ix);

A_r0  = core.prep_image(S.Img_s.RAW_R0, nzi, nxi, z, dx, cfg.bg_sigma_m, cfg.z_mute_m);
A_s4h = core.prep_image(S.Img_s.RAW_UQW1b_Unanc_100(:,3), nzi, nxi, z, dx, cfg.bg_sigma_m, cfg.z_mute_m);
A_h25 = core.prep_image(S.Img_s.DSA_UQW1b_PAnc_25(:,3), nzi, nxi, z, dx, cfg.bg_sigma_m, cfg.z_mute_m);
A_amp = core.prep_image(S.Img_s.DSA_UQW1b_PAnc_Amp(:,3), nzi, nxi, z, dx, cfg.bg_sigma_m, cfg.z_mute_m);

fig = figure('Visible', 'off', 'Position', [50, 50, 1400, 920], 'Color', 'w');
tl = tiledlayout(2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
M = S.res.methods;

panels = {
    {A_r0,  sprintf('(a) RAW R0 不补偿（4–45 Hz | F-5=%.3f | 画弧=%.1f%%）', M.RAW_R0.F5_fix, M.RAW_R0.Smile_Frac*100)}, ...
    {A_s4h, sprintf('(b) Stage 4H RAW-UQW1b 100Hz（F-5=%.3f | 画弧=%.1f%%）', M.RAW_UQW1b_Unanc_100.F5_fix(3), M.RAW_UQW1b_Unanc_100.Smile_Frac(3)*100)}, ...
    {A_h25, sprintf('(c) 合流主推 DSA-UQW1b-PAnc 25Hz（F-5=%.3f | 画弧=%.1f%%）', M.DSA_UQW1b_PAnc_25.F5_fix(3), M.DSA_UQW1b_PAnc_25.Smile_Frac(3)*100)}, ...
    {A_amp, sprintf('(d) 合流纯幅伴随 DSA-UQW1b-PAnc Amp（F-5=%.3f | 画弧=%.1f%%）', M.DSA_UQW1b_PAnc_Amp.F5_fix(3), M.DSA_UQW1b_PAnc_Amp.Smile_Frac(3)*100)}
};

for ip = 1:4
    nexttile(tl);
    Pimg = panels{ip}{1}(iz, ix);
    sc = prctile(abs(Pimg(zv >= 2400 & zv <= 3250, :)), 97, 'all');
    imagesc(xv, zv, Pimg, [-sc sc]); colormap(gca, gray_red_blue());
    hold on;
    plot(w5.sw, w5.zw, 'y-', 'LineWidth', 2.0);
    rectangle('Position', [-2500, 2600, 4000, 600], 'EdgeColor', [0.1 0.85 0.3], 'LineStyle', '--', 'LineWidth', 1.1);
    xlabel('测线坐标 x (m)'); ylabel('深度 z (m)');
    title(panels{ip}{2}, 'FontSize', 9.5);
end

exportgraphics(fig, out_png, 'Resolution', 300);
close(fig);
end

function cmap = gray_red_blue()
n = 128;
r = [linspace(0.08, 1, n), ones(1, n)]';
g = [linspace(0.22, 1, n), linspace(1, 0.14, n)]';
b = [ones(1, n),           linspace(1, 0.12, n)]';
cmap = [r, g, b];
end
