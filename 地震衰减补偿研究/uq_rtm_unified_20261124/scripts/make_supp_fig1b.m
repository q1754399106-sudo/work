function make_supp_fig1b()
% Supplementary publication figure:
%   fig1b_synth_midSNR_and_snr_discovery.png
% Shows:
%   (a) Theorem 3 in DSA domain on U1b-U4b (+2 to +5 dB realistic l30 image SNR)
%   (b) Four-factor monotonic confluence on U1b-U3b (E1 progression)
%   (c) U4b Volve full-pipeline twin well-tie recovery (RAW R0 -> RAW UQW1 drop -> DSA-UQW1b-PAnc/EnvAnc)
%   (d) Landmark Volve 4-45 Hz SNR discovery: Stage 4H df=0.50 Hz wraparound vs full df=0.25 Hz DSA B=300!

here = fileparts(mfilename('fullpath'));
root = fileparts(here);
res_dir = fullfile(root, 'results');
fig_dir = fullfile(root, 'figures');
if ~isfolder(fig_dir), mkdir(fig_dir); end

set(0, 'DefaultAxesFontName', 'Microsoft YaHei', ...
       'DefaultTextFontName', 'Microsoft YaHei', ...
       'DefaultAxesFontSize', 9.5, ...
       'DefaultLineLineWidth', 1.6);

U1 = jsondecode(fileread(fullfile(res_dir, 'synth_U1b.json')));
U2 = jsondecode(fileread(fullfile(res_dir, 'synth_U2b.json')));
U3 = jsondecode(fileread(fullfile(res_dir, 'synth_U3b.json')));
U4 = jsondecode(fileread(fullfile(res_dir, 'synth_U4b.json')));
V45 = jsondecode(fileread(fullfile(res_dir, 'volve45_unified.json')));
Us = {U1, U2, U3};

fig = figure('Visible', 'off', 'Position', [60, 60, 1320, 880], 'Color', 'w');
tl = tiledlayout(2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

% (a) Theorem 3 Noise Acc on U1b-U4b
nexttile(tl);
all4 = {U1, U2, U3, U4};
M_acc = zeros(4, 4);
for i = 1:4
    na = all4{i}.noise_acc;
    M_acc(i, :) = [na.RAW_B1_log2err, na.DSA_B1_log2err, na.RAW_B300_log2err, na.DSA_B300_log2err];
end
b = bar(M_acc);
b(1).FaceColor = [0.82 0.25 0.25]; b(2).FaceColor = [0.95 0.58 0.35];
b(3).FaceColor = [0.30 0.60 0.85]; b(4).FaceColor = [0.08 0.50 0.28];
hold on; yline(1, '--k', 'LineWidth', 1); yline(-1, '--k', 'LineWidth', 1); yline(0, ':k');
set(gca, 'XTickLabel', {'U1b 气烟囱', 'U2b 22°逆断层', 'U3b 单边拖缆', 'U4b 实测孪生'});
ylabel('log_2(iv_{est} / iv_{true})');
title('(a) 合理 SNR 档定理 3 盲估噪精度（DSA B=300 全落 \pm0.62 档内）');
legend({'RAW B=1', 'DSA B=1', 'RAW B=300', 'DSA B=300 (推荐)'}, 'Location', 'southwest', 'FontSize', 8);
grid on;

% (b) 4-factor monotonic E1 progression on U1b-U3b
nexttile(tl);
keys = {'RAW_UQW1b_100', 'RAW_UQW1_25', 'RAW_UQW1b_25', 'DSA_UQW1_25', 'DSA_UQW1b_25'};
lbls = {'RAW UQW1b(100)', 'RAW UQW1(B1)', 'RAW UQW1b(B300)', 'DSA UQW1(B1)', 'DSA-UQW1b(合流)'};
E_mat = zeros(3, numel(keys));
for i = 1:3
    for j = 1:numel(keys)
        E_mat(i, j) = Us{i}.W_deep.(keys{j}).E1_mean;
    end
end
b2 = bar(E_mat');
b2(1).FaceColor = [0.22 0.45 0.72]; b2(2).FaceColor = [0.85 0.45 0.15]; b2(3).FaceColor = [0.14 0.58 0.34];
set(gca, 'XTick', 1:numel(keys), 'XTickLabel', lbls); xtickangle(16);
ylim([0.50 0.80]);
ylabel('深部窗 5-m 平均 E1');
title('(b) U1b–U3b 四因子逐级合流 E1 单调跃升（3/3 严格递增）');
legend({'U1b (+1.7 dB)', 'U2b (+5.1 dB)', 'U3b (+2.5 dB)'}, 'Location', 'southeast', 'FontSize', 8);
grid on;

% (c) U4b Volve Twin well-tie curves
nexttile(tl);
ms = U4.m(:)'; M4 = U4.methods;
plot(ms, M4.DSA_UQW1b_EnvAnc_25.tie, '-d', 'Color', [0.10 0.45 0.80], 'LineWidth', 2.0); hold on;
plot(ms, M4.DSA_UQW1b_PAnc_25.tie,   '-o', 'Color', [0.05 0.50 0.25], 'LineWidth', 2.2, 'MarkerFaceColor', [0.05 0.50 0.25]);
plot(ms, M4.DSA_W1b_PAnc_25.tie,     '--s', 'Color', [0.85 0.48 0.10]);
plot(ms, M4.RAW_UQW1_25.tie,         '-.x', 'Color', [0.82 0.20 0.20]);
yline(M4.RAW_R0.mean, '--k', sprintf('RAW R0 (%.3f)', M4.RAW_R0.mean), 'LineWidth', 1.3);
xlabel('Q 模型乘子 m'); ylabel('5 井平均固定带相关 [4–18 cyc/km]');
title('(c) U4b Volve 全链路孪生体（未锚定掉分 \rightarrow 合流跃升 +0.083~+0.112）');
legend({'DSA-UQW1b-EnvAnc (25Hz)', 'DSA-UQW1b-PAnc (25Hz)', 'DSA-W1b-PAnc', 'RAW-UQW1 (未锚定)'}, 'Location', 'southeast', 'FontSize', 7.8);
grid on;

% (d) Volve 4-45 Hz SNR discovery: Stage 4H df=0.50 Hz wraparound vs full df=0.25 Hz DSA B=300
nexttile(tl);
f45 = V45.f(:)';
s1  = 10*log10(max(V45.snr_med_b1(:)',       1e-3));
s3r = 10*log10(max(V45.snr_med_b300_raw(:)', 1e-3));
s3d = 10*log10(max(V45.snr_med_b300_dsa(:)', 1e-3));
plot(f45, s1,  '--r', 'LineWidth', 1.5); hold on;
plot(f45, s3r, '-b',  'LineWidth', 1.7);
plot(f45, s3d, '-',   'Color', [0.05 0.50 0.25], 'LineWidth', 2.3);
yline(0, ':k', '0 dB');
xlabel('频率 f (Hz)'); ylabel('储层 [2600–3200 m] 中位成像域 SNR (dB)');
title('(d) Volve 4–45 Hz 重要机制发现：Pass-1 隔频混叠 vs 全频 DSA B=300');
legend({'Stage 4H B=1 (df=0.50 Hz, 混叠)', 'Stage 4H RAW B=300 (df=0.50 Hz)', '全频 DSA B=300 (df=0.25 Hz, 40Hz仍+4.1dB)'}, ...
       'Location', 'southwest', 'FontSize', 7.6);
xlim([4 45]); ylim([-32 22]); grid on;

out_png = fullfile(fig_dir, 'fig1b_synth_midSNR_and_snr_discovery.png');
exportgraphics(fig, out_png, 'Resolution', 300);
close(fig);
fprintf('Exported %s\n', out_png);
end
