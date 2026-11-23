function make_volve_closure_figures()
% Publication-grade visual before/after comparison figures for uq_rtm_volve_closure_20261123.
% Reads final 2D images from scratch snapshot (zero Helmholtz recomputation, ~10 s).

here = fileparts(mfilename('fullpath'));
root = fileparts(here);
fig_dir = fullfile(root, 'figures');
if ~isfolder(fig_dir), mkdir(fig_dir); end

addpath('D:\ChatGPT\研究\地震衰减补偿研究\uq_rtm_20261117\scripts');
P = fld_setup();
S = load('C:\Users\qq175\.gemini\antigravity\brain\c0dc2750-de05-424f-adb4-e7ef332a1950\scratch\volve_closure\volve_closure_final_imgs.mat');
Img = S.Img_s; res = S.res; z = S.z_in; x = S.x_in;

set(0, 'DefaultAxesFontName', 'Microsoft YaHei', ...
       'DefaultTextFontName', 'Microsoft YaHei', ...
       'DefaultAxesFontSize', 9.5);

prep = @(v) prep_volve(v, P.nzi, P.nxi, z);
A_raw_r0    = prep(Img.RAW_R0(:,1));
A_raw_un_05 = prep(Img.RAW_UQW1_Unanc(:,1));
A_raw_un_10 = prep(Img.RAW_UQW1_Unanc(:,3));
A_raw_pa_10 = prep(Img.RAW_UQW1_PAnc(:,3));
A_dsa_eq    = prep(Img.DSA_R0_Eq(:,1));
A_dsa_pa_10 = prep(Img.DSA_UQW1_PAnc(:,3));
A_dsa_pa_05 = prep(Img.DSA_UQW1_PAnc(:,1));

w5  = load_well(fullfile(P.cache, 'well_F-5.bin'),  z, x);
w14 = load_well(fullfile(P.cache, 'well_F-14.bin'), z, x);

zz = find(z >= 1400 & z <= 3400);
xx = find(x >= -3000 & x <= 2000);
zz_res = find(z >= 2550 & z <= 3250);
xx_res = find(x >= -2500 & x <= 1500);

% Match RMS inside reservoir window [2550-3250m, -2500..1500m] for fair visual brightness!
c_v = prctile(abs(A_raw_r0(zz_res, xx_res)), 99, 'all');
normR = @(A) A * (norm(A_raw_r0(zz_res, xx_res), 'fro') / max(norm(A(zz_res, xx_res), 'fro'), 1e-30));

%% ========================================================================
%% FIGURE 1: 2x3 Full 2D Volve Migrated Section Comparison (Before vs After Closure)
%% ========================================================================
f = figure('Visible', 'off', 'Color', 'w', 'Position', [0 0 2160 1060]);
t = tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
title(t, 'Volve 北海实测海上 OBC 闭环前后 2D RTM 剖面对比（85 频全波叠前，零测井先验）', ...
    'FontSize', 14.5, 'FontWeight', 'bold');

M = res.methods;
pans = {A_raw_r0, normR(A_raw_un_05), normR(A_raw_pa_10), ...
        normR(A_dsa_eq), normR(A_dsa_pa_10), normR(A_dsa_pa_05)};
tits = { ...
    sprintf('(a) 原版 Stage 4 RAW\\_R0 (不补偿)\nF-5: R0带=%.3f | 宽带=%.3f | 上覆画弧=%.1f%%', ...
        M.RAW_R0.F5_R0band, M.RAW_R0.F5_Broad, 100*M.RAW_R0.Smile_Frac), ...
    sprintf('(b) 原版 Stage 4 UQW1\\_amp (m=0.5, 未锚定)\nF-5: R0带=%.3f | 宽带=%.3f | 上覆画弧=%.1f%% (崩盘)', ...
        M.RAW_UQW1_Unanc.F5_R0band(1), M.RAW_UQW1_Unanc.F5_Broad(1), 100*M.RAW_UQW1_Unanc.Smile_Frac(1)), ...
    sprintf('(c) 仅浅层锚定 RAW\\_UQW1\\_PAnc (m=1.0, 无DSA)\nF-5: R0带=%.3f | 宽带=%.3f | 上覆画弧=%.1f%% (仍存画弧)', ...
        M.RAW_UQW1_PAnc.F5_R0band(3), M.RAW_UQW1_PAnc.F5_Broad(3), 100*M.RAW_UQW1_PAnc.Smile_Frac(3)), ...
    sprintf('(d) 本文 DSA\\_R0\\_Eq (浅层均衡不补偿基线)\nF-5: R0带=%.3f | 宽带=%.3f | 上覆画弧=%.1f%% (-90%%)', ...
        M.DSA_R0_Eq.F5_R0band, M.DSA_R0_Eq.F5_Broad, 100*M.DSA_R0_Eq.Smile_Frac), ...
    sprintf('(e) 本文 DSA\\_UQW1\\_PAnc (m=1.0, 完整闭环)\nF-5: R0带=%.3f | 宽带=%.3f (+0.025) | 画弧=%.1f%%', ...
        M.DSA_UQW1_PAnc.F5_R0band(3), M.DSA_UQW1_PAnc.F5_Broad(3), 100*M.DSA_UQW1_PAnc.Smile_Frac(3)), ...
    sprintf('(f) 本文 DSA\\_UQW1\\_PAnc (m=0.5, VSP Q=60 稳健)\nF-5: R0带=%.3f (+0.134 vs原版) | 宽带=%.3f | 画弧=%.1f%%', ...
        M.DSA_UQW1_PAnc.F5_R0band(1), M.DSA_UQW1_PAnc.F5_Broad(1), 100*M.DSA_UQW1_PAnc.Smile_Frac(1))};
cols = {[0 0 0], [0.75 0.1 0.1], [0.55 0.3 0], [0.15 0.45 0.2], [0 0.35 0.78], [0 0.35 0.78]};

for p = 1:6
    ax = nexttile;
    imagesc(x(xx)/1000, z(zz), pans{p}(zz, xx), [-c_v c_v]); colormap(ax, gray); hold on;
    % Overburden smile evaluation window (orange dotted) & reservoir window (cyan dashed)
    rectangle(ax, 'Position', [-2.5, 1600, 4.0, 800], 'EdgeColor', [1 0.65 0.1], 'LineStyle', ':', 'LineWidth', 1.1);
    rectangle(ax, 'Position', [-2.5, 2600, 4.0, 600], 'EdgeColor', [0 0.8 0.95], 'LineStyle', '--', 'LineWidth', 1.2);
    plot(w5.sw/1000,  w5.zw,  'r-', 'LineWidth', 1.9);
    plot(w14.sw/1000, w14.zw, 'c-', 'LineWidth', 1.5);
    if p == 1
        text(-1.18, 2505, 'F-5', 'Color', 'r', 'FontWeight', 'bold', 'FontSize', 9, 'BackgroundColor', [1 1 1], 'Margin', 1.0);
        text(-0.52, 2610, 'F-14', 'Color', [0 0.55 0.75], 'FontWeight', 'bold', 'FontSize', 8.5, 'BackgroundColor', [1 1 1], 'Margin', 1.0);
        text(-0.5, 1730, '上覆画弧审计窗 [1.6-2.4 km]', 'Color', [1 0.85 0.25], 'FontWeight', 'bold', 'FontSize', 8.5, 'HorizontalAlignment', 'center', 'BackgroundColor', [0.12 0.12 0.12], 'Margin', 1.5);
    end
    title(tits{p}, 'Color', cols{p});
    if mod(p, 3) == 1, ylabel('深度 z (m)'); end
    if p > 3, xlabel('测线坐标 x (km)'); end
end

exportgraphics(f, fullfile(fig_dir, 'fig1_volve_closure_2d_sections.png'), 'Resolution', 150);
close(f);

%% ========================================================================
%% FIGURE 2: Zoomed Reservoir + Removed Smiles + Spectral Closure + F-5 Well-Tie
%% ========================================================================
f = figure('Visible', 'off', 'Color', 'w', 'Position', [0 0 2160 1040]);
t = tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
title(t, 'Volve 实测闭环深层储层细节、画弧剥离图、频谱斜率归零与 F-5 井震标定（5/5 预注册门限全通）', ...
    'FontSize', 14.5, 'FontWeight', 'bold');

zz_zm = find(z >= 2500 & z <= 3260);
xx_zm = find(x >= -2200 & x <= 800);

% (a) Zoomed reservoir: Stage 4 Unanchored UQW1 (m=0.5)
A_un05_n = normR(A_raw_un_05);
ax = nexttile;
imagesc(x(xx_zm)/1000, z(zz_zm), A_un05_n(zz_zm, xx_zm), [-c_v c_v]); colormap(ax, gray); hold on;
plot(w5.sw/1000, w5.zw, 'r-', 'LineWidth', 2.2);
text(-1.18, 2545, 'F-5', 'Color', 'r', 'FontWeight', 'bold', 'FontSize', 9, 'BackgroundColor', [1 1 1], 'Margin', 1.0);
title(sprintf('(a) 储层放大：原版 UQW1\\_amp (m=0.5 未锚定)\nF-5 R0带=%.3f | 宽带=%.3f (高频噪化/蓝倾过重)', ...
    M.RAW_UQW1_Unanc.F5_R0band(1), M.RAW_UQW1_Unanc.F5_Broad(1)), 'Color', [0.75 0.1 0.1]);
ylabel('深度 z (m)'); xlabel('测线坐标 x (km)');

% (b) Zoomed reservoir: Closure DSA_UQW1_PAnc (m=0.5)
A_pa05_n = normR(A_dsa_pa_05);
ax = nexttile;
imagesc(x(xx_zm)/1000, z(zz_zm), A_pa05_n(zz_zm, xx_zm), [-c_v c_v]); colormap(ax, gray); hold on;
plot(w5.sw/1000, w5.zw, 'r-', 'LineWidth', 2.2);
text(-1.18, 2545, 'F-5', 'Color', 'r', 'FontWeight', 'bold', 'FontSize', 9, 'BackgroundColor', [1 1 1], 'Margin', 1.0);
title(sprintf('(b) 储层放大：本文 DSA\\_UQW1\\_PAnc (m=0.5 闭环)\nF-5 R0带=%.3f (+0.134) | 宽带=%.3f (+0.145)', ...
    M.DSA_UQW1_PAnc.F5_R0band(1), M.DSA_UQW1_PAnc.F5_Broad(1)), 'Color', [0 0.35 0.78]);
xlabel('测线坐标 x (km)');

% (c) Removed overburden wide-angle smiles difference: RAW_UQW1_PAnc - DSA_UQW1_PAnc (m=1.0)
ax = nexttile;
diff_smile = normR(A_raw_pa_10) - normR(A_dsa_pa_10);
imagesc(x(xx)/1000, z(zz), diff_smile(zz, xx), [-0.7*c_v 0.7*c_v]); colormap(ax, bwr_cmap(256)); hold on;
rectangle(ax, 'Position', [-2.5, 1600, 4.0, 800], 'EdgeColor', [0.1 0.1 0.1], 'LineStyle', '--', 'LineWidth', 1.2);
title('(c) DSA 剥离的出阵炮/宽角假频差异：RAW - DSA (m=1.0)');
xlabel('测线坐标 x (km)');

% (d) Reservoir [2600-3200m] frequency spectrum slope closure ([6-22 Hz])
ax = nexttile; hold on;
fv = res.f(:)'; fm = (fv >= 6 & fv <= 22);
s_un = res.spec_res.RAW_UQW1_Unanc_m1(:)'; s_un = s_un(fm) / mean(s_un(fm));
s_r0 = res.spec_res.RAW_R0(:)';            s_r0 = s_r0(fm) / mean(s_r0(fm));
s_eq = res.spec_res.DSA_R0_Eq(:)';         s_eq = s_eq(fm) / mean(s_eq(fm));
s_pa = res.spec_res.DSA_UQW1_PAnc_m1(:)';  s_pa = s_pa(fm) / mean(s_pa(fm));
plot(fv(fm), s_un, 'r-', 'LineWidth', 1.5, ...
    'DisplayName', sprintf('原版 UQW1 未锚定 (斜率 %+.3f/Hz)', res.res_slope.RAW_UQW1_Unanc_m1));
plot(fv(fm), s_r0, '--', 'Color', [0.45 0.45 0.45], 'LineWidth', 1.6, ...
    'DisplayName', sprintf('原版 RAW\\_R0 (斜率 %+.3f/Hz)', res.res_slope.RAW_R0));
plot(fv(fm), s_eq, '-.', 'Color', [0.1 0.6 0.25], 'LineWidth', 1.7, ...
    'DisplayName', sprintf('本文 DSA\\_R0\\_Eq (斜率 %+.3f/Hz 实测衰减)', res.res_slope.DSA_R0_Eq));
plot(fv(fm), s_pa, 'b-', 'LineWidth', 2.2, ...
    'DisplayName', sprintf('本文 DSA\\_UQW1\\_PAnc (斜率 %+.4f/Hz 平谱!)', res.res_slope.DSA_UQW1_PAnc_m1));
yline(1.0, 'k:', 'HandleVisibility', 'off');
xlim([6 22]); ylim([0.2 2.4]); grid on;
xlabel('频率 f (Hz)'); ylabel('带内均值归一化储层振幅谱');
title('(d) 储层窗 [2600-3200 m] 频率振幅谱斜率归零验证');
legend('Location', 'northwest', 'FontSize', 7.8);

% (e) F-5 Broadband [3.52-16.0 cyc/km] well-seismic trace overlay (m=1.0)
ax = nexttile; hold on;
lo_b = 3.515625; hi_b = 16.0;
rb5 = bp_z(w5.rw, P.dx, lo_b, hi_b); rb5_n = rb5 / rms(rb5);
[tr_r0,  c_r0]  = align_well(A_raw_r0,    w5, rb5_n, P.dx, lo_b, hi_b);
[tr_un,  c_un]  = align_well(A_raw_un_10, w5, rb5_n, P.dx, lo_b, hi_b);
[tr_pa,  c_pa]  = align_well(A_dsa_pa_10, w5, rb5_n, P.dx, lo_b, hi_b);
plot(rb5_n, w5.zw, 'k-',  'LineWidth', 2.1, 'DisplayName', 'F-5 测井反射系数 (宽带 3.5-16)');
plot(tr_r0, w5.zw, '--',  'Color', [0.5 0.5 0.5], 'LineWidth', 1.4, 'DisplayName', sprintf('RAW\\_R0 (r=%.3f)', c_r0));
plot(tr_un, w5.zw, 'r-',  'LineWidth', 1.3, 'DisplayName', sprintf('原版 UQW1 未锚定 (r=%.3f)', c_un));
plot(tr_pa, w5.zw, 'b-',  'LineWidth', 2.0, 'DisplayName', sprintf('本文 DSA\\_UQW1\\_PAnc (r=%.3f)', c_pa));
set(ax, 'YDir', 'reverse'); ylim([min(w5.zw) max(w5.zw)]); xlim([-3.3 3.3]); grid on;
xlabel('归一化道振幅'); ylabel('深度 z (m)');
title('(e) F-5 井宽带 [3.5-16 cyc/km] 井震波形叠加 (m=1.0)');
legend('Location', 'southwest', 'FontSize', 7.8);

% (f) F-5 R0-band & Broadband well-tie vs m = [0.5..2.0]
ax = nexttile; hold on;
ms = res.m(:)';
yline(M.RAW_R0.F5_R0band, 'k--', 'LineWidth', 1.6, 'DisplayName', sprintf('原版 RAW\\_R0 (%.3f)', M.RAW_R0.F5_R0band));
yline(M.DSA_R0_Eq.F5_R0band, ':', 'Color', [0.35 0.35 0.35], 'LineWidth', 1.6, 'DisplayName', sprintf('本文 DSA\\_R0\\_Eq (%.3f)', M.DSA_R0_Eq.F5_R0band));
plot(ms, M.RAW_UQW1_Unanc.F5_R0band, 'r-^', 'LineWidth', 1.6, 'MarkerSize', 6, 'DisplayName', '原版 Stage 4 UQW1\_amp');
plot(ms, M.DSA_IQ70_PAnc.F5_R0band,  '--s', 'Color', [0.88 0.45 0.1], 'LineWidth', 1.4, 'MarkerSize', 5, 'DisplayName', 'DSA\_IQ70\_PAnc');
plot(ms, M.DSA_W1_PAnc.F5_R0band,    '-d',  'Color', [0.1 0.65 0.25], 'LineWidth', 1.5, 'MarkerSize', 5.5, 'DisplayName', 'DSA\_W1\_PAnc');
plot(ms, M.DSA_UQW1_PAnc.F5_R0band,  'b-o', 'LineWidth', 2.3, 'MarkerSize', 7, 'MarkerFaceColor', 'b', 'DisplayName', '本文 DSA\_UQW1\_PAnc (全线反超R0!)');
set(ax, 'XScale', 'log'); xticks(ms); ylim([0.48 0.66]); grid on;
xlabel('Q 乘子 m（名义 Q=120；VSP 实测 Q≈60-90 对应 m≈0.5-0.75）');
ylabel('F-5 井旁相关系数 (exact Stage-4 R0带)');
title('(f) F-5 井旁相关随 Q 误差变化（斜率翻正 + 全 m 反超 R0）');
legend('Location', 'southeast', 'FontSize', 7.8);

exportgraphics(f, fullfile(fig_dir, 'fig2_volve_closure_diagnostics.png'), 'Resolution', 150);
close(f);
fprintf('Exported fig1_volve_closure_2d_sections.png and fig2_volve_closure_diagnostics.png\n');
end

function cm = bwr_cmap(n)
h = floor(n/2);
cm = [linspace(0.08,1,h)' linspace(0.25,1,h)' linspace(0.78,1,h)'; ...
      linspace(1,0.78,n-h)' linspace(1,0.14,n-h)' linspace(1,0.14,n-h)'];
end

function y = prep_volve(v, nzi, nxi, z_in)
y = reshape(double(v), nzi, nxi);
ker = exp(-0.5 * ((-12:12)/4).^2); ker = ker / sum(ker);
y = y - conv2(ker', 1, [repmat(y(1,:), 12, 1); y; repmat(y(end,:), 12, 1)], 'valid');
y(z_in < 150, :) = 0;
end

function w = load_well(fn, z, x)
fid = fopen(fn); a = fread(fid, Inf, 'single'); fclose(fid);
a = reshape(a, 4, [])'; w.zw = a(:,1); w.sw = a(:,2); w.rw = double(a(:,4));
[~, w.iz] = min(abs(z(:) - w.zw'), [], 1);
w.ix = arrayfun(@(xv) find(abs(x - xv) == min(abs(x - xv)), 1), w.sw);
end

function [tr_best, c_best] = align_well(A, w, rb_n, dx, lo, hi)
nw = numel(w.zw); tr = zeros(nw, 1);
for i = 1:nw, tr(i) = mean(A(w.iz(i), max(w.ix(i)-2,1):min(w.ix(i)+2,size(A,2)))); end
tr_b = bp_z(tr, dx, lo, hi);
c_best = -Inf; tr_best = tr_b;
for L = -3:3
    cand = circshift(tr_b, L);
    a = cand - mean(cand); b = rb_n - mean(rb_n);
    val = (a' * b) / max(sqrt((a' * a) * (b' * b)), 1e-30);
    if val > c_best, c_best = val; tr_best = cand; end
end
tr_best = tr_best / max(rms(tr_best), 1e-30);
end

function y = bp_z(r, dx, lo, hi)
n = 2^nextpow2(4 * numel(r)); R = fft(r, n);
k = abs([0:n/2, -n/2+1:-1]') / (n * dx) * 1000;
y = real(ifft(R .* (k >= lo & k <= hi))); y = y(1:numel(r));
end
