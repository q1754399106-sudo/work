function make_xmir_qrtm_figures()
% Generate 4 publication-grade 300-DPI figures for uq_rtm_xmir_qrtm_20261126
%   fig1_wp1_xmir_snell_and_volve.png  : Analytical Snell circular-arc error (Str vs Bent) + Volve xmir audit
%   fig2_synth_s1_s5_head_to_head.png  : S1-S5 E1(m) & E_amp(m) head-to-head + S3 NF-vs-Noisy crossover
%   fig3_synth_2d_sections_s2_s3.png   : 2D migrated sections on S2 (Str vs Bent vs WE-QRTM) & S3 (m=0.5 stress)
%   fig4_volve25_xmir_qrtm_showdown.png: Volve 4-25 Hz F-5/F-14 curves, spectra, 2D sections & LU-cost Pareto

here = fileparts(mfilename('fullpath'));
root = fileparts(here);
res_dir = fullfile(root, 'results');
fig_dir = fullfile(root, 'figures');
if ~isfolder(fig_dir), mkdir(fig_dir); end

set(0, 'DefaultAxesFontName', 'DejaVu Sans', 'DefaultTextFontName', 'DejaVu Sans');

%% ========================================================================
%  FIG 1: WP1 Analytical Circular-Arc xmir_bent Proof & Volve Discrepancy
%  ========================================================================
wp1 = jsondecode(fileread(fullfile(res_dir, 'wp1_physics_verify.json')));
t1a = wp1.test1a_analytical_ray;
keys_d = {'dip_m20','dip_m10','dip_p0','dip_p10','dip_p20'};
dips = [-20, -10, 0, 10, 20];
e_str_med = arrayfun(@(i) t1a.(keys_d{i}).xmir_med_err_str_m,  1:5);
e_bnt_med = arrayfun(@(i) t1a.(keys_d{i}).xmir_med_err_bent_m, 1:5);
e_str_p90 = arrayfun(@(i) t1a.(keys_d{i}).xmir_p90_err_str_m,  1:5);
e_bnt_p90 = arrayfun(@(i) t1a.(keys_d{i}).xmir_p90_err_bent_m, 1:5);
th_str    = arrayfun(@(i) t1a.(keys_d{i}).th_s_med_err_str_deg,  1:5);
th_bnt    = arrayfun(@(i) t1a.(keys_d{i}).th_s_med_err_bent_deg, 1:5);

f1 = figure('Visible','off','Position',[50 50 1350 460],'Color','w');
subplot(1,3,1);
b = bar(dips, [e_str_med(:), e_bnt_med(:)], 'grouped');
b(1).FaceColor = [0.85 0.33 0.10]; b(2).FaceColor = [0.00 0.45 0.74];
hold on;
plot(dips-1.4, e_str_p90, 'v', 'Color', [0.65 0.15 0.05], 'MarkerFaceColor', [0.85 0.33 0.10], 'MarkerSize', 6);
plot(dips+1.4, e_bnt_p90, '^', 'Color', [0.00 0.25 0.55], 'MarkerFaceColor', [0.00 0.45 0.74], 'MarkerSize', 6);
grid on; xlabel('Reflector Dip \phi (deg)'); ylabel('|x_{mir} - x_{mir,exact}| (m)');
title('(a) Analytical v(z)=1500+0.8z Specular Error');
legend({'Straight-Ray x_{mir} (Med)', 'Bent-Ray x_{mir} (Med)', 'Straight p90', 'Bent p90'}, 'Location','northwest','FontSize',8);

subplot(1,3,2);
b2 = bar(dips, [th_str(:), th_bnt(:)], 'grouped');
b2(1).FaceColor = [0.85 0.33 0.10]; b2(2).FaceColor = [0.00 0.45 0.74];
grid on; xlabel('Reflector Dip \phi (deg)'); ylabel('|\theta_s - \theta_{s,exact}| (deg)');
title('(b) Incident Ray Angle Error (Chord vs Eikonal)');
legend({'Straight Chord \theta_s', 'Eikonal \nablaT \theta_s'}, 'Location','northeast','FontSize',8);

subplot(1,3,3);
t1c = wp1.test1c_volve_xmir_audit;
cats = categorical({'Res All\n(Med)', 'Res All\n(p90)', 'Dip>=5°\n(Med)', 'Dip>=5°\n(p90)'});
cats = reordercats(cats, {'Res All\n(Med)', 'Res All\n(p90)', 'Dip>=5°\n(Med)', 'Dip>=5°\n(p90)'});
vals = [t1c.res_med_dxmir_m, t1c.res_p90_dxmir_m, t1c.res_steep5_med_dxmir_m, t1c.res_steep5_p90_dxmir_m];
bh = bar(cats, vals, 0.55, 'FaceColor', [0.20 0.55 0.45]);
for k = 1:numel(vals)
    text(k, vals(k) + 22, sprintf('%.0f m', vals(k)), 'HorizontalAlignment','center','FontSize',9,'FontWeight','bold');
end
ylim([0, max(vals)*1.18]); grid on; ylabel('|x_{mir,bent} - x_{mir,str}| (m)');
title(sprintf('(c) Volve 10m Field Audit (Flip=%.1f%%)', t1c.res_gate_flip_pct));
exportgraphics(f1, fullfile(fig_dir, 'fig1_wp1_xmir_snell_and_volve.png'), 'Resolution', 300);
close(f1);

%% ========================================================================
%  FIG 2: Synthetic S1-S5 Head-to-Head Curves & S3 Multipathing Crossover
%  ========================================================================
s_files = {'synth_S1_Chimney_Base.json', 'synth_S2_StrongGrad_SteepDip.json', ...
           'synth_S3_SharpQ_Multipath.json', 'synth_S4_NonScalar_SpatialQ.json', ...
           'synth_S5_FD_GSLS_NoInvCrime.json'};
s_titles = {'S1: Chimney Base', 'S2: StrongGrad + 24° Dip + AsymCable', ...
            'S3: Sharp Q=12/-22%v Lens', 'S4: NonScalar Spatial Q Err', ...
            'S5: Indep. 8th-Order GSLS FD'};

f2 = figure('Visible','off','Position',[30 30 1550 820],'Color','w');
plot_methods = {'RAW_R0','WE_QRTM_Unreg','WE_QRTM_ZhuLP','WE_QRTM_Reg','DSA_WE_QRTM_W1b','DSA_Str_UQW1b_25','DSA_Bent_UQW1b_25','WE_Disp_UQW1b_25'};
plot_labels  = {'RAW R0','WE-QRTM Unreg','WE-QRTM ZhuLP','WE-QRTM Reg','DSA WE-QRTM W1b','DSA Str UQW1b','DSA Bent UQW1b (Prop)','WE-Disp UQW1b (Hyb)'};
cols = [0.45 0.45 0.45; 0.85 0.15 0.15; 0.90 0.50 0.10; 0.65 0.20 0.70; 0.15 0.60 0.25; 0.30 0.70 0.90; 0.00 0.30 0.80; 0.00 0.00 0.00];
lsty = {':','--','-.','--','-','-.','-','-'};
lwid = [1.5, 1.5, 1.5, 1.6, 2.0, 1.8, 2.6, 2.0];
mk   = {'none','x','s','d','^','v','o','p'};

for idx = 1:5
    D = jsondecode(fileread(fullfile(res_dir, s_files{idx})));
    ms = D.m(:)';
    subplot(2, 3, idx); hold on;
    for im = 1:numel(plot_methods)
        k = plot_methods{im};
        y = D.W_deep.(k).E1(:)';
        if numel(y) == 1, y = repmat(y, size(ms)); end
        plot(ms, y, 'Color', cols(im,:), 'LineStyle', lsty{im}, 'LineWidth', lwid(im), 'Marker', mk{im}, 'MarkerSize', 5);
    end
    grid on; xlabel('Q Perturbation Factor m = Q_{est}/Q_{true}'); ylabel('W_{deep} E1');
    title(sprintf('(%c) %s (SNR=%+.1fdB)', 'a'+idx-1, s_titles{idx}, D.W_deep.RAW_R0.SNR_img_dB_mean), 'FontSize', 10);
    if idx == 1
        legend(plot_labels, 'Location','southeast','FontSize',7);
    end
end

% Panel (f): S3 Exact Crossover Bar Chart (Noise-Free m=1 vs Noisy m=0.5 vs Noisy 5-m Mean)
D3 = jsondecode(fileread(fullfile(res_dir, 'synth_S3_SharpQ_Multipath.json')));
subplot(2, 3, 6);
cm = {'RAW_R0','WE_QRTM_Unreg','WE_QRTM_Reg','DSA_WE_QRTM_W1b','DSA_Bent_UQW1b_25','WE_Disp_UQW1b_25'};
cl = {'R0','WE Unr','WE Reg','DSA-WE W1b','Bent UQW1b','WE-Disp UQ'};
nf_v = arrayfun(@(i) D3.NoiseFree_m1.(cm{i}).E1, 1:numel(cm));
m05_v = arrayfun(@(i) D3.W_deep.(cm{i}).E1(1), 1:numel(cm));
avg_v = arrayfun(@(i) D3.W_deep.(cm{i}).E1_mean, 1:numel(cm));
bc = bar([nf_v(:), avg_v(:), m05_v(:)], 'grouped');
bc(1).FaceColor = [0.20 0.65 0.45]; bc(2).FaceColor = [0.00 0.45 0.74]; bc(3).FaceColor = [0.85 0.33 0.10];
set(gca, 'XTickLabel', cl, 'XTickLabelRotation', 25, 'FontSize', 8);
ylabel('W_{deep} E1'); grid on;
title('(f) S3 Multipath Crossover: NF(m=1) vs Noisy');
legend({'Noise-Free m=1','Noisy 5-m Mean','Noisy m=0.5 Stress'}, 'Location','southwest','FontSize',7);
exportgraphics(f2, fullfile(fig_dir, 'fig2_synth_s1_s5_head_to_head.png'), 'Resolution', 300);
close(f2);

%% ========================================================================
%  FIG 3: 2D Migrated Sections on S2 (Str vs Bent vs WE-QRTM) & S3 (m=0.5 Stress)
%  ========================================================================
img_dir = 'D:\tmp_uqrtm_xmir_qrtm_20261126';
if isfile(fullfile(img_dir, 'imgs_S2_StrongGrad_SteepDip.mat')) && isfile(fullfile(img_dir, 'imgs_S3_SharpQ_Multipath.mat'))
    I2 = load(fullfile(img_dir, 'imgs_S2_StrongGrad_SteepDip.mat'));
    I3 = load(fullfile(img_dir, 'imgs_S3_SharpQ_Multipath.mat'));
    f3 = figure('Visible','off','Position',[40 40 1550 760],'Color','w');
    panels2 = {'ref_2d', 'WE_QRTM_Unreg_mlow', 'DSA_WE_QRTM_W1b_mlow', 'DSA_Str_UQW1b_25_m1', 'DSA_Bent_UQW1b_25_mlow'};
    t2 = {'S2 Elastic Ref', 'S2 WE-QRTM Unreg (m=0.5)', 'S2 DSA WE-QRTM W1b (m=0.5)', 'S2 DSA Str UQW1b (m=1)', 'S2 DSA Bent UQW1b (m=0.5)'};
    for ip = 1:5
        subplot(2, 5, ip);
        A = double(I2.(panels2{ip}));
        clm = prctile(abs(A(I2.z_in>=450,:)), 98.5, 'all'); if clm <= 0, clm = 1; end
        imagesc(I2.x_in, I2.z_in, A, [-clm clm]); colormap(gca, gray);
        ylim([160 980]); title(t2{ip}, 'FontSize', 9);
        if ip == 1, ylabel('Depth z (m)'); else, set(gca,'YTickLabel',[]); end
    end
    panels3 = {'ref_2d', 'WE_QRTM_Unreg_mlow', 'WE_QRTM_Reg_mlow', 'DSA_WE_QRTM_W1b_mlow', 'DSA_Bent_UQW1b_25_mlow'};
    t3 = {'S3 Elastic Ref', 'S3 WE-QRTM Unreg (m=0.5)', 'S3 WE-QRTM Reg (m=0.5)', 'S3 DSA WE-QRTM W1b (m=0.5)', 'S3 DSA Bent UQW1b (m=0.5)'};
    for ip = 1:5
        subplot(2, 5, 5 + ip);
        A = double(I3.(panels3{ip}));
        clm = prctile(abs(A(I3.z_in>=450,:)), 98.5, 'all'); if clm <= 0, clm = 1; end
        imagesc(I3.x_in, I3.z_in, A, [-clm clm]); colormap(gca, gray);
        ylim([160 980]); xlabel('x (m)'); title(t3{ip}, 'FontSize', 9);
        if ip == 1, ylabel('Depth z (m)'); else, set(gca,'YTickLabel',[]); end
    end
    exportgraphics(f3, fullfile(fig_dir, 'fig3_synth_2d_sections_s2_s3.png'), 'Resolution', 300);
    close(f3);
end

%% ========================================================================
%  FIG 4: Volve 4-25 Hz Showdown (F-5/F-14, Reservoir Spectra, 2D Sections)
%  ========================================================================
v25_p = fullfile(res_dir, 'volve25_xmir_qrtm.json');
if isfile(v25_p)
    V = jsondecode(fileread(v25_p));
    ms = V.m(:)'; fv = V.f(:)';
    f4 = figure('Visible','off','Position',[40 40 1500 840],'Color','w');

    v_keys = {'RAW_R0','DSA_Bent_R0_Eq','WE_QRTM_Unreg_PAnc','WE_QRTM_ZhuLP_PAnc','DSA_WE_QRTM_Reg_PAnc','DSA_WE_QRTM_W1b_PAnc','DSA_Str_UQW1b_PAnc_25','DSA_Bent_UQW1b_PAnc_25'};
    v_labs = {'RAW R0','DSA Bent R0 Eq','WE-QRTM Unreg PAnc','WE-QRTM ZhuLP PAnc','DSA WE-QRTM Reg PAnc','DSA WE-QRTM W1b PAnc','DSA Str UQW1b PAnc','DSA Bent UQW1b PAnc'};
    v_cols = [0.45 0.45 0.45; 0.25 0.65 0.65; 0.85 0.15 0.15; 0.90 0.50 0.10; 0.65 0.20 0.70; 0.15 0.60 0.25; 0.30 0.70 0.90; 0.00 0.30 0.80];

    subplot(2,3,1); hold on;
    for i = 1:numel(v_keys)
        y = V.methods.(v_keys{i}).F5_R0band(:)'; if numel(y)==1, y=repmat(y,size(ms)); end
        plot(ms, y, '-o', 'Color', v_cols(i,:), 'LineWidth', 1.6 + 0.8*(i==numel(v_keys)), 'MarkerSize', 5);
    end
    grid on; xlabel('m = Q_{est}/Q_{nom}'); ylabel('F-5 Well Tie (R0 band)');
    title('(a) Volve 4-25 Hz F-5 Well Tie vs m');
    legend(v_labs, 'Location','southeast','FontSize',7);

    subplot(2,3,2); hold on;
    for i = 1:numel(v_keys)
        y = V.methods.(v_keys{i}).F14_R0band(:)'; if numel(y)==1, y=repmat(y,size(ms)); end
        plot(ms, y, '-s', 'Color', v_cols(i,:), 'LineWidth', 1.6 + 0.8*(i==numel(v_keys)), 'MarkerSize', 5);
    end
    grid on; xlabel('m = Q_{est}/Q_{nom}'); ylabel('Blind Well F-14 Tie');
    title('(b) Volve 4-25 Hz Blind Well F-14 Tie');

    subplot(2,3,3); hold on;
    ref_sp = V.spec_res.DSA_Bent_R0_Eq(:)';
    for i = [2,3,5,6,7,8]
        sp = V.spec_res.(v_keys{i})(:)';
        plot(fv, sp ./ max(ref_sp, 1e-30), 'Color', v_cols(i,:), 'LineWidth', 1.7 + 0.7*(i==8));
    end
    yline(1.0, 'k--'); grid on; xlabel('Frequency (Hz)'); ylabel('Reservoir Spec Ratio vs DSA Bent R0 Eq');
    title('(c) Volve 4-25 Hz Active HF Gain (m=1)');
    legend(v_labs([2,3,5,6,7,8]), 'Location','northwest','FontSize',7);

    vimg_p = fullfile(img_dir, 'volve25_xmir_qrtm_imgs.mat');
    if isfile(vimg_p)
        IV = load(vimg_p);
        sec_k = {'RAW_R0_m1', 'WE_QRTM_Unreg_PAnc_m1', 'DSA_Bent_UQW1b_PAnc_25_m1'};
        sec_t = {'(d) Volve RAW R0 (Sm=38.4%)', '(e) Volve WE-QRTM Unreg PAnc (m=1)', '(f) Volve DSA Bent UQW1b PAnc (m=1)'};
        for is = 1:3
            subplot(2,3,3+is);
            A = double(IV.(sec_k{is}));
            clm = prctile(abs(A(IV.z_in>=1600 & IV.z_in<=3200, :)), 98.5, 'all');
            imagesc(IV.x_in, IV.z_in, A, [-clm clm]); colormap(gca, gray);
            ylim([1000 3400]); xlim([-2500 1500]);
            xlabel('Inline x (m)'); if is==1, ylabel('Depth z (m)'); end
            title(sec_t{is}, 'FontSize', 9);
        end
    end
    exportgraphics(f4, fullfile(fig_dir, 'fig4_volve25_xmir_qrtm_showdown.png'), 'Resolution', 300);
    close(f4);
end
fprintf('Generated all 4 publication-grade 300-DPI figures in %s\n', fig_dir);
end
