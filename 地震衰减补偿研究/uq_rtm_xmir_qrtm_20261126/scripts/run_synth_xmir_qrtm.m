function run_synth_xmir_qrtm(model_filter)
% WP2: Head-to-Head Showdown — Upgraded Bent-Ray xmir_bent vs Legacy xmir_str
%      vs 8 True Viscoacoustic Wave-Equation Q-RTM Variants (Models S1-S4)
%
% Models (brand-new seeds 20261501-20261504, pre-calibrated via smoke check to RAW_R0 SNR_img in [+2,+6] dB):
%   S1_Chimney_Base        : Gas chimney + gentle dip + l30 corr noise, m in {0.5, 0.7, 1.0, 1.4, 2.0}
%   S2_StrongGrad_SteepDip : Strong gradient (1550->3950 m/s) + 24 deg steep fault + asymmetric cable + chimney
%   S3_SharpQ_Multipath    : Ultra-sharp Q_min=12 & -22% velocity pull-down gas lens (multipathing stress test)
%   S4_NonScalar_SpatialQ  : Shifted (+160m) & widened (1.35x) chimney + shallow(1.3x)/deep(0.75x) opposite Q error
%
% All 16 methods evaluated for BOTH Noise-Free (m=1) Oracle AND Noisy All-m W_deep:
%   Acoustic R0 (3)         : RAW_R0, DSA_Str_R0, DSA_Bent_R0
%   True Wave-Eq Q-RTM (8)  : WE_Q_Adj, WE_Q_DispOnly, DSA_WE_DispOnly, WE_QRTM_Unreg,
%                             WE_QRTM_ZhuLP, WE_QRTM_Reg, DSA_WE_QRTM_Reg, DSA_WE_QRTM_W1b
%   Acoustic Ray UQ/Det (4) : DSA_Str_IQ70_25, DSA_Str_UQW1b_25, DSA_Bent_W1b_25, DSA_Bent_UQW1b_25
%   Wave-Eq + UQ Hybrid (1) : WE_Disp_UQW1b_25

if nargin < 1, model_filter = {}; end
if ischar(model_filter) || isstring(model_filter), model_filter = {char(model_filter)}; end
maxNumCompThreads(10);

here = fileparts(mfilename('fullpath'));
root = fileparts(here);
res_dir = fullfile(root, 'results');
if ~isfolder(res_dir), mkdir(res_dir); end

addpath('D:\ChatGPT\研究\地震衰减补偿研究\uq_rtm_20261117\scripts');
core = uq_rtm_xmir_qrtm_core();
cfg  = core.default_cfg('synth_10m');

all_models = {'S1_Chimney_Base', 'S2_StrongGrad_SteepDip', 'S3_SharpQ_Multipath', 'S4_NonScalar_SpatialQ'};
% Pre-freeze smoke-calibrated snr_db -> RAW_R0 W_deep SNR_img ~ [+3.8, +4.4] dB
snr_map = struct('S1_Chimney_Base', 0.5, 'S2_StrongGrad_SteepDip', 11.5, ...
                 'S3_SharpQ_Multipath', 1.5, 'S4_NonScalar_SpatialQ', 0.5);

if isempty(model_filter)
    run_list = all_models;
else
    run_list = model_filter;
end

for im = 1:numel(run_list)
    mname = run_list{im};
    fn = fullfile(res_dir, ['synth_' mname '.json']);
    if isfile(fn)
        fprintf('Skipping %s (already exists: %s)\n', mname, fn);
        continue;
    end
    res = run_one_synth_model(mname, snr_map.(mname), core, cfg, root);
    fid = fopen(fn, 'w'); fwrite(fid, jsonencode(res, 'PrettyPrint', true), 'char'); fclose(fid);
    fprintf('Saved %s (total %.1fs)\n\n', fn, res.elapsed_s);
end
end

function out = run_one_synth_model(mname, snr_db, core, cfg, root)
t_all = tic;
[P, v_true, v_mig, Q_true, Q_mig_base, rec_mask, W_targ, W_flank, ms, seed, q_desc] = build_model(mname, core);
ns = numel(P.ixs); nr = numel(P.ixr); N = P.nz * P.nx; nm = numel(ms);
[Zi, Xi] = ndgrid(P.in{1}, P.in{2}); lin = sub2ind([P.nz P.nx], Zi(:), Xi(:)); nin = numel(lin);
isrc = sub2ind([P.nz P.nx], P.izs*ones(1,ns), P.ixs);
irec = sub2ind([P.nz P.nx], P.izs*ones(1,nr), P.ixr);
z_in = P.z(P.in{1}); x_in = P.x(P.in{2});
[Zin, Xin] = ndgrid(z_in, x_in);
xs = P.x(P.ixs); xr = P.x(P.ixr);
xr_lo = arrayfun(@(k) min(xr(rec_mask(:,k))), 1:ns);
xr_hi = arrayfun(@(k) max(xr(rec_mask(:,k))), 1:ns);

iz_w = find(z_in >= 560 & z_in <= 950);
ix_w = find(x_in >= 180 & x_in <= 1420);
w_deep_mask = (Zin(:) >= 560 & Zin(:) <= 950) & (Xin(:) >= 180 & Xin(:) <= 1420);
ovb_z = find(z_in >= 180 & z_in <= 520);
ovb_x = find(x_in >= 180 & x_in <= 1420);

freqs = 4.0:1.0:37.0; nf = numel(freqs);
fref = 25.0; fpeak = 22.0; sig = cfg.sig; s70 = gain_s(70);
m1_idx = find(abs(ms - 1.0) < 1e-9, 1);

rng(seed, 'twister');
corr_sig = 5.0; % l30 50 m laterally correlated noise
rc = ceil(3 * corr_sig); g_nc = exp(-0.5 * ((-rc:rc)' / corr_sig).^2); g_nc = g_nc / sqrt(sum(g_nc.^2));

% Calibrate noise amplitude at 22 Hz
W22 = exp(-1.0);
b22 = zeros(N, ns); b22(sub2ind([N ns], isrc, 1:ns)) = -W22 / P.dx^2;
[Lm,Um,Pm,Qm,Rm] = lu(core.helm_op_visco(v_mig, v_mig, 22, P));
Ub = core.solve_helm_shot(Lm,Um,Pm,Qm,Rm, b22); Dbg22 = Ub(irec,:); clear Ub Lm Um Pm Qm Rm
vc22 = core.make_visco_vel(v_true, Q_true, 22, fref, P.npml, 'att', 1.0);
[Lv,Uv,Pv,Qv,Rv] = lu(core.helm_op_visco(vc22, v_true, 22, P));
Uv22 = core.solve_helm_shot(Lv,Uv,Pv,Qv,Rv, b22); Dcal = (Uv22(irec,:) - Dbg22) .* rec_mask; clear Uv22 Lv Uv Pv Qv Rv
sig_n = sqrt(mean(abs(Dcal(rec_mask)).^2)) * 10^(-snr_db/20);

% Eikonal traveltime F_T and single-leg attenuation F_ts (using Q_mig_base)
s_mig = 1.0 ./ v_mig; eps_q = 0.5;
F_T  = zeros(P.nz, P.nx, ns, 'single');
F_ts = zeros(P.nz, P.nx, ns, 'single');
for k = 1:ns
    F_T(:,:,k)  = single(eikonal_fs(s_mig, P.dx, P.izs, P.ixs(k)));
    F_ts(:,:,k) = single((eikonal_fs(s_mig.*(1+eps_q./Q_mig_base), P.dx, P.izs, P.ixs(k)) - ...
                          eikonal_fs(s_mig.*(1-eps_q./Q_mig_base), P.dx, P.izs, P.ixs(k))) / (2*eps_q));
end

rng(seed + 1,   'twister'); EPS_B1   = 2*(rand(nr, ns) > 0.5) - 1;
rng(seed + 300, 'twister'); sgn_B300 = 2*(rand(1, ns)  > 0.5) - 1;

%% ========================================================================
%  PHASE 1: Data Generation + Acoustic RTM Wavefields (1 LU / freq!)
%  ========================================================================
t_ac = tic;
D_vi_all = zeros(nr, ns, nf, 'single');
N_f_all  = zeros(nr, ns, nf, 'single');
Cs_sig_all = zeros(nin, ns, nf, 'single');
Cs_nu_all  = zeros(nin, ns, nf, 'single');
Cs_b1_all  = zeros(nin, ns, nf, 'single');
I_ref = zeros(nin, 1); I_R0_tot = zeros(nin, 1);

for jf = 1:nf
    f = freqs(jf); w = 2*pi*f;
    W_src = (f/fpeak)^2 * exp(-(f/fpeak)^2);
    b_src = zeros(N, ns); b_src(sub2ind([N ns], isrc, 1:ns)) = -W_src / P.dx^2;

    % Acoustic migration operator (ONLY 1 LU per freq for ALL acoustic UQ/Det methods!)
    [Lm,Um,Pm,Qm,Rm] = lu(core.helm_op_visco(v_mig, v_mig, f, P));
    U_bg = core.solve_helm_shot(Lm,Um,Pm,Qm,Rm, b_src); D_bg = U_bg(irec,:);
    Us_c = conj(U_bg(lin,:)) * w^2; clear U_bg

    % Elastic true reference data
    [Le,Ue,Pe,Qe,Re] = lu(core.helm_op_visco(v_true, v_true, f, P));
    U_el = core.solve_helm_shot(Le,Ue,Pe,Qe,Re, b_src); D_el = (U_el(irec,:) - D_bg) .* rec_mask; clear U_el Le Ue Pe Qe Re

    % Viscoacoustic true forward data
    vc_t = core.make_visco_vel(v_true, Q_true, f, fref, P.npml, 'att', 1.0);
    [Lv,Uv,Pv,Qv,Rv] = lu(core.helm_op_visco(vc_t, v_true, f, P));
    U_vi = core.solve_helm_shot(Lv,Uv,Pv,Qv,Rv, b_src); D_vi = (U_vi(irec,:) - D_bg) .* rec_mask; clear U_vi Lv Uv Pv Qv Rv

    % Correlated l30 noise
    N_f = conv2(g_nc, 1, sig_n * (1/sqrt(2)) * (randn(nr, ns) + 1i*randn(nr, ns)), 'same') .* rec_mask;
    D_vi_all(:,:,jf) = single(D_vi);
    N_f_all(:,:,jf)  = single(N_f);

    b_el = zeros(N,ns); b_el(irec,:)=D_el; V_el = core.solve_helm_adj(Lm,Um,Pm,Qm,Rm, b_el, P.dx);
    b_vi = zeros(N,ns); b_vi(irec,:)=D_vi; V_vi = core.solve_helm_adj(Lm,Um,Pm,Qm,Rm, b_vi, P.dx);
    b_nu = zeros(N,ns); b_nu(irec,:)=N_f;  V_nu = core.solve_helm_adj(Lm,Um,Pm,Qm,Rm, b_nu, P.dx);
    b_b1 = zeros(N,ns); b_b1(irec,:)=(D_vi + N_f) .* EPS_B1; V_b1 = core.solve_helm_adj(Lm,Um,Pm,Qm,Rm, b_b1, P.dx);
    clear Lm Um Pm Qm Rm b_el b_vi b_nu b_b1

    I_ref = I_ref + real(sum(Us_c .* V_el(lin,:), 2)); clear V_el
    cs = Us_c .* V_vi(lin,:); cn = Us_c .* V_nu(lin,:); cb = Us_c .* V_b1(lin,:);
    Cs_sig_all(:,:,jf) = single(cs); Cs_nu_all(:,:,jf) = single(cn); Cs_b1_all(:,:,jf) = single(cb);
    I_R0_tot = I_R0_tot + real(sum(cs + cn, 2));
end
time_acoustic_mig = toc(t_ac); % includes forward gen; pure mig is ~1/3 of this, recorded explicitly below

%% ========================================================================
%  PHASE 2: Blind Dip + Legacy Straight-Ray xmir_str vs Upgraded Bent-Ray xmir_bent
%  ========================================================================
R0_2d = core.prep_image(I_R0_tot, P.nzi, P.nxi, z_in, P.dx, cfg.bg_sigma_m, cfg.z_mute_m);
[phi_eff, c_st] = core.blind_dip(R0_2d, P.dx, cfg.dip_sigma_m, cfg.dip_cap_deg);

[w_str, TS_str, xm_str] = core.build_dsa_str( ...
    Xin(:), Zin(:), P.z(P.izs), P.z(P.izs), xs, xr_lo, xr_hi, phi_eff(:), ...
    F_ts, Zi(:), Xi(:), P.nz, P.nx, cfg);
[w_bent, TS_bent, db] = core.build_dsa_bent( ...
    v_mig, P, P.z(P.izs), P.z(P.izs), xs, xr_lo, xr_hi, phi_eff, ...
    F_T, F_ts, cfg);

% Geometric discrepancy audit (Str vs Bent) inside W_deep
wd_mat = repmat(w_deep_mask, 1, ns);
act_wd = wd_mat & ((w_str > 0.1) | (w_bent > 0.1)) & (~db.turned);
dxmir_wd = abs(double(db.xmir(act_wd)) - double(xm_str(act_wd)));
dts_wd   = abs(double(TS_bent(act_wd)) - double(TS_str(act_wd))) * 1e3;
xmir_stats = struct( ...
    'W_deep_med_dxmir_m', median(dxmir_wd), 'W_deep_p90_dxmir_m', prctile(dxmir_wd, 90), ...
    'W_deep_gate_flip_pct', mean((w_str(wd_mat) > 0.5) ~= (w_bent(wd_mat) > 0.5)) * 100, ...
    'W_deep_med_dts_ms', median(dts_wd), 'W_deep_p90_dts_ms', prctile(dts_wd, 90), ...
    'W_deep_turned_pct', mean(db.turned(wd_mat)) * 100);
dip_stats = struct('med_abs_phi_deg', median(abs(rad2deg(phi_eff(w_deep_mask)))), ...
                   'p90_abs_phi_deg', prctile(abs(rad2deg(phi_eff(w_deep_mask))), 90), ...
                   'med_coh', median(c_st(w_deep_mask)));

% Blind B=300 & B=1 noise priors in Str and Bent acoustic domains
I_str_s = zeros(nin, nf, 'single'); I_str_n = I_str_s; I_str_t = I_str_s; I_str_f1 = I_str_s; I_str_f3 = I_str_s;
I_bnt_s = zeros(nin, nf, 'single'); I_bnt_n = I_bnt_s; I_bnt_t = I_bnt_s; I_bnt_f1 = I_bnt_s; I_bnt_f3 = I_bnt_s;
w_str_d = double(w_str); w_bnt_d = double(w_bent);
for jf = 1:nf
    cs = double(Cs_sig_all(:,:,jf)); cn = double(Cs_nu_all(:,:,jf)); ct = cs + cn; cb = double(Cs_b1_all(:,:,jf));
    I_str_s(:,jf)=single(real(sum(w_str_d.*cs,2))); I_str_n(:,jf)=single(real(sum(w_str_d.*cn,2)));
    I_str_t(:,jf)=single(real(sum(w_str_d.*ct,2))); I_str_f1(:,jf)=single(real(sum(w_str_d.*cb,2)));
    I_str_f3(:,jf)=single(real(sum((w_str_d.*ct).*sgn_B300,2)));

    I_bnt_s(:,jf)=single(real(sum(w_bnt_d.*cs,2))); I_bnt_n(:,jf)=single(real(sum(w_bnt_d.*cn,2)));
    I_bnt_t(:,jf)=single(real(sum(w_bnt_d.*ct,2))); I_bnt_f1(:,jf)=single(real(sum(w_bnt_d.*cb,2)));
    I_bnt_f3(:,jf)=single(real(sum((w_bnt_d.*ct).*sgn_B300,2)));
end
clear Cs_b1_all

iv_str_b3 = core.band_inv_snr(I_str_t, I_str_f3, freqs, P.nzi, P.nxi, P.dx, cfg.snr_band_hz, cfg.snr_box_m, cfg.snr_cap);
iv_bnt_b1 = core.band_inv_snr(I_bnt_t, I_bnt_f1, freqs, P.nzi, P.nxi, P.dx, cfg.snr_band_hz, cfg.snr_box_m, cfg.snr_cap);
iv_bnt_b3 = core.band_inv_snr(I_bnt_t, I_bnt_f3, freqs, P.nzi, P.nxi, P.dx, cfg.snr_band_hz, cfg.snr_box_m, cfg.snr_cap);
iv_tr_bnt = core.true_inv_snr(I_bnt_s, I_bnt_n, freqs, P.nzi, P.nxi, P.dx, cfg.snr_band_hz, cfg.snr_box_m, cfg.snr_cap);

fa = (freqs >= 10.0 & freqs <= 35.0);
l2e = @(e, t) median(arrayfun(@(j) median(log2(max(double(e(w_deep_mask,j)),1e-6)./max(double(t(w_deep_mask,j)),1e-6))), find(fa)));
noise_acc = struct('DSA_B1_log2err', l2e(iv_bnt_b1, iv_tr_bnt), 'DSA_B300_log2err', l2e(iv_bnt_b3, iv_tr_bnt));

fprintf('=== WP2 %s (%s, snr_db=%+.1f dB) ===\n', mname, q_desc, snr_db);
fprintf('  xmir_bent vs xmir_str (W_deep): med|dx|=%.1fm (p90=%.1fm) | flip=%.1f%% | p90|dTS|=%.2fms | turned=%.1f%%\n', ...
    xmir_stats.W_deep_med_dxmir_m, xmir_stats.W_deep_p90_dxmir_m, xmir_stats.W_deep_gate_flip_pct, ...
    xmir_stats.W_deep_p90_dts_ms, xmir_stats.W_deep_turned_pct);
fprintf('  Noise Acc (Bent-DSA): B1=%+.2f | B300=%+.2f log2\n', noise_acc.DSA_B1_log2err, noise_acc.DSA_B300_log2err);

%% ========================================================================
%  PHASE 3: Evaluate All Acoustic / Eikonal Ray Methods (7 methods)
%  ========================================================================
mnames = { ...
    'RAW_R0', 'DSA_Str_R0', 'DSA_Bent_R0', ...
    'WE_Q_Adj', 'WE_Q_DispOnly', 'DSA_WE_DispOnly', ...
    'WE_QRTM_Unreg', 'WE_QRTM_ZhuLP', 'WE_QRTM_Reg', ...
    'DSA_WE_QRTM_Reg', 'DSA_WE_QRTM_W1b', ...
    'DSA_Str_IQ70_25', 'DSA_Str_UQW1b_25', ...
    'DSA_Bent_W1b_25', 'DSA_Bent_UQW1b_25', 'WE_Disp_UQW1b_25'};
r0_list = {'RAW_R0', 'DSA_Str_R0', 'DSA_Bent_R0'};

IS = struct(); IN = struct(); IS_nf0 = struct();
for k = 1:numel(mnames)
    nc = 1 + (nm - 1) * ~ismember(mnames{k}, r0_list);
    IS.(mnames{k}) = zeros(nin, nc);
    IN.(mnames{k}) = zeros(nin, nc);
    IS_nf0.(mnames{k}) = zeros(nin, 1);
end

t_uq_eval = tic;
xmax = 1.05 * max(max(double(TS_str(:))), max(double(TS_bent(:)))) / min(ms);
xg = linspace(0, xmax, 1201); dxg = xg(2);
TS_str_d  = double(TS_str);
TS_bent_d = double(TS_bent);

for jf = 1:nf
    f = freqs(jf);
    Cs_s = double(Cs_sig_all(:,:,jf)); Cs_n = double(Cs_nu_all(:,:,jf));
    T0_25 = uq_tab(f, xg, 0, fref);
    T1_25 = uq_tab(f, xg, sig, fref);

    IS.RAW_R0      = IS.RAW_R0      + real(sum(Cs_s, 2));          IN.RAW_R0      = IN.RAW_R0      + real(sum(Cs_n, 2));
    IS.DSA_Str_R0  = IS.DSA_Str_R0  + real(sum(w_str_d .* Cs_s, 2)); IN.DSA_Str_R0  = IN.DSA_Str_R0  + real(sum(w_str_d .* Cs_n, 2));
    IS.DSA_Bent_R0 = IS.DSA_Bent_R0 + real(sum(w_bnt_d .* Cs_s, 2)); IN.DSA_Bent_R0 = IN.DSA_Bent_R0 + real(sum(w_bnt_d .* Cs_n, 2));

    ivs3 = double(iv_str_b3(:, jf));
    ivb3 = double(iv_bnt_b3(:, jf));

    for im = 1:nm
        % Legacy Straight-Ray xmir_str methods
        u = (TS_str_d / ms(im)) / dxg; i0 = min(floor(u), numel(xg)-2); ww = u - i0; i0 = i0 + 1;
        it = @(tb) tb(i0).*(1-ww) + tb(i0+1).*ww;
        A0 = it(T0_25.A); C0 = it(T0_25.C); P0 = it(T0_25.P);
        A1 = it(T1_25.A); C1 = it(T1_25.C); P1 = it(T1_25.P);
        Hiq_s = (A0 + s70*P0.*C0) ./ (1 + s70*C0);
        Huq_s = (A1 + s70*P1.*C1 + P1.*ivs3) ./ (1 + s70*C1 + ivs3);
        IS.DSA_Str_IQ70_25(:,im)  = IS.DSA_Str_IQ70_25(:,im)  + real(sum((w_str_d .* Hiq_s) .* Cs_s, 2));
        IN.DSA_Str_IQ70_25(:,im)  = IN.DSA_Str_IQ70_25(:,im)  + real(sum((w_str_d .* Hiq_s) .* Cs_n, 2));
        IS.DSA_Str_UQW1b_25(:,im) = IS.DSA_Str_UQW1b_25(:,im) + real(sum((w_str_d .* Huq_s) .* Cs_s, 2));
        IN.DSA_Str_UQW1b_25(:,im) = IN.DSA_Str_UQW1b_25(:,im) + real(sum((w_str_d .* Huq_s) .* Cs_n, 2));
        if im == m1_idx
            Huq_s0 = (A1 + s70*P1.*C1) ./ (1 + s70*C1);
            IS_nf0.DSA_Str_UQW1b_25 = IS_nf0.DSA_Str_UQW1b_25 + real(sum((w_str_d .* Huq_s0) .* Cs_s, 2));
        end

        % Upgraded Bent-Ray xmir_bent methods
        u = (TS_bent_d / ms(im)) / dxg; i0 = min(floor(u), numel(xg)-2); ww = u - i0; i0 = i0 + 1;
        it = @(tb) tb(i0).*(1-ww) + tb(i0+1).*ww;
        A0 = it(T0_25.A); C0 = it(T0_25.C); P0 = it(T0_25.P);
        A1 = it(T1_25.A); C1 = it(T1_25.C); P1 = it(T1_25.P);
        Hw1_b = (A0 + s70*P0.*C0 + P0.*ivb3) ./ (1 + s70*C0 + ivb3);
        Huq_b = (A1 + s70*P1.*C1 + P1.*ivb3) ./ (1 + s70*C1 + ivb3);
        IS.DSA_Bent_W1b_25(:,im)   = IS.DSA_Bent_W1b_25(:,im)   + real(sum((w_bnt_d .* Hw1_b) .* Cs_s, 2));
        IN.DSA_Bent_W1b_25(:,im)   = IN.DSA_Bent_W1b_25(:,im)   + real(sum((w_bnt_d .* Hw1_b) .* Cs_n, 2));
        IS.DSA_Bent_UQW1b_25(:,im) = IS.DSA_Bent_UQW1b_25(:,im) + real(sum((w_bnt_d .* Huq_b) .* Cs_s, 2));
        IN.DSA_Bent_UQW1b_25(:,im) = IN.DSA_Bent_UQW1b_25(:,im) + real(sum((w_bnt_d .* Huq_b) .* Cs_n, 2));
        if im == m1_idx
            Hw1_b0 = (A0 + s70*P0.*C0) ./ (1 + s70*C0);
            Huq_b0 = (A1 + s70*P1.*C1) ./ (1 + s70*C1);
            IS_nf0.DSA_Bent_W1b_25   = IS_nf0.DSA_Bent_W1b_25   + real(sum((w_bnt_d .* Hw1_b0) .* Cs_s, 2));
            IS_nf0.DSA_Bent_UQW1b_25 = IS_nf0.DSA_Bent_UQW1b_25 + real(sum((w_bnt_d .* Huq_b0) .* Cs_s, 2));
        end
    end
end
clear Cs_sig_all Cs_nu_all
time_uq_eval = toc(t_uq_eval);

%% ========================================================================
%  PHASE 4: True Viscoacoustic Wave-Equation Q-RTM Showdown (9 WE methods!)
%  ========================================================================
% Precompute bent-ray specular shot interpolation weights for exact wave-equation
% receiver-leg attenuation envelope b_R(x, xmir_bent, f) from noise-free forward b_S!
[xs_sort, ord_xs] = sort(xs(:));
xm_bnt_cl = min(max(double(db.xmir), xs_sort(1)), xs_sort(end));
q_rec = interp1(xs_sort, 1:ns, xm_bnt_cl);
k0_rec = min(floor(q_rec), ns - 1);
ww_rec = q_rec - k0_rec;
row_v = repmat((1:nin)', 1, ns);
idx_r0 = sub2ind([nin ns], row_v, ord_xs(k0_rec));
idx_r1 = sub2ind([nin ns], row_v, ord_xs(k0_rec + 1));

% Zhu et al. (2014) raised-cosine frequency taper on amplitude-loss reversal [18, 28] Hz
w_lp_vec = core.cos_taper(28.0 - freqs, 0.0, 10.0); % 1 for f<=18, cosine 18..28, 0 for f>=28

t_we = tic;
we_lu_counts = struct('WE_Q_Adj', 0, 'WE_Q_DispOnly', 0, 'WE_QRTM_Unreg', 0, 'WE_QRTM_ZhuLP_extra', 0);

% For exact wave-eq noise-adaptive DSA_WE_QRTM_W1b & WE_Disp_UQW1b_25, we need B=300
% inverse SNR iv_we_b3 in the WE_Q_DispOnly Bent-DSA domain for each m!
% We compute WE_Q_DispOnly + WE_Q_Adj + WE_QRTM_Unreg per m in a clean 2-step pass:
%   Step 4a (per m): solve A_disp, A_att, A_comp, and (for transition freqs only) A_zhu;
%            store per-freq Disp/Comp/b_WE2 slices (compact single, ~34*16k*49 only if needed?
%            Even better: stream per-m inside loop -> only 1 m in RAM at a time: < 220 MB!)
xg_we = linspace(0, 0.35, 1401); dxg_we = xg_we(2); % wave-eq t*_WE grid [0, 350 ms]

for im = 1:nm
    Q_est = ms(im) * Q_mig_base;
    t_m = tic;

    % Buffers for THIS m only (compact single: nin x nf for blind B=300 SNR, plus per-freq Disp & Comp sums)
    I_wed_t  = zeros(nin, nf, 'single');
    I_wed_f3 = zeros(nin, nf, 'single');
    % Store per-freq (Cs_disp_s, Cs_disp_n, Cs_comp_s, Cs_comp_n, b_WE2) for THIS m only (~420 MB single total!)
    Cs_ds_all = zeros(nin, ns, nf, 'single'); Cs_dn_all = Cs_ds_all;
    Cs_cs_all = zeros(nin, ns, nf, 'single'); Cs_cn_all = Cs_cs_all;
    bWE2_all  = zeros(nin, ns, nf, 'single');

    for jf = 1:nf
        f = freqs(jf); w = 2*pi*f;
        W_src = (f/fpeak)^2 * exp(-(f/fpeak)^2);
        b_src = zeros(N, ns); b_src(sub2ind([N ns], isrc, 1:ns)) = -W_src / P.dx^2;
        D_vi = double(D_vi_all(:,:,jf));
        N_f  = double(N_f_all(:,:,jf));
        b_vi = zeros(N, ns); b_vi(irec,:) = D_vi;
        b_nu = zeros(N, ns); b_nu(irec,:) = N_f;

        % 1) Dispersive-only Viscoacoustic Wave Equation A_disp
        vc_d = core.make_visco_vel(v_mig, Q_est, f, fref, P.npml, 'disp', 0.0);
        [Ld,Ud,Pd,Qd,Rd] = lu(core.helm_op_visco(vc_d, v_mig, f, P));
        we_lu_counts.WE_Q_DispOnly = we_lu_counts.WE_Q_DispOnly + 1;
        Us_d_raw = core.solve_helm_shot(Ld,Ud,Pd,Qd,Rd, b_src);
        Us_d_in  = Us_d_raw(lin, :); clear Us_d_raw
        Us_d = conj(Us_d_in) * w^2;
        Vd_s = core.solve_helm_adj(Ld,Ud,Pd,Qd,Rd, b_vi, P.dx);
        Vd_n = core.solve_helm_adj(Ld,Ud,Pd,Qd,Rd, b_nu, P.dx);
        clear Ld Ud Pd Qd Rd

        cds = Us_d .* Vd_s(lin,:); cdn = Us_d .* Vd_n(lin,:); clear Vd_s Vd_n Us_d
        IS.WE_Q_DispOnly(:,im)   = IS.WE_Q_DispOnly(:,im)   + real(sum(cds, 2));
        IN.WE_Q_DispOnly(:,im)   = IN.WE_Q_DispOnly(:,im)   + real(sum(cdn, 2));
        IS.DSA_WE_DispOnly(:,im) = IS.DSA_WE_DispOnly(:,im) + real(sum(w_bnt_d .* cds, 2));
        IN.DSA_WE_DispOnly(:,im) = IN.DSA_WE_DispOnly(:,im) + real(sum(w_bnt_d .* cdn, 2));
        cdt_b = w_bnt_d .* (cds + cdn);
        I_wed_t(:,jf)  = single(real(sum(cdt_b, 2)));
        I_wed_f3(:,jf) = single(real(sum(cdt_b .* sgn_B300, 2))); clear cdt_b
        Cs_ds_all(:,:,jf) = single(cds); Cs_dn_all(:,:,jf) = single(cdn);

        % 2) Attenuating (Adjoint) Viscoacoustic Wave Equation A_att (unconditionally stable -> also gives exact b_WE2!)
        vc_a = core.make_visco_vel(v_mig, Q_est, f, fref, P.npml, 'att', 1.0);
        [La,Ua,Pa,Qa,Ra] = lu(core.helm_op_visco(vc_a, v_mig, f, P));
        we_lu_counts.WE_Q_Adj = we_lu_counts.WE_Q_Adj + 1;
        Us_a_raw = core.solve_helm_shot(La,Ua,Pa,Qa,Ra, b_src);
        Us_a_in  = Us_a_raw(lin, :); clear Us_a_raw
        Us_a = conj(Us_a_in) * w^2;
        Va_s = core.solve_helm_adj(La,Ua,Pa,Qa,Ra, b_vi, P.dx);
        Va_n = core.solve_helm_adj(La,Ua,Pa,Qa,Ra, b_nu, P.dx);
        clear La Ua Pa Qa Ra
        IS.WE_Q_Adj(:,im) = IS.WE_Q_Adj(:,im) + real(sum(Us_a .* Va_s(lin,:), 2));
        IN.WE_Q_Adj(:,im) = IN.WE_Q_Adj(:,im) + real(sum(Us_a .* Va_n(lin,:), 2));
        clear Us_a Va_s Va_n

        % Exact finite-frequency wave-equation attenuation envelope b_S = |U_att| / |U_disp| (3x3 smoothed power)
        b_S = zeros(nin, ns);
        for ks = 1:ns
            pa = movmean(movmean(reshape(abs(Us_a_in(:,ks)).^2, P.nzi, P.nxi), 3, 1), 3, 2);
            pd = movmean(movmean(reshape(abs(Us_d_in(:,ks)).^2, P.nzi, P.nxi), 3, 1), 3, 2);
            b_S(:,ks) = min(max(sqrt(pa(:) ./ max(pd(:), 1e-30)), 1e-4), 1.0);
        end
        clear Us_a_in Us_d_in
        b_R = (1 - ww_rec) .* b_S(idx_r0) + ww_rec .* b_S(idx_r1);
        b_WE2 = min(max(b_S .* b_R, 1e-5), 1.0);
        bWE2_all(:,:,jf) = single(b_WE2);

        % 3) Compensated True Viscoacoustic Wave-Equation Q-RTM A_comp (anti-attenuating interior, real absorbing PML)
        vc_c = core.make_visco_vel(v_mig, Q_est, f, fref, P.npml, 'comp', 1.0);
        [Lc,Uc,Pc,Qc,Rc] = lu(core.helm_op_visco(vc_c, v_mig, f, P));
        we_lu_counts.WE_QRTM_Unreg = we_lu_counts.WE_QRTM_Unreg + 1;
        Us_c_raw = core.solve_helm_shot(Lc,Uc,Pc,Qc,Rc, b_src);
        Us_c = conj(Us_c_raw(lin,:)) * w^2; clear Us_c_raw
        Vc_s = core.solve_helm_adj(Lc,Uc,Pc,Qc,Rc, b_vi, P.dx);
        Vc_n = core.solve_helm_adj(Lc,Uc,Pc,Qc,Rc, b_nu, P.dx);
        clear Lc Uc Pc Qc Rc
        ccs = Us_c .* Vc_s(lin,:); ccn = Us_c .* Vc_n(lin,:); clear Us_c Vc_s Vc_n
        Cs_cs_all(:,:,jf) = single(ccs); Cs_cn_all(:,:,jf) = single(ccn);

        IS.WE_QRTM_Unreg(:,im) = IS.WE_QRTM_Unreg(:,im) + real(sum(ccs, 2));
        IN.WE_QRTM_Unreg(:,im) = IN.WE_QRTM_Unreg(:,im) + real(sum(ccn, 2));

        % 4) Zhu et al. (2014) Frequency-Tapered Stabilized True Q-RTM
        wl = w_lp_vec(jf);
        if wl >= 1.0 - 1e-9
            IS.WE_QRTM_ZhuLP(:,im) = IS.WE_QRTM_ZhuLP(:,im) + real(sum(ccs, 2));
            IN.WE_QRTM_ZhuLP(:,im) = IN.WE_QRTM_ZhuLP(:,im) + real(sum(ccn, 2));
        elseif wl <= 1e-9
            IS.WE_QRTM_ZhuLP(:,im) = IS.WE_QRTM_ZhuLP(:,im) + real(sum(cds, 2));
            IN.WE_QRTM_ZhuLP(:,im) = IN.WE_QRTM_ZhuLP(:,im) + real(sum(cdn, 2));
        else
            vc_z = core.make_visco_vel(v_mig, Q_est, f, fref, P.npml, 'comp', wl);
            [Lz,Uz,Pz,Qz,Rz] = lu(core.helm_op_visco(vc_z, v_mig, f, P));
            we_lu_counts.WE_QRTM_ZhuLP_extra = we_lu_counts.WE_QRTM_ZhuLP_extra + 1;
            Us_z = conj(core.solve_helm_shot(Lz,Uz,Pz,Qz,Rz, b_src)); Us_z = Us_z(lin,:)*w^2;
            Vz_s = core.solve_helm_adj(Lz,Uz,Pz,Qz,Rz, b_vi, P.dx);
            Vz_n = core.solve_helm_adj(Lz,Uz,Pz,Qz,Rz, b_nu, P.dx);
            clear Lz Uz Pz Qz Rz
            IS.WE_QRTM_ZhuLP(:,im) = IS.WE_QRTM_ZhuLP(:,im) + real(sum(Us_z .* Vz_s(lin,:), 2));
            IN.WE_QRTM_ZhuLP(:,im) = IN.WE_QRTM_ZhuLP(:,im) + real(sum(Us_z .* Vz_n(lin,:), 2));
            clear Us_z Vz_s Vz_n
        end

        % 5) Wang / Sun-Zhu (2018) Tikhonov Gain-Stabilized True Q-RTM (RAW and +Bent-DSA)
        b2 = b_WE2.^2;
        ws_reg = b2 ./ (b2 + s70);
        crs = ws_reg .* ccs + (1 - ws_reg) .* cds;
        crn = ws_reg .* ccn + (1 - ws_reg) .* cdn;
        IS.WE_QRTM_Reg(:,im)     = IS.WE_QRTM_Reg(:,im)     + real(sum(crs, 2));
        IN.WE_QRTM_Reg(:,im)     = IN.WE_QRTM_Reg(:,im)     + real(sum(crn, 2));
        IS.DSA_WE_QRTM_Reg(:,im) = IS.DSA_WE_QRTM_Reg(:,im) + real(sum(w_bnt_d .* crs, 2));
        IN.DSA_WE_QRTM_Reg(:,im) = IN.DSA_WE_QRTM_Reg(:,im) + real(sum(w_bnt_d .* crn, 2));
    end

    % Compute blind B=300 inverse SNR in WE_Q_DispOnly Bent-DSA domain for THIS m
    iv_wed_b3 = core.band_inv_snr(I_wed_t, I_wed_f3, freqs, P.nzi, P.nxi, P.dx, cfg.snr_band_hz, cfg.snr_box_m, cfg.snr_cap);

    % 6) Noise-Adaptive Stabilized True Wave-Eq Q-RTM (DSA_WE_QRTM_W1b) & Hybrid WE_Disp_UQW1b_25
    for jf = 1:nf
        f = freqs(jf);
        cds = double(Cs_ds_all(:,:,jf)); cdn = double(Cs_dn_all(:,:,jf));
        ccs = double(Cs_cs_all(:,:,jf)); ccn = double(Cs_cn_all(:,:,jf));
        b_WE2 = double(bWE2_all(:,:,jf)); b2 = b_WE2.^2;
        iv_w = double(iv_wed_b3(:, jf));

        % DSA_WE_QRTM_W1b: exact wave-eq noise-adaptive counterpart of W1b!
        ws_w1b = b2 ./ (b2 + s70 + b2 .* iv_w);
        IS.DSA_WE_QRTM_W1b(:,im) = IS.DSA_WE_QRTM_W1b(:,im) + real(sum(w_bnt_d .* (ws_w1b .* ccs + (1 - ws_w1b) .* cds), 2));
        IN.DSA_WE_QRTM_W1b(:,im) = IN.DSA_WE_QRTM_W1b(:,im) + real(sum(w_bnt_d .* (ws_w1b .* ccn + (1 - ws_w1b) .* cdn), 2));

        % WE_Disp_UQW1b_25: exact dispersive wave-eq + finite-freq wave-eq t*_WE inside Bayesian UQW1b amp kernel!
        % (fref = f -> k = 0 -> zero double-dispersion; pure Bayesian MMSE amplitude risk on dispersive wavefield!)
        ts_we = min(max(-log(b_WE2) / (pi * f), 0.0), xg_we(end) - 1e-6);
        T1_amp = uq_tab(f, xg_we, sig, f);
        u = ts_we / dxg_we; i0 = min(floor(u), numel(xg_we)-2); ww = u - i0; i0 = i0 + 1;
        A1a = T1_amp.A(i0).*(1-ww) + T1_amp.A(i0+1).*ww;
        C1a = T1_amp.C(i0).*(1-ww) + T1_amp.C(i0+1).*ww;
        H_hyb = (A1a + s70*C1a + iv_w) ./ (1 + s70*C1a + iv_w);
        IS.WE_Disp_UQW1b_25(:,im) = IS.WE_Disp_UQW1b_25(:,im) + real(sum((w_bnt_d .* H_hyb) .* cds, 2));
        IN.WE_Disp_UQW1b_25(:,im) = IN.WE_Disp_UQW1b_25(:,im) + real(sum((w_bnt_d .* H_hyb) .* cdn, 2));
        if im == m1_idx
            H_hyb0 = (A1a + s70*C1a) ./ (1 + s70*C1a);
            IS_nf0.WE_Disp_UQW1b_25 = IS_nf0.WE_Disp_UQW1b_25 + real(sum((w_bnt_d .* H_hyb0) .* cds, 2));
            ws_reg0 = b2 ./ (b2 + s70);
            IS_nf0.DSA_WE_QRTM_W1b = IS_nf0.DSA_WE_QRTM_W1b + real(sum(w_bnt_d .* (ws_reg0 .* ccs + (1 - ws_reg0) .* cds), 2));
        end
    end
    clear Cs_ds_all Cs_dn_all Cs_cs_all Cs_cn_all bWE2_all
    fprintf('  [WE-QRTM m=%.2f (%d/%d)] done in %.1fs\n', ms(im), im, nm, toc(t_m));
end
time_we_total = toc(t_we);

%% ========================================================================
%  PHASE 5: Score ALL 16 Methods (Noise-Free m=1 Oracle + Noisy All-m W_deep)
%  ========================================================================
ref_2d = core.prep_image(I_ref, P.nzi, P.nxi, z_in, P.dx, cfg.bg_sigma_m, cfg.z_mute_m);
ref_hf = core.highpass_kz(ref_2d, P.dx, 10.0);
rms_w = @(A, W) sqrt(mean(A(z_in>=W{2}(1)&z_in<=W{2}(2), x_in>=W{1}(1)&x_in<=W{1}(2)).^2, 'all'));
ref_rcf = log10(rms_w(ref_2d, W_targ) / max(rms_w(ref_2d, W_flank), 1e-30));

out = struct('model', mname, 'q_err_desc', q_desc, 'seed', seed, 'snr_db', snr_db, 'm', ms, ...
             'dip_stats', dip_stats, 'xmir_stats', xmir_stats, 'noise_acc', noise_acc, ...
             'compute_cost', struct('acoustic_LU_total', nf, 'we_qrtm_LU_per_m', nf, ...
                 'we_qrtm_LU_total_all_m', nf * nm, 'we_lu_counts', we_lu_counts, ...
                 'time_acoustic_plus_fwd_s', time_acoustic_mig, 'time_uq_eval_s', time_uq_eval, ...
                 'time_we_qrtm_all_s', time_we_total), ...
             'NoiseFree_m1', struct(), 'W_deep', struct(), 'elapsed_s', 0);

fprintf('\n  --- %s Results (Noise-Free m=1 Oracle | Noisy All-m W_deep) ---\n', mname);
save_imgs = struct('ref_2d', single(ref_2d), 'z_in', z_in, 'x_in', x_in, 'v_true', single(v_true(P.in{1}, P.in{2})), ...
                   'Q_true', single(Q_true(P.in{1}, P.in{2})), 'Q_mig_base', single(Q_mig_base(P.in{1}, P.in{2})));

for k = 1:numel(mnames)
    key = mnames{k}; Sm = IS.(key); Nm = IN.(key); nc = size(Sm, 2);
    c_m1 = min(m1_idx, nc);

    % Noise-free m=1 oracle score (true iv=0 for noise-adaptive methods, pure signal S2!)
    if ismember(key, {'DSA_Str_UQW1b_25','DSA_Bent_W1b_25','DSA_Bent_UQW1b_25','WE_Disp_UQW1b_25','DSA_WE_QRTM_W1b'})
        S_m1 = core.prep_image(IS_nf0.(key), P.nzi, P.nxi, z_in, P.dx, cfg.bg_sigma_m, cfg.z_mute_m);
    else
        S_m1 = core.prep_image(Sm(:, c_m1), P.nzi, P.nxi, z_in, P.dx, cfg.bg_sigma_m, cfg.z_mute_m);
    end
    S_m1_hf = core.highpass_kz(S_m1, P.dx, 10.0);
    nf_e1 = core.score_e1(S_m1, ref_2d, iz_w, ix_w);
    nf_hf = core.score_e1(S_m1_hf, ref_hf, iz_w, ix_w);
    nf_ea = abs(log10(rms_w(S_m1, W_targ) / max(rms_w(S_m1, W_flank), 1e-30)) - ref_rcf);
    out.NoiseFree_m1.(key) = struct('E1', nf_e1, 'E1_HF', nf_hf, 'E_amp', nf_ea);

    % Noisy all-m W_deep metrics
    e1 = zeros(1, nc); e1h = e1; eamp = e1; snr = e1; smf = e1;
    for c = 1:nc
        S2 = core.prep_image(Sm(:,c), P.nzi, P.nxi, z_in, P.dx, cfg.bg_sigma_m, cfg.z_mute_m);
        N2 = core.prep_image(Nm(:,c), P.nzi, P.nxi, z_in, P.dx, cfg.bg_sigma_m, cfg.z_mute_m);
        T2 = S2 + N2; Th = core.highpass_kz(T2, P.dx, 10.0);
        e1(c)   = core.score_e1(T2, ref_2d, iz_w, ix_w);
        e1h(c)  = core.score_e1(Th, ref_hf, iz_w, ix_w);
        snr(c)  = 10*log10(sum(S2(iz_w,ix_w).^2,'all') / max(sum(N2(iz_w,ix_w).^2,'all'), 1e-30));
        eamp(c) = abs(log10(rms_w(T2, W_targ) / max(rms_w(T2, W_flank), 1e-30)) - ref_rcf);
        smf(c)  = core.steep_smile_frac(T2(ovb_z, ovb_x), P.dx, 35.0, 4.0, 35.0);
        if c == c_m1 && ismember(key, {'RAW_R0','WE_QRTM_Unreg','WE_QRTM_Reg','DSA_WE_QRTM_W1b','DSA_Str_UQW1b_25','DSA_Bent_UQW1b_25','WE_Disp_UQW1b_25'})
            save_imgs.([key '_m1']) = single(T2);
        end
        if c == 1 && nc > 1 && ismember(key, {'WE_QRTM_Unreg','WE_QRTM_Reg','DSA_WE_QRTM_W1b','DSA_Bent_UQW1b_25','WE_Disp_UQW1b_25'})
            save_imgs.([key '_mlow']) = single(T2);
        end
    end
    out.W_deep.(key) = struct('E1', e1, 'E1_mean', mean(e1), 'E1_HF', e1h, 'E1_HF_mean', mean(e1h), ...
                              'E_amp', eamp, 'E_amp_mean', mean(eamp), ...
                              'SNR_img_dB', snr, 'SNR_img_dB_mean', mean(snr), ...
                              'Smile_Frac', smf, 'Smile_Frac_mean', mean(smf));
    fprintf('  %-18s | NF(m=1): E1=%.4f Eamp=%.3f | Noisy: E1=%.4f (m_lo=%.3f, m=1=%.3f) HF=%.4f Eamp=%.3f SNR=%+5.2fdB Sm=%.1f%%\n', ...
        key, nf_e1, nf_ea, mean(e1), e1(1), e1(c_m1), mean(e1h), mean(eamp), mean(snr), mean(smf)*100);
end

out.elapsed_s = toc(t_all);
img_dir = fullfile('D:\tmp_uqrtm_xmir_qrtm_20261126');
if ~isfolder(img_dir), mkdir(img_dir); end
save(fullfile(img_dir, ['imgs_' mname '.mat']), '-struct', 'save_imgs', '-v7.3');
end

function [P, v_true, v_mig, Q_true, Q_mig_base, rec_mask, W_targ, W_flank, ms, seed, q_desc] = build_model(mname, core)
P.dx = 10.0; P.npml = 18; P.nzi = 100; P.nxi = 160;
P.nz = P.nzi + 2*P.npml; P.nx = P.nxi + 2*P.npml;
P.z = ((1:P.nz)' - P.npml - 1) * P.dx; P.x = ((1:P.nx) - P.npml - 1) * P.dx;
P.in = {P.npml + (1:P.nzi), P.npml + (1:P.nxi)}; P.izs = P.npml + 3;
P.ixs = P.npml + 1 + (8:3:152); P.ixr = P.npml + 1 + (0:1:P.nxi-1);
ns = numel(P.ixs); nr = numel(P.ixr);
[Zfull, Xfull] = ndgrid(P.z, P.x);
xs = P.x(P.ixs); xr = P.x(P.ixr);
rec_mask = true(nr, ns);
ms = [0.5, 0.7, 1.0, 1.4, 2.0];

switch mname
    case 'S1_Chimney_Base'
        seed = 20261501; q_desc = 'scalar_m_chimney_base';
        z1 = 230.0 - 14.0 * exp(-((Xfull - 770.0)/330.0).^2);
        z2 = 625.0 - 36.0 * exp(-((Xfull - 790.0)/300.0).^2) + 0.08*(Xfull - 800.0);
        z3 = 760.0 - 26.0 * exp(-((Xfull - 790.0)/330.0).^2) + 0.06*(Xfull - 800.0);
        z4 = 890.0 - 14.0 * exp(-((Xfull - 800.0)/370.0).^2);
        v_true = 1850.0 * ones(P.nz, P.nx); v_true(Zfull>=z1)=2160; v_true(Zfull>=z2)=2580; v_true(Zfull>=z3)=2940; v_true(Zfull>=z4)=3320;
        invQ_t = (1/90)*ones(P.nz, P.nx); invQ_t(Zfull < 80) = 1/1000;
        invQ_t = invQ_t + (1/14 - 1/90)*exp(-(((Xfull - 790)/210).^2 + ((Zfull - 415)/140).^2).^1.5);
        Q_true = 1 ./ invQ_t; Q_mig_base = Q_true;
        W_targ = {[580 1000], [560 950]}; W_flank = {[200 480], [560 950]};
    case 'S2_StrongGrad_SteepDip'
        seed = 20261502; q_desc = 'strong_grad_24deg_dip_asymmetric_cable';
        v_bg = 1550.0 + 2.1 * max(Zfull, 0.0) + 0.35 * (Xfull - 800.0) .* (max(Zfull,0)/1000.0);
        Zg = Zfull;
        hw = (Xfull < 840.0 - (Zfull - 180.0) / tand(62.0)) & (Zfull > 180.0); Zg(hw) = Zg(hw) + 35.0;
        z1 = 220.0 - 0.06*(Xfull - 800);
        z2 = 640.0 + 0.44*(Xfull - 800);
        z3 = 775.0 - 0.36*(Xfull - 800);
        z4 = 895.0 + 0.15*(Xfull - 800);
        v_true = v_bg;
        v_true(Zg>=z1) = v_true(Zg>=z1) + 220.0;
        v_true(Zg>=z2) = v_true(Zg>=z2) + 340.0;
        v_true(Zg>=z3) = v_true(Zg>=z3) + 320.0;
        v_true(Zg>=z4) = v_true(Zg>=z4) + 300.0;
        v_true = min(max(v_true, 1500.0), 4200.0);
        invQ_t = (1/85)*ones(P.nz, P.nx); invQ_t(Zfull < 80) = 1/1000;
        invQ_t = invQ_t + (1/15 - 1/85)*exp(-(((Xfull - 880)/230).^2 + ((Zfull - 380)/120).^2).^1.5);
        Q_true = 1 ./ invQ_t; Q_mig_base = Q_true;
        for k = 1:ns
            rec_mask(:, k) = (xr(:) >= xs(k) - 120.0) & (xr(:) <= xs(k) + 760.0);
        end
        W_targ = {[700 1120], [560 950]}; W_flank = {[220 540], [560 950]};
    case 'S3_SharpQ_Multipath'
        seed = 20261503; q_desc = 'sharp_Q12_and_minus22pct_velocity_lens_multipath';
        z1 = 220.0 - 10.0 * exp(-((Xfull - 800.0)/350.0).^2);
        z2 = 640.0 - 20.0 * exp(-((Xfull - 800.0)/320.0).^2) + 0.05*(Xfull - 800.0);
        z3 = 770.0 - 15.0 * exp(-((Xfull - 800.0)/340.0).^2) + 0.04*(Xfull - 800.0);
        z4 = 895.0 * ones(P.nz, P.nx);
        v_true = 1900.0 * ones(P.nz, P.nx); v_true(Zfull>=z1)=2250; v_true(Zfull>=z2)=2700; v_true(Zfull>=z3)=3050; v_true(Zfull>=z4)=3400;
        lens = exp(-(((Xfull - 800.0)/130.0).^2 + ((Zfull - 410.0)/85.0).^2).^2.0);
        v_true = v_true .* (1.0 - 0.22 * lens);
        invQ_t = (1/95)*ones(P.nz, P.nx); invQ_t(Zfull < 80) = 1/1000;
        invQ_t = invQ_t + (1/12 - 1/95) * lens;
        Q_true = 1 ./ invQ_t; Q_mig_base = Q_true;
        W_targ = {[620 980], [560 950]}; W_flank = {[200 480], [560 950]};
    case 'S4_NonScalar_SpatialQ'
        seed = 20261504; q_desc = 'nonscalar_shifted160m_widened1.35x_plus_shallow1.3_deep0.75';
        ms = [0.6, 1.0, 1.4];
        z1 = 230.0 - 14.0 * exp(-((Xfull - 770.0)/330.0).^2);
        z2 = 625.0 - 36.0 * exp(-((Xfull - 790.0)/300.0).^2) + 0.08*(Xfull - 800.0);
        z3 = 760.0 - 26.0 * exp(-((Xfull - 790.0)/330.0).^2) + 0.06*(Xfull - 800.0);
        z4 = 890.0 - 14.0 * exp(-((Xfull - 800.0)/370.0).^2);
        v_true = 1850.0 * ones(P.nz, P.nx); v_true(Zfull>=z1)=2160; v_true(Zfull>=z2)=2580; v_true(Zfull>=z3)=2940; v_true(Zfull>=z4)=3320;
        invQ_t = (1/90)*ones(P.nz, P.nx); invQ_t(Zfull < 80) = 1/1000;
        invQ_t = invQ_t + (1/14 - 1/90)*exp(-(((Xfull - 790)/210).^2 + ((Zfull - 415)/140).^2).^1.5);
        Q_true = 1 ./ invQ_t;
        invQ_m = (1/90)*ones(P.nz, P.nx); invQ_m(Zfull < 80) = 1/1000;
        invQ_m = invQ_m + (1/14 - 1/90)*exp(-(((Xfull - 950)/(210*1.35)).^2 + ((Zfull - 415)/140).^2).^1.5);
        Q_mig_base = 1 ./ invQ_m;
        Q_mig_base(Zfull >= 80 & Zfull < 480) = 1.30 * Q_mig_base(Zfull >= 80 & Zfull < 480);
        Q_mig_base(Zfull >= 480)              = 0.75 * Q_mig_base(Zfull >= 480);
        W_targ = {[580 1000], [560 950]}; W_flank = {[200 480], [560 950]};
end
v_mig = core.smooth2d(v_true, 2.5);
end
