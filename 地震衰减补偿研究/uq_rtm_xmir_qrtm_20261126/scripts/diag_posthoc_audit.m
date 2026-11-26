function diag_posthoc_audit()
% POST-HOC ROOT-CAUSE DIAGNOSTICS & SENSITIVITY AUDIT (Explicitly labeled Post-Hoc, NOT pre-registered!)
%   Diag A: Volve 4-25 Hz root cause of DSA_Bent F-5 drop (0.5704 vs Str 0.6222) & F-14 jump (+0.275 vs +0.224)
%           -> Factor-by-factor decomposition of w_bent = t_arr * t_inc * t_up * t_s_srf * t_r_srf at wells F-5 & F-14
%           -> Plus ablation of inc_max_deg (55 deg vs 70 deg vs off) on Volve cached acoustic wavefield?
%              We have saved volve25_geom_str_bent.mat & F_T, so we decompose the exact weight factors at F-5 & F-14!
%   Diag B: S2_StrongGrad_SteepDip root cause of DSA zero-noise -0.07 drop & B300 +1.90 log2 bias
%           -> Structure-tensor dip bias on opposing +24/-20 deg layers & active fold n_eff collapse under asymmetric cable!
%   Diag C: Post-hoc Tighter Gain-Cap Sensitivity Sweep (10 / 15 / 20 / 30 / 70 dB) for True Wave-Equation Q-RTM
%           on S1_Chimney_Base & S3_SharpQ_Multipath (testing whether tighter fixed gain caps rescue WE-QRTM at m=0.5!)

maxNumCompThreads(10);
here = fileparts(mfilename('fullpath'));
root = fileparts(here);
res_dir = fullfile(root, 'results');
addpath('D:\ChatGPT\研究\地震衰减补偿研究\uq_rtm_20261117\scripts');
core = uq_rtm_xmir_qrtm_core();

out = struct();

%% ========================================================================
%  DIAG A: Volve 4-25 Hz Factor-by-Factor Decomposition of w_bent at Wells F-5 & F-14
%  ========================================================================
fprintf('=== POST-HOC DIAG A: Volve 4-25 Hz DSA_Bent F-5 Drop vs F-14 Jump Root Cause ===\n');
Pf = fld_setup();
cfg_f = core.default_cfg('field_obc');
z_f = Pf.z(Pf.in{1}); x_f = Pf.x(Pf.in{2});
[Zfi, Xfi] = ndgrid(z_f, x_f);
w5  = core.load_well(fullfile(Pf.cache, 'well_F-5.bin'),  z_f, x_f);
w14 = core.load_well(fullfile(Pf.cache, 'well_F-14.bin'), z_f, x_f);

S1 = load(fullfile(Pf.cache, 'state1.mat'), 'I');
R0_f = core.prep_image(double(sum(S1.I, 2)), Pf.nzi, Pf.nxi, z_f, Pf.dx, cfg_f.bg_sigma_m, cfg_f.z_mute_m);
clear S1
[phi_f, c_st_f] = core.blind_dip(R0_f, Pf.dx, cfg_f.dip_sigma_m, cfg_f.dip_cap_deg);

% Recompute F_T (fast: ~20s) to inspect exact per-factor values at F-5 and F-14 well trajectories!
s_f = 1.0 ./ Pf.vmig; ns = Pf.nshot;
F_Tf = zeros(Pf.nz, Pf.nx, ns, 'single');
for k = 1:ns
    F_Tf(:,:,k) = single(eikonal_fs(s_f, Pf.dx, Pf.src_iz(k), Pf.src_ix(k)));
end
TF = load(fullfile(Pf.cache, 'tfields.mat'), 'F');
rx = Pf.x(Pf.rec_ix);
xr_lo = arrayfun(@(s) min(rx(Pf.shot_of == s)), 1:ns);
xr_hi = arrayfun(@(s) max(rx(Pf.shot_of == s)), 1:ns);
[Zi_f, Xi_f] = ndgrid(Pf.in{1}, Pf.in{2});

[w_str, ~, xm_str, ths_str, thr_str] = core.build_dsa_str(Xfi(:), Zfi(:), Pf.z(Pf.src_iz), 90.0, Pf.src_x, xr_lo, xr_hi, phi_f(:), TF.F, Zi_f(:), Xi_f(:), Pf.nz, Pf.nx, cfg_f);
[w_bnt, ~, db] = core.build_dsa_bent(Pf.vmig, Pf, Pf.z(Pf.src_iz), 90.0, Pf.src_x, xr_lo, xr_hi, phi_f, F_Tf, TF.F, cfg_f);
clear F_Tf TF

% Extract well trajectory interior linear indices for F-5 and F-14
idx_f5  = sub2ind([Pf.nzi Pf.nxi], w5.iz(:),  w5.ix(:));
idx_f14 = sub2ind([Pf.nzi Pf.nxi], w14.iz(:), w14.ix(:));
v_in_f = double(Pf.vmig(Pf.in{1}, Pf.in{2}));

diagA = struct();
wlist = {{'F5', idx_f5, w5}, {'F14', idx_f14, w14}};
for iw = 1:2
    wn = wlist{iw}{1}; idw = wlist{iw}{2}; ww = wlist{iw}{3};
    nw = numel(idw);
    % Decompose w_bent factors along well trajectory
    t_arr_b = zeros(nw, ns); t_inc_b = t_arr_b; t_up_b = t_arr_b; t_ss_b = t_arr_b; t_rs_b = t_arr_b;
    t_arr_s = zeros(nw, ns); t_gs_s  = t_arr_s; t_gr_s = t_arr_s;
    phi_w = rad2deg(phi_f(idw));
    for k = 1:ns
        xm_b = double(db.xmir(idw, k)); tr_b = db.turned(idw, k);
        d_in_b = min(xm_b - xr_lo(k), xr_hi(k) - xm_b); d_in_b(tr_b) = -1e6;
        t_arr_b(:,k) = core.cos_taper(d_in_b, cfg_f.arr_taper_m(1), cfg_f.arr_taper_m(2));
        th_s_b = double(db.th_s_deg(idw, k));
        th_r_b = double(db.th_r_deg(idw, k));
        th_rs  = double(db.th_r_surf_deg(idw, k));
        th_inc = th_s_b + phi_w;
        t_inc_b(:,k) = core.cos_taper(cfg_f.inc_max_deg - abs(th_inc), 0, cfg_f.theta_tap_deg);
        t_up_b(:,k)  = core.cos_taper(84.0 - abs(th_r_b), 0, 8.0);
        t_rs_b(:,k)  = core.cos_taper(cfg_f.surf_max_deg - abs(th_rs), 0, cfg_f.theta_tap_deg);
        % Source surf angle
        v_src = double(Pf.vmig(Pf.src_iz(k), Pf.src_ix(k)));
        th_ss = rad2deg(asin(min(max(sin(deg2rad(th_s_b)) ./ v_in_f(idw) * v_src, -1), 1)));
        t_ss_b(:,k)  = core.cos_taper(cfg_f.surf_max_deg - abs(th_ss), 0, cfg_f.theta_tap_deg);

        xm_s = double(xm_str(idw, k));
        d_in_s = min(xm_s - xr_lo(k), xr_hi(k) - xm_s);
        t_arr_s(:,k) = core.cos_taper(d_in_s, cfg_f.arr_taper_m(1), cfg_f.arr_taper_m(2));
        t_gs_s(:,k)  = core.cos_taper(cfg_f.theta_max_deg - abs(double(ths_str(idw, k))), 0, cfg_f.theta_tap_deg);
        t_gr_s(:,k)  = core.cos_taper(cfg_f.theta_max_deg - abs(double(thr_str(idw, k))), 0, cfg_f.theta_tap_deg);
    end
    fold_str_w  = median(sum(double(w_str(idw,:)), 2));
    fold_bnt_w  = median(sum(double(w_bnt(idw,:)), 2));
    fold_arr_s  = median(sum(t_arr_s, 2));
    fold_arr_b  = median(sum(t_arr_b, 2));
    fold_arr_inc = median(sum(t_arr_b .* t_inc_b, 2));
    fold_arr_srf = median(sum(t_arr_b .* t_ss_b .* t_rs_b, 2));
    % Shots killed by t_inc alone (t_arr_b > 0.5 & t_inc_b < 0.5)
    killed_by_inc  = mean(sum((t_arr_b > 0.5) & (t_inc_b < 0.5), 2));
    killed_by_srf  = mean(sum((t_arr_b > 0.5) & ((t_ss_b .* t_rs_b) < 0.5), 2));
    killed_by_turn = mean(sum(db.turned(idw,:), 2));
    med_v_w = median(v_in_f(idw));
    med_phi_w = median(phi_w);
    med_ths_eik_active = median(abs(double(db.th_s_deg(idw, :)) .* (t_arr_b > 0.5) ./ max(t_arr_b > 0.5, eps)), 'all', 'omitnan');
    act_mask = (t_arr_b > 0.5);
    th_eik_vals = abs(double(db.th_s_deg(idw, :)));
    th_str_vals = abs(double(ths_str(idw, :)));
    diagA.(wn) = struct('x_med_m', median(ww.sw), 'z_range_m', [min(ww.zw), max(ww.zw)], ...
        'v_med_m_s', med_v_w, 'phi_med_deg', med_phi_w, ...
        'fold_str_total', fold_str_w, 'fold_bent_total', fold_bnt_w, ...
        'fold_str_arr_only', fold_arr_s, 'fold_bent_arr_only', fold_arr_b, ...
        'fold_bent_arr_x_inc55', fold_arr_inc, 'fold_bent_arr_x_surf40', fold_arr_srf, ...
        'shots_killed_by_inc55', killed_by_inc, 'shots_killed_by_surf40', killed_by_srf, ...
        'shots_turned', killed_by_turn, ...
        'p50_p90_ths_eik_in_arr_deg', [prctile(th_eik_vals(act_mask), 50), prctile(th_eik_vals(act_mask), 90)], ...
        'p50_p90_ths_str_in_arr_deg', [prctile(th_str_vals(act_mask), 50), prctile(th_str_vals(act_mask), 90)]);
    fprintf('  [%s] x=%.0fm z=[%.0f,%.0f]m v_med=%.0fm/s phi_med=%+.1f deg:\n', wn, median(ww.sw), min(ww.zw), max(ww.zw), med_v_w, med_phi_w);
    fprintf('       Fold: Str_total=%.1f (arr=%.1f) vs Bent_total=%.1f (arr=%.1f, arr*inc55=%.1f, arr*surf40=%.1f)\n', ...
        fold_str_w, fold_arr_s, fold_bnt_w, fold_arr_b, fold_arr_inc, fold_arr_srf);
    fprintf('       Active shots killed by inc55=%.1f | killed by surf40=%.1f | turned=%.1f\n', ...
        killed_by_inc, killed_by_srf, killed_by_turn);
    fprintf('       In-array |th_s|: Str p50/p90=%.1f/%.1f deg vs Bent(eik) p50/p90=%.1f/%.1f deg\n', ...
        prctile(th_str_vals(act_mask),50), prctile(th_str_vals(act_mask),90), ...
        prctile(th_eik_vals(act_mask),50), prctile(th_eik_vals(act_mask),90));
end
out.diagA_volve_well_factor_decomp = diagA;

%% ========================================================================
%  DIAG B: S2_StrongGrad_SteepDip Root Cause (Dip Bias + Asymmetric Fold Collapse)
%  ========================================================================
fprintf('\n=== POST-HOC DIAG B: S2_StrongGrad_SteepDip DSA Zero-Noise Drop & B300 +1.90 log2 Root Cause ===\n');
I2 = load('D:\tmp_uqrtm_xmir_qrtm_20261126\imgs_S2_StrongGrad_SteepDip.mat');
cfg_s = core.default_cfg('synth_10m');
[phi_s2, c_s2] = core.blind_dip(double(I2.RAW_R0_m1), 10.0, cfg_s.dip_sigma_m, cfg_s.dip_cap_deg);
[phi_s2_nf, ~] = core.blind_dip(double(I2.ref_2d),    10.0, cfg_s.dip_sigma_m, cfg_s.dip_cap_deg);
z_s2 = I2.z_in; x_s2 = I2.x_in; [Z2, X2] = ndgrid(z_s2, x_s2);
% True dip of z2 = 640 + 0.44*(x-800) (+23.75 deg) and z3 = 775 - 0.36*(x-800) (-19.80 deg)
z2_line = 640.0 + 0.44*(X2 - 800.0);
z3_line = 775.0 - 0.36*(X2 - 800.0);
m_z2 = abs(Z2 - z2_line) <= 25.0 & X2 >= 300 & X2 <= 1300;
m_z3 = abs(Z2 - z3_line) <= 25.0 & X2 >= 300 & X2 <= 1300;
diagB = struct( ...
    'z2_true_dip_deg', rad2deg(atan(0.44)), ...
    'z2_blind_phi_noisy_med_deg', median(rad2deg(phi_s2(m_z2))), ...
    'z2_blind_phi_nf_med_deg',    median(rad2deg(phi_s2_nf(m_z2))), ...
    'z3_true_dip_deg', rad2deg(atan(-0.36)), ...
    'z3_blind_phi_noisy_med_deg', median(rad2deg(phi_s2(m_z3))), ...
    'z3_blind_phi_nf_med_deg',    median(rad2deg(phi_s2_nf(m_z3))), ...
    'W_deep_med_coh_noisy', median(c_s2(Z2>=560 & Z2<=950 & X2>=180 & X2<=1420)));
fprintf('  S2 Layer z2 (true=%+.1f deg): blind phi_eff noisy=%+.1f deg, NF=%+.1f deg\n', ...
    diagB.z2_true_dip_deg, diagB.z2_blind_phi_noisy_med_deg, diagB.z2_blind_phi_nf_med_deg);
fprintf('  S2 Layer z3 (true=%+.1f deg): blind phi_eff noisy=%+.1f deg, NF=%+.1f deg\n', ...
    diagB.z3_true_dip_deg, diagB.z3_blind_phi_noisy_med_deg, diagB.z3_blind_phi_nf_med_deg);
fprintf('  S2 W_deep median structure-tensor coherence (SNR=+1.06 dB): %.3f\n', diagB.W_deep_med_coh_noisy);
out.diagB_s2_dip_bias = diagB;

%% ========================================================================
%  DIAG C: Post-Hoc Tighter Gain-Cap Sweep (10 / 15 / 20 / 30 / 70 dB) for True WE-QRTM on S1
%  ========================================================================
fprintf('\n=== POST-HOC DIAG C: True WE-QRTM Gain-Cap Sensitivity Sweep (10/15/20/30/70 dB) on S1 ===\n');
% Fast direct check on S1 (m = 0.5, 1.0, 2.0) to test whether tighter gain-cap rescues deterministic WE-QRTM!
P.dx = 10.0; P.npml = 18; P.nzi = 100; P.nxi = 160;
P.nz = P.nzi + 2*P.npml; P.nx = P.nxi + 2*P.npml;
P.z = ((1:P.nz)' - P.npml - 1) * P.dx; P.x = ((1:P.nx) - P.npml - 1) * P.dx;
P.in = {P.npml + (1:P.nzi), P.npml + (1:P.nxi)}; P.izs = P.npml + 3;
P.ixs = P.npml + 1 + (8:3:152); P.ixr = P.npml + 1 + (0:1:P.nxi-1);
ns = numel(P.ixs); nr = numel(P.ixr); N = P.nz * P.nx;
[Zi, Xi] = ndgrid(P.in{1}, P.in{2}); lin = sub2ind([P.nz P.nx], Zi(:), Xi(:)); nin = numel(lin);
isrc = sub2ind([P.nz P.nx], P.izs*ones(1,ns), P.ixs);
irec = sub2ind([P.nz P.nx], P.izs*ones(1,nr), P.ixr);
[Zfull, Xfull] = ndgrid(P.z, P.x);
z_in = P.z(P.in{1}); x_in = P.x(P.in{2}); xs = P.x(P.ixs); xr = P.x(P.ixr);

z1 = 230.0 - 14.0 * exp(-((Xfull - 770.0)/330.0).^2);
z2 = 625.0 - 36.0 * exp(-((Xfull - 790.0)/300.0).^2) + 0.08*(Xfull - 800.0);
z3 = 760.0 - 26.0 * exp(-((Xfull - 790.0)/330.0).^2) + 0.06*(Xfull - 800.0);
z4 = 890.0 - 14.0 * exp(-((Xfull - 800.0)/370.0).^2);
v_true = 1850.0 * ones(P.nz, P.nx); v_true(Zfull>=z1)=2160; v_true(Zfull>=z2)=2580; v_true(Zfull>=z3)=2940; v_true(Zfull>=z4)=3320;
v_mig  = core.smooth2d(v_true, 2.5);
invQ_t = (1/90)*ones(P.nz, P.nx); invQ_t(Zfull < 80) = 1/1000;
invQ_t = invQ_t + (1/14 - 1/90)*exp(-(((Xfull - 790)/210).^2 + ((Zfull - 415)/140).^2).^1.5);
Q_true = 1 ./ invQ_t;
W_targ = {[580 1000], [560 950]}; W_flank = {[200 480], [560 950]};
iz_w = find(z_in >= 560 & z_in <= 950); ix_w = find(x_in >= 180 & x_in <= 1420);

% Exact same S1 seed & noise as WP2!
rng(20261501, 'twister');
corr_sig = 5.0; rc = ceil(3*corr_sig); g_nc = exp(-0.5*((-rc:rc)'/corr_sig).^2); g_nc = g_nc / sqrt(sum(g_nc.^2));
fref = 25.0; fpeak = 22.0; freqs = 4.0:1.0:37.0; nf = numel(freqs);
W22 = exp(-1.0); b22 = zeros(N,ns); b22(sub2ind([N ns],isrc,1:ns))=-W22/P.dx^2;
[Lm,Um,Pm,Qm,Rm] = lu(core.helm_op_visco(v_mig, v_mig, 22, P)); Ub = core.solve_helm_shot(Lm,Um,Pm,Qm,Rm, b22); Dbg22 = Ub(irec,:); clear Ub Lm Um Pm Qm Rm
vc22 = core.make_visco_vel(v_true, Q_true, 22, fref, P.npml, 'att', 1.0);
[Lv,Uv,Pv,Qv,Rv] = lu(core.helm_op_visco(vc22, v_true, 22, P)); Uv22 = core.solve_helm_shot(Lv,Uv,Pv,Qv,Rv, b22); Dcal = Uv22(irec,:)-Dbg22; clear Uv22 Lv Uv Pv Qv Rv
sig_n = sqrt(mean(abs(Dcal(:)).^2)) * 10^(-0.5/20);

D_vi_all = zeros(nr, ns, nf, 'single'); N_f_all = D_vi_all; I_ref = zeros(nin,1); I_r0 = I_ref;
for jf = 1:nf
    f = freqs(jf); w = 2*pi*f; W_src = (f/fpeak)^2*exp(-(f/fpeak)^2);
    b_src = zeros(N,ns); b_src(sub2ind([N ns],isrc,1:ns))=-W_src/P.dx^2;
    [Lm,Um,Pm,Qm,Rm] = lu(core.helm_op_visco(v_mig, v_mig, f, P)); U_bg = core.solve_helm_shot(Lm,Um,Pm,Qm,Rm, b_src); D_bg = U_bg(irec,:); Us = conj(U_bg(lin,:))*w^2; clear U_bg
    [Le,Ue,Pe,Qe,Re] = lu(core.helm_op_visco(v_true, v_true, f, P)); U_el = core.solve_helm_shot(Le,Ue,Pe,Qe,Re, b_src); D_el = U_el(irec,:)-D_bg; clear U_el Le Ue Pe Qe Re
    vc_t = core.make_visco_vel(v_true, Q_true, f, fref, P.npml, 'att', 1.0);
    [Lv,Uv,Pv,Qv,Rv] = lu(core.helm_op_visco(vc_t, v_true, f, P)); U_vi = core.solve_helm_shot(Lv,Uv,Pv,Qv,Rv, b_src); D_vi = U_vi(irec,:)-D_bg; clear U_vi Lv Uv Pv Qv Rv
    N_f = conv2(g_nc, 1, sig_n*(1/sqrt(2))*(randn(nr,ns)+1i*randn(nr,ns)), 'same');
    D_vi_all(:,:,jf) = single(D_vi); N_f_all(:,:,jf) = single(N_f);
    be = zeros(N,ns); be(irec,:)=D_el; Ve = core.solve_helm_adj(Lm,Um,Pm,Qm,Rm, be, P.dx);
    bt = zeros(N,ns); bt(irec,:)=D_vi+N_f; Vt = core.solve_helm_adj(Lm,Um,Pm,Qm,Rm, bt, P.dx); clear Lm Um Pm Qm Rm
    I_ref = I_ref + real(sum(Us .* Ve(lin,:), 2));
    I_r0  = I_r0  + real(sum(Us .* Vt(lin,:), 2));
end
% Rebuild w_bent
s_mig = 1.0 ./ v_mig; F_T = zeros(P.nz,P.nx,ns,'single'); F_ts = F_T;
for k = 1:ns, F_T(:,:,k)=single(eikonal_fs(s_mig,P.dx,P.izs,P.ixs(k))); end
R0_2d = core.prep_image(I_r0, P.nzi, P.nxi, z_in, P.dx, 50, 160);
[phi_e, ~] = core.blind_dip(R0_2d, P.dx, 40, 35);
[w_bent, ~, db_s1] = core.build_dsa_bent(v_mig, P, P.z(P.izs), P.z(P.izs), xs, repmat(min(xr),1,ns), repmat(max(xr),1,ns), phi_e, F_T, F_ts, cfg_s);
w_bnt_d = double(w_bent);
[xs_s, ord_s] = sort(xs(:)); q_r = interp1(xs_s, 1:ns, min(max(double(db_s1.xmir), xs_s(1)), xs_s(end)));
k0_r = min(floor(q_r), ns-1); ww_r = q_r - k0_r; row_v = repmat((1:nin)',1,ns);
idx_r0 = sub2ind([nin ns], row_v, ord_s(k0_r)); idx_r1 = sub2ind([nin ns], row_v, ord_s(k0_r+1));

ref_2d = core.prep_image(I_ref, P.nzi, P.nxi, z_in, P.dx, 50, 160);
rms_w = @(A, W) sqrt(mean(A(z_in>=W{2}(1)&z_in<=W{2}(2), x_in>=W{1}(1)&x_in<=W{1}(2)).^2, 'all'));
ref_rcf = log10(rms_w(ref_2d, W_targ) / max(rms_w(ref_2d, W_flank), 1e-30));

g_caps = [6, 10, 15, 20, 30, 70]; % dB gain caps for Wang/Sun-Zhu WE-QRTM!
s_caps = arrayfun(@(g) gain_s(g), g_caps);
ms_c = [0.5, 1.0, 2.0];
IT_raw = zeros(nin, numel(g_caps), numel(ms_c));
IT_dsa = zeros(nin, numel(g_caps), numel(ms_c));

for im = 1:numel(ms_c)
    Q_est = ms_c(im) * Q_true;
    for jf = 1:nf
        f = freqs(jf); w = 2*pi*f; W_src = (f/fpeak)^2*exp(-(f/fpeak)^2);
        b_src = zeros(N,ns); b_src(sub2ind([N ns],isrc,1:ns))=-W_src/P.dx^2;
        b_tot = zeros(N,ns); b_tot(irec,:) = double(D_vi_all(:,:,jf) + N_f_all(:,:,jf));

        vc_d = core.make_visco_vel(v_mig, Q_est, f, fref, P.npml, 'disp', 0.0);
        [Ld,Ud,Pd,Qd,Rd] = lu(core.helm_op_visco(vc_d, v_mig, f, P));
        Usd_in = core.solve_helm_shot(Ld,Ud,Pd,Qd,Rd, b_src); Usd_in = Usd_in(lin,:);
        Vd = core.solve_helm_adj(Ld,Ud,Pd,Qd,Rd, b_tot, P.dx); clear Ld Ud Pd Qd Rd
        cd_tot = real((conj(Usd_in)*w^2) .* Vd(lin,:)); clear Vd

        vc_a = core.make_visco_vel(v_mig, Q_est, f, fref, P.npml, 'att', 1.0);
        [La,Ua,Pa,Qa,Ra] = lu(core.helm_op_visco(vc_a, v_mig, f, P));
        Usa_in = core.solve_helm_shot(La,Ua,Pa,Qa,Ra, b_src); Usa_in = Usa_in(lin,:); clear La Ua Pa Qa Ra
        b_S = zeros(nin, ns);
        for ks = 1:ns
            pa = movmean(movmean(reshape(abs(Usa_in(:,ks)).^2, P.nzi, P.nxi), 3, 1), 3, 2);
            pd = movmean(movmean(reshape(abs(Usd_in(:,ks)).^2, P.nzi, P.nxi), 3, 1), 3, 2);
            b_S(:,ks) = min(max(sqrt(pa(:) ./ max(pd(:), 1e-30)), 1e-4), 1.0);
        end
        clear Usa_in Usd_in
        b_R = (1 - ww_r).*b_S(idx_r0) + ww_r.*b_S(idx_r1);
        b2 = (min(max(b_S .* b_R, 1e-5), 1.0)).^2;

        vc_c = core.make_visco_vel(v_mig, Q_est, f, fref, P.npml, 'comp', 1.0);
        [Lc,Uc,Pc,Qc,Rc] = lu(core.helm_op_visco(vc_c, v_mig, f, P));
        Usc = conj(core.solve_helm_shot(Lc,Uc,Pc,Qc,Rc, b_src)); Usc = Usc(lin,:)*w^2;
        Vc = core.solve_helm_adj(Lc,Uc,Pc,Qc,Rc, b_tot, P.dx); clear Lc Uc Pc Qc Rc
        cc_tot = real(Usc .* Vc(lin,:)); clear Usc Vc

        for ig = 1:numel(g_caps)
            ws_g = b2 ./ (b2 + s_caps(ig));
            c_g = ws_g .* cc_tot + (1 - ws_g) .* cd_tot;
            IT_raw(:, ig, im) = IT_raw(:, ig, im) + sum(c_g, 2);
            IT_dsa(:, ig, im) = IT_dsa(:, ig, im) + sum(w_bnt_d .* c_g, 2);
        end
    end
end

diagC = struct('g_caps_dB', g_caps, 'm', ms_c, 'WE_QRTM_Reg_RAW', struct(), 'DSA_WE_QRTM_Reg', struct());
for ig = 1:numel(g_caps)
    gk = sprintf('cap_%02ddB', g_caps(ig));
    e1_r = zeros(1,3); ea_r = e1_r; e1_d = e1_r; ea_d = e1_r;
    for im = 1:3
        Ar = core.prep_image(IT_raw(:,ig,im), P.nzi, P.nxi, z_in, P.dx, 50, 160);
        Ad = core.prep_image(IT_dsa(:,ig,im), P.nzi, P.nxi, z_in, P.dx, 50, 160);
        e1_r(im) = core.score_e1(Ar, ref_2d, iz_w, ix_w);
        ea_r(im) = abs(log10(rms_w(Ar, W_targ)/max(rms_w(Ar, W_flank),1e-30)) - ref_rcf);
        e1_d(im) = core.score_e1(Ad, ref_2d, iz_w, ix_w);
        ea_d(im) = abs(log10(rms_w(Ad, W_targ)/max(rms_w(Ad, W_flank),1e-30)) - ref_rcf);
    end
    diagC.WE_QRTM_Reg_RAW.(gk) = struct('E1', e1_r, 'E1_mean', mean(e1_r), 'E_amp', ea_r, 'E_amp_mean', mean(ea_r));
    diagC.DSA_WE_QRTM_Reg.(gk) = struct('E1', e1_d, 'E1_mean', mean(e1_d), 'E_amp', ea_d, 'E_amp_mean', mean(ea_d));
    fprintf('  GainCap=%2d dB (s=%.2e) | RAW WE-Reg: E1=%.4f (m0.5=%.3f m1=%.3f m2=%.3f) Eamp=%.3f | DSA WE-Reg: E1=%.4f (m0.5=%.3f m1=%.3f m2=%.3f) Eamp=%.3f\n', ...
        g_caps(ig), s_caps(ig), mean(e1_r), e1_r(1), e1_r(2), e1_r(3), mean(ea_r), mean(e1_d), e1_d(1), e1_d(2), e1_d(3), mean(ea_d));
end
out.diagC_we_qrtm_gaincap_sweep_S1 = diagC;

fn = fullfile(res_dir, 'posthoc_diagnostics.json');
fid = fopen(fn, 'w'); fwrite(fid, jsonencode(out, 'PrettyPrint', true), 'char'); fclose(fid);
fprintf('\nSaved %s\n', fn);
end
