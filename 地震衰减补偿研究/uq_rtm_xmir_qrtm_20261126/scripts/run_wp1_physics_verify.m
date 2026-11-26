function run_wp1_physics_verify()
% WP1: Pre-Registration Physics & Numerical Unit Verification Suite
%   Test 1A: Exact analytical circular-arc ray test in v(z) = v0 + g*z (xmir_str vs xmir_bent)
%   Test 1B: Exact noise-free Viscoacoustic Wave-Equation Q-RTM recovery & PML stability proof
%   Test 1C: Volve 10m real velocity field xmir_str vs xmir_bent geometric discrepancy audit

maxNumCompThreads(12);
here = fileparts(mfilename('fullpath'));
root = fileparts(here);
res_dir = fullfile(root, 'results');
if ~isfolder(res_dir), mkdir(res_dir); end

addpath('D:\ChatGPT\研究\地震衰减补偿研究\uq_rtm_20261117\scripts');
core = uq_rtm_xmir_qrtm_core();

out = struct();

%% ========================================================================
%  TEST 1A: Exact Analytical Circular-Arc Snell Ray Test in v(z) = 1500 + 0.8*z
%  ========================================================================
fprintf('=== Test 1A: Exact Circular-Arc Snell Ray Verification (v(z) = 1500 + 0.8*z) ===\n');
P.dx = 10.0; P.npml = 18; P.nzi = 200; P.nxi = 300;
P.nz = P.nzi + 2*P.npml; P.nx = P.nxi + 2*P.npml;
P.z = ((1:P.nz)' - P.npml - 1) * P.dx;
P.x = ((1:P.nx)  - P.npml - 1) * P.dx;
P.in = {P.npml + (1:P.nzi), P.npml + (1:P.nxi)};
v0 = 1500.0; g = 0.8;
[Zfull, ~] = ndgrid(P.z, P.x);
v_lin = v0 + g * max(Zfull, 0.0);
s_lin = 1.0 ./ v_lin;

z_src = 0.0; iz_src = P.npml + 1;
ix_src = P.npml + 101; xs_val = P.x(ix_src); % x_s = 1000 m
F_T = eikonal_fs(s_lin, P.dx, iz_src, ix_src);
F_ts = zeros(P.nz, P.nx, 1, 'single');

z_in = P.z(P.in{1}); x_in = P.x(P.in{2});
[Zin, Xin] = ndgrid(z_in, x_in);
[Zi_f, Xi_f] = ndgrid(P.in{1}, P.in{2});

% Evaluation window: z in [500, 1800] m, x - xs in [200, 1400] m (non-zero offset, deep enough)
eval_mask = (Zin >= 500 & Zin <= 1800) & (Xin - xs_val >= 200 & Xin - xs_val <= 1400);

% Exact analytical circular-arc incident angle th_s_exact at (x, z) from (xs_val, 0):
% Circle center lies on z_c = -v0/g. Passing through (xs, 0) and (x, z):
%   x_c = 0.5*(xs + x) + (z / (2*(x - xs))) * (z + 2*v0/g)
%   R   = hypot(xs - x_c, 0 - (-v0/g))
%   sin(th_s_exact) = (z + v0/g) / R,  cos(th_s_exact) = (x_c - x) / R
zc = -v0 / g;
dx_m = Xin - xs_val;
xc_arc = 0.5 * (xs_val + Xin) + (Zin ./ max(2 * dx_m, 1e-6)) .* (Zin - 2*zc);
R_arc  = hypot(xs_val - xc_arc, -zc);
th_s_exact = atan2(Zin - zc, xc_arc - Xin); % exact local downgoing angle from +z!
% Check: at z=0, (0 - zc)/(xc_arc - xs) = tan(th0) = (v0/g)/(xc-xs) -> sin(th0)=(v0/g)/R. Exact!

cfg_t = core.default_cfg('synth_10m');
dips_deg = [-20, -10, 0, 10, 20];
t1a = struct();
for id = 1:numel(dips_deg)
    phi_d = dips_deg(id);
    phi_2d = deg2rad(phi_d) * ones(P.nzi, P.nxi);

    % Exact analytical reflected emergence xmir_exact:
    th_r_exact = th_s_exact + 2 * phi_2d;
    v_loc = v0 + g * Zin;
    px_ex = sin(th_r_exact) ./ v_loc;
    % Upgoing circular arc from z to 0: xmir = x + (1/(g*px)) * (cos(th_surf) - cos(th_r_exact))
    sin_surf = px_ex * v0;
    valid_ex = eval_mask & (abs(th_r_exact) < deg2rad(80)) & (abs(sin_surf) < 0.98) & (abs(px_ex) > 1e-9);
    th_surf_ex = asin(min(max(sin_surf, -1), 1));
    xmir_exact = Xin + (cos(th_surf_ex) - cos(th_r_exact)) ./ (g * max(abs(px_ex), 1e-12)) .* sign(px_ex);
    zero_px = abs(px_ex) <= 1e-9;
    xmir_exact(zero_px) = Xin(zero_px);

    % Legacy straight-ray xmir_str
    [~, ~, xm_str, ths_str] = core.build_dsa_str(Xin(:), Zin(:), z_src, z_src, xs_val, ...
        -1e5, 1e5, phi_2d(:), F_ts, Zi_f(:), Xi_f(:), P.nz, P.nx, cfg_t);
    xm_str = reshape(double(xm_str), P.nzi, P.nxi);
    ths_str = reshape(double(ths_str), P.nzi, P.nxi);

    % Upgraded bent-ray xmir_bent
    [~, ~, d_bent] = core.build_dsa_bent(v_lin, P, z_src, z_src, xs_val, ...
        -1e5, 1e5, phi_2d, F_T, F_ts, cfg_t);
    xm_bent = reshape(double(d_bent.xmir), P.nzi, P.nxi);
    ths_bent = reshape(double(d_bent.th_s_deg), P.nzi, P.nxi);

    err_str  = abs(xm_str(valid_ex)  - xmir_exact(valid_ex));
    err_bent = abs(xm_bent(valid_ex) - xmir_exact(valid_ex));
    err_ths_str  = abs(ths_str(valid_ex)  - rad2deg(th_s_exact(valid_ex)));
    err_ths_bent = abs(ths_bent(valid_ex) - rad2deg(th_s_exact(valid_ex)));

    red_pct = (1 - median(err_bent) / max(median(err_str), 1e-9)) * 100;
    fprintf('  dip=%+3d deg (N=%5d): th_s med_err Str=%.2f vs Bent=%.2f deg | xmir med_err Str=%6.1f m (p90=%6.1f m) vs Bent=%5.1f m (p90=%5.1f m) | red=%.1f%%\n', ...
        phi_d, nnz(valid_ex), median(err_ths_str), median(err_ths_bent), ...
        median(err_str), prctile(err_str, 90), median(err_bent), prctile(err_bent, 90), red_pct);

    key = sprintf('dip_%s%d', ternary(phi_d<0,'m','p'), abs(phi_d));
    t1a.(key) = struct('dip_deg', phi_d, 'n_pts', nnz(valid_ex), ...
        'th_s_med_err_str_deg', median(err_ths_str), 'th_s_med_err_bent_deg', median(err_ths_bent), ...
        'xmir_med_err_str_m', median(err_str), 'xmir_p90_err_str_m', prctile(err_str, 90), ...
        'xmir_med_err_bent_m', median(err_bent), 'xmir_p90_err_bent_m', prctile(err_bent, 90), ...
        'error_reduction_pct', red_pct);
end
out.test1a_analytical_ray = t1a;

%% ========================================================================
%  TEST 1B: Exact Noise-Free Viscoacoustic Wave-Equation Q-RTM Verification
%  ========================================================================
fprintf('\n=== Test 1B: Noise-Free True Viscoacoustic Wave-Equation Q-RTM Physics Proof ===\n');
P2.dx = 10.0; P2.npml = 18; P2.nzi = 100; P2.nxi = 160;
P2.nz = P2.nzi + 2*P2.npml; P2.nx = P2.nxi + 2*P2.npml;
P2.z = ((1:P2.nz)' - P2.npml - 1) * P2.dx; P2.x = ((1:P2.nx) - P2.npml - 1) * P2.dx;
P2.in = {P2.npml + (1:P2.nzi), P2.npml + (1:P2.nxi)};
P2.izs = P2.npml + 3; P2.ixs = P2.npml + 1 + (8:3:152); P2.ixr = P2.npml + 1 + (0:1:P2.nxi-1);
ns = numel(P2.ixs); nr = numel(P2.ixr); N = P2.nz * P2.nx;
[Zi2, Xi2] = ndgrid(P2.in{1}, P2.in{2}); lin2 = sub2ind([P2.nz P2.nx], Zi2(:), Xi2(:)); nin2 = numel(lin2);
isrc2 = sub2ind([P2.nz P2.nx], P2.izs*ones(1,ns), P2.ixs);
irec2 = sub2ind([P2.nz P2.nx], P2.izs*ones(1,nr), P2.ixr);
[Zf2, Xf2] = ndgrid(P2.z, P2.x);
z_in2 = P2.z(P2.in{1}); x_in2 = P2.x(P2.in{2});

% Heterogeneous gas-chimney model (same geometry as U6)
z1 = 230.0 - 14.0 * exp(-((Xf2 - 770.0)/330.0).^2);
z2 = 625.0 - 36.0 * exp(-((Xf2 - 790.0)/300.0).^2) + 0.08*(Xf2 - 800.0);
z3 = 760.0 - 26.0 * exp(-((Xf2 - 790.0)/330.0).^2) + 0.06*(Xf2 - 800.0);
z4 = 890.0 - 14.0 * exp(-((Xf2 - 800.0)/370.0).^2);
v_true = 1850.0 * ones(P2.nz, P2.nx); v_true(Zf2>=z1)=2160; v_true(Zf2>=z2)=2580; v_true(Zf2>=z3)=2940; v_true(Zf2>=z4)=3320;
v_mig  = core.smooth2d(v_true, 2.5);
invQ_t = (1/90)*ones(P2.nz, P2.nx); invQ_t(Zf2 < 80) = 1/1000;
invQ_t = invQ_t + (1/14 - 1/90)*exp(-(((Xf2 - 790)/210).^2 + ((Zf2 - 415)/140).^2).^1.5);
Q_true = 1 ./ invQ_t;

freqs = 6.0:2.0:36.0; nf = numel(freqs); fref = 25.0; fpeak = 22.0;
I_ref = zeros(nin2, 1); I_r0 = zeros(nin2, 1);
I_adj = zeros(nin2, 1); I_disp = zeros(nin2, 1); I_comp = zeros(nin2, 1);
max_comp_ratio = zeros(1, nf);
for jf = 1:nf
    f = freqs(jf); w = 2*pi*f;
    W_src = (f/fpeak)^2 * exp(-(f/fpeak)^2);
    b_src = zeros(N, ns); b_src(sub2ind([N ns], isrc2, 1:ns)) = -W_src / P2.dx^2;

    % Background & elastic data
    [Lm,Um,Pm,Qm,Rm] = lu(core.helm_op_visco(v_mig, v_mig, f, P2));
    U_bg = core.solve_helm_shot(Lm,Um,Pm,Qm,Rm, b_src); D_bg = U_bg(irec2,:); Us_m = conj(U_bg(lin2,:))*w^2; clear U_bg
    [Le,Ue,Pe,Qe,Re] = lu(core.helm_op_visco(v_true, v_true, f, P2));
    U_el = core.solve_helm_shot(Le,Ue,Pe,Qe,Re, b_src); D_el = U_el(irec2,:) - D_bg; clear U_el Le Ue Pe Qe Re

    % Attenuated data (true viscoacoustic forward)
    vc_att_true = core.make_visco_vel(v_true, Q_true, f, fref, P2.npml, 'att', 1.0);
    [Lv,Uv,Pv,Qv,Rv] = lu(core.helm_op_visco(vc_att_true, v_true, f, P2));
    U_vi = core.solve_helm_shot(Lv,Uv,Pv,Qv,Rv, b_src); D_vi = U_vi(irec2,:) - D_bg; clear U_vi Lv Uv Pv Qv Rv

    b_el = zeros(N, ns); b_el(irec2,:) = D_el; V_el = core.solve_helm_adj(Lm,Um,Pm,Qm,Rm, b_el, P2.dx);
    b_vi = zeros(N, ns); b_vi(irec2,:) = D_vi; V_r0 = core.solve_helm_adj(Lm,Um,Pm,Qm,Rm, b_vi, P2.dx);
    clear Lm Um Pm Qm Rm
    I_ref = I_ref + real(sum(Us_m .* V_el(lin2,:), 2));
    I_r0  = I_r0  + real(sum(Us_m .* V_r0(lin2,:), 2)); clear V_el V_r0

    % True Viscoacoustic Wave-Equation Adjoint RTM (exact Q_true)
    vc_att = core.make_visco_vel(v_mig, Q_true, f, fref, P2.npml, 'att', 1.0);
    [La,Ua,Pa,Qa,Ra] = lu(core.helm_op_visco(vc_att, v_mig, f, P2));
    Us_a = conj(core.solve_helm_shot(La,Ua,Pa,Qa,Ra, b_src)); Us_a = Us_a(lin2,:)*w^2;
    V_a  = core.solve_helm_adj(La,Ua,Pa,Qa,Ra, b_vi, P2.dx); clear La Ua Pa Qa Ra
    I_adj = I_adj + real(sum(Us_a .* V_a(lin2,:), 2)); clear Us_a V_a

    % True Viscoacoustic Dispersion-Only Wave-Equation RTM (exact Q_true)
    vc_disp = core.make_visco_vel(v_mig, Q_true, f, fref, P2.npml, 'disp', 0.0);
    [Ld,Ud,Pd,Qd,Rd] = lu(core.helm_op_visco(vc_disp, v_mig, f, P2));
    Us_d = conj(core.solve_helm_shot(Ld,Ud,Pd,Qd,Rd, b_src)); Us_d = Us_d(lin2,:)*w^2;
    V_d  = core.solve_helm_adj(Ld,Ud,Pd,Qd,Rd, b_vi, P2.dx); clear Ld Ud Pd Qd Rd
    I_disp = I_disp + real(sum(Us_d .* V_d(lin2,:), 2)); clear Us_d V_d

    % True Viscoacoustic Compensated Q-RTM (Unregularized, exact Q_true)
    vc_comp = core.make_visco_vel(v_mig, Q_true, f, fref, P2.npml, 'comp', 1.0);
    [Lc,Uc,Pc,Qc,Rc] = lu(core.helm_op_visco(vc_comp, v_mig, f, P2));
    Us_c = core.solve_helm_shot(Lc,Uc,Pc,Qc,Rc, b_src);
    max_comp_ratio(jf) = max(abs(Us_c(lin2,:)), [], 'all') / max(abs(Us_m(:)/w^2), [], 'all');
    Us_c = conj(Us_c(lin2,:))*w^2;
    V_c  = core.solve_helm_adj(Lc,Uc,Pc,Qc,Rc, b_vi, P2.dx); clear Lc Uc Pc Qc Rc
    I_comp = I_comp + real(sum(Us_c .* V_c(lin2,:), 2)); clear Us_c V_c
end

iz_w = find(z_in2>=560 & z_in2<=950); ix_w = find(x_in2>=180 & x_in2<=1420);
W_targ = {[580 1000], [560 950]}; W_flank = {[200 480], [560 950]};
rms_w = @(A, W) sqrt(mean(A(z_in2>=W{2}(1)&z_in2<=W{2}(2), x_in2>=W{1}(1)&x_in2<=W{1}(2)).^2, 'all'));

R_ref  = core.prep_image(I_ref,  P2.nzi, P2.nxi, z_in2, P2.dx, 50, 160);
R_r0   = core.prep_image(I_r0,   P2.nzi, P2.nxi, z_in2, P2.dx, 50, 160);
R_adj  = core.prep_image(I_adj,  P2.nzi, P2.nxi, z_in2, P2.dx, 50, 160);
R_disp = core.prep_image(I_disp, P2.nzi, P2.nxi, z_in2, P2.dx, 50, 160);
R_comp = core.prep_image(I_comp, P2.nzi, P2.nxi, z_in2, P2.dx, 50, 160);
ref_rcf = log10(rms_w(R_ref, W_targ) / max(rms_w(R_ref, W_flank), 1e-30));

eval_nf = @(A) struct('E1', core.score_e1(A, R_ref, iz_w, ix_w), ...
    'E_amp', abs(log10(rms_w(A, W_targ)/max(rms_w(A, W_flank),1e-30)) - ref_rcf), ...
    'rel_L2_scaled', norm(A(iz_w,ix_w)/norm(A(iz_w,ix_w),'fro') - R_ref(iz_w,ix_w)/norm(R_ref(iz_w,ix_w),'fro'), 'fro'));

t1b = struct('RAW_R0', eval_nf(R_r0), 'WE_Q_Adj', eval_nf(R_adj), ...
             'WE_Q_DispOnly', eval_nf(R_disp), 'WE_QRTM_Unreg', eval_nf(R_comp), ...
             'max_comp_wave_ratio_f36', max_comp_ratio(end));
fprintf('  [Noise-Free m=1] RAW_R0        : E1=%.4f | E_amp=%.4f | relL2=%.4f\n', t1b.RAW_R0.E1, t1b.RAW_R0.E_amp, t1b.RAW_R0.rel_L2_scaled);
fprintf('  [Noise-Free m=1] WE_Q_Adj      : E1=%.4f | E_amp=%.4f | relL2=%.4f\n', t1b.WE_Q_Adj.E1, t1b.WE_Q_Adj.E_amp, t1b.WE_Q_Adj.rel_L2_scaled);
fprintf('  [Noise-Free m=1] WE_Q_DispOnly : E1=%.4f | E_amp=%.4f | relL2=%.4f\n', t1b.WE_Q_DispOnly.E1, t1b.WE_Q_DispOnly.E_amp, t1b.WE_Q_DispOnly.rel_L2_scaled);
fprintf('  [Noise-Free m=1] WE_QRTM_Unreg : E1=%.4f | E_amp=%.4f | relL2=%.4f | max_wave_ratio@36Hz=%.2f\n', ...
    t1b.WE_QRTM_Unreg.E1, t1b.WE_QRTM_Unreg.E_amp, t1b.WE_QRTM_Unreg.rel_L2_scaled, t1b.max_comp_wave_ratio_f36);

% Also test stability of A_comp at extreme m=0.5 (Q_min = 7) at f=36 Hz
vc_m05 = core.make_visco_vel(v_mig, 0.5*Q_true, 36.0, fref, P2.npml, 'comp', 1.0);
[Lc,Uc,Pc,Qc,Rc] = lu(core.helm_op_visco(vc_m05, v_mig, 36.0, P2));
U_m05 = core.solve_helm_shot(Lc,Uc,Pc,Qc,Rc, b_src);
t1b.m05_f36_finite = all(isfinite(U_m05(:)));
t1b.m05_f36_max_amp = max(abs(U_m05(lin2,:)), [], 'all');
fprintf('  [Stability m=0.5 Q_min=7 @ 36 Hz] finite=%d | max|U_comp|=%.3e\n', t1b.m05_f36_finite, t1b.m05_f36_max_amp);
out.test1b_we_qrtm_physics = t1b;

%% ========================================================================
%  TEST 1C: Volve 10m Real Velocity Field xmir_str vs xmir_bent Audit
%  ========================================================================
fprintf('\n=== Test 1C: Volve 10m Real Velocity Field xmir_str vs xmir_bent Audit ===\n');
Pf = fld_setup();
cfg_f = core.default_cfg('field_obc');
z_f = Pf.z(Pf.in{1}); x_f = Pf.x(Pf.in{2});
[Zfi, Xfi] = ndgrid(z_f, x_f);
[Zif, Xif] = ndgrid(Pf.in{1}, Pf.in{2});

% Compute eikonal traveltimes F_T on Volve 10m (fast: 124 shots)
t_eik = tic;
s_f = 1.0 ./ Pf.vmig;
F_Tf = zeros(Pf.nz, Pf.nx, Pf.nshot, 'single');
for k = 1:Pf.nshot
    F_Tf(:, :, k) = single(eikonal_fs(s_f, Pf.dx, Pf.src_iz(k), Pf.src_ix(k)));
end
fprintf('  Volve 10m 124-shot eikonal F_T computed in %.1fs\n', toc(t_eik));

% Load Volve cached R0 & tfields.mat (read-only)
S1 = load(fullfile(Pf.cache, 'state1.mat'), 'I');
R0_f = core.prep_image(double(sum(S1.I, 2)), Pf.nzi, Pf.nxi, z_f, Pf.dx, cfg_f.bg_sigma_m, cfg_f.z_mute_m);
clear S1
[phi_f, ~] = core.blind_dip(R0_f, Pf.dx, cfg_f.dip_sigma_m, cfg_f.dip_cap_deg);
TF = load(fullfile(Pf.cache, 'tfields.mat'), 'F');
rx = Pf.x(Pf.rec_ix);
xr_lo = arrayfun(@(s) min(rx(Pf.shot_of == s)), 1:Pf.nshot);
xr_hi = arrayfun(@(s) max(rx(Pf.shot_of == s)), 1:Pf.nshot);

t_str = tic;
[w_str_f, TS_str_f, xm_str_f, ths_str_f] = core.build_dsa_str( ...
    Xfi(:), Zfi(:), Pf.z(Pf.src_iz), 90.0, Pf.src_x, xr_lo, xr_hi, phi_f(:), ...
    TF.F, Zif(:), Xif(:), Pf.nz, Pf.nx, cfg_f);
dt_str = toc(t_str);

t_bnt = tic;
[w_bnt_f, TS_bnt_f, db_f] = core.build_dsa_bent( ...
    Pf.vmig, Pf, Pf.z(Pf.src_iz), 90.0, Pf.src_x, xr_lo, xr_hi, phi_f, ...
    F_Tf, TF.F, cfg_f);
dt_bnt = toc(t_bnt);

res_m = (Zfi(:) >= 2600 & Zfi(:) <= 3200) & (Xfi(:) >= -2500 & Xfi(:) <= 1500);
res_mat = repmat(res_m, 1, Pf.nshot);
act_either = res_mat & ((w_str_f > 0.1) | (w_bnt_f > 0.1)) & (~db_f.turned);
dx_mir = abs(double(db_f.xmir(act_either)) - double(xm_str_f(act_either)));
dts_ms = abs(double(TS_bnt_f(act_either)) - double(TS_str_f(act_either))) * 1e3;
flip_rate = mean((w_str_f(res_mat) > 0.5) ~= (w_bnt_f(res_mat) > 0.5)) * 100;
fold_str  = median(sum(double(w_str_f(res_m, :)), 2));
fold_bnt  = median(sum(double(w_bnt_f(res_m, :)), 2));

% Also check dip > 5 deg subset in reservoir
steep_m = res_m & (abs(rad2deg(phi_f(:))) >= 5.0);
steep_act = repmat(steep_m, 1, Pf.nshot) & ((w_str_f > 0.1) | (w_bnt_f > 0.1)) & (~db_f.turned);
dx_mir_steep = abs(double(db_f.xmir(steep_act)) - double(xm_str_f(steep_act)));

t1c = struct('dt_str_s', dt_str, 'dt_bent_s', dt_bnt, ...
    'vmig_seabed_m_s', median(Pf.vmig(Pf.npml+11, Pf.in{2})), ...
    'vmig_res_m_s', median(Pf.vmig(Pf.in{1}(z_f>=2600 & z_f<=3200), Pf.in{2}), 'all'), ...
    'res_med_dxmir_m', median(dx_mir), 'res_p90_dxmir_m', prctile(dx_mir, 90), ...
    'res_steep5_med_dxmir_m', median(dx_mir_steep), 'res_steep5_p90_dxmir_m', prctile(dx_mir_steep, 90), ...
    'res_med_dts_ms', median(dts_ms), 'res_p90_dts_ms', prctile(dts_ms, 90), ...
    'res_gate_flip_pct', flip_rate, 'res_med_fold_str', fold_str, 'res_med_fold_bent', fold_bnt, ...
    'res_med_ths_str_deg', median(abs(double(ths_str_f(act_either)))), ...
    'res_med_ths_eik_deg', median(abs(double(db_f.th_s_deg(act_either)))), ...
    'res_med_thr_surf_deg', median(abs(double(db_f.th_r_surf_deg(act_either)))));
fprintf('  Volve timing: build_dsa_str=%.1fs | build_dsa_bent=%.1fs\n', dt_str, dt_bnt);
fprintf('  Volve vmig: seabed=%.0f m/s -> reservoir=%.0f m/s\n', t1c.vmig_seabed_m_s, t1c.vmig_res_m_s);
fprintf('  Volve Res |xmir_bent - xmir_str|: med=%.1f m, p90=%.1f m (dip>=5 deg: med=%.1f m, p90=%.1f m)\n', ...
    t1c.res_med_dxmir_m, t1c.res_p90_dxmir_m, t1c.res_steep5_med_dxmir_m, t1c.res_steep5_p90_dxmir_m);
fprintf('  Volve Res gate flip=%.1f%% | fold Str=%.1f vs Bent=%.1f | |dTS| p90=%.2f ms\n', ...
    t1c.res_gate_flip_pct, t1c.res_med_fold_str, t1c.res_med_fold_bent, t1c.res_p90_dts_ms);
fprintf('  Volve Res med angles: th_s_str=%.1f deg | th_s_eik=%.1f deg | th_r_surf=%.1f deg\n', ...
    t1c.res_med_ths_str_deg, t1c.res_med_ths_eik_deg, t1c.res_med_thr_surf_deg);
out.test1c_volve_xmir_audit = t1c;

fn = fullfile(res_dir, 'wp1_physics_verify.json');
fid = fopen(fn, 'w'); fwrite(fid, jsonencode(out, 'PrettyPrint', true), 'char'); fclose(fid);
fprintf('\nSaved %s\n', fn);
end

function s = ternary(cond, a, b)
if cond, s = a; else, s = b; end
end
