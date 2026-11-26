function run_closure_synth_u6_u9(smoke_snr)
% WP4: Pre-Registered New Synthetic Suite U6-U9 (Seeds 20261401-20261404)
%   U6: Gas chimney + gentle dip + l30 correlated noise (Theorem 3 + E_amp + E1 + SA vs DSA)
%   U7: 22 deg steep reverse fault + l30 correlated noise (SA vs DSA physical tradeoff)
%   U8: Non-scalar spatial Q error #1: Chimney laterally shifted +180m & width x1.4, plus m in {0.6, 1.0, 1.4}
%   U9: Non-scalar spatial Q error #2: Shallow/deep opposite Q error (z<480m: Q_mig=1.5*Q_true; z>=480m: Q_mig=0.6*Q_true) + l30
%
% Usage:
%   run_closure_synth_u6_u9(true)   % 3-freq pre-freeze SNR calibration smoke check (MUST verify RAW_R0 SNR_img in [+2, +6] dB)
%   run_closure_synth_u6_u9(false)  % Full 34-freq production run (4:1:37 Hz)

if nargin < 1, smoke_snr = false; end
maxNumCompThreads(8);
here = fileparts(mfilename('fullpath'));
root = fileparts(here);
res_dir = fullfile(root, 'results');
if ~isfolder(res_dir), mkdir(res_dir); end

wt_root = fileparts(root);
addpath(fullfile(wt_root, 'uq_rtm_20261117', 'scripts'));
core = uq_rtm_unified_core();
cfg  = core.default_cfg('synth_10m');

% Pre-calibrated snr_db per model (verified via multi-freq smoke check so RAW_R0 W_deep SNR_img in [+2.0, +6.0] dB!)
snr_cfg = struct('U6', 0.5, 'U7', -3.0, 'U8', 0.5, 'U9', 0.5);

models = {'U6', 'U7', 'U8', 'U9'};
if smoke_snr
    fprintf('=== WP4 3-Freq Pre-Freeze Image SNR Calibration Smoke Check ===\n');
    all_ok = true;
    for im = 1:numel(models)
        mname = models{im};
        snr_img = run_one_model(mname, snr_cfg.(mname), core, cfg, true);
        in_rng = (snr_img >= 2.0 && snr_img <= 6.0);
        if ~in_rng, all_ok = false; end
        fprintf('  [%s] snr_db=%+.1f dB -> RAW_R0 W_deep SNR_img = %+.2f dB [%s]\n', ...
            mname, snr_cfg.(mname), snr_img, ternary(in_rng, 'PASS [+2,+6]', 'OUT OF RANGE'));
    end
    fprintf('=== Smoke SNR Check Overall: %s ===\n', ternary(all_ok, 'ALL 4 PASS', 'NEEDS RECALIBRATION'));
    return;
end

for im = 1:numel(models)
    mname = models{im};
    fn = fullfile(res_dir, ['synth_' mname '.json']);
    if isfile(fn)
        fprintf('Skipping %s (already exists: %s)\n', mname, fn);
        continue;
    end
    res = run_one_model(mname, snr_cfg.(mname), core, cfg, false);
    fid = fopen(fn, 'w'); fwrite(fid, jsonencode(res, 'PrettyPrint', true), 'char'); fclose(fid);
    fprintf('Saved %s (%.1fs)\n\n', fn, res.elapsed_s);
end
end

function s = ternary(cond, a, b)
if cond, s = a; else, s = b; end
end

function out = run_one_model(mname, snr_db, core, cfg, smoke_snr)
t0_all = tic;
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
z_in = P.z(P.in{1}); x_in = P.x(P.in{2}); [Zin, Xin] = ndgrid(z_in, x_in);
xs = P.x(P.ixs); xr = P.x(P.ixr);

rec_mask = true(nr, ns);
corr_noise_sig = 5.0; % l30 laterally correlated noise (50 m)
ms = [0.5, 0.7, 1.0, 1.4, 2.0];

switch mname
    case 'U6'
        seed = 20261401;
        z1 = 230.0 - 14.0 * exp(-((Xfull - 770.0)/330.0).^2);
        z2 = 625.0 - 36.0 * exp(-((Xfull - 790.0)/300.0).^2) + 0.08*(Xfull - 800.0);
        z3 = 760.0 - 26.0 * exp(-((Xfull - 790.0)/330.0).^2) + 0.06*(Xfull - 800.0);
        z4 = 890.0 - 14.0 * exp(-((Xfull - 800.0)/370.0).^2);
        v_true = 1850.0 * ones(P.nz, P.nx); v_true(Zfull>=z1)=2160; v_true(Zfull>=z2)=2580; v_true(Zfull>=z3)=2940; v_true(Zfull>=z4)=3320;
        invQ_t = (1/90)*ones(P.nz, P.nx); invQ_t(Zfull < 80) = 1/1000;
        invQ_t = invQ_t + (1/14 - 1/90)*exp(-(((Xfull - 790)/210).^2 + ((Zfull - 415)/140).^2).^1.5);
        Q_true = 1./invQ_t;
        Q_mig_base = Q_true;
        W_targ = {[580 1000], [560 950]}; W_flank = {[200 480], [560 950]};
        q_err_desc = 'scalar_m';
    case 'U7'
        seed = 20261402;
        Zg = Zfull;
        hw = (Xfull < 820.0 - (Zfull - 200.0) / tand(63.0)) & (Zfull > 200.0); Zg(hw) = Zg(hw) + 32.0;
        z1 = 235.0 - 0.05*(Xfull-800); z2 = 655.0 + 0.40*(Xfull-800); z3 = 785.0 - 0.29*(Xfull-800); z4 = 900.0 - 0.10*(Xfull-800);
        v_true = 1860.0 * ones(P.nz, P.nx); v_true(Zg>=z1)=2170; v_true(Zg>=z2)=2630; v_true(Zg>=z3)=2990; v_true(Zg>=z4)=3360;
        invQ_t = (1/85)*ones(P.nz, P.nx); invQ_t(Zfull < 80) = 1/1000;
        invQ_t = invQ_t + (1/15 - 1/85)*exp(-(((Xfull - 960)/250).^2 + ((Zfull - 365)/100).^2).^1.5);
        Q_true = 1./invQ_t;
        Q_mig_base = Q_true;
        W_targ = {[780 1200], [560 950]}; W_flank = {[220 560], [560 950]};
        q_err_desc = 'scalar_m';
    case 'U8'
        seed = 20261403;
        ms = [0.6, 1.0, 1.4]; % Per Task Spec WP4 U8
        z1 = 230.0 - 14.0 * exp(-((Xfull - 770.0)/330.0).^2);
        z2 = 625.0 - 36.0 * exp(-((Xfull - 790.0)/300.0).^2) + 0.08*(Xfull - 800.0);
        z3 = 760.0 - 26.0 * exp(-((Xfull - 790.0)/330.0).^2) + 0.06*(Xfull - 800.0);
        z4 = 890.0 - 14.0 * exp(-((Xfull - 800.0)/370.0).^2);
        v_true = 1850.0 * ones(P.nz, P.nx); v_true(Zfull>=z1)=2160; v_true(Zfull>=z2)=2580; v_true(Zfull>=z3)=2940; v_true(Zfull>=z4)=3320;
        % True chimney: x=790m, rx=210m
        invQ_t = (1/90)*ones(P.nz, P.nx); invQ_t(Zfull < 80) = 1/1000;
        invQ_t = invQ_t + (1/14 - 1/90)*exp(-(((Xfull - 790)/210).^2 + ((Zfull - 415)/140).^2).^1.5);
        Q_true = 1./invQ_t;
        % Non-scalar spatial Q error #1: laterally shifted +180m (x=970m) + width x1.4 (rx=294m)
        invQ_m = (1/90)*ones(P.nz, P.nx); invQ_m(Zfull < 80) = 1/1000;
        invQ_m = invQ_m + (1/14 - 1/90)*exp(-(((Xfull - 970)/(210*1.4)).^2 + ((Zfull - 415)/140).^2).^1.5);
        Q_mig_base = 1./invQ_m;
        W_targ = {[580 1000], [560 950]}; W_flank = {[200 480], [560 950]};
        q_err_desc = 'nonscalar_shifted_chimney_x970_rx294_plus_m';
    case 'U9'
        seed = 20261404;
        z1 = 230.0 - 14.0 * exp(-((Xfull - 770.0)/330.0).^2);
        z2 = 625.0 - 36.0 * exp(-((Xfull - 790.0)/300.0).^2) + 0.08*(Xfull - 800.0);
        z3 = 760.0 - 26.0 * exp(-((Xfull - 790.0)/330.0).^2) + 0.06*(Xfull - 800.0);
        z4 = 890.0 - 14.0 * exp(-((Xfull - 800.0)/370.0).^2);
        v_true = 1850.0 * ones(P.nz, P.nx); v_true(Zfull>=z1)=2160; v_true(Zfull>=z2)=2580; v_true(Zfull>=z3)=2940; v_true(Zfull>=z4)=3320;
        invQ_t = (1/90)*ones(P.nz, P.nx); invQ_t(Zfull < 80) = 1/1000;
        invQ_t = invQ_t + (1/14 - 1/90)*exp(-(((Xfull - 790)/210).^2 + ((Zfull - 415)/140).^2).^1.5);
        Q_true = 1./invQ_t;
        % Non-scalar spatial Q error #2: shallow (z<480m) 1.5*Q_true (undercomp), deep (z>=480m) 0.6*Q_true (overcomp)
        Q_mig_base = Q_true;
        Q_mig_base(Zfull <  480.0 & Zfull >= 80.0) = 1.5 * Q_true(Zfull <  480.0 & Zfull >= 80.0);
        Q_mig_base(Zfull >= 480.0)                 = 0.6 * Q_true(Zfull >= 480.0);
        W_targ = {[580 1000], [560 950]}; W_flank = {[200 480], [560 950]};
        q_err_desc = 'nonscalar_opposite_shallow1.5_deep0.6';
end

nm = numel(ms);
rng(seed, 'twister');
v_mig = core.smooth2d(v_true, 2.5);

if smoke_snr
    freqs = 6.0:3.0:33.0; % 10-freq representative bandwidth check (~40s total)
else
    freqs = 4.0:1.0:37.0;
end
nf = numel(freqs); fref_phys = 25.0; fpeak = 22.0; sig = cfg.sig; s70 = gain_s(70);

% Calibrate noise amplitude at 22 Hz
W22 = (22/fpeak)^2 * exp(-(22/fpeak)^2);
b22 = zeros(N, ns); b22(sub2ind([N ns], isrc, 1:ns)) = -W22 / P.dx^2;
[Lm,Um,Pm,Qm,Rm] = lu(helm_op(v_mig, 22, P)); Ub = Qm*(Um\(Lm\(Pm*(Rm\b22)))); Dbg22 = Ub(irec,:); clear Ub Lm Um Pm Qm Rm
vc22 = v_true ./ (1 - (1./(pi*Q_true))*log(22/fref_phys) - 1i./(2*Q_true));
[Lv,Uv,Pv,Qv,Rv] = lu(helm_op(vc22, 22, P)); Uv22 = Qv*(Uv\(Lv\(Pv*(Rv\b22)))); Dcal = (Uv22(irec,:) - Dbg22) .* rec_mask; clear Uv22 Lv Uv Pv Qv Rv
sig_n = sqrt(mean(abs(Dcal(rec_mask)).^2)) * 10^(-snr_db/20);
rc = ceil(3 * corr_noise_sig); g_nc = exp(-0.5 * ((-rc:rc)' / corr_noise_sig).^2); g_nc = g_nc / sqrt(sum(g_nc.^2));

iz_w = find(z_in>=560 & z_in<=950); ix_w = find(x_in>=180 & x_in<=1420);

if smoke_snr
    % Quick 3-freq RAW_R0 image SNR check in W_deep
    I_s = zeros(nin, 1); I_n = zeros(nin, 1);
    for jf = 1:nf
        f = freqs(jf); w = 2*pi*f;
        W_src = (f/fpeak)^2 * exp(-(f/fpeak)^2);
        b_src = zeros(N, ns); b_src(sub2ind([N ns], isrc, 1:ns)) = -W_src / P.dx^2;
        [Lm,Um,Pm,Qm,Rm] = lu(helm_op(v_mig, f, P));
        U_bg = Qm*(Um\(Lm\(Pm*(Rm\b_src)))); D_bg = U_bg(irec,:); Us_c = conj(U_bg(lin,:)) * w^2; clear U_bg
        vc = v_true ./ (1 - (1./(pi*Q_true))*log(max(f,1e-6)/fref_phys) - 1i./(2*Q_true));
        [Lv,Uv,Pv,Qv,Rv] = lu(helm_op(vc, f, P));
        U_vi = Qv*(Uv\(Lv\(Pv*(Rv\b_src)))); D_vi = (U_vi(irec,:) - D_bg) .* rec_mask; clear U_vi Lv Uv Pv Qv Rv b_src
        N_f = conv2(g_nc, 1, sig_n * (1/sqrt(2)) * (randn(nr, ns) + 1i*randn(nr, ns)), 'same') .* rec_mask;
        V_vi = extract_bp(D_vi, Lm, Um, Pm, Qm, Rm, irec, lin, N, ns, P.dx);
        V_nu = extract_bp(N_f,  Lm, Um, Pm, Qm, Rm, irec, lin, N, ns, P.dx);
        clear Lm Um Pm Qm Rm
        I_s = I_s + real(sum(Us_c .* V_vi, 2));
        I_n = I_n + real(sum(Us_c .* V_nu, 2));
    end
    S2 = core.prep_image(I_s, P.nzi, P.nxi, z_in, P.dx, cfg.bg_sigma_m, cfg.z_mute_m);
    N2 = core.prep_image(I_n, P.nzi, P.nxi, z_in, P.dx, cfg.bg_sigma_m, cfg.z_mute_m);
    out = 10*log10(sum(S2(iz_w,ix_w).^2,'all') / max(sum(N2(iz_w,ix_w).^2,'all'), 1e-30));
    return;
end

% Eikonal t* computed using Q_mig_base (captures non-scalar spatial Q error in U8 & U9!)
s_mig = 1.0 ./ v_mig; eps_q = 0.5;
F.T = zeros(P.nz, P.nx, ns, 'single'); F.ts = F.T;
for k = 1:ns
    F.T(:,:,k)  = eikonal_fs(s_mig, P.dx, P.izs, P.ixs(k));
    F.ts(:,:,k) = (eikonal_fs(s_mig.*(1+eps_q./Q_mig_base), P.dx, P.izs, P.ixs(k)) - ...
                   eikonal_fs(s_mig.*(1-eps_q./Q_mig_base), P.dx, P.izs, P.ixs(k))) / (2*eps_q);
end
TS_flat = zeros(nin, ns); Fts_in = zeros(nin, ns);
for k = 1:ns
    [~, tt] = two_way(F, k, P); TS_flat(:, k) = double(tt(lin));
    tk = F.ts(:,:,k); Fts_in(:, k) = double(tk(lin));
end
clear F
xr_lo = arrayfun(@(k) min(xr(rec_mask(:,k))), 1:ns);
xr_hi = arrayfun(@(k) max(xr(rec_mask(:,k))), 1:ns);
w_deep_mask = (Zin(:) >= 560 & Zin(:) <= 950) & (Xin(:) >= 180 & Xin(:) <= 1420);

rng(seed + 1, 'twister');   EPS_B1   = 2*(rand(nr, ns) > 0.5) - 1;
rng(seed + 300, 'twister'); sgn_B300 = 2*(rand(1, ns)  > 0.5) - 1;

Cs_sig_all = zeros(nin, ns, nf, 'single'); Cs_nu_all = Cs_sig_all; Cs_b1_all = Cs_sig_all;
I_ref = zeros(nin, 1); I_R0_tot = zeros(nin, 1);
for jf = 1:nf
    f = freqs(jf); w = 2*pi*f;
    W_src = (f/fpeak)^2 * exp(-(f/fpeak)^2);
    b_src = zeros(N, ns); b_src(sub2ind([N ns], isrc, 1:ns)) = -W_src / P.dx^2;
    [Lm,Um,Pm,Qm,Rm] = lu(helm_op(v_mig, f, P));
    U_bg = Qm*(Um\(Lm\(Pm*(Rm\b_src)))); D_bg = U_bg(irec,:); Us_c = conj(U_bg(lin,:)) * w^2; clear U_bg
    [Le,Ue,Pe,Qe,Re] = lu(helm_op(v_true, f, P));
    U_el = Qe*(Ue\(Le\(Pe*(Re\b_src)))); D_el = (U_el(irec,:) - D_bg) .* rec_mask; clear U_el Le Ue Pe Qe Re
    vc = v_true ./ (1 - (1./(pi*Q_true))*log(max(f,1e-6)/fref_phys) - 1i./(2*Q_true));
    [Lv,Uv,Pv,Qv,Rv] = lu(helm_op(vc, f, P));
    U_vi = Qv*(Uv\(Lv\(Pv*(Rv\b_src)))); D_vi = (U_vi(irec,:) - D_bg) .* rec_mask; clear U_vi Lv Uv Pv Qv Rv b_src
    N_f = conv2(g_nc, 1, sig_n * (1/sqrt(2)) * (randn(nr, ns) + 1i*randn(nr, ns)), 'same') .* rec_mask;
    V_el = extract_bp(D_el, Lm, Um, Pm, Qm, Rm, irec, lin, N, ns, P.dx);
    V_vi = extract_bp(D_vi, Lm, Um, Pm, Qm, Rm, irec, lin, N, ns, P.dx);
    V_nu = extract_bp(N_f,  Lm, Um, Pm, Qm, Rm, irec, lin, N, ns, P.dx);
    V_b1 = extract_bp((D_vi + N_f) .* EPS_B1, Lm, Um, Pm, Qm, Rm, irec, lin, N, ns, P.dx);
    clear Lm Um Pm Qm Rm
    I_ref = I_ref + real(sum(Us_c .* V_el, 2)); clear V_el
    cs = Us_c .* V_vi; cn = Us_c .* V_nu; cb = Us_c .* V_b1;
    Cs_sig_all(:,:,jf) = single(cs); Cs_nu_all(:,:,jf) = single(cn); Cs_b1_all(:,:,jf) = single(cb);
    I_R0_tot = I_R0_tot + real(sum(cs + cn, 2));
end

R0_2d = core.prep_image(I_R0_tot, P.nzi, P.nxi, z_in, P.dx, cfg.bg_sigma_m, cfg.z_mute_m);
[phi_eff, c_st] = core.blind_dip(R0_2d, P.dx, cfg.dip_sigma_m, cfg.dip_cap_deg);
[w_dsa, TS_dsa] = core.build_dsa_in(Xin(:), Zin(:), 20.0, xs, xr_lo, xr_hi, phi_eff(:), Fts_in, cfg);
[w_sa,  ~     ] = core.build_dsa_in(Xin(:), Zin(:), 20.0, xs, xr_lo, xr_hi, zeros(nin,1), Fts_in, cfg);
dip_stats = struct('med_abs_phi_deg', median(abs(rad2deg(phi_eff(w_deep_mask)))), ...
                   'p90_abs_phi_deg', prctile(abs(rad2deg(phi_eff(w_deep_mask))), 90), ...
                   'med_coh', median(c_st(w_deep_mask)));

I_raw_s = zeros(nin, nf, 'single'); I_raw_n = I_raw_s; I_raw_t = I_raw_s; I_raw_f1 = I_raw_s; I_raw_f3 = I_raw_s;
I_dsa_s = zeros(nin, nf, 'single'); I_dsa_n = I_dsa_s; I_dsa_t = I_dsa_s; I_dsa_f1 = I_dsa_s; I_dsa_f3 = I_dsa_s;
I_sa_t  = zeros(nin, nf, 'single'); I_sa_f3 = I_sa_t;
for jf = 1:nf
    cs = double(Cs_sig_all(:,:,jf)); cn = double(Cs_nu_all(:,:,jf)); ct = cs + cn; cb = double(Cs_b1_all(:,:,jf));
    I_raw_s(:,jf)=single(real(sum(cs,2))); I_raw_n(:,jf)=single(real(sum(cn,2)));
    I_raw_t(:,jf)=single(real(sum(ct,2))); I_raw_f1(:,jf)=single(real(sum(cb,2))); I_raw_f3(:,jf)=single(real(sum(ct.*sgn_B300,2)));
    I_dsa_s(:,jf)=single(real(sum(w_dsa.*cs,2))); I_dsa_n(:,jf)=single(real(sum(w_dsa.*cn,2)));
    I_dsa_t(:,jf)=single(real(sum(w_dsa.*ct,2))); I_dsa_f1(:,jf)=single(real(sum(w_dsa.*cb,2))); I_dsa_f3(:,jf)=single(real(sum((w_dsa.*ct).*sgn_B300,2)));
    I_sa_t(:,jf) =single(real(sum(w_sa.*ct,2)));  I_sa_f3(:,jf)=single(real(sum((w_sa.*ct).*sgn_B300,2)));
end
clear Cs_b1_all

iv_raw_b1 = core.band_inv_snr(I_raw_t, I_raw_f1, freqs, P.nzi, P.nxi, P.dx, cfg.snr_band_hz, cfg.snr_box_m, cfg.snr_cap);
iv_raw_b3 = core.band_inv_snr(I_raw_t, I_raw_f3, freqs, P.nzi, P.nxi, P.dx, cfg.snr_band_hz, cfg.snr_box_m, cfg.snr_cap);
iv_dsa_b1 = core.band_inv_snr(I_dsa_t, I_dsa_f1, freqs, P.nzi, P.nxi, P.dx, cfg.snr_band_hz, cfg.snr_box_m, cfg.snr_cap);
iv_dsa_b3 = core.band_inv_snr(I_dsa_t, I_dsa_f3, freqs, P.nzi, P.nxi, P.dx, cfg.snr_band_hz, cfg.snr_box_m, cfg.snr_cap);
iv_sa_b3  = core.band_inv_snr(I_sa_t,  I_sa_f3,  freqs, P.nzi, P.nxi, P.dx, cfg.snr_band_hz, cfg.snr_box_m, cfg.snr_cap);
iv_tr_raw = true_inv_snr(I_raw_s, I_raw_n, freqs, P.nzi, P.nxi, P.dx, cfg.snr_band_hz, cfg.snr_box_m, cfg.snr_cap);
iv_tr_dsa = true_inv_snr(I_dsa_s, I_dsa_n, freqs, P.nzi, P.nxi, P.dx, cfg.snr_band_hz, cfg.snr_box_m, cfg.snr_cap);

fa = (freqs >= 10.0 & freqs <= 35.0);
l2e = @(e, t) median(arrayfun(@(j) median(log2(max(double(e(w_deep_mask,j)),1e-6)./max(double(t(w_deep_mask,j)),1e-6))), find(fa)));
noise_acc = struct('RAW_B1_log2err', l2e(iv_raw_b1, iv_tr_raw), 'RAW_B300_log2err', l2e(iv_raw_b3, iv_tr_raw), ...
                   'DSA_B1_log2err', l2e(iv_dsa_b1, iv_tr_dsa), 'DSA_B300_log2err', l2e(iv_dsa_b3, iv_tr_dsa));
fprintf('=== WP4 %s (%s, snr_db=%+.1f dB) | Noise Acc: RAW B1=%+.2f B300=%+.2f | DSA B1=%+.2f B300=%+.2f ===\n', ...
    mname, q_err_desc, snr_db, noise_acc.RAW_B1_log2err, noise_acc.RAW_B300_log2err, noise_acc.DSA_B1_log2err, noise_acc.DSA_B300_log2err);

mnames = {'RAW_R0','DSA_R0','RAW_UQW1_25','RAW_UQW1b_25','RAW_UQW1b_100','SA_UQW1b_25','DSA_PH_25','DSA_IQ70_25','DSA_W1b_25','DSA_UQW1_25','DSA_UQW1b_25'};
IS = struct(); IN = struct();
for k = 1:numel(mnames), nc = 1+(nm-1)*~endsWith(mnames{k},'_R0'); IS.(mnames{k})=zeros(nin,nc); IN.(mnames{k})=zeros(nin,nc); end
xmax = 1.05 * max(max(TS_flat(:)), max(TS_dsa(:))) / min(ms); xg = linspace(0, xmax, 1201); dxg = xg(2);

for jf = 1:nf
    f = freqs(jf); Cs_s = double(Cs_sig_all(:,:,jf)); Cs_n = double(Cs_nu_all(:,:,jf));
    T0_25 = uq_tab(f, xg, 0, 25); T1_25 = uq_tab(f, xg, sig, 25); T1_100 = uq_tab(f, xg, sig, 100);
    IS.RAW_R0 = IS.RAW_R0 + real(sum(Cs_s,2)); IN.RAW_R0 = IN.RAW_R0 + real(sum(Cs_n,2));
    IS.DSA_R0 = IS.DSA_R0 + real(sum(w_dsa.*Cs_s,2)); IN.DSA_R0 = IN.DSA_R0 + real(sum(w_dsa.*Cs_n,2));
    ivr1=double(iv_raw_b1(:,jf)); ivr3=double(iv_raw_b3(:,jf)); ivd1=double(iv_dsa_b1(:,jf)); ivd3=double(iv_dsa_b3(:,jf)); ivs3=double(iv_sa_b3(:,jf));
    for im = 1:nm
        u=(TS_flat/ms(im))/dxg; i0=min(floor(u),numel(xg)-2); ww=u-i0; i0=i0+1; it=@(tb) tb(i0).*(1-ww)+tb(i0+1).*ww;
        A1=it(T1_25.A); C1=it(T1_25.C); P1=it(T1_25.P);
        H1=(A1+s70*P1.*C1+P1.*ivr1)./(1+s70*C1+ivr1); H3=(A1+s70*P1.*C1+P1.*ivr3)./(1+s70*C1+ivr3); Hs=(A1+s70*P1.*C1+P1.*ivs3)./(1+s70*C1+ivs3);
        IS.RAW_UQW1_25(:,im)=IS.RAW_UQW1_25(:,im)+real(sum(H1.*Cs_s,2)); IN.RAW_UQW1_25(:,im)=IN.RAW_UQW1_25(:,im)+real(sum(H1.*Cs_n,2));
        IS.RAW_UQW1b_25(:,im)=IS.RAW_UQW1b_25(:,im)+real(sum(H3.*Cs_s,2)); IN.RAW_UQW1b_25(:,im)=IN.RAW_UQW1b_25(:,im)+real(sum(H3.*Cs_n,2));
        IS.SA_UQW1b_25(:,im)=IS.SA_UQW1b_25(:,im)+real(sum((w_sa.*Hs).*Cs_s,2)); IN.SA_UQW1b_25(:,im)=IN.SA_UQW1b_25(:,im)+real(sum((w_sa.*Hs).*Cs_n,2));
        A1h=it(T1_100.A); C1h=it(T1_100.C); P1h=it(T1_100.P); Hh=(A1h+s70*P1h.*C1h+P1h.*ivr3)./(1+s70*C1h+ivr3);
        IS.RAW_UQW1b_100(:,im)=IS.RAW_UQW1b_100(:,im)+real(sum(Hh.*Cs_s,2)); IN.RAW_UQW1b_100(:,im)=IN.RAW_UQW1b_100(:,im)+real(sum(Hh.*Cs_n,2));

        u=(TS_dsa/ms(im))/dxg; i0=min(floor(u),numel(xg)-2); ww=u-i0; i0=i0+1; it=@(tb) tb(i0).*(1-ww)+tb(i0+1).*ww;
        A0=it(T0_25.A); C0=it(T0_25.C); P0=it(T0_25.P); A1=it(T1_25.A); C1=it(T1_25.C); P1=it(T1_25.P);
        Hiq=(A0+s70*P0.*C0)./(1+s70*C0); Hw1b=(A0+s70*P0.*C0+P0.*ivd3)./(1+s70*C0+ivd3);
        Huq1=(A1+s70*P1.*C1+P1.*ivd1)./(1+s70*C1+ivd1); Huqb=(A1+s70*P1.*C1+P1.*ivd3)./(1+s70*C1+ivd3);
        IS.DSA_PH_25(:,im)=IS.DSA_PH_25(:,im)+real(sum((w_dsa.*P0).*Cs_s,2)); IN.DSA_PH_25(:,im)=IN.DSA_PH_25(:,im)+real(sum((w_dsa.*P0).*Cs_n,2));
        IS.DSA_IQ70_25(:,im)=IS.DSA_IQ70_25(:,im)+real(sum((w_dsa.*Hiq).*Cs_s,2)); IN.DSA_IQ70_25(:,im)=IN.DSA_IQ70_25(:,im)+real(sum((w_dsa.*Hiq).*Cs_n,2));
        IS.DSA_W1b_25(:,im)=IS.DSA_W1b_25(:,im)+real(sum((w_dsa.*Hw1b).*Cs_s,2)); IN.DSA_W1b_25(:,im)=IN.DSA_W1b_25(:,im)+real(sum((w_dsa.*Hw1b).*Cs_n,2));
        IS.DSA_UQW1_25(:,im)=IS.DSA_UQW1_25(:,im)+real(sum((w_dsa.*Huq1).*Cs_s,2)); IN.DSA_UQW1_25(:,im)=IN.DSA_UQW1_25(:,im)+real(sum((w_dsa.*Huq1).*Cs_n,2));
        IS.DSA_UQW1b_25(:,im)=IS.DSA_UQW1b_25(:,im)+real(sum((w_dsa.*Huqb).*Cs_s,2)); IN.DSA_UQW1b_25(:,im)=IN.DSA_UQW1b_25(:,im)+real(sum((w_dsa.*Huqb).*Cs_n,2));
    end
end

ref_2d = core.prep_image(I_ref, P.nzi, P.nxi, z_in, P.dx, cfg.bg_sigma_m, cfg.z_mute_m);
ref_hf = core.highpass_kz(ref_2d, P.dx, 10.0);
rms_w = @(A, W) sqrt(mean(A(z_in>=W{2}(1)&z_in<=W{2}(2), x_in>=W{1}(1)&x_in<=W{1}(2)).^2, 'all'));
ref_rcf = log10(rms_w(ref_2d, W_targ) / max(rms_w(ref_2d, W_flank), 1e-30));

out = struct('model', mname, 'q_err_desc', q_err_desc, 'seed', seed, 'snr_db', snr_db, 'm', ms, ...
             'dip_stats', dip_stats, 'noise_acc', noise_acc, 'W_deep', struct(), 'elapsed_s', 0);
for k = 1:numel(mnames)
    key = mnames{k}; Sm = IS.(key); Nm = IN.(key); nc = size(Sm,2);
    e1=zeros(1,nc); e1h=e1; eamp=e1; snr=e1;
    for c=1:nc
        S2=core.prep_image(Sm(:,c),P.nzi,P.nxi,z_in,P.dx,cfg.bg_sigma_m,cfg.z_mute_m);
        N2=core.prep_image(Nm(:,c),P.nzi,P.nxi,z_in,P.dx,cfg.bg_sigma_m,cfg.z_mute_m);
        T2=S2+N2; Th=core.highpass_kz(T2,P.dx,10.0);
        e1(c)=core.score_e1(T2,ref_2d,iz_w,ix_w); e1h(c)=core.score_e1(Th,ref_hf,iz_w,ix_w);
        snr(c)=10*log10(sum(S2(iz_w,ix_w).^2,'all')/max(sum(N2(iz_w,ix_w).^2,'all'),1e-30));
        eamp(c)=abs(log10(rms_w(T2,W_targ)/max(rms_w(T2,W_flank),1e-30))-ref_rcf);
    end
    out.W_deep.(key)=struct('E1',e1,'E1_mean',mean(e1),'E1_HF',e1h,'E1_HF_mean',mean(e1h),'E_amp',eamp,'E_amp_mean',mean(eamp),'SNR_img_dB',snr,'SNR_img_dB_mean',mean(snr));
    fprintf('  %-15s | E1=%.4f (min=%.3f) | HF=%.4f | Eamp=%.3f | SNR=%+.2f dB\n', key, mean(e1), min(e1), mean(e1h), mean(eamp), mean(snr));
end
out.elapsed_s = toc(t0_all);
end

function iv = true_inv_snr(Is_mat, In_mat, freqs, nzi, nxi, dx, band_hz, box_m, cap)
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
