function precheck_calib_and_fd()
% Pre-Freeze Calibration & Compatibility Check (NO production scores!)
%   1) Calibrate snr_db for S1-S4 so RAW_R0 W_deep SNR_img in [+2.5, +5.5] dB
%   2) Verify fd2d_visco.m <-> helm_op_visco exact wavelet/polarity/phase sign match + timing
%   3) Measure steady-state Volve 10m viscoacoustic A_comp LU timing & stability

maxNumCompThreads(12);
addpath('D:\ChatGPT\研究\地震衰减补偿研究\uq_rtm_20261117\scripts');
addpath(fileparts(mfilename('fullpath')));
core = uq_rtm_xmir_qrtm_core();
cfg  = core.default_cfg('synth_10m');

%% 1) Check S1-S4 10-freq SNR calibration
fprintf('=== 1. Pre-Freeze SNR Calibration Smoke Check (S1-S4) ===\n');
models = {'S1_Chimney_Base', 'S2_StrongGrad_SteepDip', 'S3_SharpQ_Multipath', 'S4_NonScalar_SpatialQ'};
snr_try = [0.5, -1.5, 0.5, 0.5];
for im = 1:numel(models)
    [P, v_true, v_mig, Q_true, ~, rec_mask, iz_w, ix_w] = build_synth_model(models{im}, core);
    ns = numel(P.ixs); nr = numel(P.ixr); N = P.nz * P.nx;
    [Zi, Xi] = ndgrid(P.in{1}, P.in{2}); lin = sub2ind([P.nz P.nx], Zi(:), Xi(:)); nin = numel(lin);
    isrc = sub2ind([P.nz P.nx], P.izs*ones(1,ns), P.ixs);
    irec = sub2ind([P.nz P.nx], P.izs*ones(1,nr), P.ixr);
    z_in = P.z(P.in{1});

    rng(20261500 + im, 'twister');
    fref = 25.0; fpeak = 22.0; corr_sig = 5.0;
    rc = ceil(3*corr_sig); g_nc = exp(-0.5*((-rc:rc)'/corr_sig).^2); g_nc = g_nc / sqrt(sum(g_nc.^2));

    % Calibrate at 22 Hz
    W22 = 1.0 * exp(-1.0);
    b22 = zeros(N, ns); b22(sub2ind([N ns], isrc, 1:ns)) = -W22 / P.dx^2;
    [Lm,Um,Pm,Qm,Rm] = lu(core.helm_op_visco(v_mig, v_mig, 22, P));
    Ub = core.solve_helm_shot(Lm,Um,Pm,Qm,Rm, b22); Dbg22 = Ub(irec,:); clear Ub Lm Um Pm Qm Rm
    vc22 = core.make_visco_vel(v_true, Q_true, 22, fref, P.npml, 'att', 1.0);
    [Lv,Uv,Pv,Qv,Rv] = lu(core.helm_op_visco(vc22, v_true, 22, P));
    Uv22 = core.solve_helm_shot(Lv,Uv,Pv,Qv,Rv, b22); Dcal = (Uv22(irec,:) - Dbg22) .* rec_mask; clear Uv22 Lv Uv Pv Qv Rv
    sig_n = sqrt(mean(abs(Dcal(rec_mask)).^2)) * 10^(-snr_try(im)/20);

    freqs = 6.0:3.0:33.0;
    I_s = zeros(nin, 1); I_n = zeros(nin, 1);
    for jf = 1:numel(freqs)
        f = freqs(jf); w = 2*pi*f;
        W_src = (f/fpeak)^2 * exp(-(f/fpeak)^2);
        b_src = zeros(N, ns); b_src(sub2ind([N ns], isrc, 1:ns)) = -W_src / P.dx^2;
        [Lm,Um,Pm,Qm,Rm] = lu(core.helm_op_visco(v_mig, v_mig, f, P));
        U_bg = core.solve_helm_shot(Lm,Um,Pm,Qm,Rm, b_src); D_bg = U_bg(irec,:); Us_c = conj(U_bg(lin,:))*w^2; clear U_bg
        vc = core.make_visco_vel(v_true, Q_true, f, fref, P.npml, 'att', 1.0);
        [Lv,Uv,Pv,Qv,Rv] = lu(core.helm_op_visco(vc, v_true, f, P));
        U_vi = core.solve_helm_shot(Lv,Uv,Pv,Qv,Rv, b_src); D_vi = (U_vi(irec,:) - D_bg) .* rec_mask; clear U_vi Lv Uv Pv Qv Rv
        N_f = conv2(g_nc, 1, sig_n * (1/sqrt(2)) * (randn(nr, ns) + 1i*randn(nr, ns)), 'same') .* rec_mask;
        bs = zeros(N,ns); bs(irec,:)=D_vi; Vs = core.solve_helm_adj(Lm,Um,Pm,Qm,Rm, bs, P.dx);
        bn = zeros(N,ns); bn(irec,:)=N_f;  Vn = core.solve_helm_adj(Lm,Um,Pm,Qm,Rm, bn, P.dx);
        clear Lm Um Pm Qm Rm
        I_s = I_s + real(sum(Us_c .* Vs(lin,:), 2));
        I_n = I_n + real(sum(Us_c .* Vn(lin,:), 2));
    end
    S2 = core.prep_image(I_s, P.nzi, P.nxi, z_in, P.dx, cfg.bg_sigma_m, cfg.z_mute_m);
    N2 = core.prep_image(I_n, P.nzi, P.nxi, z_in, P.dx, cfg.bg_sigma_m, cfg.z_mute_m);
    snr_img = 10*log10(sum(S2(iz_w,ix_w).^2,'all') / max(sum(N2(iz_w,ix_w).^2,'all'), 1e-30));
    % Note: 34-freq df=1 Hz SNR_img is ~+3.0 to +3.5 dB higher than 10-freq df=3 Hz!
    fprintf('  [%-24s] snr_db=%+.1f dB -> 10-freq SNR_img = %+.2f dB (est 34-freq ~ %+.2f dB)\n', ...
        models{im}, snr_try(im), snr_img, snr_img + 3.2);
end

%% 2) Check fd2d_visco.m <-> helm_op_visco phase/polarity match & speed
fprintf('\n=== 2. Time-Domain GSLS FD (fd2d_visco.m) <-> Helmholtz Phase/Polarity Verification ===\n');
[P5, v_true5, v_mig5, Q_true5, ~, ~, iz_w5, ix_w5] = build_synth_model('S1_Chimney_Base', core);
ns5 = numel(P5.ixs); nr5 = numel(P5.ixr); N5 = P5.nz * P5.nx;
[Zi5, Xi5] = ndgrid(P5.in{1}, P5.in{2}); lin5 = sub2ind([P5.nz P5.nx], Zi5(:), Xi5(:));
isrc5 = sub2ind([P5.nz P5.nx], P5.izs*ones(1,ns5), P5.ixs);
irec5 = sub2ind([P5.nz P5.nx], P5.izs*ones(1,nr5), P5.ixr);
z_in5 = P5.z(P5.in{1});

opt.dx = P5.dx; opt.dt = 0.001; opt.nt = 1500; opt.fref = 25.0; opt.npml = P5.npml;
opt.rec = [P5.izs * ones(nr5, 1), P5.ixr(:)];
t_vec = (0:opt.nt-1)' * opt.dt; t0 = 0.08; fpeak = 22.0;
tau_r = pi * fpeak * (t_vec - t0);
q_rik = (1 - 2*tau_r.^2) .* exp(-tau_r.^2); % zero-phase Ricker centered at t0
% Note: fd2d_visco adds dp/dt += q(t), so pressure source spectrum in Helmholtz A U = -F is i*w*Q(w)/c^2 or -i*w*Q(w)!
% Let's check exact frequency-domain match on 5 test shots!
test_sh = [9, 17, 25, 33, 41];
M_el = struct('v', v_true5, 'rho', ones(P5.nz, P5.nx), 'Q', inf(P5.nz, P5.nx));
M_bg = struct('v', v_mig5,  'rho', ones(P5.nz, P5.nx), 'Q', inf(P5.nz, P5.nx));
M_vi = struct('v', v_true5, 'rho', ones(P5.nz, P5.nx), 'Q', Q_true5);

t_fd = tic;
D_fd_el = zeros(opt.nt, nr5, numel(test_sh));
D_fd_vi = zeros(opt.nt, nr5, numel(test_sh));
for k = 1:numel(test_sh)
    opt.src = [P5.izs, P5.ixs(test_sh(k))];
    r_bg = fd2d_visco(M_bg, q_rik, opt);
    r_el = fd2d_visco(M_el, q_rik, opt);
    r_vi = fd2d_visco(M_vi, q_rik, opt);
    D_fd_el(:,:,k) = r_el - r_bg;
    D_fd_vi(:,:,k) = r_vi - r_bg;
end
dt_fd5 = toc(t_fd);
fprintf('  fd2d_visco 5 shots x 3 runs (15 FD sims): %.2fs (est 49 shots x 3 = %.1fs)\n', dt_fd5, dt_fd5 * 49/5);

% FFT of time-domain data and source wavelet
nfft = 2000; % df = 1/(2000*0.001) = 0.5 Hz -> exact integer freqs 4:1:37 Hz!
df_fd = 1 / (nfft * opt.dt);
Q_fft = fft(q_rik, nfft);
freqs_fd = 4.0:2.0:36.0;
I_fd_plus = zeros(numel(lin5), 1); I_fd_minus = zeros(numel(lin5), 1);
I_fd_iw   = zeros(numel(lin5), 1); I_helm_5   = zeros(numel(lin5), 1);
ns_t = numel(test_sh); isrc_t = isrc5(test_sh);
for jf = 1:numel(freqs_fd)
    f = freqs_fd(jf); w = 2*pi*f;
    idx_f = round(f / df_fd) + 1;
    % Shift t=0 to t0: multiply by exp(+i*w*t0)
    %Wait: fd2d_visco stores rec(n,:) after step n (time n*dt), whereas q_rik(n) at (n-1)*dt
    W_q = Q_fft(idx_f) * exp(1i * w * t0);
    D_el_f = squeeze(fft(D_fd_el, nfft, 1)); D_el_f = squeeze(D_el_f(idx_f, :, :)) * exp(1i * w * t0); % nr x ns_t

    [Lm,Um,Pm,Qm,Rm] = lu(core.helm_op_visco(v_mig5, v_mig5, f, P5));
    % In fd2d_visco: dp/dt = -K div v + q(t) -> d^2 p/dt^2 - c^2 lap p = dq/dt -> (-w^2 - c^2 lap) P = i*w*Q -> (lap + w^2/c^2) P = -(i*w/c^2)*Q
    % Let's compare source weights: W_src = i*w*W_q / v_src^2 vs real W_q!
    v_src = v_mig5(P5.izs, P5.ixs(test_sh(1)));
    b_fd = zeros(N5, ns_t); b_fd(sub2ind([N5 ns_t], isrc_t, 1:ns_t)) = -(1i*w * W_q / v_src^2) / P5.dx^2;
    U_fd = core.solve_helm_shot(Lm,Um,Pm,Qm,Rm, b_fd);
    Us_fd = conj(U_fd(lin5,:)) * w^2;
    b_r = zeros(N5, ns_t); b_r(irec5,:) = D_el_f;
    V_fd = core.solve_helm_adj(Lm,Um,Pm,Qm,Rm, b_r, P5.dx);
    I_fd_iw = I_fd_iw + real(sum(Us_fd .* V_fd(lin5,:), 2));

    % Also check Helmholtz-generated elastic reference for the same 5 shots!
    [Le,Ue,Pe,Qe,Re] = lu(core.helm_op_visco(v_true5, v_true5, f, P5));
    U_eh = core.solve_helm_shot(Le,Ue,Pe,Qe,Re, b_fd);
    D_eh = U_eh(irec5,:) - U_fd(irec5,:);
    b_eh = zeros(N5, ns_t); b_eh(irec5,:) = D_eh;
    V_eh = core.solve_helm_adj(Lm,Um,Pm,Qm,Rm, b_eh, P5.dx);
    I_helm_5 = I_helm_5 + real(sum(Us_fd .* V_eh(lin5,:), 2));
    clear Lm Um Pm Qm Rm Le Ue Pe Qe Re
end
R_fd = core.prep_image(I_fd_iw,  P5.nzi, P5.nxi, z_in5, P5.dx, 50, 160);
R_hm = core.prep_image(I_helm_5, P5.nzi, P5.nxi, z_in5, P5.dx, 50, 160);
cc_pos = core.cc(R_fd(iz_w5, ix_w5), R_hm(iz_w5, ix_w5));
e1_fd  = core.score_e1(R_fd, R_hm, iz_w5, ix_w5);
fprintf('  Elastic FD vs Helmholtz 5-shot image correlation: cc(L=0)=%+.4f | E1=%.4f\n', cc_pos, e1_fd);

%% 3) Measure steady-state Volve 10m viscoacoustic A_comp timing & stability
fprintf('\n=== 3. Volve 10m Viscoacoustic A_comp Steady-State Timing & Stability ===\n');
Pf = fld_setup();
ns_f = Pf.nshot; N_f = Pf.nz * Pf.nx;
isrc_f = sub2ind([Pf.nz Pf.nx], Pf.src_iz, Pf.src_ix);
irec_f = sub2ind([Pf.nz Pf.nx], Pf.rec_iz, Pf.rec_ix);
col_f  = sub2ind([N_f ns_f], irec_f, Pf.shot_of);
for f_test = [10.0, 18.0, 24.0]
    [~, jf] = min(abs(Pf.f - f_test));
    t_v = tic;
    vc_f = core.make_visco_vel(Pf.vmig, 0.5 * Pf.Qm, f_test, 25.0, Pf.npml, 'comp', 1.0);
    [L,U,Pp,Qq,R] = lu(core.helm_op_visco(vc_f, Pf.vmig, f_test, Pf));
    bs = zeros(N_f, ns_f); bs(sub2ind([N_f ns_f], isrc_f, (1:ns_f)')) = -Pf.W(jf) / Pf.dx^2;
    Us = core.solve_helm_shot(L,U,Pp,Qq,R, bs);
    br = zeros(N_f, ns_f); br(col_f) = Pf.D(:, jf);
    Vr = core.solve_helm_adj(L,U,Pp,Qq,R, br, Pf.dx);
    dt_v = toc(t_v);
    fprintf('  Volve 10m m=0.5 (Q=60) f=%4.1f Hz: %.2fs | finite=%d | max|Us|=%.2e max|Vr|=%.2e\n', ...
        f_test, dt_v, all(isfinite(Us(:))) && all(isfinite(Vr(:))), max(abs(Us(:))), max(abs(Vr(:))));
    clear L U Pp Qq R Us Vr bs br
end
end

function [P, v_true, v_mig, Q_true, Q_mig_base, rec_mask, iz_w, ix_w] = build_synth_model(mname, core)
P.dx = 10.0; P.npml = 18; P.nzi = 100; P.nxi = 160;
P.nz = P.nzi + 2*P.npml; P.nx = P.nxi + 2*P.npml;
P.z = ((1:P.nz)' - P.npml - 1) * P.dx; P.x = ((1:P.nx) - P.npml - 1) * P.dx;
P.in = {P.npml + (1:P.nzi), P.npml + (1:P.nxi)}; P.izs = P.npml + 3;
P.ixs = P.npml + 1 + (8:3:152); P.ixr = P.npml + 1 + (0:1:P.nxi-1);
ns = numel(P.ixs); nr = numel(P.ixr);
[Zfull, Xfull] = ndgrid(P.z, P.x);
z_in = P.z(P.in{1}); x_in = P.x(P.in{2});
xs = P.x(P.ixs); xr = P.x(P.ixr);
rec_mask = true(nr, ns);

switch mname
    case 'S1_Chimney_Base'
        z1 = 230.0 - 14.0 * exp(-((Xfull - 770.0)/330.0).^2);
        z2 = 625.0 - 36.0 * exp(-((Xfull - 790.0)/300.0).^2) + 0.08*(Xfull - 800.0);
        z3 = 760.0 - 26.0 * exp(-((Xfull - 790.0)/330.0).^2) + 0.06*(Xfull - 800.0);
        z4 = 890.0 - 14.0 * exp(-((Xfull - 800.0)/370.0).^2);
        v_true = 1850.0 * ones(P.nz, P.nx); v_true(Zfull>=z1)=2160; v_true(Zfull>=z2)=2580; v_true(Zfull>=z3)=2940; v_true(Zfull>=z4)=3320;
        invQ_t = (1/90)*ones(P.nz, P.nx); invQ_t(Zfull < 80) = 1/1000;
        invQ_t = invQ_t + (1/14 - 1/90)*exp(-(((Xfull - 790)/210).^2 + ((Zfull - 415)/140).^2).^1.5);
        Q_true = 1 ./ invQ_t; Q_mig_base = Q_true;
    case 'S2_StrongGrad_SteepDip'
        % Strong marine velocity gradient (1550 -> 3950 m/s) + lateral velocity wedge + 24 deg steep faulted dip
        % + Asymmetric marine towed-streamer array (off-end/split right cable xr in [xs-100, xs+750] m)!
        v_bg = 1550.0 + 2.1 * max(Zfull, 0.0) + 0.35 * (Xfull - 800.0) .* (max(Zfull,0)/1000.0);
        Zg = Zfull;
        hw = (Xfull < 840.0 - (Zfull - 180.0) / tand(62.0)) & (Zfull > 180.0); Zg(hw) = Zg(hw) + 35.0;
        z1 = 220.0 - 0.06*(Xfull - 800);
        z2 = 640.0 + 0.44*(Xfull - 800); % ~23.7 deg steep dip!
        z3 = 775.0 - 0.36*(Xfull - 800); % ~-19.8 deg opposite steep dip!
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
        % Asymmetric finite marine cable: each shot only records xr in [xs(k) - 120, xs(k) + 760] m
        for k = 1:ns
            rec_mask(:, k) = (xr(:) >= xs(k) - 120.0) & (xr(:) <= xs(k) + 760.0);
        end
    case 'S3_SharpQ_Multipath'
        % Ultra-sharp low-Q gas lens (Q_min = 12) + strong -22% velocity pull-down lens (severe multipathing/healing!)
        z1 = 220.0 - 10.0 * exp(-((Xfull - 800.0)/350.0).^2);
        z2 = 640.0 - 20.0 * exp(-((Xfull - 800.0)/320.0).^2) + 0.05*(Xfull - 800.0);
        z3 = 770.0 - 15.0 * exp(-((Xfull - 800.0)/340.0).^2) + 0.04*(Xfull - 800.0);
        z4 = 895.0 * ones(P.nz, P.nx);
        v_true = 1900.0 * ones(P.nz, P.nx); v_true(Zfull>=z1)=2250; v_true(Zfull>=z2)=2700; v_true(Zfull>=z3)=3050; v_true(Zfull>=z4)=3400;
        % Sharp compact low-velocity + ultra-low-Q gas lens (width rx=130m, rz=85m -> strong wavefield multipathing!)
        lens = exp(-(((Xfull - 800.0)/130.0).^2 + ((Zfull - 410.0)/85.0).^2).^2.0);
        v_true = v_true .* (1.0 - 0.22 * lens);
        invQ_t = (1/95)*ones(P.nz, P.nx); invQ_t(Zfull < 80) = 1/1000;
        invQ_t = invQ_t + (1/12 - 1/95) * lens;
        Q_true = 1 ./ invQ_t; Q_mig_base = Q_true;
    case 'S4_NonScalar_SpatialQ'
        z1 = 230.0 - 14.0 * exp(-((Xfull - 770.0)/330.0).^2);
        z2 = 625.0 - 36.0 * exp(-((Xfull - 790.0)/300.0).^2) + 0.08*(Xfull - 800.0);
        z3 = 760.0 - 26.0 * exp(-((Xfull - 790.0)/330.0).^2) + 0.06*(Xfull - 800.0);
        z4 = 890.0 - 14.0 * exp(-((Xfull - 800.0)/370.0).^2);
        v_true = 1850.0 * ones(P.nz, P.nx); v_true(Zfull>=z1)=2160; v_true(Zfull>=z2)=2580; v_true(Zfull>=z3)=2940; v_true(Zfull>=z4)=3320;
        invQ_t = (1/90)*ones(P.nz, P.nx); invQ_t(Zfull < 80) = 1/1000;
        invQ_t = invQ_t + (1/14 - 1/90)*exp(-(((Xfull - 790)/210).^2 + ((Zfull - 415)/140).^2).^1.5);
        Q_true = 1 ./ invQ_t;
        % Combined spatial error: shifted +160m, widened 1.35x, AND shallow 1.3x / deep 0.75x bias
        invQ_m = (1/90)*ones(P.nz, P.nx); invQ_m(Zfull < 80) = 1/1000;
        invQ_m = invQ_m + (1/14 - 1/90)*exp(-(((Xfull - 950)/(210*1.35)).^2 + ((Zfull - 415)/140).^2).^1.5);
        Q_mig_base = 1 ./ invQ_m;
        Q_mig_base(Zfull >= 80 & Zfull < 480) = 1.30 * Q_mig_base(Zfull >= 80 & Zfull < 480);
        Q_mig_base(Zfull >= 480)              = 0.75 * Q_mig_base(Zfull >= 480);
end
v_mig = core.smooth2d(v_true, 2.5);
iz_w = find(z_in >= 560 & z_in <= 950);
ix_w = find(x_in >= 180 & x_in <= 1420);
end
