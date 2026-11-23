function run_unified_synth(models_req, smoke)
% Tier 1 Pre-registered Unified Multi-Physics Synthetic Confluence Suite (U1-U5)
% Plan: uq_rtm_unified_20261124/PLAN_ZH.md
%
% Simultaneously couples ALL FOUR physical & statistical challenges under in-band
% migration velocity convention v_band(25 Hz) with brand-new independent seeds (20261301-20261305):
%   U1: Gas chimney + dipping reflectors + laterally correlated noise l30 (60 m) + v_band(25 Hz)
%   U2: Reverse fault + steep opposing dips (+22/-16 deg) + l30 (60 m) + v_band(25 Hz)
%   U3: Volve-like one-sided trailing streamer (~65% out-of-array) + chimney + l30 + v_band(25 Hz)
%   U4: Exact Volve field twin (unknown Ricker -> flat W_mig=1 + (2*pi*f)^2 + one-sided streamer
%       + l30 + v_band(25 Hz) + 5 finely layered synthetic wells)
%   U5: High-SNR (15 dB) + -3% velocity error + lateral Q gradient safety boundary test
%
% Usage:
%   run_unified_synth({'U1'}, true)   % 3-freq smoke check
%   run_unified_synth([], false)      % full U1-U5 production suite

if nargin < 1 || isempty(models_req), models_req = {'U1', 'U2', 'U3', 'U4', 'U5'}; end
if ischar(models_req) || isstring(models_req), models_req = {char(models_req)}; end
if nargin < 2, smoke = false; end

maxNumCompThreads(8);
here = fileparts(mfilename('fullpath'));
root = fileparts(here);
res_dir = fullfile(root, 'results');
if ~isfolder(res_dir), mkdir(res_dir); end

% Add uq_rtm_20261117 engine helpers from same worktree (read-only)
wt_root = fileparts(root);
addpath(fullfile(wt_root, 'uq_rtm_20261117', 'scripts'));
core = uq_rtm_unified_core();
cfg  = core.default_cfg('synth_10m');

for imod = 1:numel(models_req)
    mname = models_req{imod};
    tag = ''; if smoke, tag = '_smoke'; end
    fn = fullfile(res_dir, ['synth_' mname tag '.json']);
    if isfile(fn) && ~smoke
        fprintf('Skipping %s (already exists: %s)\n', mname, fn);
        continue;
    end
    if strcmp(mname, 'U4')
        res = run_u4_volve_twin(core, cfg, smoke);
    else
        res = run_u_model(mname, core, cfg, smoke);
    end
    fid = fopen(fn, 'w'); fwrite(fid, jsonencode(res, 'PrettyPrint', true), 'char'); fclose(fid);
    fprintf('Saved %s (%.1fs)\n\n', fn, res.elapsed_s);
end
end

%% ========================================================================
%% U1, U2, U3, U5: Full-wave viscoelastic Helmholtz under v_band(25 Hz)
%% ========================================================================
function res = run_u_model(mname, core, cfg, smoke)
t0_all = tic;
P.dx = 10.0; P.npml = 18; P.nzi = 100; P.nxi = 160;
P.nz = P.nzi + 2*P.npml; P.nx = P.nxi + 2*P.npml;
P.z = ((1:P.nz)' - P.npml - 1) * P.dx;
P.x = ((1:P.nx)  - P.npml - 1) * P.dx;
P.in = {P.npml + (1:P.nzi), P.npml + (1:P.nxi)};
P.izs = P.npml + 3;
P.ixs = P.npml + 1 + (8:3:152);       % 49 shots, 30 m spacing (x = 80..1520 m, matching est_block 51 shots)
P.ixr = P.npml + 1 + (0:1:P.nxi-1);   % 160 receivers, 10 m spacing (x = 0..1590 m)
ns = numel(P.ixs); nr = numel(P.ixr); N = P.nz * P.nx;
[Zi, Xi] = ndgrid(P.in{1}, P.in{2}); lin = sub2ind([P.nz P.nx], Zi(:), Xi(:)); nin = numel(lin);
isrc = sub2ind([P.nz P.nx], P.izs*ones(1,ns), P.ixs);
irec = sub2ind([P.nz P.nx], P.izs*ones(1,nr), P.ixr);
[Zfull, Xfull] = ndgrid(P.z, P.x);
z_in = P.z(P.in{1}); x_in = P.x(P.in{2});
[Zin, Xin] = ndgrid(z_in, x_in);
xs = P.x(P.ixs); xr = P.x(P.ixr);

rec_mask = true(nr, ns);
vscale_mig = 1.0;
corr_noise_sig = 5.0;  % 5 traces = 50 m laterally correlated noise (exact est_block l30 scale!)
snr_db = -18.0;        % prestack l30 calibration -> realistic ~2-4 dB image SNR in W_deep!

switch mname
    case 'U1'
        % Gas chimney + anticline/dipping reflectors + l30 (60 m) + v_band(25 Hz)
        seed = 20261301;
        z1 = 230.0 - 15.0 * exp(-((Xfull - 760.0)/330.0).^2);
        z2 = 625.0 - 38.0 * exp(-((Xfull - 800.0)/300.0).^2) + 0.08*(Xfull - 800.0);
        z3 = 760.0 - 28.0 * exp(-((Xfull - 800.0)/330.0).^2) + 0.06*(Xfull - 800.0);
        z4 = 890.0 - 14.0 * exp(-((Xfull - 800.0)/370.0).^2);
        v_true = 1850.0 * ones(P.nz, P.nx);
        v_true(Zfull >= z1) = 2160.0;
        v_true(Zfull >= z2) = 2580.0;
        v_true(Zfull >= z3) = 2940.0;
        v_true(Zfull >= z4) = 3320.0;
        invQ = (1/90)*ones(P.nz, P.nx); invQ(Zfull < 80) = 1/1000;
        chim = exp(-(((Xfull - 790)/210).^2 + ((Zfull - 415)/140).^2).^1.5);
        invQ = invQ + (1/14 - 1/90)*chim;
        Q_true = 1./invQ;
        W_targ = {[580 1000], [560 950]}; W_flank = {[200 480], [560 950]};

    case 'U2'
        % Reverse fault + opposing steep dips (+22/-16 deg) + low-Q wedge + l30 (60 m) + v_band(25 Hz)
        seed = 20261302;
        Zg = Zfull;
        hw = (Xfull < 820.0 - (Zfull - 200.0) / tand(63.0)) & (Zfull > 200.0);
        Zg(hw) = Zg(hw) + 32.0;
        z1 = 235.0 - 0.05 * (Xfull - 800.0);
        z2 = 655.0 + 0.40 * (Xfull - 800.0);   % +21.8 deg
        z3 = 785.0 - 0.29 * (Xfull - 800.0);   % -16.2 deg
        z4 = 900.0 - 0.10 * (Xfull - 800.0);
        v_true = 1860.0 * ones(P.nz, P.nx);
        v_true(Zg >= z1) = 2170.0;
        v_true(Zg >= z2) = 2630.0;
        v_true(Zg >= z3) = 2990.0;
        v_true(Zg >= z4) = 3360.0;
        invQ = (1/85)*ones(P.nz, P.nx); invQ(Zfull < 80) = 1/1000;
        wedge = exp(-(((Xfull - 960)/250).^2 + ((Zfull - 365)/100).^2).^1.5);
        invQ = invQ + (1/15 - 1/85)*wedge;
        Q_true = 1./invQ;
        W_targ = {[780 1200], [560 950]}; W_flank = {[220 560], [560 950]};

    case 'U3'
        % Volve-like one-sided marine streamer (~65% out-of-array) + chimney + 11 deg dip + l30 + v_band(25 Hz)
        seed = 20261303;
        for s = 1:ns
            rec_mask(:, s) = (xr >= xs(s) - 240.0) & (xr <= xs(s) + 860.0);
        end
        z1 = 240.0 + 0.04 * (Xfull - 800.0);
        z2 = 640.0 + 0.20 * (Xfull - 800.0);   % +11.3 deg
        z3 = 770.0 + 0.17 * (Xfull - 800.0);   % +9.6 deg
        z4 = 890.0 - 0.04 * (Xfull - 800.0);
        v_true = 1840.0 * ones(P.nz, P.nx);
        v_true(Zfull >= z1) = 2140.0;
        v_true(Zfull >= z2) = 2560.0;
        v_true(Zfull >= z3) = 2910.0;
        v_true(Zfull >= z4) = 3280.0;
        invQ = (1/90)*ones(P.nz, P.nx); invQ(Zfull < 80) = 1/1000;
        chim = exp(-(((Xfull - 780)/195).^2 + ((Zfull - 415)/145).^2).^1.5);
        invQ = invQ + (1/14 - 1/90)*chim;
        Q_true = 1./invQ;
        W_targ = {[580 980], [560 950]}; W_flank = {[1040 1380], [560 950]};

    case 'U5'
        % High-SNR (clean) + -3% velocity error + lateral Q gradient + folds safety boundary test
        seed = 20261305;
        vscale_mig = 0.97;
        corr_noise_sig = 0.0;  % clean independent noise (~18 dB image SNR)
        snr_db = 0.0;
        z1 = 235.0 + 14.0 * sin(2*pi*Xfull/880.0);
        z2 = 635.0 + 22.0 * sin(2*pi*(Xfull - 180)/840.0);
        z3 = 765.0 + 18.0 * cos(2*pi*(Xfull - 140)/920.0);
        z4 = 885.0 - 12.0 * sin(2*pi*Xfull/1080.0);
        v_true = 1850.0 * ones(P.nz, P.nx);
        v_true(Zfull >= z1) = 2150.0;
        v_true(Zfull >= z2) = 2560.0;
        v_true(Zfull >= z3) = 2920.0;
        v_true(Zfull >= z4) = 3300.0;
        Q_true = 110.0 - 82.0 * min(max(Xfull / 1600.0, 0), 1);
        Q_true(Zfull < 80) = 1000.0;
        W_targ = {[980 1400], [560 950]}; W_flank = {[200 600], [560 950]};

    otherwise
        error('Unknown model %s', mname);
end

rng(seed, 'twister');
v_mig_exact = core.smooth2d(v_true, 2.5);
v_mig = vscale_mig * v_mig_exact;

% Single-leg eikonal t* fields (reused for BOTH flat & dip-aware two-way t*)
s_mig = 1.0 ./ v_mig; eps_q = 0.5;
F.T = zeros(P.nz, P.nx, ns, 'single'); F.ts = F.T;
for k = 1:ns
    F.T(:,:,k)  = eikonal_fs(s_mig, P.dx, P.izs, P.ixs(k));
    F.ts(:,:,k) = (eikonal_fs(s_mig.*(1+eps_q./Q_true), P.dx, P.izs, P.ixs(k)) - ...
                   eikonal_fs(s_mig.*(1-eps_q./Q_true), P.dx, P.izs, P.ixs(k))) / (2*eps_q);
end
TS_flat = zeros(nin, ns); Fts_in = zeros(nin, ns);
for k = 1:ns
    [~, tt] = two_way(F, k, P); TS_flat(:, k) = double(tt(lin));
    tk = F.ts(:,:,k); Fts_in(:, k) = double(tk(lin));
end
clear F

xr_lo = zeros(1, ns); xr_hi = zeros(1, ns);
for k = 1:ns
    live_xr = xr(rec_mask(:, k));
    xr_lo(k) = min(live_xr); xr_hi(k) = max(live_xr);
end
w_deep_mask = (Zin(:) >= 560 & Zin(:) <= 950) & (Xin(:) >= 180 & Xin(:) <= 1420);
out_pairs = 0; tot_pairs = 0;
for k = 1:ns
    xm = 2*Xin(w_deep_mask) - xs(k);
    out_pairs = out_pairs + sum(xm < xr_lo(k) | xm > xr_hi(k));
    tot_pairs = tot_pairs + numel(xm);
end
out_frac = out_pairs / max(tot_pairs, 1);

if smoke, freqs = [10.0, 22.0, 34.0]; else, freqs = 4.0:1.0:37.0; end
nf = numel(freqs);
fref_phys = 25.0;  % In-band velocity convention v_band(25 Hz)!
fpeak = 22.0; sig = cfg.sig; s70 = gain_s(70);
ms = [0.5, 0.7, 1.0, 1.4, 2.0]; nm = numel(ms);

% Theorem 3 polarity signs: B=1 (per receiver-shot) and B=300 (one sign per shot)
rng(seed + 1, 'twister');   EPS_B1   = 2*(rand(nr, ns) > 0.5) - 1;
rng(seed + 300, 'twister'); sgn_B300 = 2*(rand(1, ns)  > 0.5) - 1;

fprintf('=== Tier-1 Model %s (seed=%d, out=%.1f%%, vscale=%.2f, corr=%.1f tr, snr=%.1f dB) ===\n', ...
    mname, seed, 100*out_frac, vscale_mig, corr_noise_sig, snr_db);

% Calibrate prestack noise at f=22 Hz
W22 = (22/fpeak)^2 * exp(-(22/fpeak)^2);
b22 = zeros(N, ns); b22(sub2ind([N ns], isrc, 1:ns)) = -W22 / P.dx^2;
[Lm,Um,Pm,Qm,Rm] = lu(helm_op(v_mig_exact, 22, P)); Ub = Qm*(Um\(Lm\(Pm*(Rm\b22)))); Dbg22 = Ub(irec,:); clear Ub Lm Um Pm Qm Rm
vc22 = v_true ./ (1 - (1./(pi*Q_true))*log(22/fref_phys) - 1i./(2*Q_true));
[Lv,Uv,Pv,Qv,Rv] = lu(helm_op(vc22, 22, P)); Uv22 = Qv*(Uv\(Lv\(Pv*(Rv\b22)))); Dcal = (Uv22(irec,:) - Dbg22) .* rec_mask; clear Uv22 Lv Uv Pv Qv Rv
sig_n = sqrt(mean(abs(Dcal(rec_mask)).^2)) * 10^(-snr_db/20);

if corr_noise_sig > 0
    rc = ceil(3 * corr_noise_sig);
    g_nc = exp(-0.5 * ((-rc:rc)' / corr_noise_sig).^2);
    g_nc = g_nc / sqrt(sum(g_nc.^2));
end

% Pass 1: Wavefield solves (store Cs_sig, Cs_nu, and exact B=1 flipped Cs_b1 in single)
Cs_sig_all = zeros(nin, ns, nf, 'single');
Cs_nu_all  = zeros(nin, ns, nf, 'single');
Cs_b1_all  = zeros(nin, ns, nf, 'single');
I_ref      = zeros(nin, 1);
I_R0_tot   = zeros(nin, 1);

for jf = 1:nf
    f = freqs(jf); w = 2*pi*f;
    W_src = (f/fpeak)^2 * exp(-(f/fpeak)^2);
    b_src = zeros(N, ns); b_src(sub2ind([N ns], isrc, 1:ns)) = -W_src / P.dx^2;

    [Lb,Ub_,Pb,Qb,Rb] = lu(helm_op(v_mig_exact, f, P));
    U_bg_ex = Qb*(Ub_\(Lb\(Pb*(Rb\b_src)))); D_bg = U_bg_ex(irec,:);
    if vscale_mig == 1.0
        Lm = Lb; Um = Ub_; Pm = Pb; Qm = Qb; Rm = Rb;
        Us_c = conj(U_bg_ex(lin,:)) * w^2; clear U_bg_ex
    else
        clear U_bg_ex Lb Ub_ Pb Qb Rb
        [Lm,Um,Pm,Qm,Rm] = lu(helm_op(v_mig, f, P));
        U_mig = Qm*(Um\(Lm\(Pm*(Rm\b_src))));
        Us_c = conj(U_mig(lin,:)) * w^2; clear U_mig
    end

    [Le,Ue,Pe,Qe,Re] = lu(helm_op(v_true, f, P));
    U_el = Qe*(Ue\(Le\(Pe*(Re\b_src)))); D_el = (U_el(irec,:) - D_bg) .* rec_mask; clear U_el Le Ue Pe Qe Re

    vc = v_true ./ (1 - (1./(pi*Q_true))*log(max(f,1e-6)/fref_phys) - 1i./(2*Q_true));
    [Lv,Uv,Pv,Qv,Rv] = lu(helm_op(vc, f, P));
    U_vi = Qv*(Uv\(Lv\(Pv*(Rv\b_src)))); D_vi = (U_vi(irec,:) - D_bg) .* rec_mask; clear U_vi Lv Uv Pv Qv Rv b_src

    N_f = sig_n * (1/sqrt(2)) * (randn(nr, ns) + 1i*randn(nr, ns));
    if corr_noise_sig > 0, N_f = conv2(g_nc, 1, N_f, 'same'); end
    N_f = N_f .* rec_mask;

    bp = @(D) extract_bp(D, Lm, Um, Pm, Qm, Rm, irec, lin, N, ns, P.dx);
    V_el = bp(D_el);
    V_vi = bp(D_vi);
    V_nu = bp(N_f);
    V_b1 = bp((D_vi + N_f) .* EPS_B1);
    clear Lm Um Pm Qm Rm

    I_ref = I_ref + real(sum(Us_c .* V_el, 2)); clear V_el
    cs = Us_c .* V_vi;
    cn = Us_c .* V_nu;
    cb = Us_c .* V_b1;
    Cs_sig_all(:,:,jf) = single(cs);
    Cs_nu_all(:,:,jf)  = single(cn);
    Cs_b1_all(:,:,jf)  = single(cb);
    I_R0_tot = I_R0_tot + real(sum(cs + cn, 2));
end

% Blind dip estimation from noisy R0 -> DSA mask & dip-consistent t*
R0_2d = core.prep_image(I_R0_tot, P.nzi, P.nxi, z_in, P.dx, cfg.bg_sigma_m, cfg.z_mute_m);
[phi_eff, c_st] = core.blind_dip(R0_2d, P.dx, cfg.dip_sigma_m, cfg.dip_cap_deg);
[w_dsa, TS_dsa] = core.build_dsa_in(Xin(:), Zin(:), 20.0, xs, xr_lo, xr_hi, phi_eff(:), Fts_in, cfg);
[w_sa,  ~     ] = core.build_dsa_in(Xin(:), Zin(:), 20.0, xs, xr_lo, xr_hi, zeros(nin,1), Fts_in, cfg);
dip_stats = struct('med_abs_phi_deg', median(abs(rad2deg(phi_eff(w_deep_mask)))), ...
                   'p90_abs_phi_deg', prctile(abs(rad2deg(phi_eff(w_deep_mask))), 90), ...
                   'med_coh', median(c_st(w_deep_mask)));

% Pass 1.5: Exact Theorem 3 Band-Summed (+-5 Hz, 160 m box) Noise Estimation:
% Compare RAW B=1, RAW B=300, DSA B=1, DSA B=300 against EXACT ORACLE iv_true!
I_raw_s = zeros(nin, nf, 'single'); I_raw_n = I_raw_s; I_raw_t = I_raw_s; I_raw_f1 = I_raw_s; I_raw_f3 = I_raw_s;
I_dsa_s = zeros(nin, nf, 'single'); I_dsa_n = I_dsa_s; I_dsa_t = I_dsa_s; I_dsa_f1 = I_dsa_s; I_dsa_f3 = I_dsa_s;
I_sa_t  = zeros(nin, nf, 'single'); I_sa_f3 = I_sa_t;
for jf = 1:nf
    cs = double(Cs_sig_all(:,:,jf)); cn = double(Cs_nu_all(:,:,jf)); ct = cs + cn; cb = double(Cs_b1_all(:,:,jf));
    I_raw_s(:,jf)  = single(real(sum(cs, 2)));
    I_raw_n(:,jf)  = single(real(sum(cn, 2)));
    I_raw_t(:,jf)  = single(real(sum(ct, 2)));
    I_raw_f1(:,jf) = single(real(sum(cb, 2)));
    I_raw_f3(:,jf) = single(real(sum(ct .* sgn_B300, 2)));

    I_dsa_s(:,jf)  = single(real(sum(w_dsa .* cs, 2)));
    I_dsa_n(:,jf)  = single(real(sum(w_dsa .* cn, 2)));
    I_dsa_t(:,jf)  = single(real(sum(w_dsa .* ct, 2)));
    I_dsa_f1(:,jf) = single(real(sum(w_dsa .* cb, 2)));
    I_dsa_f3(:,jf) = single(real(sum((w_dsa .* ct) .* sgn_B300, 2)));

    I_sa_t(:,jf)   = single(real(sum(w_sa  .* ct, 2)));
    I_sa_f3(:,jf)  = single(real(sum((w_sa .* ct) .* sgn_B300, 2)));
end
clear Cs_b1_all

iv_raw_b1   = core.band_inv_snr(I_raw_t, I_raw_f1, freqs, P.nzi, P.nxi, P.dx, cfg.snr_band_hz, cfg.snr_box_m, cfg.snr_cap);
iv_raw_b300 = core.band_inv_snr(I_raw_t, I_raw_f3, freqs, P.nzi, P.nxi, P.dx, cfg.snr_band_hz, cfg.snr_box_m, cfg.snr_cap);
iv_dsa_b1   = core.band_inv_snr(I_dsa_t, I_dsa_f1, freqs, P.nzi, P.nxi, P.dx, cfg.snr_band_hz, cfg.snr_box_m, cfg.snr_cap);
iv_dsa_b300 = core.band_inv_snr(I_dsa_t, I_dsa_f3, freqs, P.nzi, P.nxi, P.dx, cfg.snr_band_hz, cfg.snr_box_m, cfg.snr_cap);
iv_sa_b300  = core.band_inv_snr(I_sa_t,  I_sa_f3,  freqs, P.nzi, P.nxi, P.dx, cfg.snr_band_hz, cfg.snr_box_m, cfg.snr_cap);

% Exact oracle true iv (same +-5 Hz band sum & spatial box on pure signal vs pure noise!)
iv_true_raw = true_inv_snr(I_raw_s, I_raw_n, freqs, P.nzi, P.nxi, P.dx, cfg.snr_band_hz, cfg.snr_box_m, cfg.snr_cap);
iv_true_dsa = true_inv_snr(I_dsa_s, I_dsa_n, freqs, P.nzi, P.nxi, P.dx, cfg.snr_band_hz, cfg.snr_box_m, cfg.snr_cap);

f_acc_mask = (freqs >= 10.0 & freqs <= 35.0);
if ~any(f_acc_mask), f_acc_mask = true(size(freqs)); end
log2_err = @(est, tru) median(arrayfun(@(j) median(log2(max(double(est(w_deep_mask,j)),1e-6) ./ max(double(tru(w_deep_mask,j)),1e-6))), find(f_acc_mask)));
noise_acc = struct( ...
    'RAW_B1_log2err',   log2_err(iv_raw_b1,   iv_true_raw), ...
    'RAW_B300_log2err', log2_err(iv_raw_b300, iv_true_raw), ...
    'DSA_B1_log2err',   log2_err(iv_dsa_b1,   iv_true_dsa), ...
    'DSA_B300_log2err', log2_err(iv_dsa_b300, iv_true_dsa));
fprintf('  Theorem 3 Noise Acc log2(est/true) [10-35 Hz]: RAW B1=%+.2f, B300=%+.2f | DSA B1=%+.2f, B300=%+.2f\n', ...
    noise_acc.RAW_B1_log2err, noise_acc.RAW_B300_log2err, noise_acc.DSA_B1_log2err, noise_acc.DSA_B300_log2err);

% Pass 2: Full Factorial Evaluation (B=1 vs B=300, RAW vs SA_flat vs DSA, f_ref=25 vs 100 Hz, IQ70/W1b/PH/R0)
mnames = {'RAW_R0', 'DSA_R0', 'RAW_UQW1_25', 'RAW_UQW1b_25', 'RAW_UQW1b_100', ...
          'SA_UQW1b_25', 'DSA_PH_25', 'DSA_IQ70_25', 'DSA_W1b_25', 'DSA_UQW1_25', 'DSA_UQW1b_25'};
IS = struct(); IN = struct();
for k = 1:numel(mnames)
    nc = 1 + (nm-1)*~endsWith(mnames{k}, '_R0');
    IS.(mnames{k}) = zeros(nin, nc);
    IN.(mnames{k}) = zeros(nin, nc);
end

xmax = 1.05 * max(max(TS_flat(:)), max(TS_dsa(:))) / min(ms);
xg = linspace(0, xmax, 1201); dxg = xg(2);

for jf = 1:nf
    f = freqs(jf);
    Cs_s = double(Cs_sig_all(:,:,jf));
    Cs_n = double(Cs_nu_all(:,:,jf));
    T0_25  = uq_tab(f, xg, 0,   25.0);
    T1_25  = uq_tab(f, xg, sig, 25.0);
    T1_100 = uq_tab(f, xg, sig, 100.0);

    IS.RAW_R0 = IS.RAW_R0 + real(sum(Cs_s, 2));
    IN.RAW_R0 = IN.RAW_R0 + real(sum(Cs_n, 2));
    IS.DSA_R0 = IS.DSA_R0 + real(sum(w_dsa .* Cs_s, 2));
    IN.DSA_R0 = IN.DSA_R0 + real(sum(w_dsa .* Cs_n, 2));

    ivr1 = double(iv_raw_b1(:,jf));
    ivr3 = double(iv_raw_b300(:,jf));
    ivd1 = double(iv_dsa_b1(:,jf));
    ivd3 = double(iv_dsa_b300(:,jf));
    ivs3 = double(iv_sa_b300(:,jf));

    for im = 1:nm
        % RAW & SA_flat (using TS_flat)
        u = (TS_flat / ms(im)) / dxg; i0 = min(floor(u), numel(xg)-2); ww = u - i0; i0 = i0 + 1;
        it = @(tb) tb(i0).*(1-ww) + tb(i0+1).*ww;
        A1 = it(T1_25.A); C1 = it(T1_25.C); P1 = it(T1_25.P);
        H_r1 = (A1 + s70*P1.*C1 + P1.*ivr1) ./ (1 + s70*C1 + ivr1);
        H_r3 = (A1 + s70*P1.*C1 + P1.*ivr3) ./ (1 + s70*C1 + ivr3);
        H_s3 = (A1 + s70*P1.*C1 + P1.*ivs3) ./ (1 + s70*C1 + ivs3);
        IS.RAW_UQW1_25(:,im)  = IS.RAW_UQW1_25(:,im)  + real(sum(H_r1 .* Cs_s, 2));
        IN.RAW_UQW1_25(:,im)  = IN.RAW_UQW1_25(:,im)  + real(sum(H_r1 .* Cs_n, 2));
        IS.RAW_UQW1b_25(:,im) = IS.RAW_UQW1b_25(:,im) + real(sum(H_r3 .* Cs_s, 2));
        IN.RAW_UQW1b_25(:,im) = IN.RAW_UQW1b_25(:,im) + real(sum(H_r3 .* Cs_n, 2));
        IS.SA_UQW1b_25(:,im)  = IS.SA_UQW1b_25(:,im)  + real(sum((w_sa.*H_s3) .* Cs_s, 2));
        IN.SA_UQW1b_25(:,im)  = IN.SA_UQW1b_25(:,im)  + real(sum((w_sa.*H_s3) .* Cs_n, 2));

        A1h = it(T1_100.A); C1h = it(T1_100.C); P1h = it(T1_100.P);
        H_100 = (A1h + s70*P1h.*C1h + P1h.*ivr3) ./ (1 + s70*C1h + ivr3);
        IS.RAW_UQW1b_100(:,im) = IS.RAW_UQW1b_100(:,im) + real(sum(H_100 .* Cs_s, 2));
        IN.RAW_UQW1b_100(:,im) = IN.RAW_UQW1b_100(:,im) + real(sum(H_100 .* Cs_n, 2));

        % DSA (using dip-consistent TS_dsa)
        u = (TS_dsa / ms(im)) / dxg; i0 = min(floor(u), numel(xg)-2); ww = u - i0; i0 = i0 + 1;
        it = @(tb) tb(i0).*(1-ww) + tb(i0+1).*ww;
        A0 = it(T0_25.A); C0 = it(T0_25.C); P0 = it(T0_25.P);
        A1 = it(T1_25.A); C1 = it(T1_25.C); P1 = it(T1_25.P);
        H_ph  = P0;
        H_iq  = (A0 + s70*P0.*C0) ./ (1 + s70*C0);
        H_w1b = (A0 + s70*P0.*C0 + P0.*ivd3) ./ (1 + s70*C0 + ivd3);
        H_uq1 = (A1 + s70*P1.*C1 + P1.*ivd1) ./ (1 + s70*C1 + ivd1);
        H_uqb = (A1 + s70*P1.*C1 + P1.*ivd3) ./ (1 + s70*C1 + ivd3);

        IS.DSA_PH_25(:,im)    = IS.DSA_PH_25(:,im)    + real(sum((w_dsa.*H_ph)  .* Cs_s, 2));
        IN.DSA_PH_25(:,im)    = IN.DSA_PH_25(:,im)    + real(sum((w_dsa.*H_ph)  .* Cs_n, 2));
        IS.DSA_IQ70_25(:,im)  = IS.DSA_IQ70_25(:,im)  + real(sum((w_dsa.*H_iq)  .* Cs_s, 2));
        IN.DSA_IQ70_25(:,im)  = IN.DSA_IQ70_25(:,im)  + real(sum((w_dsa.*H_iq)  .* Cs_n, 2));
        IS.DSA_W1b_25(:,im)   = IS.DSA_W1b_25(:,im)   + real(sum((w_dsa.*H_w1b) .* Cs_s, 2));
        IN.DSA_W1b_25(:,im)   = IN.DSA_W1b_25(:,im)   + real(sum((w_dsa.*H_w1b) .* Cs_n, 2));
        IS.DSA_UQW1_25(:,im)  = IS.DSA_UQW1_25(:,im)  + real(sum((w_dsa.*H_uq1) .* Cs_s, 2));
        IN.DSA_UQW1_25(:,im)  = IN.DSA_UQW1_25(:,im)  + real(sum((w_dsa.*H_uq1) .* Cs_n, 2));
        IS.DSA_UQW1b_25(:,im) = IS.DSA_UQW1b_25(:,im) + real(sum((w_dsa.*H_uqb) .* Cs_s, 2));
        IN.DSA_UQW1b_25(:,im) = IN.DSA_UQW1b_25(:,im) + real(sum((w_dsa.*H_uqb) .* Cs_n, 2));
    end
end

% Score against lossless elastic reference I_ref
ref_2d = core.prep_image(I_ref, P.nzi, P.nxi, z_in, P.dx, cfg.bg_sigma_m, cfg.z_mute_m);
ref_hf = core.highpass_kz(ref_2d, P.dx, 10.0);
W_deep = {[180 1420], [560 950]};
rms_w = @(A, W) sqrt(mean(A(z_in>=W{2}(1)&z_in<=W{2}(2), x_in>=W{1}(1)&x_in<=W{1}(2)).^2, 'all'));
ref_rcf = log10(rms_w(ref_2d, W_targ) / max(rms_w(ref_2d, W_flank), 1e-30));

res = struct('model', mname, 'seed', seed, 'smoke', smoke, 'n_freq', nf, 'm', ms, ...
             'out_of_array_frac', out_frac, 'vscale_mig', vscale_mig, ...
             'corr_noise_sig_tr', corr_noise_sig, 'snr_db', snr_db, ...
             'dip_stats', dip_stats, 'noise_acc', noise_acc, 'elapsed_s', toc(t0_all));

for wn = {'W_deep', 'W_targ'}
    wname = wn{1}; W = W_deep; if strcmp(wname, 'W_targ'), W = W_targ; end
    iz_w = find(z_in >= W{2}(1) & z_in <= W{2}(2));
    ix_w = find(x_in >= W{1}(1) & x_in <= W{1}(2));
    if strcmp(wname, 'W_deep'), fprintf('  [%s W_deep] dip med=%.1f, p90=%.1f deg\n', mname, dip_stats.med_abs_phi_deg, dip_stats.p90_abs_phi_deg); end
    for k = 1:numel(mnames)
        key = mnames{k};
        Sm = IS.(key); Nm = IN.(key); nc = size(Sm, 2);
        e1 = zeros(1,nc); e1h = zeros(1,nc); eamp = zeros(1,nc); snr = zeros(1,nc);
        for c = 1:nc
            S2 = core.prep_image(Sm(:,c), P.nzi, P.nxi, z_in, P.dx, cfg.bg_sigma_m, cfg.z_mute_m);
            N2 = core.prep_image(Nm(:,c), P.nzi, P.nxi, z_in, P.dx, cfg.bg_sigma_m, cfg.z_mute_m);
            T2 = S2 + N2; Th = core.highpass_kz(T2, P.dx, 10.0);
            e1(c)   = core.score_e1(T2, ref_2d, iz_w, ix_w);
            e1h(c)  = core.score_e1(Th, ref_hf, iz_w, ix_w);
            snr(c)  = 10*log10(sum(S2(iz_w,ix_w).^2,'all') / max(sum(N2(iz_w,ix_w).^2,'all'), 1e-30));
            eamp(c) = abs(log10(rms_w(T2, W_targ)/max(rms_w(T2, W_flank),1e-30)) - ref_rcf);
        end
        res.(wname).(key) = struct('E1', e1, 'E1_mean', mean(e1), 'E1_HF', e1h, 'E1_HF_mean', mean(e1h), ...
                                   'E_amp', eamp, 'E_amp_mean', mean(eamp), 'SNR_img_dB', snr, 'SNR_img_dB_mean', mean(snr));
        if strcmp(wname, 'W_deep')
            fprintf('    %-15s | E1=%.4f (m0.5=%.4f,m1=%.4f,m2=%.4f) | HF=%.4f | Eamp=%.3f | SNR=%.2f dB\n', ...
                key, mean(e1), e1(1), e1(min(3,nc)), e1(end), mean(e1h), mean(eamp), mean(snr));
        end
    end
end
res.elapsed_s = toc(t0_all);
end

%% ========================================================================
%% U4: Exact Volve OBC Full-Pipeline Twin (Unknown wavelet + flat W_mig=1
%%     + (2*pi*f)^2 + one-sided streamer + l30 + v_band(25 Hz) + 5 wells)
%% ========================================================================
function res = run_u4_volve_twin(core, cfg, smoke)
t0_all = tic;
seed = 20261304; rng(seed, 'twister');
P.dx = 10.0; P.npml = 18; P.nzi = 100; P.nxi = 160;
P.nz = P.nzi + 2*P.npml; P.nx = P.nxi + 2*P.npml;
P.z = ((1:P.nz)' - P.npml - 1) * P.dx;
P.x = ((1:P.nx)  - P.npml - 1) * P.dx;
P.in = {P.npml + (1:P.nzi), P.npml + (1:P.nxi)};
P.izs = P.npml + 3;
P.ixs = P.npml + 1 + (8:3:152);
P.ixr = P.npml + 1 + (0:1:P.nxi-1);
ns = numel(P.ixs); nr = numel(P.ixr); N = P.nz * P.nx;
[Zi, Xi] = ndgrid(P.in{1}, P.in{2}); lin = sub2ind([P.nz P.nx], Zi(:), Xi(:)); nin = numel(lin);
isrc = sub2ind([P.nz P.nx], P.izs*ones(1,ns), P.ixs);
irec = sub2ind([P.nz P.nx], P.izs*ones(1,nr), P.ixr);
[Zfull, Xfull] = ndgrid(P.z, P.x);
z_in = P.z(P.in{1}); x_in = P.x(P.in{2});
[Zin, Xin] = ndgrid(z_in, x_in);
xs = P.x(P.ixs); xr = P.x(P.ixr);

% One-sided trailing streamer (~60% out-of-array)
rec_mask = false(nr, ns); xr_lo = zeros(1, ns); xr_hi = zeros(1, ns);
for s = 1:ns
    rec_mask(:, s) = (xr >= xs(s) - 240.0) & (xr <= xs(s) + 860.0);
    xr_lo(s) = min(xr(rec_mask(:, s))); xr_hi(s) = max(xr(rec_mask(:, s)));
end

% Finely layered earth with gentle 6-deg dip (well-log reflectivity style)
v_bg = 1800.0 + 1.45 * max(Zfull, 0);
z_dip = Zfull + 0.11 * (Xfull - 800.0);
dv_layer = zeros(size(Zfull));
z_edges = [175, 225, 285, 345, 415, 485, 555, 590, 630, 665, 705, 745, 785, 825, 860, 900, 945];
for il = 1:numel(z_edges)-1
    step = max(min(0.055 * randn(), 0.09), -0.09);
    dv_layer(z_dip >= z_edges(il) & z_dip < z_edges(il+1)) = step;
end
v_true = v_bg .* (1 + dv_layer);
v_mig  = core.smooth2d(v_true, 3.0);

Q_true = 55.0 * ones(P.nz, P.nx);
Q_true(Zfull < 80) = 1000.0;
Q_true(Zfull >= 520) = 35.0;

s_mig = 1.0 ./ v_mig; eps_q = 0.5;
F.T = zeros(P.nz, P.nx, ns, 'single'); F.ts = F.T;
for k = 1:ns
    F.T(:,:,k)  = eikonal_fs(s_mig, P.dx, P.izs, P.ixs(k));
    F.ts(:,:,k) = (eikonal_fs(s_mig.*(1+eps_q./Q_true), P.dx, P.izs, P.ixs(k)) - ...
                   eikonal_fs(s_mig.*(1-eps_q./Q_true), P.dx, P.izs, P.ixs(k))) / (2*eps_q);
end
TS_flat = zeros(nin, ns); Fts_in = zeros(nin, ns);
for k = 1:ns
    [~, tt] = two_way(F, k, P); TS_flat(:, k) = double(tt(lin));
    tk = F.ts(:,:,k); Fts_in(:, k) = double(tk(lin));
end
clear F

if smoke, freqs = [10.0, 22.0, 34.0]; else, freqs = 4.0:1.0:37.0; end
nf = numel(freqs);
fref = 25.0; fpeak_true = 20.0; sig = cfg.sig; s70 = gain_s(70);
ms = [0.5, 0.7, 1.0, 1.4, 2.0]; nm = numel(ms);
corr_noise_sig = 5.0; snr_db = -14.0;
rc = ceil(3 * corr_noise_sig);
g_nc = exp(-0.5 * ((-rc:rc)' / corr_noise_sig).^2); g_nc = g_nc / sqrt(sum(g_nc.^2));

rng(seed + 1, 'twister');   EPS_B1   = 2*(rand(nr, ns) > 0.5) - 1;
rng(seed + 300, 'twister'); sgn_B300 = 2*(rand(1, ns)  > 0.5) - 1;

% Flat boxcar migration source wavelet W_mig(f) (exact Volve fld_run.m convention!)
W_mig = ones(1, nf);
a1 = freqs < 6; W_mig(a1) = 0.5 * (1 - cos(pi*(freqs(a1)-4)/2));
a2 = freqs > 34; W_mig(a2) = 0.5 * (1 + cos(pi*(freqs(a2)-34)/3));

% Calibrate 3 dB l30 noise at f=20 Hz
W20 = (20/fpeak_true)^2 * exp(-(20/fpeak_true)^2);
b20 = zeros(N, ns); b20(sub2ind([N ns], isrc, 1:ns)) = -W20 / P.dx^2;
[Lm,Um,Pm,Qm,Rm] = lu(helm_op(v_mig, 20, P)); Ub = Qm*(Um\(Lm\(Pm*(Rm\b20)))); Dbg20 = Ub(irec,:); clear Ub Lm Um Pm Qm Rm
vc20 = v_true ./ (1 - (1./(pi*Q_true))*log(20/fref) - 1i./(2*Q_true));
[Lv,Uv,Pv,Qv,Rv] = lu(helm_op(vc20, 20, P)); Uv20 = Qv*(Uv\(Lv\(Pv*(Rv\b20)))); Dcal = (Uv20(irec,:) - Dbg20) .* rec_mask; clear Uv20 Lv Uv Pv Qv Rv
sig_field = sqrt(mean(abs(Dcal(rec_mask)).^2)) * 10^(-snr_db/20);

Cs_sig_all = zeros(nin, ns, nf, 'single');
Cs_nu_all  = zeros(nin, ns, nf, 'single');
Cs_b1_all  = zeros(nin, ns, nf, 'single');
I_RAW_R0   = zeros(nin, 1);

for jf = 1:nf
    f = freqs(jf); w = 2*pi*f;
    W_true = (f/fpeak_true)^2 * exp(-(f/fpeak_true)^2);
    b_true = zeros(N, ns); b_true(sub2ind([N ns], isrc, 1:ns)) = -W_true / P.dx^2;
    b_mig  = zeros(N, ns); b_mig(sub2ind([N ns], isrc, 1:ns))  = -W_mig(jf) / P.dx^2;

    [Lm,Um,Pm,Qm,Rm] = lu(helm_op(v_mig, f, P));
    U_bg = Qm*(Um\(Lm\(Pm*(Rm\b_true)))); D_bg = U_bg(irec,:); clear U_bg
    U_ms = Qm*(Um\(Lm\(Pm*(Rm\b_mig))));  Us_c = conj(U_ms(lin,:)) * (w^2); clear U_ms

    vc = v_true ./ (1 - (1./(pi*Q_true))*log(max(f,1e-6)/fref) - 1i./(2*Q_true));
    [Lv,Uv,Pv,Qv,Rv] = lu(helm_op(vc, f, P));
    U_vi = Qv*(Uv\(Lv\(Pv*(Rv\b_true)))); D_vi = (U_vi(irec,:) - D_bg) .* rec_mask; clear U_vi Lv Uv Pv Qv Rv

    N_f = sig_field * (1/sqrt(2)) * (randn(nr, ns) + 1i*randn(nr, ns));
    N_f = conv2(g_nc, 1, N_f, 'same') .* rec_mask;

    bp = @(D) extract_bp(D, Lm, Um, Pm, Qm, Rm, irec, lin, N, ns, P.dx);
    V_s = bp(D_vi); V_n = bp(N_f); V_b1 = bp((D_vi + N_f) .* EPS_B1);
    clear Lm Um Pm Qm Rm
    cs = Us_c .* V_s; cn = Us_c .* V_n; cb = Us_c .* V_b1;
    Cs_sig_all(:,:,jf) = single(cs);
    Cs_nu_all(:,:,jf)  = single(cn);
    Cs_b1_all(:,:,jf)  = single(cb);
    I_RAW_R0 = I_RAW_R0 + real(sum(cs + cn, 2));
end

% Blind dip -> DSA mask & dip-consistent t*
R0_2d = core.prep_image(I_RAW_R0, P.nzi, P.nxi, z_in, P.dx, cfg.bg_sigma_m, cfg.z_mute_m);
[phi_eff, ~] = core.blind_dip(R0_2d, P.dx, cfg.dip_sigma_m, cfg.dip_cap_deg);
[w_dsa, TS_dsa] = core.build_dsa_in(Xin(:), Zin(:), 20.0, xs, xr_lo, xr_hi, phi_eff(:), Fts_in, cfg);

% Band-summed B=1 & B=300 noise estimates (RAW & DSA) + oracle check
I_raw_s = zeros(nin, nf, 'single'); I_raw_n = I_raw_s; I_raw_t = I_raw_s; I_raw_f1 = I_raw_s; I_raw_f3 = I_raw_s;
I_dsa_s = zeros(nin, nf, 'single'); I_dsa_n = I_dsa_s; I_dsa_t = I_dsa_s; I_dsa_f1 = I_dsa_s; I_dsa_f3 = I_dsa_s;
for jf = 1:nf
    cs = double(Cs_sig_all(:,:,jf)); cn = double(Cs_nu_all(:,:,jf)); ct = cs + cn; cb = double(Cs_b1_all(:,:,jf));
    I_raw_s(:,jf) = single(real(sum(cs, 2))); I_raw_n(:,jf) = single(real(sum(cn, 2)));
    I_raw_t(:,jf) = single(real(sum(ct, 2))); I_raw_f1(:,jf) = single(real(sum(cb, 2)));
    I_raw_f3(:,jf) = single(real(sum(ct .* sgn_B300, 2)));

    I_dsa_s(:,jf) = single(real(sum(w_dsa .* cs, 2))); I_dsa_n(:,jf) = single(real(sum(w_dsa .* cn, 2)));
    I_dsa_t(:,jf) = single(real(sum(w_dsa .* ct, 2))); I_dsa_f1(:,jf) = single(real(sum(w_dsa .* cb, 2)));
    I_dsa_f3(:,jf) = single(real(sum((w_dsa .* ct) .* sgn_B300, 2)));
end
clear Cs_b1_all

iv_raw_b1   = core.band_inv_snr(I_raw_t, I_raw_f1, freqs, P.nzi, P.nxi, P.dx, cfg.snr_band_hz, cfg.snr_box_m, cfg.snr_cap);
iv_raw_b300 = core.band_inv_snr(I_raw_t, I_raw_f3, freqs, P.nzi, P.nxi, P.dx, cfg.snr_band_hz, cfg.snr_box_m, cfg.snr_cap);
iv_dsa_b1   = core.band_inv_snr(I_dsa_t, I_dsa_f1, freqs, P.nzi, P.nxi, P.dx, cfg.snr_band_hz, cfg.snr_box_m, cfg.snr_cap);
iv_dsa_b300 = core.band_inv_snr(I_dsa_t, I_dsa_f3, freqs, P.nzi, P.nxi, P.dx, cfg.snr_band_hz, cfg.snr_box_m, cfg.snr_cap);
iv_true_raw = true_inv_snr(I_raw_s, I_raw_n, freqs, P.nzi, P.nxi, P.dx, cfg.snr_band_hz, cfg.snr_box_m, cfg.snr_cap);
iv_true_dsa = true_inv_snr(I_dsa_s, I_dsa_n, freqs, P.nzi, P.nxi, P.dx, cfg.snr_band_hz, cfg.snr_box_m, cfg.snr_cap);

w_res_mask = (Zin(:) >= 540 & Zin(:) <= 940) & (Xin(:) >= 240 & Xin(:) <= 1360);
f_acc = (freqs >= 10.0 & freqs <= 32.0); if ~any(f_acc), f_acc = true(size(freqs)); end
log2_err = @(est, tru) median(arrayfun(@(j) median(log2(max(double(est(w_res_mask,j)),1e-6) ./ max(double(tru(w_res_mask,j)),1e-6))), find(f_acc)));
noise_acc = struct('RAW_B1_log2err', log2_err(iv_raw_b1, iv_true_raw), ...
                   'RAW_B300_log2err', log2_err(iv_raw_b300, iv_true_raw), ...
                   'DSA_B1_log2err', log2_err(iv_dsa_b1, iv_true_dsa), ...
                   'DSA_B300_log2err', log2_err(iv_dsa_b300, iv_true_dsa));

% Blind Shallow Anchor [190, 320] m -> PSAE & broad EnvAnc & relative dt*
anc_mask = (Zin(:) >= cfg.anc_z_m(1) & Zin(:) <= cfg.anc_z_m(2)) & (Xin(:) >= cfg.anc_x_m(1) & Xin(:) <= cfg.anc_x_m(2));
S_anc_raw = arrayfun(@(j) rms(double(I_raw_t(anc_mask, j))), 1:nf);
S_anc_dsa = arrayfun(@(j) rms(double(I_dsa_t(anc_mask, j))), 1:nf);
f_mid = 20.0;
[W_eq_raw, beta_anc_raw, ~] = core.fit_psae(freqs, S_anc_raw, [6.0, 32.0], f_mid);
[W_eq_dsa, beta_anc_dsa, ~] = core.fit_psae(freqs, S_anc_dsa, [6.0, 32.0], f_mid);
W_env_dsa = core.fit_env_anc(freqs, S_anc_dsa, W_mig, 6.0, f_mid);

tstar_anc_dsa  = median(TS_dsa(anc_mask, :), 'all');
tstar_anc_flat = median(TS_flat(anc_mask, :), 'all');
TS_dsa_rel  = max(TS_dsa  - tstar_anc_dsa,  0.0);
TS_flat_rel = max(TS_flat - tstar_anc_flat, 0.0);

Img = struct('RAW_R0', zeros(nin,1), 'DSA_R0', zeros(nin,1), 'DSA_R0_Eq', zeros(nin,1), ...
             'RAW_UQW1_25', zeros(nin,nm), 'RAW_UQW1b_25', zeros(nin,nm), ...
             'DSA_UQW1b_Unanc_25', zeros(nin,nm), 'RAW_UQW1b_PAnc_25', zeros(nin,nm), ...
             'DSA_IQ70_PAnc_25', zeros(nin,nm), 'DSA_W1b_PAnc_25', zeros(nin,nm), ...
             'DSA_UQW1_PAnc_25', zeros(nin,nm), 'DSA_UQW1b_PAnc_25', zeros(nin,nm), ...
             'DSA_UQW1b_EnvAnc_25', zeros(nin,nm));

xmax = 1.05 * max(max(TS_flat(:)), max(TS_dsa(:))) / min(ms);
xg = linspace(0, xmax, 1201); dxg = xg(2);
spec_res_eq = zeros(1, nf); spec_res_un = zeros(1, nf); spec_res_hero = zeros(1, nf);

for jf = 1:nf
    f = freqs(jf);
    Cs_t = double(Cs_sig_all(:,:,jf)) + double(Cs_nu_all(:,:,jf));
    wr = W_eq_raw(jf); wd = W_eq_dsa(jf); we = W_env_dsa(jf);

    raw_f = real(sum(Cs_t, 2)); dsa_f = real(sum(w_dsa .* Cs_t, 2));
    Img.RAW_R0    = Img.RAW_R0    + raw_f;
    Img.DSA_R0    = Img.DSA_R0    + dsa_f;
    Img.DSA_R0_Eq = Img.DSA_R0_Eq + wd * dsa_f;
    spec_res_eq(jf) = rms(wd * dsa_f(w_res_mask));

    T0 = uq_tab(f, xg, 0,   fref);
    T1 = uq_tab(f, xg, sig, fref);
    ivr1 = double(iv_raw_b1(:,jf)); ivr3 = double(iv_raw_b300(:,jf));
    ivd1 = double(iv_dsa_b1(:,jf)); ivd3 = double(iv_dsa_b300(:,jf));

    for im = 1:nm
        % Unanchored RAW (TS_flat)
        u = (TS_flat / ms(im)) / dxg; i0 = min(floor(u), numel(xg)-2); ww = u - i0; i0 = i0 + 1;
        it = @(tb) tb(i0).*(1-ww) + tb(i0+1).*ww;
        A1 = it(T1.A); C1 = it(T1.C); P1 = it(T1.P);
        f_r1 = real(sum(((A1 + s70*P1.*C1 + P1.*ivr1)./(1 + s70*C1 + ivr1)) .* Cs_t, 2));
        f_r3 = real(sum(((A1 + s70*P1.*C1 + P1.*ivr3)./(1 + s70*C1 + ivr3)) .* Cs_t, 2));
        Img.RAW_UQW1_25(:,im)  = Img.RAW_UQW1_25(:,im)  + f_r1;
        Img.RAW_UQW1b_25(:,im) = Img.RAW_UQW1b_25(:,im) + f_r3;
        if im == 3, spec_res_un(jf) = rms(f_r3(w_res_mask)); end

        % Unanchored DSA (TS_dsa)
        u = (TS_dsa / ms(im)) / dxg; i0 = min(floor(u), numel(xg)-2); ww = u - i0; i0 = i0 + 1;
        it = @(tb) tb(i0).*(1-ww) + tb(i0+1).*ww;
        A1 = it(T1.A); C1 = it(T1.C); P1 = it(T1.P);
        Img.DSA_UQW1b_Unanc_25(:,im) = Img.DSA_UQW1b_Unanc_25(:,im) + real(sum((w_dsa .* ((A1 + s70*P1.*C1 + P1.*ivd3)./(1 + s70*C1 + ivd3))) .* Cs_t, 2));

        % Anchored RAW (TS_flat_rel)
        u = (TS_flat_rel / ms(im)) / dxg; i0 = min(floor(u), numel(xg)-2); ww = u - i0; i0 = i0 + 1;
        it = @(tb) tb(i0).*(1-ww) + tb(i0+1).*ww;
        A1 = it(T1.A); C1 = it(T1.C); P1 = it(T1.P);
        Img.RAW_UQW1b_PAnc_25(:,im) = Img.RAW_UQW1b_PAnc_25(:,im) + wr * real(sum(((A1 + s70*P1.*C1 + P1.*ivr3)./(1 + s70*C1 + ivr3)) .* Cs_t, 2));

        % Anchored DSA (TS_dsa_rel)
        u = (TS_dsa_rel / ms(im)) / dxg; i0 = min(floor(u), numel(xg)-2); ww = u - i0; i0 = i0 + 1;
        it = @(tb) tb(i0).*(1-ww) + tb(i0+1).*ww;
        A0 = it(T0.A); C0 = it(T0.C); P0 = it(T0.P);
        A1 = it(T1.A); C1 = it(T1.C); P1 = it(T1.P);
        f_iq  = real(sum((w_dsa .* ((A0 + s70*P0.*C0)./(1 + s70*C0))) .* Cs_t, 2));
        f_w1b = real(sum((w_dsa .* ((A0 + s70*P0.*C0 + P0.*ivd3)./(1 + s70*C0 + ivd3))) .* Cs_t, 2));
        f_uq1 = real(sum((w_dsa .* ((A1 + s70*P1.*C1 + P1.*ivd1)./(1 + s70*C1 + ivd1))) .* Cs_t, 2));
        f_uqb = real(sum((w_dsa .* ((A1 + s70*P1.*C1 + P1.*ivd3)./(1 + s70*C1 + ivd3))) .* Cs_t, 2));

        Img.DSA_IQ70_PAnc_25(:,im)    = Img.DSA_IQ70_PAnc_25(:,im)    + wd * f_iq;
        Img.DSA_W1b_PAnc_25(:,im)     = Img.DSA_W1b_PAnc_25(:,im)     + wd * f_w1b;
        Img.DSA_UQW1_PAnc_25(:,im)    = Img.DSA_UQW1_PAnc_25(:,im)    + wd * f_uq1;
        Img.DSA_UQW1b_PAnc_25(:,im)   = Img.DSA_UQW1b_PAnc_25(:,im)   + wd * f_uqb;
        Img.DSA_UQW1b_EnvAnc_25(:,im) = Img.DSA_UQW1b_EnvAnc_25(:,im) + we * f_uqb;
        if im == 3, spec_res_hero(jf) = rms(wd * f_uqb(w_res_mask)); end
    end
end

% Evaluate exact 5-well flat-bandpass [4-18 cyc/km] well tie
well_xs = [400, 600, 800, 1000, 1200];
iz_res = find(z_in >= 540 & z_in <= 940);
v_in = v_true(P.in{1}, P.in{2});
rw_map = zeros(size(v_in)); rw_map(2:end, :) = diff(log(v_in), 1, 1);
lo_k = 4.0; hi_k = 18.0;
eval_w = @(Ivec) score_wells5(core, core.prep_image(Ivec, P.nzi, P.nxi, z_in, P.dx, cfg.bg_sigma_m, cfg.z_mute_m), rw_map, x_in, iz_res, well_xs, P.dx, lo_k, hi_k);

fm = (freqs >= 6.0 & freqs <= 32.0);
sl = @(s) feval(@(p) p(1), polyfit(freqs(fm), log(max(s(fm), 1e-30)), 1));
res = struct('model', 'U4', 'seed', seed, 'smoke', smoke, 'm', ms, ...
             'tstar_anc_dsa_ms', tstar_anc_dsa*1e3, 'beta_anc_raw', beta_anc_raw, 'beta_anc_dsa', beta_anc_dsa, ...
             'noise_acc', noise_acc, 'band_cyckm', [lo_k hi_k], ...
             'res_slope', struct('DSA_R0_Eq', sl(spec_res_eq), 'RAW_UQW1b_Unanc_m1', sl(spec_res_un), 'DSA_UQW1b_PAnc_m1', sl(spec_res_hero)), ...
             'methods', struct(), 'elapsed_s', 0);

fprintf('=== Tier-1 U4 Volve-Exact Twin (l30 + one-sided + flat W_mig + v_band=25 Hz) ===\n');
fprintf('  Noise Acc log2(est/true): RAW B1=%+.2f, B300=%+.2f | DSA B1=%+.2f, B300=%+.2f\n', ...
    noise_acc.RAW_B1_log2err, noise_acc.RAW_B300_log2err, noise_acc.DSA_B1_log2err, noise_acc.DSA_B300_log2err);
mkeys = fieldnames(Img)';
for ik = 1:numel(mkeys)
    k = mkeys{ik}; M = Img.(k); nc = size(M, 2);
    ts_ = arrayfun(@(c) eval_w(M(:, c)), 1:nc);
    res.methods.(k) = struct('tie', ts_, 'mean', mean(ts_));
    if nc == 1
        fprintf('  %-20s : %.4f\n', k, ts_(1));
    else
        fprintf('  %-20s : mean=%.4f (m0.5=%.4f, m1.0=%.4f, m2.0=%.4f)\n', k, mean(ts_), ts_(1), ts_(3), ts_(5));
    end
end
res.elapsed_s = toc(t0_all);
end

function c_mean = score_wells5(core, A, rw_map, x_in, iz_res, well_xs, dx, lo_k, hi_k)
cs = zeros(1, numel(well_xs));
for iw = 1:numel(well_xs)
    [~, jx] = min(abs(x_in - well_xs(iw)));
    tr = mean(A(iz_res, max(jx-1,1):min(jx+1,numel(x_in))), 2);
    rw = rw_map(iz_res, jx);
    tr_b = core.bp_z(tr, dx, lo_k, hi_k);
    rw_b = core.bp_z(rw, dx, lo_k, hi_k);
    best = -Inf;
    for L = -3:3
        best = max(best, core.cc(circshift(tr_b, L), rw_b));
    end
    cs(iw) = best;
end
c_mean = mean(cs);
end

function iv = true_inv_snr(Is_mat, In_mat, freqs, nzi, nxi, dx, band_hz, box_m, cap)
[nin, nf] = size(Is_mat); bw = max(1, round(box_m / dx));
iv = zeros(nin, nf, 'single');
for j = 1:nf
    bnd = abs(freqs - freqs(j)) <= band_hz + 1e-9;
    ps = reshape(sum(double(Is_mat(:,bnd)), 2), nzi, nxi).^2; ps = movmean(movmean(ps, bw, 1), bw, 2);
    pn = reshape(sum(double(In_mat(:,bnd)), 2), nzi, nxi).^2; pn = movmean(movmean(pn, bw, 1), bw, 2);
    iv(:, j) = single(min(pn(:) ./ max(ps(:), realmin('single')), cap));
end
end

function V = extract_bp(Dmat, L, U, Pp, Qq, R, irec, lin, N, ns, dx)
b = zeros(N, ns); b(irec, :) = Dmat;
Vf = -(R \ (Pp' * (L' \ (U' \ (Qq' * b))))) / dx^2;
V = Vf(lin, :);
end
