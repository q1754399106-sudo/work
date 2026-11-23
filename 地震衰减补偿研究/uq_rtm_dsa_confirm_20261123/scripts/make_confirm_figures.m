function make_confirm_figures(mode)
% Publication-grade 2D RTM imaging comparison figures for uq_rtm_dsa_confirm_20261123.
% Never modifies pre-registered scripts/results/hashes.
% Reconstructs bit-exact noisy (3 dB) migrated images for H1-H4 (m=0.5, 1.0) and H5_VolveAnc,
% verifies exact metric agreement against results/*.json, caches lightweight 2D image .mat files
% outside git (scratch dir) for fast restyling, and exports 5 PNG comparison figures to figures/.
%
% Usage:
%   make_confirm_figures()           % run everything (compute if not cached + plot all 5 figures)
%   make_confirm_figures('plotonly') % re-plot all 5 figures from cached scratch .mat in ~5 s

if nargin < 1 || isempty(mode), mode = 'all'; end
maxNumCompThreads(4);

here = fileparts(mfilename('fullpath'));
root = fileparts(here);
fig_dir = fullfile(root, 'figures');
res_dir = fullfile(root, 'results');
if ~isfolder(fig_dir), mkdir(fig_dir); end

scratch_dir = 'C:\Users\qq175\.gemini\antigravity\brain\c0dc2750-de05-424f-adb4-e7ef332a1950\scratch\confirm_imgs';
if ~isfolder(scratch_dir), mkdir(scratch_dir); end

addpath('D:\ChatGPT\研究\地震衰减补偿研究\uq_rtm_20261117\scripts');
set(0, 'DefaultAxesFontName', 'Microsoft YaHei', ...
       'DefaultTextFontName', 'Microsoft YaHei', ...
       'DefaultAxesFontSize', 9.5);

% 1. Reconstruct or load H1-H4 and H5_VolveAnc exact images
if ~strcmpi(mode, 'plotonly')
    for mcell = {'H2', 'H3', 'H1', 'H4'}
        mname = mcell{1};
        cfile = fullfile(scratch_dir, [mname '_imgs.mat']);
        if isfile(cfile)
            fprintf('Cache hit: %s\n', cfile);
        else
            D = compute_h_model_images(mname, res_dir);
            save(cfile, '-struct', 'D', '-v7.3');
            fprintf('Saved scratch cache: %s\n', cfile);
        end
    end
    cfile5 = fullfile(scratch_dir, 'H5_VolveAnc_imgs.mat');
    if isfile(cfile5)
        fprintf('Cache hit: %s\n', cfile5);
    else
        D5 = compute_h5_images(res_dir);
        save(cfile5, '-struct', 'D5', '-v7.3');
        fprintf('Saved scratch cache: %s\n', cfile5);
    end
end

% 2. Render all 5 figures
fprintf('\n--- Rendering Figures to %s ---\n', fig_dir);
D2 = load(fullfile(scratch_dir, 'H2_imgs.mat'));
plot_fig1_H2(D2, fullfile(fig_dir, 'fig1_H2_steepdip_compare.png'));
fprintf('Exported fig1_H2_steepdip_compare.png\n');

D3 = load(fullfile(scratch_dir, 'H3_imgs.mat'));
plot_fig2_H3(D3, fullfile(fig_dir, 'fig2_H3_streamer_chimney_compare.png'));
fprintf('Exported fig2_H3_streamer_chimney_compare.png\n');

D1 = load(fullfile(scratch_dir, 'H1_imgs.mat'));
D4 = load(fullfile(scratch_dir, 'H4_imgs.mat'));
plot_fig3_H1_H4(D1, D4, fullfile(fig_dir, 'fig3_H1_H4_syncline_l30_compare.png'));
fprintf('Exported fig3_H1_H4_syncline_l30_compare.png\n');

D5 = load(fullfile(scratch_dir, 'H5_VolveAnc_imgs.mat'));
plot_fig4_H5(D5, fullfile(fig_dir, 'fig4_H5_volve_anc_welltie.png'));
fprintf('Exported fig4_H5_volve_anc_welltie.png\n');

plot_fig5_volve_real(fullfile(fig_dir, 'fig5_volve_real_stage4_compare.png'), scratch_dir);
fprintf('Exported fig5_volve_real_stage4_compare.png\n');
fprintf('ALL 5 FIGURES COMPLETED.\n');
end

%% ========================================================================
%% EXACT DETERMINISTIC RECONSTRUCTION FOR H1-H4 (Matches dsa_confirm_bench.m)
%% ========================================================================
function D = compute_h_model_images(mname, res_dir)
t0 = tic;
P.dx = 10.0; P.npml = 18; P.nzi = 100; P.nxi = 160;
P.nz = P.nzi + 2*P.npml; P.nx = P.nxi + 2*P.npml;
P.z = ((1:P.nz)' - P.npml - 1) * P.dx;
P.x = ((1:P.nx)  - P.npml - 1) * P.dx;
P.in = {P.npml + (1:P.nzi), P.npml + (1:P.nxi)};
P.izs = P.npml + 3;
P.ixs = P.npml + 1 + (8:6:152);
P.ixr = P.npml + 1 + (0:2:P.nxi-1);
ns = numel(P.ixs); nr = numel(P.ixr); N = P.nz * P.nx;
[Zi, Xi] = ndgrid(P.in{1}, P.in{2}); lin = sub2ind([P.nz P.nx], Zi(:), Xi(:)); nin = numel(lin);
isrc = sub2ind([P.nz P.nx], P.izs*ones(1,ns), P.ixs);
irec = sub2ind([P.nz P.nx], P.izs*ones(1,nr), P.ixr);
[Zfull, Xfull] = ndgrid(P.z, P.x);
z_in = P.z(P.in{1}); x_in = P.x(P.in{2});
[Zin, Xin] = ndgrid(z_in, x_in);
xs = P.x(P.ixs); xr = P.x(P.ixr);

rec_mask = true(nr, ns); vscale_mig = 1.0; corr_noise_sig = 0.0;
switch mname
    case 'H1'
        seed = 20261201;
        z1 = 230.0 + 15.0 * exp(-((Xfull - 780.0)/350.0).^2);
        z2 = 620.0 + 40.0 * exp(-((Xfull - 820.0)/320.0).^2);
        z3 = 755.0 + 30.0 * exp(-((Xfull - 820.0)/340.0).^2);
        z4 = 885.0 + 15.0 * exp(-((Xfull - 820.0)/380.0).^2);
        v_true = 1840.0 * ones(P.nz, P.nx);
        v_true(Zfull >= z1) = 2140.0; v_true(Zfull >= z2) = 2540.0;
        v_true(Zfull >= z3) = 2900.0; v_true(Zfull >= z4) = 3300.0;
        invQ = (1/90)*ones(P.nz, P.nx); invQ(Zfull < 80) = 1/1000;
        L1 = exp(-(((Xfull - 560)/180).^2 + ((Zfull - 410)/105).^2).^1.5);
        L2 = exp(-(((Xfull - 1060)/170).^2 + ((Zfull - 460)/110).^2).^1.5);
        invQ = invQ + (1/14 - 1/90)*L1 + (1/16 - 1/90)*L2;
        Q_true = 1./invQ;
    case 'H2'
        seed = 20261202;
        Zg = Zfull;
        hw = (Xfull < 840.0 - (Zfull - 200.0) / tand(64.0)) & (Zfull > 200.0);
        Zg(hw) = Zg(hw) + 30.0;
        z1 = 235.0 - 0.06 * (Xfull - 800.0);
        z2 = 660.0 + 0.40 * (Xfull - 800.0);
        z3 = 790.0 - 0.29 * (Xfull - 800.0);
        z4 = 900.0 - 0.12 * (Xfull - 800.0);
        v_true = 1860.0 * ones(P.nz, P.nx);
        v_true(Zg >= z1) = 2160.0; v_true(Zg >= z2) = 2620.0;
        v_true(Zg >= z3) = 2980.0; v_true(Zg >= z4) = 3350.0;
        invQ = (1/85)*ones(P.nz, P.nx); invQ(Zfull < 80) = 1/1000;
        wedge = exp(-(((Xfull - 980)/260).^2 + ((Zfull - 360)/95).^2).^1.5);
        invQ = invQ + (1/16 - 1/85)*wedge;
        Q_true = 1./invQ;
    case 'H3'
        seed = 20261203;
        for s = 1:ns, rec_mask(:, s) = (xr >= xs(s) - 860.0) & (xr <= xs(s) + 240.0); end
        z1 = 240.0 - 0.05 * (Xfull - 800.0);
        z2 = 640.0 - 0.21 * (Xfull - 800.0);
        z3 = 770.0 - 0.18 * (Xfull - 800.0);
        z4 = 890.0 + 0.04 * (Xfull - 800.0);
        v_true = 1830.0 * ones(P.nz, P.nx);
        v_true(Zfull >= z1) = 2130.0; v_true(Zfull >= z2) = 2550.0;
        v_true(Zfull >= z3) = 2900.0; v_true(Zfull >= z4) = 3260.0;
        invQ = (1/90)*ones(P.nz, P.nx); invQ(Zfull < 80) = 1/1000;
        chim = exp(-(((Xfull - 760)/195).^2 + ((Zfull - 420)/150).^2).^1.5);
        invQ = invQ + (1/15 - 1/90)*chim;
        Q_true = 1./invQ;
    case 'H4'
        seed = 20261204;
        vscale_mig = 0.97; corr_noise_sig = 3.0;
        z1 = 235.0 + 12.0 * sin(2*pi*Xfull/900.0);
        z2 = 635.0 + 22.0 * sin(2*pi*(Xfull - 200)/850.0);
        z3 = 765.0 + 18.0 * cos(2*pi*(Xfull - 150)/950.0);
        z4 = 885.0 - 12.0 * sin(2*pi*Xfull/1100.0);
        v_true = 1850.0 * ones(P.nz, P.nx);
        v_true(Zfull >= z1) = 2150.0; v_true(Zfull >= z2) = 2560.0;
        v_true(Zfull >= z3) = 2920.0; v_true(Zfull >= z4) = 3300.0;
        Q_true = 110.0 - 82.0 * min(max(Xfull / 1600.0, 0), 1);
        Q_true(Zfull < 80) = 1000.0;
end

rng(seed, 'twister');
v_mig_exact = smooth2d(v_true, 2.5);
v_mig = vscale_mig * v_mig_exact;

s_mig = 1.0 ./ v_mig; eps_q = 0.5;
F.T = zeros(P.nz, P.nx, ns, 'single'); F.ts = F.T;
for k = 1:ns
    F.T(:,:,k) = eikonal_fs(s_mig, P.dx, P.izs, P.ixs(k));
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
w_sa_flat = build_dsa_mask(Xin(:), Zin(:), xs, xr_lo, xr_hi, zeros(nin, 1));

freqs = 4.0:1.0:37.0; nf = numel(freqs);
fref = 25.0; fpeak = 22.0; sig = 0.447; s70 = gain_s(70);
ms = [0.5, 1.0]; nm = numel(ms);
flip_sgn = 2 * (rand(1, ns) > 0.5) - 1;

W22 = (22/fpeak)^2 * exp(-(22/fpeak)^2);
b22 = zeros(N, ns); b22(sub2ind([N ns], isrc, 1:ns)) = -W22 / P.dx^2;
[Lm,Um,Pm,Qm,Rm] = lu(helm_op(v_mig_exact, 22, P)); Ub = Qm*(Um\(Lm\(Pm*(Rm\b22)))); Dbg22 = Ub(irec,:); clear Ub Lm Um Pm Qm Rm
vc22 = v_true ./ (1 - (1./(pi*Q_true))*log(22/fref) - 1i./(2*Q_true));
[Lv,Uv,Pv,Qv,Rv] = lu(helm_op(vc22, 22, P)); Uv22 = Qv*(Uv\(Lv\(Pv*(Rv\b22)))); Dcal = (Uv22(irec,:) - Dbg22) .* rec_mask; clear Uv22 Lv Uv Pv Qv Rv
d_rms = sqrt(mean(abs(Dcal(rec_mask)).^2));
sig_noisy = d_rms * 10^(-3/20);

if corr_noise_sig > 0
    rc = ceil(3 * corr_noise_sig);
    g_nc = exp(-0.5 * ((-rc:rc)' / corr_noise_sig).^2);
    g_nc = g_nc / sqrt(sum(g_nc.^2));
end

fprintf('Computing exact noisy images for %s...\n', mname);
Cs_sig_all = zeros(nin, ns, nf, 'single');
Cs_nu_all  = zeros(nin, ns, nf, 'single');
I_ref = zeros(nin, 1); I_R0_sig = zeros(nin, 1); I_R0_nu = zeros(nin, 1);

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

    vc = v_true ./ (1 - (1./(pi*Q_true))*log(max(f,1e-6)/fref) - 1i./(2*Q_true));
    [Lv,Uv,Pv,Qv,Rv] = lu(helm_op(vc, f, P));
    U_vi = Qv*(Uv\(Lv\(Pv*(Rv\b_src)))); D_vi = (U_vi(irec,:) - D_bg) .* rec_mask; clear U_vi Lv Uv Pv Qv Rv b_src

    N_unit = (1/sqrt(2)) * (randn(nr, ns) + 1i*randn(nr, ns));
    if corr_noise_sig > 0, N_unit = conv2(g_nc, 1, N_unit, 'same'); end
    N_unit = N_unit .* rec_mask;

    bp = @(D) extract_bp(D, Lm, Um, Pm, Qm, Rm, irec, lin, N, ns, P.dx);
    V_el = bp(D_el); V_vi = bp(D_vi); V_nu = bp(N_unit); clear Lm Um Pm Qm Rm

    I_ref = I_ref + real(sum(Us_c .* V_el, 2)); clear V_el
    cs = Us_c .* V_vi; cn = Us_c .* V_nu;
    Cs_sig_all(:,:,jf) = single(cs);
    Cs_nu_all(:,:,jf)  = single(cn);
    I_R0_sig = I_R0_sig + real(sum(cs, 2));
    I_R0_nu  = I_R0_nu  + real(sum(cn, 2));
    if mod(jf, 17) == 0, fprintf('  %s Pass-1 %d/%d (%.1fs)\n', mname, jf, nf, toc(t0)); end
end

% Pass 2: noisy regime, ms = [0.5, 1.0] (identical xmax=min(ms)=0.5 -> bit-exact!)
R0_tot_2d = prep_img(I_R0_sig + sig_noisy * I_R0_nu, P.nzi, P.nxi, z_in);
[phi_eff, c_st] = blind_structure_dip(R0_tot_2d, P.dx);
[w_dsa, TS_dsa] = build_dsa_mask_and_tstar(Xin(:), Zin(:), xs, xr_lo, xr_hi, phi_eff(:), Fts_in);
W_list  = {ones(nin, ns), w_sa_flat, w_dsa};
TS_list = {TS_flat,       TS_flat,   TS_dsa};
modes = {'RAW', 'SA_flat', 'DSA'};
ops   = {'R0', 'IQ70', 'W1', 'UQW1'};

IS = struct(); IN = struct();
for imd = 1:3
    for io = 1:4
        key = [modes{imd} '_' ops{io}];
        nc = 1 + (nm-1)*~strcmp(ops{io}, 'R0');
        IS.(key) = zeros(nin, nc); IN.(key) = zeros(nin, nc);
    end
end

xmax = 1.05 * max(max(TS_flat(:)), max(TS_dsa(:))) / 0.5;
xg = linspace(0, xmax, 1201); dxg = xg(2);

for jf = 1:nf
    f = freqs(jf);
    Cs_s = double(Cs_sig_all(:,:,jf));
    Cs_n = sig_noisy * double(Cs_nu_all(:,:,jf));
    Cs_t = Cs_s + Cs_n;
    T0 = uq_tab(f, xg, 0, fref);
    T1 = uq_tab(f, xg, sig, fref);

    for imd = 1:3
        Wm = W_list{imd}; TSm = TS_list{imd}; mstr = modes{imd};
        IS.([mstr '_R0']) = IS.([mstr '_R0']) + real(sum(Wm .* Cs_s, 2));
        IN.([mstr '_R0']) = IN.([mstr '_R0']) + real(sum(Wm .* Cs_n, 2));

        stk_t = abs(sum(Wm .* Cs_t, 2)).^2;
        stk_f = abs(sum(Wm .* Cs_t .* flip_sgn, 2)).^2;
        pt = smooth2d(reshape(stk_t, P.nzi, P.nxi), 5.0);
        pf = smooth2d(reshape(stk_f, P.nzi, P.nxi), 5.0);
        iv = min(max(pf(:) ./ max(pt(:) - pf(:), 0.05*pf(:)), 1e-3), 1e3);

        for im = 1:nm
            u = (TSm / ms(im)) / dxg; i0 = min(floor(u), numel(xg)-2); ww = u - i0; i0 = i0 + 1;
            it = @(tb) tb(i0).*(1-ww) + tb(i0+1).*ww;
            A0 = it(T0.A); C0 = it(T0.C); P0 = it(T0.P);
            A1 = it(T1.A); C1 = it(T1.C); P1 = it(T1.P);
            H_iq = (A0 + s70*P0.*C0) ./ (1 + s70*C0);
            H_w1 = (A0 + s70*P0.*C0 + P0.*iv) ./ (1 + s70*C0 + iv);
            H_uq = (A1 + s70*P1.*C1 + P1.*iv) ./ (1 + s70*C1 + iv);
            IS.([mstr '_IQ70'])(:,im) = IS.([mstr '_IQ70'])(:,im) + real(sum((Wm.*H_iq).*Cs_s, 2));
            IN.([mstr '_IQ70'])(:,im) = IN.([mstr '_IQ70'])(:,im) + real(sum((Wm.*H_iq).*Cs_n, 2));
            IS.([mstr '_W1'])(:,im)   = IS.([mstr '_W1'])(:,im)   + real(sum((Wm.*H_w1).*Cs_s, 2));
            IN.([mstr '_W1'])(:,im)   = IN.([mstr '_W1'])(:,im)   + real(sum((Wm.*H_w1).*Cs_n, 2));
            IS.([mstr '_UQW1'])(:,im) = IS.([mstr '_UQW1'])(:,im) + real(sum((Wm.*H_uq).*Cs_s, 2));
            IN.([mstr '_UQW1'])(:,im) = IN.([mstr '_UQW1'])(:,im) + real(sum((Wm.*H_uq).*Cs_n, 2));
        end
    end
end

% Prep 2D images and verify EXACT agreement with archived results/<mname>.json!
ref_2d = prep_img(I_ref, P.nzi, P.nxi, z_in);
ref_hf = highpass_kz(ref_2d, P.dx, 10.0);
iz_w = find(z_in >= 560 & z_in <= 950);
ix_w = find(x_in >= 180 & x_in <= 1420);

D = struct();
D.model = mname; D.z = z_in; D.x = x_in; D.dx = P.dx;
D.Q_in = Q_true(P.in{1}, P.in{2});
D.v_in = v_true(P.in{1}, P.in{2});
D.phi_deg = rad2deg(phi_eff);
D.c_st = c_st;
D.I_ref = single(ref_2d);
D.iz_w = iz_w; D.ix_w = ix_w;

J = jsondecode(fileread(fullfile(res_dir, [mname '.json'])));
max_err = 0;
keys_save = {'RAW_R0', 'RAW_UQW1', 'SA_flat_UQW1', 'DSA_R0', 'DSA_IQ70', 'DSA_W1', 'DSA_UQW1'};
for ik = 1:numel(keys_save)
    k = keys_save{ik};
    Sm = IS.(k); Nm = IN.(k); nc = size(Sm, 2);
    imgs = zeros(P.nzi, P.nxi, nc, 'single');
    e1 = zeros(1, nc); hf = zeros(1, nc); snr = zeros(1, nc);
    for c = 1:nc
        S2 = prep_img(Sm(:,c), P.nzi, P.nxi, z_in);
        N2 = prep_img(Nm(:,c), P.nzi, P.nxi, z_in);
        T2 = S2 + N2;
        imgs(:,:,c) = single(T2);
        e1(c)  = score_e1(T2, ref_2d, iz_w, ix_w);
        hf(c)  = score_e1(highpass_kz(T2, P.dx, 10.0), ref_hf, iz_w, ix_w);
        snr(c) = 10*log10(sum(S2(iz_w,ix_w).^2,'all') / max(sum(N2(iz_w,ix_w).^2,'all'), 1e-30));
    end
    D.img.(k) = imgs;
    D.met.(k) = struct('E1', e1, 'HF', hf, 'SNR', snr);
    % Verify against archived JSON (m=0.5 is index 1, m=1.0 is index 3 in JSON)
    ref_e1 = J.noisy.W_deep.(k).E1;
    if nc == 1
        err = abs(e1(1) - ref_e1(1));
    else
        err = max(abs(e1(1) - ref_e1(1)), abs(e1(2) - ref_e1(3)));
    end
    max_err = max(max_err, err);
end
fprintf('  [%s] Exact reproducibility check vs %s.json: max |dE1| = %.2e (elapsed %.1fs)\n', ...
    mname, mname, max_err, toc(t0));
if max_err > 1e-4
    warning('Reproducibility mismatch > 1e-4 on %s!', mname);
end
end

%% ========================================================================
%% EXACT DETERMINISTIC RECONSTRUCTION FOR H5_VolveAnc
%% ========================================================================
function D = compute_h5_images(res_dir)
t0 = tic;
seed = 20261205; rng(seed, 'twister');
P.dx = 10.0; P.npml = 18; P.nzi = 100; P.nxi = 160;
P.nz = P.nzi + 2*P.npml; P.nx = P.nxi + 2*P.npml;
P.z = ((1:P.nz)' - P.npml - 1) * P.dx;
P.x = ((1:P.nx)  - P.npml - 1) * P.dx;
P.in = {P.npml + (1:P.nzi), P.npml + (1:P.nxi)};
P.izs = P.npml + 3;
P.ixs = P.npml + 1 + (8:6:152);
P.ixr = P.npml + 1 + (0:2:P.nxi-1);
ns = numel(P.ixs); nr = numel(P.ixr); N = P.nz * P.nx;
[Zi, Xi] = ndgrid(P.in{1}, P.in{2}); lin = sub2ind([P.nz P.nx], Zi(:), Xi(:)); nin = numel(lin);
isrc = sub2ind([P.nz P.nx], P.izs*ones(1,ns), P.ixs);
irec = sub2ind([P.nz P.nx], P.izs*ones(1,nr), P.ixr);
[Zfull, Xfull] = ndgrid(P.z, P.x);
z_in = P.z(P.in{1}); x_in = P.x(P.in{2});
[Zin, Xin] = ndgrid(z_in, x_in);
xs = P.x(P.ixs); xr = P.x(P.ixr);

rec_mask = false(nr, ns); xr_lo = zeros(1, ns); xr_hi = zeros(1, ns);
for s = 1:ns
    rec_mask(:, s) = (xr >= xs(s) - 240.0) & (xr <= xs(s) + 860.0);
    xr_lo(s) = min(xr(rec_mask(:,s))); xr_hi(s) = max(xr(rec_mask(:,s)));
end

v_bg = 1800.0 + 1.45 * max(Zfull, 0);
z_dip = Zfull + 0.10 * (Xfull - 800.0);
dv_layer = zeros(size(Zfull));
z_edges = [180, 230, 290, 350, 420, 490, 560, 595, 635, 670, 710, 750, 790, 830, 865, 905, 945];
for il = 1:numel(z_edges)-1
    step = max(min(0.05 * randn(), 0.09), -0.09);
    dv_layer(z_dip >= z_edges(il) & z_dip < z_edges(il+1)) = step;
end
v_true = v_bg .* (1 + dv_layer);
v_mig  = smooth2d(v_true, 3.0);

Q_true = 55.0 * ones(P.nz, P.nx); Q_true(Zfull < 80) = 1000.0; Q_true(Zfull >= 520) = 35.0;
s_mig = 1.0 ./ v_mig; eps_q = 0.5;
F.T = zeros(P.nz, P.nx, ns, 'single'); F.ts = F.T;
for k = 1:ns
    F.T(:,:,k) = eikonal_fs(s_mig, P.dx, P.izs, P.ixs(k));
    F.ts(:,:,k) = (eikonal_fs(s_mig.*(1+eps_q./Q_true), P.dx, P.izs, P.ixs(k)) - ...
                   eikonal_fs(s_mig.*(1-eps_q./Q_true), P.dx, P.izs, P.ixs(k))) / (2*eps_q);
end
TS_flat = zeros(nin, ns); Fts_in = zeros(nin, ns);
for k = 1:ns
    [~, tt] = two_way(F, k, P); TS_flat(:, k) = double(tt(lin));
    tk = F.ts(:,:,k); Fts_in(:, k) = double(tk(lin));
end
clear F

freqs = 4.0:1.0:37.0; nf = numel(freqs);
fref = 25.0; fpeak_true = 20.0; sig = 0.447; s70 = gain_s(70);
ms = [0.5, 0.7, 1.0, 1.4, 2.0]; nm = numel(ms);
flip_sgn = 2 * (rand(1, ns) > 0.5) - 1;

W_mig = ones(1, nf);
a1 = freqs < 6; W_mig(a1) = 0.5 * (1 - cos(pi*(freqs(a1)-4)/2));
a2 = freqs > 34; W_mig(a2) = 0.5 * (1 + cos(pi*(freqs(a2)-34)/3));

W20 = (20/fpeak_true)^2 * exp(-(20/fpeak_true)^2);
b20 = zeros(N, ns); b20(sub2ind([N ns], isrc, 1:ns)) = -W20 / P.dx^2;
[Lm,Um,Pm,Qm,Rm] = lu(helm_op(v_mig, 20, P)); Ub = Qm*(Um\(Lm\(Pm*(Rm\b20)))); Dbg20 = Ub(irec,:); clear Ub Lm Um Pm Qm Rm
vc20 = v_true ./ (1 - (1./(pi*Q_true))*log(20/fref) - 1i./(2*Q_true));
[Lv,Uv,Pv,Qv,Rv] = lu(helm_op(vc20, 20, P)); Uv20 = Qv*(Uv\(Lv\(Pv*(Rv\b20)))); Dcal = (Uv20(irec,:) - Dbg20) .* rec_mask; clear Uv20 Lv Uv Pv Qv Rv
sig_field = sqrt(mean(abs(Dcal(rec_mask)).^2)) * 10^(-6/20);

fprintf('Computing exact images for H5_VolveAnc...\n');
Cs_tot_all = zeros(nin, ns, nf, 'single'); I_RAW_R0 = zeros(nin, 1);
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

    N_f = sig_field * (1/sqrt(2)) * (randn(nr, ns) + 1i*randn(nr, ns)) .* rec_mask;
    V_tot = extract_bp(D_vi + N_f, Lm, Um, Pm, Qm, Rm, irec, lin, N, ns, P.dx); clear Lm Um Pm Qm Rm
    ct = Us_c .* V_tot;
    Cs_tot_all(:,:,jf) = single(ct);
    I_RAW_R0 = I_RAW_R0 + real(sum(ct, 2));
end

R0_2d = prep_img(I_RAW_R0, P.nzi, P.nxi, z_in);
[phi_eff, ~] = blind_structure_dip(R0_2d, P.dx);
[w_dsa, TS_dsa] = build_dsa_mask_and_tstar(Xin(:), Zin(:), xs, xr_lo, xr_hi, phi_eff(:), Fts_in);

anc_mask = (Zin(:) >= 190 & Zin(:) <= 320) & (Xin(:) >= 240 & Xin(:) <= 1360);
S_anc = zeros(1, nf);
for jf = 1:nf
    img_f = real(sum(w_dsa .* double(Cs_tot_all(:,:,jf)), 2));
    S_anc(jf) = rms(img_f(anc_mask));
end
S_anc_s = movmean(S_anc, 3);
W_shape = W_mig ./ max(S_anc_s, 0.10 * max(S_anc_s));
W_shape = W_shape / median(W_shape(freqs >= 10 & freqs <= 26));
tstar_anc = median(TS_dsa(anc_mask, :), 'all');
TS_rel = max(TS_dsa - tstar_anc, 0.0);

I_DSA_R0 = zeros(nin, 1);
I_RAW_UQW1 = zeros(nin, nm); I_DSA_UQW1 = zeros(nin, nm);
I_DSA_Anc = zeros(nin, nm);  I_DSA_Anc_W1 = zeros(nin, nm);

xmax = 1.05 * max(max(TS_flat(:)), max(TS_dsa(:))) / min(ms);
xg = linspace(0, xmax, 1201); dxg = xg(2);

for jf = 1:nf
    f = freqs(jf); Cs_t = double(Cs_tot_all(:,:,jf));
    I_DSA_R0 = I_DSA_R0 + real(sum(w_dsa .* Cs_t, 2));

    pt_r = smooth2d(reshape(abs(sum(Cs_t, 2)).^2, P.nzi, P.nxi), 5.0);
    pf_r = smooth2d(reshape(abs(sum(Cs_t .* flip_sgn, 2)).^2, P.nzi, P.nxi), 5.0);
    iv_r = min(max(pf_r(:) ./ max(pt_r(:) - pf_r(:), 0.05*pf_r(:)), 1e-3), 1e3);

    pt_d = smooth2d(reshape(abs(sum(w_dsa .* Cs_t, 2)).^2, P.nzi, P.nxi), 5.0);
    pf_d = smooth2d(reshape(abs(sum(w_dsa .* Cs_t .* flip_sgn, 2)).^2, P.nzi, P.nxi), 5.0);
    iv_d = min(max(pf_d(:) ./ max(pt_d(:) - pf_d(:), 0.05*pf_d(:)), 1e-3), 1e3);

    T0 = uq_tab(f, xg, 0, fref); T1 = uq_tab(f, xg, sig, fref);
    for im = 1:nm
        u = (TS_flat / ms(im)) / dxg; i0 = min(floor(u), numel(xg)-2); ww = u - i0; i0 = i0 + 1;
        it = @(tb) tb(i0).*(1-ww) + tb(i0+1).*ww;
        A1 = it(T1.A); C1 = it(T1.C); P1 = it(T1.P);
        Huq = (A1 + s70*P1.*C1 + P1.*iv_r) ./ (1 + s70*C1 + iv_r);
        I_RAW_UQW1(:, im) = I_RAW_UQW1(:, im) + real(sum(Huq .* Cs_t, 2));

        u = (TS_dsa / ms(im)) / dxg; i0 = min(floor(u), numel(xg)-2); ww = u - i0; i0 = i0 + 1;
        it = @(tb) tb(i0).*(1-ww) + tb(i0+1).*ww;
        A1 = it(T1.A); C1 = it(T1.C); P1 = it(T1.P);
        Huq = (A1 + s70*P1.*C1 + P1.*iv_d) ./ (1 + s70*C1 + iv_d);
        I_DSA_UQW1(:, im) = I_DSA_UQW1(:, im) + real(sum((w_dsa .* Huq) .* Cs_t, 2));

        u = (TS_rel / ms(im)) / dxg; i0 = min(floor(u), numel(xg)-2); ww = u - i0; i0 = i0 + 1;
        it = @(tb) tb(i0).*(1-ww) + tb(i0+1).*ww;
        A0 = it(T0.A); C0 = it(T0.C); P0 = it(T0.P);
        A1 = it(T1.A); C1 = it(T1.C); P1 = it(T1.P);
        Hw1 = (A0 + s70*P0.*C0 + P0.*iv_d) ./ (1 + s70*C0 + iv_d);
        Huq = (A1 + s70*P1.*C1 + P1.*iv_d) ./ (1 + s70*C1 + iv_d);
        I_DSA_Anc(:, im)    = I_DSA_Anc(:, im)    + W_shape(jf) * real(sum((w_dsa .* Huq) .* Cs_t, 2));
        I_DSA_Anc_W1(:, im) = I_DSA_Anc_W1(:, im) + W_shape(jf) * real(sum((w_dsa .* Hw1) .* Cs_t, 2));
    end
end

well_xs = [400, 600, 800, 1000, 1200];
iz_res = find(z_in >= 540 & z_in <= 940);
v_in = v_true(P.in{1}, P.in{2});
rw_map = zeros(size(v_in)); rw_map(2:end, :) = diff(log(v_in), 1, 1);
lo_k = 4.0; hi_k = 18.0;

% Bandpass full reflectivity map for visual reference panel
rw_bp_2d = zeros(size(rw_map));
for ix = 1:P.nxi, rw_bp_2d(:, ix) = bp_z(rw_map(:, ix), P.dx, lo_k, hi_k); end
rw_bp_2d(z_in < 160, :) = 0;

D = struct();
D.z = z_in; D.x = x_in; D.dx = P.dx; D.m = ms; D.well_xs = well_xs;
D.iz_res = iz_res; D.lo_k = lo_k; D.hi_k = hi_k;
D.rw_map = single(rw_map); D.rw_bp_2d = single(rw_bp_2d);
D.RAW_R0 = single(prep_img(I_RAW_R0, P.nzi, P.nxi, z_in));
D.DSA_R0 = single(prep_img(I_DSA_R0, P.nzi, P.nxi, z_in));
D.RAW_UQW1_m1 = single(prep_img(I_RAW_UQW1(:,3), P.nzi, P.nxi, z_in));
D.DSA_UQW1_m1 = single(prep_img(I_DSA_UQW1(:,3), P.nzi, P.nxi, z_in));
D.DSA_Anc_m1  = single(prep_img(I_DSA_Anc(:,3),  P.nzi, P.nxi, z_in));
D.DSA_Anc_m05 = single(prep_img(I_DSA_Anc(:,1),  P.nzi, P.nxi, z_in));
D.RAW_UQW1_m05= single(prep_img(I_RAW_UQW1(:,1), P.nzi, P.nxi, z_in));

eval_w = @(Ivec) score_wells(prep_img(Ivec, P.nzi, P.nxi, z_in), rw_map, x_in, iz_res, well_xs, P.dx, lo_k, hi_k);
D.tie.RAW_R0 = eval_w(I_RAW_R0);
D.tie.DSA_R0 = eval_w(I_DSA_R0);
D.tie.RAW_UQW1 = arrayfun(@(im) eval_w(I_RAW_UQW1(:,im)), 1:nm);
D.tie.DSA_UQW1 = arrayfun(@(im) eval_w(I_DSA_UQW1(:,im)), 1:nm);
D.tie.DSA_W1_Anc = arrayfun(@(im) eval_w(I_DSA_Anc_W1(:,im)), 1:nm);
D.tie.DSA_UQW1_Anc = arrayfun(@(im) eval_w(I_DSA_Anc(:,im)), 1:nm);

J = jsondecode(fileread(fullfile(res_dir, 'H5_VolveAnc.json')));
err = max([abs(D.tie.RAW_R0 - J.RAW_R0), ...
           max(abs(D.tie.RAW_UQW1(:) - J.RAW_UQW1.tie(:))), ...
           max(abs(D.tie.DSA_UQW1_Anc(:) - J.DSA_UQW1_Anc.tie(:)))]);
fprintf('  [H5_VolveAnc] Exact reproducibility check vs H5_VolveAnc.json: max |dtie| = %.2e (elapsed %.1fs)\n', err, toc(t0));
end

%% ========================================================================
%% FIGURE 1: H2 Reverse Fault + Steep Dip (+22 deg / -16 deg)
%% ========================================================================
function plot_fig1_H2(D, out_png)
z = D.z; x = D.x; izp = find(z >= 160 & z <= 980);
ref = double(D.I_ref);
c = prctile(abs(ref(D.iz_w, D.ix_w)), 99, 'all');
normW = @(A) double(A) * (norm(ref(D.iz_w, D.ix_w), 'fro') / max(norm(double(A(D.iz_w, D.ix_w)), 'fro'), 1e-30));

raw_r0  = normW(D.img.RAW_R0(:,:,1));
raw_m1  = normW(D.img.RAW_UQW1(:,:,2));
saf_m1  = normW(D.img.SA_flat_UQW1(:,:,2));
dsa_m1  = normW(D.img.DSA_UQW1(:,:,2));
iq_m05  = normW(D.img.DSA_IQ70(:,:,1));
w1_m05  = normW(D.img.DSA_W1(:,:,1));
dsa_m05 = normW(D.img.DSA_UQW1(:,:,1));

f = figure('Visible', 'off', 'Color', 'w', 'Position', [0 0 2180 940]);
t = tiledlayout(2, 5, 'TileSpacing', 'compact', 'Padding', 'compact');
title(t, 'H2 逆断层 + 陡倾角 (+22°/-16°) 独立验证模型成像对比（3 dB 强噪声）', ...
    'FontSize', 14, 'FontWeight', 'bold');

% Row 1: m=1.0
ax = nexttile; imagesc(x, z(izp), ref(izp,:), [-c c]); colormap(ax, gray); hold on;
contour(x, z(izp), double(D.Q_in(izp,:)), [35 35], '--', 'LineColor', [1 0.65 0], 'LineWidth', 1.5);
draw_box(ax, [180 1420], [560 950], [0 0.78 0.95]);
text(980, 360, '低Q楔体(Q=16)', 'Color', [1 0.85 0.15], 'FontWeight', 'bold', 'FontSize', 8.5, 'HorizontalAlignment', 'center', 'BackgroundColor', [0.12 0.12 0.12], 'Margin', 1.5);
text(520, 620, '+22°', 'Color', [0.3 1 0.5], 'FontWeight', 'bold', 'FontSize', 9, 'BackgroundColor', [0.12 0.12 0.12], 'Margin', 1.5);
text(1230, 690, '-16°', 'Color', [0.3 1 0.5], 'FontWeight', 'bold', 'FontSize', 9, 'BackgroundColor', [0.12 0.12 0.12], 'Margin', 1.5);
title('(a) 弹性参考真值 I_{ref}'); ylabel('深度 z (m)');

ax = nexttile; imagesc(x, z(izp), raw_r0(izp,:), [-c c]); colormap(ax, gray); hold on;
draw_box(ax, [180 1420], [560 950], [0 0.78 0.95]);
title(sprintf('(b) RAW\\_R0 (不补偿)\nE1=%.3f | HF=%.3f', D.met.RAW_R0.E1(1), D.met.RAW_R0.HF(1)));

ax = nexttile; imagesc(x, z(izp), raw_m1(izp,:), [-c c]); colormap(ax, gray); hold on;
draw_box(ax, [180 1420], [560 950], [0 0.78 0.95]);
title(sprintf('(c) RAW\\_UQW1 (m=1.0, 全炮画弧)\nE1=%.3f | HF=%.3f', D.met.RAW_UQW1.E1(2), D.met.RAW_UQW1.HF(2)));

ax = nexttile; imagesc(x, z(izp), saf_m1(izp,:), [-c c]); colormap(ax, gray); hold on;
draw_box(ax, [180 1420], [560 950], [0 0.78 0.95]);
title(sprintf('(d) SA\\_flat\\_UQW1 (m=1.0, 水平镜面)\nE1=%.3f | HF=%.3f (伤陡倾翼)', D.met.SA_flat_UQW1.E1(2), D.met.SA_flat_UQW1.HF(2)));

ax = nexttile; imagesc(x, z(izp), dsa_m1(izp,:), [-c c]); colormap(ax, gray); hold on;
draw_box(ax, [180 1420], [560 950], [0 0.78 0.95]);
de1_m1 = D.met.DSA_UQW1.E1(2) - D.met.RAW_UQW1.E1(2);
title(sprintf('(e) DSA\\_UQW1 (m=1.0, 本文盲倾角)\nE1=%.3f (%+.3f) | HF=%.3f', D.met.DSA_UQW1.E1(2), de1_m1, D.met.DSA_UQW1.HF(2)), 'Color', [0 0.35 0.75]);

% Row 2: m=0.5 severe over-compensation + blind dip + rescued steep-dip difference
ax = nexttile; imagesc(x, z(izp), iq_m05(izp,:), [-c c]); colormap(ax, gray); hold on;
draw_box(ax, [180 1420], [560 950], [0 0.78 0.95]);
title(sprintf('(f) DSA\\_IQ70 (m=0.5 过补崩盘)\nE1=%.3f | SNR=%.1f dB', D.met.DSA_IQ70.E1(1), D.met.DSA_IQ70.SNR(1)), 'Color', [0.75 0.1 0.1]);
ylabel('深度 z (m)'); xlabel('水平位置 x (m)');

ax = nexttile; imagesc(x, z(izp), w1_m05(izp,:), [-c c]); colormap(ax, gray); hold on;
draw_box(ax, [180 1420], [560 950], [0 0.78 0.95]);
title(sprintf('(g) DSA\\_W1 (m=0.5 单点Wiener)\nE1=%.3f | SNR=%.1f dB', D.met.DSA_W1.E1(1), D.met.DSA_W1.SNR(1)));
xlabel('水平位置 x (m)');

ax = nexttile; imagesc(x, z(izp), dsa_m05(izp,:), [-c c]); colormap(ax, gray); hold on;
draw_box(ax, [180 1420], [560 950], [0 0.78 0.95]);
de1_m05 = D.met.DSA_UQW1.E1(1) - D.met.DSA_W1.E1(1);
title(sprintf('(h) DSA\\_UQW1 (m=0.5 抗过补稳健)\nE1=%.3f (%+.3f vs W1)', D.met.DSA_UQW1.E1(1), de1_m05), 'Color', [0 0.35 0.75]);
xlabel('水平位置 x (m)');

ax = nexttile; imagesc(x, z(izp), D.phi_deg(izp,:), [-25 25]); colormap(ax, bwr_cmap(256)); hold on;
draw_box(ax, [180 1420], [560 950], [0.1 0.1 0.1]);
cb = colorbar(ax); cb.Label.String = '倾角 (°)';
title('(i) 盲估局部倾角 \phi_{eff}(x,z)');
xlabel('水平位置 x (m)');

ax = nexttile;
diff_steep = dsa_m1 - saf_m1;
imagesc(x, z(izp), diff_steep(izp,:), [-0.55*c 0.55*c]); colormap(ax, bwr_cmap(256)); hold on;
draw_box(ax, [180 1420], [560 950], [0.1 0.1 0.1]);
title('(j) 拯救陡倾信号：DSA - SA_{flat} (m=1.0)');
xlabel('水平位置 x (m)');

exportgraphics(f, out_png, 'Resolution', 150); close(f);
end

%% ========================================================================
%% FIGURE 2: H3 Volve-like One-Sided Trailing Streamer + Q=15 Chimney
%% ========================================================================
function plot_fig2_H3(D, out_png)
z = D.z; x = D.x; izp = find(z >= 160 & z <= 980);
ref = double(D.I_ref);
c = prctile(abs(ref(D.iz_w, D.ix_w)), 99, 'all');
normW = @(A) double(A) * (norm(ref(D.iz_w, D.ix_w), 'fro') / max(norm(double(A(D.iz_w, D.ix_w)), 'fro'), 1e-30));

raw_r0  = normW(D.img.RAW_R0(:,:,1));
raw_m1  = normW(D.img.RAW_UQW1(:,:,2));
saf_m1  = normW(D.img.SA_flat_UQW1(:,:,2));
dsa_m1  = normW(D.img.DSA_UQW1(:,:,2));
iq_m05  = normW(D.img.DSA_IQ70(:,:,1));
w1_m05  = normW(D.img.DSA_W1(:,:,1));
dsa_m05 = normW(D.img.DSA_UQW1(:,:,1));

f = figure('Visible', 'off', 'Color', 'w', 'Position', [0 0 2180 940]);
t = tiledlayout(2, 5, 'TileSpacing', 'compact', 'Padding', 'compact');
title(t, 'H3 类Volve单边拖缆 + Q=15 气烟囱 + 12°倾斜储层成像对比（66.9%出阵炮，3 dB 强噪声）', ...
    'FontSize', 14, 'FontWeight', 'bold');

% Row 1: m=1.0
ax = nexttile; imagesc(x, z(izp), ref(izp,:), [-c c]); colormap(ax, gray); hold on;
contour(x, z(izp), double(D.Q_in(izp,:)), [35 35], '--', 'LineColor', [1 0.65 0], 'LineWidth', 1.5);
draw_box(ax, [180 1420], [560 950], [0 0.78 0.95]);
text(760, 410, 'Q=15 气烟囱', 'Color', [1 0.85 0.15], 'FontWeight', 'bold', 'FontSize', 9, 'HorizontalAlignment', 'center', 'BackgroundColor', [0.12 0.12 0.12], 'Margin', 1.5);
title('(a) 弹性参考真值 I_{ref}'); ylabel('深度 z (m)');

ax = nexttile; imagesc(x, z(izp), raw_r0(izp,:), [-c c]); colormap(ax, gray); hold on;
draw_box(ax, [180 1420], [560 950], [0 0.78 0.95]);
title(sprintf('(b) RAW\\_R0 (不补偿：烟囱下衰减)\nE1=%.3f | SNR=%.2f dB', D.met.RAW_R0.E1(1), D.met.RAW_R0.SNR(1)));

ax = nexttile; imagesc(x, z(izp), raw_m1(izp,:), [-c c]); colormap(ax, gray); hold on;
draw_box(ax, [180 1420], [560 950], [0 0.78 0.95]);
title(sprintf('(c) RAW\\_UQW1 (m=1.0, 单边强噪声)\nE1=%.3f | SNR=%.2f dB', D.met.RAW_UQW1.E1(2), D.met.RAW_UQW1.SNR(2)));

ax = nexttile; imagesc(x, z(izp), saf_m1(izp,:), [-c c]); colormap(ax, gray); hold on;
draw_box(ax, [180 1420], [560 950], [0 0.78 0.95]);
title(sprintf('(d) SA\\_flat\\_UQW1 (m=1.0, 水平镜面)\nE1=%.3f | SNR=%.2f dB', D.met.SA_flat_UQW1.E1(2), D.met.SA_flat_UQW1.SNR(2)));

ax = nexttile; imagesc(x, z(izp), dsa_m1(izp,:), [-c c]); colormap(ax, gray); hold on;
draw_box(ax, [180 1420], [560 950], [0 0.78 0.95]);
dsnr_m1 = D.met.DSA_UQW1.SNR(2) - D.met.RAW_UQW1.SNR(2);
title(sprintf('(e) DSA\\_UQW1 (m=1.0, 本文)\nE1=%.3f | SNR=%.2f dB (%+.2f dB!)', D.met.DSA_UQW1.E1(2), D.met.DSA_UQW1.SNR(2), dsnr_m1), 'Color', [0 0.35 0.75]);

% Row 2: m=0.5 + Residual Error Comparison (RAW vs DSA)
ax = nexttile; imagesc(x, z(izp), iq_m05(izp,:), [-c c]); colormap(ax, gray); hold on;
draw_box(ax, [180 1420], [560 950], [0 0.78 0.95]);
title(sprintf('(f) DSA\\_IQ70 (m=0.5 过补崩盘)\nE1=%.3f | SNR=%.1f dB', D.met.DSA_IQ70.E1(1), D.met.DSA_IQ70.SNR(1)), 'Color', [0.75 0.1 0.1]);
ylabel('深度 z (m)'); xlabel('水平位置 x (m)');

ax = nexttile; imagesc(x, z(izp), w1_m05(izp,:), [-c c]); colormap(ax, gray); hold on;
draw_box(ax, [180 1420], [560 950], [0 0.78 0.95]);
title(sprintf('(g) DSA\\_W1 (m=0.5 烟囱下噪化)\nE1=%.3f | SNR=%.1f dB', D.met.DSA_W1.E1(1), D.met.DSA_W1.SNR(1)));
xlabel('水平位置 x (m)');

ax = nexttile; imagesc(x, z(izp), dsa_m05(izp,:), [-c c]); colormap(ax, gray); hold on;
draw_box(ax, [180 1420], [560 950], [0 0.78 0.95]);
de1_h3 = D.met.DSA_UQW1.E1(1) - D.met.DSA_W1.E1(1);
title(sprintf('(h) DSA\\_UQW1 (m=0.5 抗过补稳健)\nE1=%.3f (%+.3f vs W1)', D.met.DSA_UQW1.E1(1), de1_h3), 'Color', [0 0.35 0.75]);
xlabel('水平位置 x (m)');

ax = nexttile; imagesc(x, z(izp), (raw_m1(izp,:) - ref(izp,:)), [-0.6*c 0.6*c]); colormap(ax, bwr_cmap(256)); hold on;
draw_box(ax, [180 1420], [560 950], [0.1 0.1 0.1]);
title('(i) 成像误差：RAW_{UQW1} - I_{ref} (m=1.0)');
xlabel('水平位置 x (m)');

ax = nexttile; imagesc(x, z(izp), (dsa_m1(izp,:) - ref(izp,:)), [-0.6*c 0.6*c]); colormap(ax, bwr_cmap(256)); hold on;
draw_box(ax, [180 1420], [560 950], [0.1 0.1 0.1]);
title('(j) 成像误差：DSA_{UQW1} - I_{ref} (显著净空)', 'Color', [0 0.35 0.75]);
xlabel('水平位置 x (m)');

exportgraphics(f, out_png, 'Resolution', 150); close(f);
end

%% ========================================================================
%% FIGURE 3: H1 Syncline Twin Lenses & H4 Correlated Noise + -3% Vel Error
%% ========================================================================
function plot_fig3_H1_H4(D1, D4, out_png)
f = figure('Visible', 'off', 'Color', 'w', 'Position', [0 0 2180 940]);
t = tiledlayout(2, 5, 'TileSpacing', 'compact', 'Padding', 'compact');
title(t, 'H1 向斜双低Q透镜 & H4 60 m横向相关噪声 + -3%速度误差 成像对比（3 dB 强噪声）', ...
    'FontSize', 14, 'FontWeight', 'bold');

Ds = {D1, D4};
row_tags = {'H1 向斜+双低Q透镜', 'H4 60m相关噪+-3%速偏'};
for ir = 1:2
    D = Ds{ir};
    z = D.z; x = D.x; izp = find(z >= 160 & z <= 980);
    ref = double(D.I_ref);
    c = prctile(abs(ref(D.iz_w, D.ix_w)), 99, 'all');
    normW = @(A) double(A) * (norm(ref(D.iz_w, D.ix_w), 'fro') / max(norm(double(A(D.iz_w, D.ix_w)), 'fro'), 1e-30));

    raw_m1  = normW(D.img.RAW_UQW1(:,:,2));
    dsa_m1  = normW(D.img.DSA_UQW1(:,:,2));
    w1_m05  = normW(D.img.DSA_W1(:,:,1));
    dsa_m05 = normW(D.img.DSA_UQW1(:,:,1));

    ax = nexttile; imagesc(x, z(izp), ref(izp,:), [-c c]); colormap(ax, gray); hold on;
    if ir == 1
        contour(x, z(izp), double(D.Q_in(izp,:)), [35 35], '--', 'LineColor', [1 0.65 0], 'LineWidth', 1.4);
        text(560, 405, 'Q=14', 'Color', [1 0.85 0.15], 'FontWeight', 'bold', 'FontSize', 8.5, 'HorizontalAlignment', 'center', 'BackgroundColor', [0.12 0.12 0.12], 'Margin', 1.5);
        text(1060, 455, 'Q=16', 'Color', [1 0.85 0.15], 'FontWeight', 'bold', 'FontSize', 8.5, 'HorizontalAlignment', 'center', 'BackgroundColor', [0.12 0.12 0.12], 'Margin', 1.5);
    else
        text(800, 320, '横向 Q 梯度: 左 Q=110 \rightarrow 右 Q=28', 'Color', [1 0.85 0.2], 'FontWeight', 'bold', 'FontSize', 8.5, 'HorizontalAlignment', 'center', 'BackgroundColor', [0.12 0.12 0.12], 'Margin', 1.5);
    end
    draw_box(ax, [180 1420], [560 950], [0 0.78 0.95]);
    title(sprintf('%s\n弹性参考真值 I_{ref}', row_tags{ir})); ylabel('深度 z (m)');
    if ir == 2, xlabel('水平位置 x (m)'); end

    ax = nexttile; imagesc(x, z(izp), raw_m1(izp,:), [-c c]); colormap(ax, gray); hold on;
    draw_box(ax, [180 1420], [560 950], [0 0.78 0.95]);
    title(sprintf('RAW\\_UQW1 (m=1.0)\nE1=%.3f | SNR=%.2f dB', D.met.RAW_UQW1.E1(2), D.met.RAW_UQW1.SNR(2)));
    if ir == 2, xlabel('水平位置 x (m)'); end

    ax = nexttile; imagesc(x, z(izp), dsa_m1(izp,:), [-c c]); colormap(ax, gray); hold on;
    draw_box(ax, [180 1420], [560 950], [0 0.78 0.95]);
    de1 = D.met.DSA_UQW1.E1(2) - D.met.RAW_UQW1.E1(2);
    dsnr = D.met.DSA_UQW1.SNR(2) - D.met.RAW_UQW1.SNR(2);
    title(sprintf('DSA\\_UQW1 (m=1.0, 本文)\nE1=%.3f (%+.3f) | SNR %+.2f dB', D.met.DSA_UQW1.E1(2), de1, dsnr), 'Color', [0 0.35 0.75]);
    if ir == 2, xlabel('水平位置 x (m)'); end

    ax = nexttile; imagesc(x, z(izp), w1_m05(izp,:), [-c c]); colormap(ax, gray); hold on;
    draw_box(ax, [180 1420], [560 950], [0 0.78 0.95]);
    title(sprintf('DSA\\_W1 (m=0.5 过补噪化)\nE1=%.3f | SNR=%.1f dB', D.met.DSA_W1.E1(1), D.met.DSA_W1.SNR(1)), 'Color', [0.75 0.1 0.1]);
    if ir == 2, xlabel('水平位置 x (m)'); end

    ax = nexttile; imagesc(x, z(izp), dsa_m05(izp,:), [-c c]); colormap(ax, gray); hold on;
    draw_box(ax, [180 1420], [560 950], [0 0.78 0.95]);
    de_uq = D.met.DSA_UQW1.E1(1) - D.met.DSA_W1.E1(1);
    title(sprintf('DSA\\_UQW1 (m=0.5 抗过补)\nE1=%.3f (%+.3f vs W1)', D.met.DSA_UQW1.E1(1), de_uq), 'Color', [0 0.35 0.75]);
    if ir == 2, xlabel('水平位置 x (m)'); end
end

exportgraphics(f, out_png, 'Resolution', 150); close(f);
end

%% ========================================================================
%% FIGURE 4: H5_VolveAnc Volve Flat-W Replica (2D Sections + Spectra + Well Ties)
%% ========================================================================
function plot_fig4_H5(D, out_png)
z = D.z; x = D.x; izp = find(z >= 160 & z <= 980);
iz_res = D.iz_res; ix_res = find(x >= 240 & x <= 1360);

rw_bp = double(D.rw_bp_2d);
raw_r0 = double(D.RAW_R0);
raw_uq = double(D.RAW_UQW1_m1);
dsa_anc = double(D.DSA_Anc_m1);

% Match RMS in reservoir window for fair visual comparison
c_ref = prctile(abs(raw_r0(iz_res, ix_res)), 99, 'all');
normR = @(A) A * (norm(raw_r0(iz_res, ix_res), 'fro') / max(norm(A(iz_res, ix_res), 'fro'), 1e-30));
rw_vis  = normR(rw_bp);
uq_vis  = normR(raw_uq);
anc_vis = normR(dsa_anc);

f = figure('Visible', 'off', 'Color', 'w', 'Position', [0 0 2050 960]);
t = tiledlayout(2, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
title(t, 'H5\_VolveAnc 类Volve平谱迁移（W_{mig}=1, (2\pif)^2）成像剖面与 5 井井震标定对比', ...
    'FontSize', 14, 'FontWeight', 'bold');

% Row 1: 4 2D sections
imgs_r1 = {rw_vis, raw_r0, uq_vis, anc_vis};
tits_r1 = { ...
    '(a) 真实反射系数剖面 r_w (4-18 cyc/km)', ...
    sprintf('(b) RAW\\_R0 (不补偿)\n5井平均相关 = %.4f', D.tie.RAW_R0), ...
    sprintf('(c) RAW\\_UQW1 (m=1.0, 地表t*未锚定)\n5井平均相关 = %.4f (过蓝倾噪化!)', D.tie.RAW_UQW1(3)), ...
    sprintf('(d) DSA\\_UQW1\\_Anc (m=1.0, 浅层锚定)\n5井平均相关 = %.4f (+0.085 挽回!)', D.tie.DSA_UQW1_Anc(3))};
cols_r1 = {[0 0 0], [0 0 0], [0.75 0.1 0.1], [0 0.35 0.75]};

for p = 1:4
    ax = nexttile;
    imagesc(x, z(izp), imgs_r1{p}(izp,:), [-c_ref c_ref]); colormap(ax, gray); hold on;
    draw_box(ax, [240 1360], [190 320], [0.2 0.9 0.3]);  % shallow anchor box
    draw_box(ax, [240 1360], [540 940], [0 0.78 0.95]);  % reservoir window
    for wx = D.well_xs
        plot([wx wx], [540 940], ':', 'Color', [1 0.85 0.2], 'LineWidth', 1.0);
    end
    if p == 1
        text(800, 255, '浅层锚定窗 [190-320 m]', 'Color', [0.3 1 0.4], 'FontWeight', 'bold', 'FontSize', 8, 'HorizontalAlignment', 'center', 'BackgroundColor', [0.12 0.12 0.12], 'Margin', 1.5);
        text(800, 580, '5口验证井 (400..1200 m)', 'Color', [1 0.9 0.3], 'FontWeight', 'bold', 'FontSize', 8, 'HorizontalAlignment', 'center', 'BackgroundColor', [0.12 0.12 0.12], 'Margin', 1.5);
    end
    title(tits_r1{p}, 'Color', cols_r1{p});
    if p == 1, ylabel('深度 z (m)'); end
    xlabel('水平位置 x (m)');
end

% Panel (e): Vertical wavenumber amplitude spectrum in reservoir [540, 940] m
ax = nexttile; hold on;
nfft = 256; kz = (0:nfft/2) / (nfft * D.dx) * 1000;
win = hann(numel(iz_res));
spec_fn = @(A) mean(abs(fft(A(iz_res, ix_res) .* win, nfft, 1)), 2);
s_rw  = spec_fn(double(D.rw_map)); s_rw = s_rw(1:nfft/2+1);
s_r0  = spec_fn(raw_r0);           s_r0 = s_r0(1:nfft/2+1);
s_uq  = spec_fn(raw_uq);           s_uq = s_uq(1:nfft/2+1);
s_anc = spec_fn(dsa_anc);          s_anc = s_anc(1:nfft/2+1);
km = (kz >= D.lo_k & kz <= D.hi_k);
patch([D.lo_k D.hi_k D.hi_k D.lo_k], [0 0 2.7 2.7], [0.92 0.96 1.0], 'EdgeColor', 'none', 'HandleVisibility', 'off');
plot(kz, s_rw  / mean(s_rw(km)),  'k-',  'LineWidth', 2.0, 'DisplayName', '真实 r_w 谱');
plot(kz, s_r0  / mean(s_r0(km)),  '--',  'Color', [0.45 0.45 0.45], 'LineWidth', 1.6, 'DisplayName', 'RAW\_R0');
plot(kz, s_uq  / mean(s_uq(km)),  'r-',  'LineWidth', 1.8, 'DisplayName', 'RAW\_UQW1 (过蓝倾)');
plot(kz, s_anc / mean(s_anc(km)), 'b-',  'LineWidth', 2.0, 'DisplayName', 'DSA\_UQW1\_Anc (均衡)');
xlim([0 30]); ylim([0 2.65]); grid on;
xlabel('深度波数 k_z (周/km)'); ylabel('带内均值归一化振幅谱');
title('(e) 储层窗 [540-940 m] 深度波数谱对比');
legend('Location', 'northeast', 'FontSize', 7.8);

% Panel (f) & (g): Well-trace wiggle overlay at x=600 m and x=1000 m
z_res = z(iz_res);
well_show = [600, 1000];
sub_lab = {'(f)', '(g)'};
for iw = 1:2
    wx_t = well_show(iw);
    [~, jx] = min(abs(x - wx_t));
    rw_tr = bp_z(double(D.rw_map(iz_res, jx)), D.dx, D.lo_k, D.hi_k);
    rw_tr = rw_tr / rms(rw_tr);
    [tr_r0, c_r0]   = best_aligned_trace(raw_r0,  jx, iz_res, rw_tr, D.dx, D.lo_k, D.hi_k);
    [tr_uq, c_uq]   = best_aligned_trace(raw_uq,  jx, iz_res, rw_tr, D.dx, D.lo_k, D.hi_k);
    [tr_anc, c_anc] = best_aligned_trace(dsa_anc, jx, iz_res, rw_tr, D.dx, D.lo_k, D.hi_k);

    ax = nexttile; hold on;
    plot(rw_tr,  z_res, 'k-',  'LineWidth', 2.0, 'DisplayName', '真实 r_w (带通)');
    plot(tr_r0,  z_res, '--',  'Color', [0.5 0.5 0.5], 'LineWidth', 1.3, 'DisplayName', sprintf('RAW\\_R0 (r=%.3f)', c_r0));
    plot(tr_uq,  z_res, 'r-',  'LineWidth', 1.3, 'DisplayName', sprintf('RAW\\_UQW1 (r=%.3f)', c_uq));
    plot(tr_anc, z_res, 'b-',  'LineWidth', 1.9, 'DisplayName', sprintf('DSA\\_UQW1\\_Anc (r=%.3f)', c_anc));
    set(ax, 'YDir', 'reverse'); ylim([540 940]); xlim([-3.2 3.2]); grid on;
    xlabel('归一化道振幅'); ylabel('深度 z (m)');
    title(sprintf('%s 井震波形叠加 (x = %d m, m=1.0)', sub_lab{iw}, wx_t));
    legend('Location', 'southwest', 'FontSize', 7.8);
end

% Panel (h): 5-well mean correlation vs m
ax = nexttile; hold on;
ms = D.m;
yline(D.tie.RAW_R0, 'k--', 'LineWidth', 1.6, 'DisplayName', sprintf('RAW\\_R0 (%.3f)', D.tie.RAW_R0));
plot(ms, D.tie.RAW_UQW1,    'r-^', 'LineWidth', 1.6, 'MarkerSize', 6, 'DisplayName', sprintf('RAW\\_UQW1 (均%.3f)', mean(D.tie.RAW_UQW1)));
plot(ms, D.tie.DSA_UQW1,    '--s', 'Color', [0.9 0.45 0], 'LineWidth', 1.4, 'MarkerSize', 5, 'DisplayName', sprintf('DSA\\_UQW1 (均%.3f)', mean(D.tie.DSA_UQW1)));
plot(ms, D.tie.DSA_W1_Anc,  '-d', 'Color', [0.1 0.65 0.25], 'LineWidth', 1.6, 'MarkerSize', 6, 'DisplayName', sprintf('DSA\\_W1\\_Anc (均%.3f)', mean(D.tie.DSA_W1_Anc)));
plot(ms, D.tie.DSA_UQW1_Anc,'b-o', 'LineWidth', 2.3, 'MarkerSize', 7, 'MarkerFaceColor', 'b', 'DisplayName', sprintf('DSA\\_UQW1\\_Anc (均%.3f)', mean(D.tie.DSA_UQW1_Anc)));
set(ax, 'XScale', 'log'); xticks(ms); ylim([0.52 0.75]); grid on;
xlabel('Q 乘子 m (m<1 过补, m>1 欠补)'); ylabel('5井平均井震相关系数');
title('(h) 5井平均井震标定随 Q 误差变化');
legend('Location', 'southeast', 'FontSize', 7.8);

exportgraphics(f, out_png, 'Resolution', 150); close(f);
end

%% ========================================================================
%% FIGURE 5: Real Volve North Sea Stage 4 Cached Comparison + Exploratory Shallow-Anchor 1D Test
%% ========================================================================
function plot_fig5_volve_real(out_png, scratch_dir)
P = fld_setup();
S1 = load(fullfile(P.cache, 'state1.mat'));
T  = load(fullfile(P.cache, 'ts.mat'));
Iv = load(fullfile(P.cache, 'inv.mat'));
Ia = load(fullfile(P.cache, 'images_amp.mat'));

z = P.z(P.in{1}); x = P.x(P.in{2});
ker = exp(-0.5*((-12:12)/4).^2); ker = ker/sum(ker);
hp = @(A) zero_top(A - conv2(ker', 1, [repmat(A(1,:),12,1); A; repmat(A(end,:),12,1)], 'valid'), z, 150);

% Load F-5 & F-14 well logs
w5 = load_volve_well(fullfile(P.cache, 'well_F-5.bin'), z, x);
w14 = load_volve_well(fullfile(P.cache, 'well_F-14.bin'), z, x);
lo_v = 3.52; hi_v = 12.89;
rb5  = bp_z(w5.rw,  P.dx, lo_v, hi_v);
rb14 = bp_z(w14.rw, P.dx, lo_v, hi_v);

% Baseline cached R0 and unanchored Stage-4 UQW1_amp / W1_amp (m = [0.5, 0.7, 1.0, 1.4, 2.0])
I_R0_2d = hp(reshape(double(sum(S1.I, 2)), P.nzi, P.nxi));
m_idx = [1 3 5 7 9]; m_vals = Ia.meta.m(m_idx);
nm = numel(m_vals);
cc5_r0  = best_cc_well(I_R0_2d, w5,  rb5,  P.dx, lo_v, hi_v);
cc14_r0 = best_cc_well(I_R0_2d, w14, rb14, P.dx, lo_v, hi_v);

cc5_uq_un  = zeros(1, nm); cc5_w1_un  = zeros(1, nm);
cc14_uq_un = zeros(1, nm);
for k = 1:nm
    im = m_idx(k);
    Au = hp(reshape(double(Ia.UQW1(:, im)), P.nzi, P.nxi));
    Aw = hp(reshape(double(Ia.W1(:, im)),   P.nzi, P.nxi));
    cc5_uq_un(k)  = best_cc_well(Au, w5,  rb5,  P.dx, lo_v, hi_v);
    cc5_w1_un(k)  = best_cc_well(Aw, w5,  rb5,  P.dx, lo_v, hi_v);
    cc14_uq_un(k) = best_cc_well(Au, w14, rb14, P.dx, lo_v, hi_v);
    if k == 2, I_UQ_m07 = Au; end
    if k == 3, I_UQ_m10 = Au; end
end

% Exploratory post-hoc 1D Shallow-Anchored UQW1_Anc on real Volve cached S1.I + T.TS + Iv.iv!
wx_anc = (x >= -2500 & x <= 1500);
wz_anc = (z >= 1000 & z <= 1600);
wz_res = (z >= 2600 & z <= 3200);
ts_med = reshape(median(double(T.TS), 2), P.nzi, P.nxi);
tstar_anc = median(ts_med(wz_anc, wx_anc), 'all');
ts_rel_vec = max(ts_med(:) - tstar_anc, 0.0);

f_vec = P.f(:)'; nf = numel(f_vec);
S_anc = zeros(1, nf); S_res = zeros(1, nf);
for jf = 1:nf
    Aj = reshape(double(S1.I(:, jf)), P.nzi, P.nxi);
    S_anc(jf) = rms(Aj(wz_anc, wx_anc), 'all');
    S_res(jf) = rms(Aj(wz_res, wx_anc), 'all');
end
S_anc_s = movmean(S_anc, 7);
% Band-limited shallow-anchor shaping (gentle cosine taper 4-6 & 22-25 Hz matching Volve usable signal band)
W_taper = ones(1, nf);
a1 = f_vec < 6;  W_taper(a1) = 0.5 * (1 - cos(pi*max(f_vec(a1)-4,0)/2));
a2 = f_vec > 22; W_taper(a2) = 0.5 * (1 + cos(pi*min(f_vec(a2)-22,3)/3));
W_shape = W_taper ./ max(S_anc_s, 0.10 * max(S_anc_s));
W_shape = W_shape / median(W_shape(f_vec >= 8 & f_vec <= 20));

s70 = gain_s(70);
cc5_anc_uq = zeros(1, nm); cc5_anc_w1 = zeros(1, nm);
cc14_anc_uq = zeros(1, nm);
for k = 1:nm
    m = m_vals(k);
    xg = linspace(0, 1.02 * max(ts_rel_vec) / m + 1e-6, 601);
    I_auq = zeros(numel(ts_rel_vec), 1);
    I_aw1 = zeros(numel(ts_rel_vec), 1);
    for jf = 1:nf
        f = P.f(jf);
        if W_taper(jf) < 1e-4, continue; end
        T0 = uq_tab(f, xg, 0, f); T1 = uq_tab(f, xg, P.sig, f);
        A0 = interp1(xg, real(T0.A), ts_rel_vec/m); C0 = interp1(xg, real(T0.C), ts_rel_vec/m);
        A1 = interp1(xg, real(T1.A), ts_rel_vec/m); C1 = interp1(xg, real(T1.C), ts_rel_vec/m);
        iv_f = double(Iv.iv(:, jf));
        Hw1 = (A0 + s70*C0 + iv_f) ./ (1 + s70*C0 + iv_f);
        Huq = (A1 + s70*C1 + iv_f) ./ (1 + s70*C1 + iv_f);
        If = double(S1.I(:, jf));
        I_aw1 = I_aw1 + W_shape(jf) * (Hw1 .* If);
        I_auq = I_auq + W_shape(jf) * (Huq .* If);
    end
    A_auq = hp(reshape(I_auq, P.nzi, P.nxi));
    A_aw1 = hp(reshape(I_aw1, P.nzi, P.nxi));
    cc5_anc_uq(k)  = best_cc_well(A_auq, w5,  rb5,  P.dx, lo_v, hi_v);
    cc5_anc_w1(k)  = best_cc_well(A_aw1, w5,  rb5,  P.dx, lo_v, hi_v);
    cc14_anc_uq(k) = best_cc_well(A_auq, w14, rb14, P.dx, lo_v, hi_v);
    if k == 3, I_Anc_m10 = A_auq; end
end

fprintf('\n=== Exploratory Real Volve Stage-4 Post-Stack Shallow-Anchor Check (read-only cache) ===\n');
fprintf('  F-5  : R0=%.4f | Unanc UQW1(m=0.5..2)=[%.3f %.3f %.3f %.3f %.3f] | Anc_1D UQW1=[%.3f %.3f %.3f %.3f %.3f] | Anc_1D W1=[%.3f %.3f %.3f %.3f %.3f]\n', ...
    cc5_r0, cc5_uq_un, cc5_anc_uq, cc5_anc_w1);
fprintf('  F-14 : R0=%.4f | Unanc UQW1(m=0.5..2)=[%.3f %.3f %.3f %.3f %.3f] | Anc_1D UQW1=[%.3f %.3f %.3f %.3f %.3f]\n', ...
    cc14_r0, cc14_uq_un, cc14_anc_uq);
save(fullfile(scratch_dir, 'volve_real_anc_summary.mat'), 'm_vals', 'cc5_r0', 'cc5_uq_un', 'cc5_w1_un', 'cc5_anc_uq', 'cc5_anc_w1', 'cc14_r0', 'cc14_uq_un', 'cc14_anc_uq');

% Plot 2x3 Figure
zz = find(z >= 1500 & z <= 3400);
xx = find(x >= -3000 & x <= 2000);
c_v = prctile(abs(I_R0_2d(zz, xx)), 99, 'all');
normV = @(A) A * (norm(I_R0_2d(zz, xx), 'fro') / max(norm(A(zz, xx), 'fro'), 1e-30));

f = figure('Visible', 'off', 'Color', 'w', 'Position', [0 0 2000 980]);
t = tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
title(t, 'Volve 北海实测海上地震数据（Stage 4 缓存）：原版未锚定 vs 浅层锚定（探索性）成像与 F-5 井对比', ...
    'FontSize', 14, 'FontWeight', 'bold');

pans = {I_R0_2d, normV(I_UQ_m10), normV(I_Anc_m10)};
ptits = { ...
    sprintf('(a) Volve 实测 R0 (不补偿)\nF-5=%.3f | F-14=%.3f', cc5_r0, cc14_r0), ...
    sprintf('(b) 原版 UQW1\\_amp (m=1.0, 地表t*未锚定)\nF-5=%.3f | F-14=%.3f (高频蓝倾过重)', cc5_uq_un(3), cc14_uq_un(3)), ...
    sprintf('(c) 浅层锚定 UQW1\\_Anc\\_1D (m=1.0, 探索性)\nF-5=%.3f | F-14=%.3f (浅层1000-1600m谱均衡)', cc5_anc_uq(3), cc14_anc_uq(3))};
pcols = {[0 0 0], [0.75 0.1 0.1], [0 0.35 0.75]};

for p = 1:3
    ax = nexttile;
    imagesc(x(xx)/1000, z(zz), pans{p}(zz, xx), [-c_v c_v]); colormap(ax, gray); hold on;
    plot(w5.sw/1000,  w5.zw,  'r-', 'LineWidth', 1.8);
    plot(w14.sw/1000, w14.zw, 'c-', 'LineWidth', 1.8);
    draw_box(ax, [-2.5 1.5], [2600 3200], [1 0.85 0.1]);
    if p == 1
        text(-1.18, 2500, 'F-5', 'Color', 'r', 'FontWeight', 'bold', 'FontSize', 9, 'BackgroundColor', [1 1 1], 'Margin', 1.2);
        text(mean(w14.sw)/1000, min(w14.zw)-70, 'F-14', 'Color', [0 0.55 0.75], 'FontWeight', 'bold', 'FontSize', 9, 'BackgroundColor', [1 1 1], 'Margin', 1.2);
    end
    title(ptits{p}, 'Color', pcols{p});
    if p == 1, ylabel('深度 z (m)'); end
    xlabel('测线坐标 x (km)');
end

% Panel (d): Shallow vs Reservoir R0 frequency spectrum (Empirical Q_eff = 77.1 & blue tilt)
ax = nexttile; hold on;
fm = f_vec >= 4 & f_vec <= 25;
plot(f_vec(fm), S_anc(fm) / max(S_anc(fm)), 'Color', [0.1 0.65 0.25], 'LineWidth', 2.0, 'DisplayName', '浅层 1000-1600 m (受(2\pif)^2强蓝倾)');
plot(f_vec(fm), S_res(fm) / max(S_res(fm)), 'k-', 'LineWidth', 2.0, 'DisplayName', '储层 2600-3200 m (高频实测衰减)');
plot(f_vec(fm), (S_res(fm)./S_anc(fm)) / max(S_res(fm)./S_anc(fm)), 'b--', 'LineWidth', 1.8, 'DisplayName', '谱比 S_{储层}/S_{浅层} (\Delta t^*=10.9ms, Q_{eff}=77)');
grid on; xlim([4 25]); ylim([0 1.15]);
xlabel('频率 f (Hz)'); ylabel('归一化 RMS 振幅谱');
title('(d) Volve 实测浅层蓝倾与深层真实衰减谱比');
legend('Location', 'southwest', 'FontSize', 8);

% Panel (e): F-5 well trace overlay (2570-3170 m)
ax = nexttile; hold on;
rb5_n = rb5 / rms(rb5);
[t5_r0,  c5_r0]  = align_volve_well(I_R0_2d,   w5, rb5_n, P.dx, lo_v, hi_v);
[t5_uq,  c5_uq]  = align_volve_well(I_UQ_m10,  w5, rb5_n, P.dx, lo_v, hi_v);
[t5_anc, c5_anc] = align_volve_well(I_Anc_m10, w5, rb5_n, P.dx, lo_v, hi_v);
plot(rb5_n,  w5.zw, 'k-',  'LineWidth', 2.0, 'DisplayName', 'F-5 测井反射系数 r_w');
plot(t5_r0,  w5.zw, '--',  'Color', [0.5 0.5 0.5], 'LineWidth', 1.4, 'DisplayName', sprintf('R0 (r=%.3f)', c5_r0));
plot(t5_uq,  w5.zw, 'r-',  'LineWidth', 1.3, 'DisplayName', sprintf('原版UQW1 未锚定 (r=%.3f)', c5_uq));
plot(t5_anc, w5.zw, 'b-',  'LineWidth', 1.9, 'DisplayName', sprintf('浅层锚定UQW1\\_1D (r=%.3f)', c5_anc));
set(ax, 'YDir', 'reverse'); ylim([min(w5.zw) max(w5.zw)]); xlim([-3.3 3.3]); grid on;
xlabel('归一化道振幅'); ylabel('深度 z (m)');
title('(e) F-5 井储层段井震波形对比 (m=1.0)');
legend('Location', 'southwest', 'FontSize', 7.8);

% Panel (f): F-5 tie vs m
ax = nexttile; hold on;
yline(cc5_r0, 'k--', 'LineWidth', 1.6, 'DisplayName', sprintf('R0 (%.3f)', cc5_r0));
plot(m_vals, cc5_w1_un,  '--s', 'Color', [0.85 0.45 0.1], 'LineWidth', 1.4, 'MarkerSize', 5, 'DisplayName', '原版 W1 (未锚定)');
plot(m_vals, cc5_uq_un,  'r-^', 'LineWidth', 1.7, 'MarkerSize', 6, 'DisplayName', '原版 UQW1 (未锚定)');
plot(m_vals, cc5_anc_w1, '-d', 'Color', [0.1 0.65 0.25], 'LineWidth', 1.6, 'MarkerSize', 6, 'DisplayName', '浅层锚定 W1\_1D (探索)');
plot(m_vals, cc5_anc_uq, 'b-o', 'LineWidth', 2.2, 'MarkerSize', 7, 'MarkerFaceColor', 'b', 'DisplayName', '浅层锚定 UQW1\_1D (探索)');
set(ax, 'XScale', 'log'); xticks(m_vals); ylim([0.30 0.76]); grid on;
xlabel('Q 乘子 m (名义 Q=120)'); ylabel('F-5 井旁相关系数');
title('(f) Volve F-5 井旁相关随 Q 误差变化');
legend('Location', 'southeast', 'FontSize', 7.8);

exportgraphics(f, out_png, 'Resolution', 150); close(f);
end

%% ========================================================================
%% HELPERS
%% ========================================================================
function draw_box(ax, xr, zr, col)
rectangle(ax, 'Position', [xr(1), zr(1), xr(2)-xr(1), zr(2)-zr(1)], ...
    'EdgeColor', col, 'LineStyle', '--', 'LineWidth', 1.1);
end

function cm = bwr_cmap(n)
h = floor(n/2);
r1 = linspace(0.08, 1, h)'; g1 = linspace(0.25, 1, h)'; b1 = linspace(0.78, 1, h)';
r2 = linspace(1, 0.78, n-h)'; g2 = linspace(1, 0.14, n-h)'; b2 = linspace(1, 0.14, n-h)';
cm = [r1 g1 b1; r2 g2 b2];
end

function A = zero_top(A, z, zcut), A(z < zcut, :) = 0; end

function w = load_volve_well(fn, z, x)
fid = fopen(fn); a = fread(fid, Inf, 'single'); fclose(fid);
a = reshape(a, 4, [])'; w.zw = a(:,1); w.sw = a(:,2); w.rw = double(a(:,4));
[~, w.iz] = min(abs(z(:) - w.zw'), [], 1);
w.ix = arrayfun(@(xv) find(abs(x-xv)==min(abs(x-xv)), 1), w.sw);
end

function c = best_cc_well(A, w, rb, dx, lo, hi)
nw = numel(w.zw); tr = zeros(nw, 1);
for i = 1:nw, tr(i) = mean(A(w.iz(i), max(w.ix(i)-2,1):min(w.ix(i)+2,size(A,2)))); end
tr_b = bp_z(tr, dx, lo, hi);
c = -Inf;
for L = -3:3, c = max(c, cc(circshift(tr_b, L), rb)); end
end

function [tr_best, c_best] = align_volve_well(A, w, rb_n, dx, lo, hi)
nw = numel(w.zw); tr = zeros(nw, 1);
for i = 1:nw, tr(i) = mean(A(w.iz(i), max(w.ix(i)-2,1):min(w.ix(i)+2,size(A,2)))); end
tr_b = bp_z(tr, dx, lo, hi);
c_best = -Inf; tr_best = tr_b;
for L = -3:3
    cand = circshift(tr_b, L);
    val = cc(cand, rb_n);
    if val > c_best, c_best = val; tr_best = cand; end
end
tr_best = tr_best / max(rms(tr_best), 1e-30);
end

function [tr_best, c_best] = best_aligned_trace(A, jx, iz_res, rw_n, dx, lo, hi)
tr = mean(A(iz_res, max(jx-1,1):min(jx+1,size(A,2))), 2);
tr_b = bp_z(tr, dx, lo, hi);
c_best = -Inf; tr_best = tr_b;
for L = -3:3
    cand = circshift(tr_b, L);
    val = cc(cand, rw_n);
    if val > c_best, c_best = val; tr_best = cand; end
end
tr_best = tr_best / max(rms(tr_best), 1e-30);
end

function c_mean = score_wells(A, rw_map, x_in, iz_res, well_xs, dx, lo_k, hi_k)
cs = zeros(1, numel(well_xs));
for iw = 1:numel(well_xs)
    [~, jx] = min(abs(x_in - well_xs(iw)));
    tr = mean(A(iz_res, max(jx-1,1):min(jx+1,numel(x_in))), 2);
    rw = rw_map(iz_res, jx);
    tr_b = bp_z(tr, dx, lo_k, hi_k);
    rw_b = bp_z(rw, dx, lo_k, hi_k);
    best = -Inf;
    for L = -3:3, best = max(best, cc(circshift(tr_b, L), rw_b)); end
    cs(iw) = best;
end
c_mean = mean(cs);
end

function [phi_eff, c_st] = blind_structure_dip(R0_2d, dx)
[gx, gz] = gradient(R0_2d, dx);
Jxx = smooth2d(gx.^2, 4.0); Jzz = smooth2d(gz.^2, 4.0); Jxz = smooth2d(gx.*gz, 4.0);
phi_raw = -0.5 * atan2(2*Jxz, Jzz - Jxx);
phi_raw = max(min(phi_raw, deg2rad(30)), -deg2rad(30));
c_st = sqrt((Jzz - Jxx).^2 + 4*Jxz.^2) ./ max(Jzz + Jxx, 1e-30);
phi_eff = c_st .* phi_raw;
end

function w_sa = build_dsa_mask(xv, zv, xs, xr_lo, xr_hi, phi_v)
nin = numel(xv); ns = numel(xs); w_sa = zeros(nin, ns); dz = max(zv - 20.0, 10.0);
for k = 1:ns
    th_s = atan2(xv - xs(k), dz);
    th_r = max(min(th_s + 2*phi_v, deg2rad(78)), -deg2rad(78));
    xmir = xv + dz .* tan(th_r);
    d_in = min(xmir - xr_lo(k), xr_hi(k) - xmir);
    w_sa(:, k) = max(cos_taper(d_in, -80.0, 120.0) .* cos_taper(80.0 - abs(rad2deg(th_s)), 0, 10) .* cos_taper(80.0 - abs(rad2deg(th_r)), 0, 10), 1e-4);
end
end

function [w_dsa, TS_dsa] = build_dsa_mask_and_tstar(xv, zv, xs, xr_lo, xr_hi, phi_v, Fts_in)
nin = numel(xv); ns = numel(xs); w_dsa = zeros(nin, ns); TS_dsa = zeros(nin, ns);
dz = max(zv - 20.0, 10.0); row_idx = (1:nin)';
for k = 1:ns
    th_s = atan2(xv - xs(k), dz);
    th_r = max(min(th_s + 2*phi_v, deg2rad(78)), -deg2rad(78));
    xmir = xv + dz .* tan(th_r);
    d_in = min(xmir - xr_lo(k), xr_hi(k) - xmir);
    w_dsa(:, k) = max(cos_taper(d_in, -80.0, 120.0) .* cos_taper(80.0 - abs(rad2deg(th_s)), 0, 10) .* cos_taper(80.0 - abs(rad2deg(th_r)), 0, 10), 1e-4);
    q = min(max((xmir - xs(1)) / (xs(2) - xs(1)) + 1, 1), ns);
    k0 = floor(q); k0(k0 == ns) = ns - 1; ww = q - k0;
    i0 = sub2ind([nin ns], row_idx, k0); i1 = sub2ind([nin ns], row_idx, k0 + 1);
    TS_dsa(:, k) = Fts_in(:, k) + (1 - ww).*Fts_in(i0) + ww.*Fts_in(i1);
end
end

function y = bp_z(r, dx, lo, hi)
n = 2^nextpow2(4 * numel(r)); R = fft(r, n);
k = abs([0:n/2, -n/2+1:-1]') / (n * dx) * 1000;
y = real(ifft(R .* (k >= lo & k <= hi))); y = y(1:numel(r));
end
function t = cos_taper(x, x0, x1), u = min(max((x - x0)/max(x1 - x0, 1e-12), 0), 1); t = 0.5 - 0.5*cos(pi*u); end
function V = extract_bp(Dmat, L, U, Pp, Qq, R, irec, lin, N, ns, dx)
b = zeros(N, ns); b(irec, :) = Dmat;
Vf = -(R \ (Pp' * (L' \ (U' \ (Qq' * b))))) / dx^2; V = Vf(lin, :);
end
function B = smooth2d(A, sig)
r = max(1, ceil(3*sig)); g = exp(-0.5*((-r:r)/sig).^2); g = g/sum(g);
Ap = [repmat(A(1,:), r, 1); A; repmat(A(end,:), r, 1)];
B1 = conv2(g', 1, Ap, 'valid');
B1p = [repmat(B1(:,1), 1, r), B1, repmat(B1(:,end), 1, r)];
B = conv2(1, g, B1p, 'valid');
end
function y = prep_img(v, nzi, nxi, z_in), y = reshape(double(v), nzi, nxi); y = y - smooth2d(y, 5.0); y(z_in < 160, :) = 0; end
function yh = highpass_kz(y, dx, kz_cut)
n = 2^nextpow2(2*size(y,1)); Y = fft(y, n, 1);
kz = abs([0:n/2, -n/2+1:-1]')/(n*dx)*1000; g = 0.5 - 0.5*cos(pi*min(max((kz-(kz_cut-2))/4,0),1));
yh = real(ifft(Y.*g, n, 1)); yh = yh(1:size(y,1), :);
end
function s = score_e1(a, r, iz, ix)
rr = r(iz,ix); c0 = cc(a(iz,ix), rr); cl = c0; nz = size(a,1);
for L = [-2 -1 1 2], cl = max(cl, cc(a(min(max(iz+L,1),nz), ix), rr)); end
s = 0.5*(c0 + cl);
end
function c = cc(a, b), a = a(:)-mean(a(:)); b = b(:)-mean(b(:)); c = (a'*b)/max(sqrt((a'*a)*(b'*b)),1e-30); end
