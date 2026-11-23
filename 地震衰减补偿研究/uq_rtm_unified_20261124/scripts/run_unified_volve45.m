function run_unified_volve45(smoke)
% Tier 3 Pre-registered Real Volve North Sea OBC 4-45 Hz Broadband Stage 4H Confluence
% Plan: uq_rtm_unified_20261124/PLAN_ZH.md
%
% Strictly read-only on Stage 4H cache:
%   D:\ChatGPT\研究\.claude\worktrees\study-research-docs-plan-c4d4c6\地震衰减补偿研究\uq_rtm_20261117\cache\fld45
% Reuses Stage 4H cached inv.mat (iv [B=1] & iv3 [B=300]), tfields.mat, ts.mat, images.mat (R0).
% 100% blind to F-5 / F-14 well logs throughout all migration, dip, DSA, and PSAE steps!
% Resumable per frequency via lightweight slices in D:\tmp_uqrtm_unified_20261124\slices_fld45.
%
% Usage:
%   run_unified_volve45(true)   % 3-freq smoke check (~5 min)
%   run_unified_volve45(false)  % full 165-freq broadband production run (~4.5-5 hr)

if nargin < 1, smoke = false; end
maxNumCompThreads(12);
t_all = tic;

here = fileparts(mfilename('fullpath'));
root = fileparts(here);
res_dir = fullfile(root, 'results');
if ~isfolder(res_dir), mkdir(res_dir); end

wt_root = fileparts(root);
addpath(fullfile(wt_root, 'uq_rtm_20261117', 'scripts'));
core = uq_rtm_unified_core();
cfg  = core.default_cfg('field_obc');

% Strictly read-only Stage 4H cache paths
fld45_cache = 'D:\ChatGPT\研究\.claude\worktrees\study-research-docs-plan-c4d4c6\地震衰减补偿研究\uq_rtm_20261117\cache\fld45';
well_dir    = 'D:\ChatGPT\研究\.claude\worktrees\study-research-docs-plan-c4d4c6\地震衰减补偿研究\uq_rtm_20261117\cache\fld';

% Heavy intermediates strictly outside git worktree on D:
scratch_dir = 'D:\tmp_uqrtm_unified_20261124';
if smoke
    slice_dir = fullfile(scratch_dir, 'slices_fld45_smoke');
else
    slice_dir = fullfile(scratch_dir, 'slices_fld45');
end
if ~isfolder(slice_dir), mkdir(slice_dir); end

% Load 6.25 m grid & 4-45 Hz data (read-only)
P = load_fld_ro(fld45_cache);
ns = P.nshot; N = P.nz * P.nx; nf_all = numel(P.f);
[Zi, Xi] = ndgrid(P.in{1}, P.in{2});
lin = sub2ind([P.nz P.nx], Zi(:), Xi(:)); nin = numel(lin);
irec = sub2ind([P.nz P.nx], P.rec_iz, P.rec_ix);
isrc = sub2ind([P.nz P.nx], P.src_iz, P.src_ix);
col  = sub2ind([N ns], irec, P.shot_of);
z_in = P.z(P.in{1}); x_in = P.x(P.in{2});
[Zin, Xin] = ndgrid(z_in, x_in);

wres_mask = (Zin >= cfg.res_z_m(1) & Zin <= cfg.res_z_m(2)) & (Xin >= cfg.res_x_m(1) & Xin <= cfg.res_x_m(2));
anc_mask  = (Zin >= cfg.anc_z_m(1) & Zin <= cfg.anc_z_m(2)) & (Xin >= cfg.anc_x_m(1) & Xin <= cfg.anc_x_m(2));

% 1. Precompute/load lean DSA geometry (avoid keeping 340 MB tfields.mat in RAM during LU!)
geom_file = fullfile(scratch_dir, 'fld45_dsa_geom.mat');
S_im = load(fullfile(fld45_cache, 'images.mat'), 'R0', 'UQW1b');
R0_stage4h = double(S_im.R0);
UQW1b_stage4h_m1 = double(S_im.UQW1b(:, 3));
clear S_im

if isfile(geom_file)
    G = load(geom_file);
    w_dsa = G.w_dsa; TS_dsa = G.TS_dsa; phi_eff = G.phi_eff; c_st = G.c_st;
    tstar_anc_dsa = G.tstar_anc_dsa; tstar_anc_flat = G.tstar_anc_flat;
    clear G
else
    R0_2d = core.prep_image(R0_stage4h, P.nzi, P.nxi, z_in, P.dx, cfg.bg_sigma_m, cfg.z_mute_m);
    [phi_eff, c_st] = core.blind_dip(R0_2d, P.dx, cfg.dip_sigma_m, cfg.dip_cap_deg);
    rx = P.x(P.rec_ix);
    xr_lo = arrayfun(@(s) min(rx(P.shot_of == s)), 1:ns);
    xr_hi = arrayfun(@(s) max(rx(P.shot_of == s)), 1:ns);
    TF = load(fullfile(fld45_cache, 'tfields.mat'), 'F');
    [w_dsa, TS_dsa] = core.build_dsa_full(Xin(:), Zin(:), P.z(P.src_iz), 90.0, P.src_x, ...
        xr_lo, xr_hi, phi_eff(:), TF.F, Zi(:), Xi(:), P.nz, P.nx, cfg);
    clear TF
    TS_flat_S = load(fullfile(fld45_cache, 'ts.mat'), 'TS');
    tstar_anc_flat = median(double(TS_flat_S.TS(anc_mask(:), :)), 'all');
    clear TS_flat_S
    tstar_anc_dsa = median(double(TS_dsa(anc_mask(:), :)), 'all');
    save(geom_file, 'w_dsa', 'TS_dsa', 'phi_eff', 'c_st', 'tstar_anc_dsa', 'tstar_anc_flat', '-v7.3');
end

TS_flat_S = load(fullfile(fld45_cache, 'ts.mat'), 'TS');
TS_flat = TS_flat_S.TS; clear TS_flat_S

fprintf('=== Tier-3 Volve 4-45 Hz Broadband Confluence (smoke=%d, 6.25 m grid) ===\n', smoke);
fprintf('Blind dip reservoir: med=%.2f deg, p90=%.2f deg, med_coh=%.3f\n', ...
    median(abs(rad2deg(phi_eff(wres_mask)))), prctile(abs(rad2deg(phi_eff(wres_mask))), 90), median(c_st(wres_mask)));
fprintf('Shallow anchor [1000-1600m] med t*: DSA=%.2f ms, flat=%.2f ms | Res fold med=%.1f / %d\n', ...
    tstar_anc_dsa*1e3, tstar_anc_flat*1e3, median(sum(double(w_dsa(wres_mask(:), :)), 2)), ns);

% Load Stage 4H cached B=1 (iv) and B=300 (iv3) noise priors (exact apples-to-apples!)
IV_S = load(fullfile(fld45_cache, 'inv.mat'), 'iv', 'iv3');
INV1 = IV_S.iv; INV3 = IV_S.iv3; clear IV_S

% Same B=300 shot signs as fld_run_h.m (for diagnostic DSA-domain B=300 SNR audit)
[~, sgn_b300] = core.make_block_signs(P.shot_of, P.rec_ix, ns, 300, cfg.flip_seed);
sgn_b300 = sgn_b300(:).';

% Stage 4H exact m-grid: [0.6, 0.85, 1.0, 1.2, 1.4]
ms = [0.6, 0.85, 1.0, 1.2, 1.4]; nm = numel(ms);
s70 = single(gain_s(70));
tanc_d_s = single(tstar_anc_dsa);
tanc_f_s = single(tstar_anc_flat);

xmax_full = 1.05 * double(max(max(TS_flat(:)), max(TS_dsa(:)))) / min(ms);
xg_full   = single(linspace(0, xmax_full, 1201)); dxg_full = single(xg_full(2));
xmax_rel  = 1.05 * max(double(max(TS_flat(:)) - tanc_f_s), double(max(TS_dsa(:)) - tanc_d_s)) / min(ms) + 1e-6;
xg_rel    = single(linspace(0, xmax_rel, 1201));  dxg_rel  = single(xg_rel(2));
ntab_m2   = single(1201 - 2);

if smoke
    jf_list = [17, 85, 145];  % 8.0, 25.0, 40.0 Hz
else
    jf_list = 1:nf_all;
end
nf_run = numel(jf_list); f_run = P.f(jf_list);

% 2. Single resumable frequency loop over 165 frequencies (6.25 m grid)
for idx_f = 1:nf_run
    jf = jf_list(idx_f); f = P.f(jf);
    sfile = fullfile(slice_dir, sprintf('slice_%03d.mat', jf));
    if isfile(sfile)
        if mod(idx_f, 15) == 0 || idx_f == nf_run
            fprintf('  [Cache hit] f=%5.2f Hz (%3d/%d)\n', f, idx_f, nf_run);
        end
        continue;
    end
    t0 = tic;
    [L, U, Pp, Qq, R] = lu(helm_op(P.vmig, f, P));
    b = zeros(N, ns); b(sub2ind([N ns], isrc, (1:ns)')) = -P.W(jf) / P.dx^2;
    Us = Qq * (U \ (L \ (Pp * (R \ b))));
    Us = conj(single(Us(lin, :))) * single((2*pi*f)^2);
    b = zeros(N, ns); b(col) = P.D(:, jf);
    V = -(R \ (Pp' * (L' \ (U' \ (Qq' * b))))) / P.dx^2;
    V_in = single(V(lin, :));
    clear L U Pp Qq R b V
    t_helm = toc(t0);

    iv1 = INV1(:, jf);
    iv3 = INV3(:, jf);

    % Operator lookup tables (precomputed for 100% pure-real fast inner loop: Re(H*Cs) = (Pr*Cr - Pi*Ci) + (Dr*Cr - Di*Ci)/(E + iv))
    T1_f100 = uq_tab(f, xg_full, cfg.sig, 100.0);
    P1r_f100 = single(real(T1_f100.P)); P1i_f100 = single(imag(T1_f100.P));
    D1r_f100 = single(real(T1_f100.A - T1_f100.P)); D1i_f100 = single(imag(T1_f100.A - T1_f100.P));
    E1_f100  = single(1 + double(s70)*real(T1_f100.C));

    T1_f25  = uq_tab(f, xg_full, cfg.sig, 25.0);
    P1r_f25 = single(real(T1_f25.P)); P1i_f25 = single(imag(T1_f25.P));
    D1r_f25 = single(real(T1_f25.A - T1_f25.P)); D1i_f25 = single(imag(T1_f25.A - T1_f25.P));
    E1_f25  = single(1 + double(s70)*real(T1_f25.C));

    T0_r25  = uq_tab(f, xg_rel, 0,       25.0);
    P0r_r25 = single(real(T0_r25.P)); P0i_r25 = single(imag(T0_r25.P));
    D0r_r25 = single(real(T0_r25.A - T0_r25.P)); D0i_r25 = single(imag(T0_r25.A - T0_r25.P));
    E0_r25  = single(1 + double(s70)*real(T0_r25.C));
    IQ_r25  = (T0_r25.A + double(s70)*T0_r25.P.*T0_r25.C) ./ (1 + double(s70)*T0_r25.C);
    IQr_r25 = single(real(IQ_r25)); IQi_r25 = single(imag(IQ_r25));

    T1_r25  = uq_tab(f, xg_rel, cfg.sig, 25.0);
    P1r_r25 = single(real(T1_r25.P)); P1i_r25 = single(imag(T1_r25.P));
    D1r_r25 = single(real(T1_r25.A - T1_r25.P)); D1i_r25 = single(imag(T1_r25.A - T1_r25.P));
    E1_r25  = single(1 + double(s70)*real(T1_r25.C));

    T1_ramp = uq_tab(f, xg_rel, cfg.sig, f);
    D1_ramp = single(real(T1_ramp.A) - 1);
    E1_ramp = single(1 + double(s70)*real(T1_ramp.C));

    T1_r100 = uq_tab(f, xg_rel, cfg.sig, 100.0);
    P1r_r100 = single(real(T1_r100.P)); P1i_r100 = single(imag(T1_r100.P));
    D1r_r100 = single(real(T1_r100.A - T1_r100.P)); D1i_r100 = single(imag(T1_r100.A - T1_r100.P));
    E1_r100  = single(1 + double(s70)*real(T1_r100.C));

    RAW_R0_f              = zeros(nin, 1,  'single');
    DSA_R0_f              = zeros(nin, 1,  'single');
    DSA_R0_flip300_f      = zeros(nin, 1,  'single');
    RAW_UQW1b_Unanc_100_f = zeros(nin, nm, 'single');
    DSA_UQW1b_Unanc_100_f = zeros(nin, nm, 'single');
    DSA_UQW1b_Unanc_25_f  = zeros(nin, nm, 'single');
    RAW_UQW1b_Rel_25_f    = zeros(nin, nm, 'single');
    DSA_PH_Rel_25_f       = zeros(nin, nm, 'single');
    DSA_IQ70_Rel_25_f     = zeros(nin, nm, 'single');
    DSA_W1b_Rel_25_f      = zeros(nin, nm, 'single');
    DSA_UQW1_Rel_25_f     = zeros(nin, nm, 'single');
    DSA_UQW1b_Rel_25_f    = zeros(nin, nm, 'single');
    DSA_UQW1b_Rel_Amp_f   = zeros(nin, nm, 'single');
    DSA_UQW1b_Rel_100_f   = zeros(nin, nm, 'single');

    for c = 1:ceil(ns/31)
        sh = (c-1)*31+1 : min(c*31, ns);
        Cs   = Us(:, sh) .* V_in(:, sh);
        Cr   = real(Cs); Ci = imag(Cs); clear Cs
        w_sh = w_dsa(:, sh);
        Cr_d = w_sh .* Cr; Ci_d = w_sh .* Ci;

        RAW_R0_f         = RAW_R0_f         + sum(Cr, 2);
        DSA_R0_f         = DSA_R0_f         + sum(Cr_d, 2);
        DSA_R0_flip300_f = DSA_R0_flip300_f + sum(Cr_d .* sgn_b300(sh), 2);

        ts_f  = TS_flat(:, sh);
        ts_d  = TS_dsa(:, sh);
        ts_fr = max(ts_f - tanc_f_s, single(0));
        ts_dr = max(ts_d - tanc_d_s, single(0));

        for k = 1:nm
            invm = single(1.0 / ms(k));

            % 1) Exact Stage 4H RAW_UQW1b_Unanc_100
            [idx, ww] = core.tab_idx(ts_f * invm, dxg_full, ntab_m2); ow = 1 - ww;
            pc = (P1r_f100(idx).*ow + P1r_f100(idx+1).*ww).*Cr - (P1i_f100(idx).*ow + P1i_f100(idx+1).*ww).*Ci;
            dc = (D1r_f100(idx).*ow + D1r_f100(idx+1).*ww).*Cr - (D1i_f100(idx).*ow + D1i_f100(idx+1).*ww).*Ci;
            E  = E1_f100(idx).*ow + E1_f100(idx+1).*ww;
            RAW_UQW1b_Unanc_100_f(:,k) = RAW_UQW1b_Unanc_100_f(:,k) + sum(pc + dc ./ (E + iv3), 2);

            % 2-3) DSA_UQW1b_Unanc_100 & DSA_UQW1b_Unanc_25 (shared TS_dsa lookup!)
            [idx, ww] = core.tab_idx(ts_d * invm, dxg_full, ntab_m2); ow = 1 - ww;
            pc = (P1r_f100(idx).*ow + P1r_f100(idx+1).*ww).*Cr_d - (P1i_f100(idx).*ow + P1i_f100(idx+1).*ww).*Ci_d;
            dc = (D1r_f100(idx).*ow + D1r_f100(idx+1).*ww).*Cr_d - (D1i_f100(idx).*ow + D1i_f100(idx+1).*ww).*Ci_d;
            E  = E1_f100(idx).*ow + E1_f100(idx+1).*ww;
            DSA_UQW1b_Unanc_100_f(:,k) = DSA_UQW1b_Unanc_100_f(:,k) + sum(pc + dc ./ (E + iv3), 2);
            pc = (P1r_f25(idx).*ow + P1r_f25(idx+1).*ww).*Cr_d - (P1i_f25(idx).*ow + P1i_f25(idx+1).*ww).*Ci_d;
            dc = (D1r_f25(idx).*ow + D1r_f25(idx+1).*ww).*Cr_d - (D1i_f25(idx).*ow + D1i_f25(idx+1).*ww).*Ci_d;
            E  = E1_f25(idx).*ow + E1_f25(idx+1).*ww;
            DSA_UQW1b_Unanc_25_f(:,k)  = DSA_UQW1b_Unanc_25_f(:,k)  + sum(pc + dc ./ (E + iv3), 2);

            % 4) RAW_UQW1b_Rel_25
            [idx, ww] = core.tab_idx(ts_fr * invm, dxg_rel, ntab_m2); ow = 1 - ww;
            pc = (P1r_r25(idx).*ow + P1r_r25(idx+1).*ww).*Cr - (P1i_r25(idx).*ow + P1i_r25(idx+1).*ww).*Ci;
            dc = (D1r_r25(idx).*ow + D1r_r25(idx+1).*ww).*Cr - (D1i_r25(idx).*ow + D1i_r25(idx+1).*ww).*Ci;
            E  = E1_r25(idx).*ow + E1_r25(idx+1).*ww;
            RAW_UQW1b_Rel_25_f(:,k) = RAW_UQW1b_Rel_25_f(:,k) + sum(pc + dc ./ (E + iv3), 2);

            % 5-11) Shared idx for TS_dsa_rel: 100% pure real vectorized arithmetic!
            [idx, ww] = core.tab_idx(ts_dr * invm, dxg_rel, ntab_m2); ow = 1 - ww;
            p0c = (P0r_r25(idx).*ow + P0r_r25(idx+1).*ww).*Cr_d - (P0i_r25(idx).*ow + P0i_r25(idx+1).*ww).*Ci_d;
            d0c = (D0r_r25(idx).*ow + D0r_r25(idx+1).*ww).*Cr_d - (D0i_r25(idx).*ow + D0i_r25(idx+1).*ww).*Ci_d;
            E0  = E0_r25(idx).*ow + E0_r25(idx+1).*ww;
            iqc = (IQr_r25(idx).*ow + IQr_r25(idx+1).*ww).*Cr_d - (IQi_r25(idx).*ow + IQi_r25(idx+1).*ww).*Ci_d;

            p1c = (P1r_r25(idx).*ow + P1r_r25(idx+1).*ww).*Cr_d - (P1i_r25(idx).*ow + P1i_r25(idx+1).*ww).*Ci_d;
            d1c = (D1r_r25(idx).*ow + D1r_r25(idx+1).*ww).*Cr_d - (D1i_r25(idx).*ow + D1i_r25(idx+1).*ww).*Ci_d;
            E1  = E1_r25(idx).*ow + E1_r25(idx+1).*ww;

            DSA_PH_Rel_25_f(:,k)    = DSA_PH_Rel_25_f(:,k)    + sum(p0c, 2);
            DSA_IQ70_Rel_25_f(:,k)  = DSA_IQ70_Rel_25_f(:,k)  + sum(iqc, 2);
            DSA_W1b_Rel_25_f(:,k)   = DSA_W1b_Rel_25_f(:,k)   + sum(p0c + d0c ./ (E0 + iv3), 2);
            DSA_UQW1_Rel_25_f(:,k)  = DSA_UQW1_Rel_25_f(:,k)  + sum(p1c + d1c ./ (E1 + iv1), 2);
            DSA_UQW1b_Rel_25_f(:,k) = DSA_UQW1b_Rel_25_f(:,k) + sum(p1c + d1c ./ (E1 + iv3), 2);

            Da = D1_ramp(idx).*ow + D1_ramp(idx+1).*ww; Ea = E1_ramp(idx).*ow + E1_ramp(idx+1).*ww;
            DSA_UQW1b_Rel_Amp_f(:,k) = DSA_UQW1b_Rel_Amp_f(:,k) + sum((1 + Da ./ (Ea + iv3)) .* Cr_d, 2);

            phc = (P1r_r100(idx).*ow + P1r_r100(idx+1).*ww).*Cr_d - (P1i_r100(idx).*ow + P1i_r100(idx+1).*ww).*Ci_d;
            dhc = (D1r_r100(idx).*ow + D1r_r100(idx+1).*ww).*Cr_d - (D1i_r100(idx).*ow + D1i_r100(idx+1).*ww).*Ci_d;
            Eh  = E1_r100(idx).*ow + E1_r100(idx+1).*ww;
            DSA_UQW1b_Rel_100_f(:,k) = DSA_UQW1b_Rel_100_f(:,k) + sum(phc + dhc ./ (Eh + iv3), 2);
        end
    end
    clear Us V_in Cr Ci Cr_d Ci_d

    save(sfile, 'f', 'RAW_R0_f', 'DSA_R0_f', 'DSA_R0_flip300_f', ...
        'RAW_UQW1b_Unanc_100_f', 'DSA_UQW1b_Unanc_100_f', 'DSA_UQW1b_Unanc_25_f', ...
        'RAW_UQW1b_Rel_25_f', 'DSA_PH_Rel_25_f', 'DSA_IQ70_Rel_25_f', ...
        'DSA_W1b_Rel_25_f', 'DSA_UQW1_Rel_25_f', 'DSA_UQW1b_Rel_25_f', ...
        'DSA_UQW1b_Rel_Amp_f', 'DSA_UQW1b_Rel_100_f', '-v6');
    fprintf('  f=%5.2f Hz (%3d/%d): total=%.1fs (helm=%.1fs, op+io=%.1fs)\n', ...
        f, idx_f, nf_run, toc(t0), t_helm, toc(t0) - t_helm);
end

clear w_dsa TS_flat TS_dsa

% 3. Blind Shallow Anchor PSAE [1000-1600m] over [6, 38] Hz + Broad 8-Hz EnvAnc + DSA-B300 SNR audit
S_anc_raw = zeros(1, nf_run); S_anc_dsa = zeros(1, nf_run);
I_dsa_all = zeros(nin, nf_run, 'single'); I_dsaf_all = I_dsa_all;
am = anc_mask(:);
for idx_f = 1:nf_run
    jf = jf_list(idx_f);
    D = load(fullfile(slice_dir, sprintf('slice_%03d.mat', jf)), 'RAW_R0_f', 'DSA_R0_f', 'DSA_R0_flip300_f');
    S_anc_raw(idx_f) = rms(double(D.RAW_R0_f(am)));
    S_anc_dsa(idx_f) = rms(double(D.DSA_R0_f(am)));
    I_dsa_all(:, idx_f)  = D.DSA_R0_f;
    I_dsaf_all(:, idx_f) = D.DSA_R0_flip300_f;
end
if ~smoke
    INV3_DSA = core.band_inv_snr(I_dsa_all, I_dsaf_all, f_run, P.nzi, P.nxi, P.dx, cfg.snr_band_hz, cfg.snr_box_m, cfg.snr_cap);
    snr_med_b300_dsa = 1 ./ median(double(INV3_DSA(wres_mask(:), :)), 1);
    clear INV3_DSA
else
    snr_med_b300_dsa = zeros(1, nf_run);
end
clear I_dsa_all I_dsaf_all

f_mid = 25.0;
[W_eq_raw, beta_anc_raw, ~] = core.fit_psae(f_run, S_anc_raw, [6.0, 38.0], f_mid);
[W_eq_dsa, beta_anc_dsa, ~] = core.fit_psae(f_run, S_anc_dsa, [6.0, 38.0], f_mid);
W_env_dsa = core.fit_env_anc(f_run, S_anc_dsa, P.W(jf_list), 8.0, f_mid);

fprintf('\n=== Blind 4-45 Hz PSAE Shallow [1000-1600m] Slopes ([6-38 Hz]) ===\n');
fprintf('  RAW beta_anc = %+.5f /Hz (eff_t* = %.2f ms)\n', beta_anc_raw, -beta_anc_raw/pi*1e3);
fprintf('  DSA beta_anc = %+.5f /Hz (eff_t* = %.2f ms)\n', beta_anc_dsa, -beta_anc_dsa/pi*1e3);

% 4. Stack final broadband 2D images
mkeys = {'RAW_R0', 'DSA_R0', 'DSA_R0_Eq', ...
         'RAW_UQW1b_Unanc_100', 'DSA_UQW1b_Unanc_100', 'DSA_UQW1b_Unanc_25', 'RAW_UQW1b_PAnc_25', ...
         'DSA_PH_PAnc_25', 'DSA_IQ70_PAnc_25', 'DSA_W1b_PAnc_25', 'DSA_UQW1_PAnc_25', ...
         'DSA_UQW1b_PAnc_25', 'DSA_UQW1b_PAnc_Amp', 'DSA_UQW1b_PAnc_100', ...
         'DSA_UQW1b_EnvAnc_25', 'DSA_UQW1b_EnvAnc_Amp'};
Img = struct();
for ik = 1:numel(mkeys)
    nc = 1 + (nm-1)*~ismember(mkeys{ik}, {'RAW_R0','DSA_R0','DSA_R0_Eq'});
    Img.(mkeys{ik}) = zeros(nin, nc);
end
spec_res = struct();
for ik = 1:numel(mkeys), spec_res.(mkeys{ik}) = zeros(1, nf_run); end

rm = wres_mask(:);
for idx_f = 1:nf_run
    jf = jf_list(idx_f);
    D = load(fullfile(slice_dir, sprintf('slice_%03d.mat', jf)));
    wr = W_eq_raw(idx_f); wd = W_eq_dsa(idx_f); we = W_env_dsa(idx_f);

    Img.RAW_R0               = Img.RAW_R0               + double(D.RAW_R0_f);
    Img.DSA_R0               = Img.DSA_R0               + double(D.DSA_R0_f);
    Img.DSA_R0_Eq            = Img.DSA_R0_Eq            + wd * double(D.DSA_R0_f);
    Img.RAW_UQW1b_Unanc_100  = Img.RAW_UQW1b_Unanc_100  + double(D.RAW_UQW1b_Unanc_100_f);
    Img.DSA_UQW1b_Unanc_100  = Img.DSA_UQW1b_Unanc_100  + double(D.DSA_UQW1b_Unanc_100_f);
    Img.DSA_UQW1b_Unanc_25   = Img.DSA_UQW1b_Unanc_25   + double(D.DSA_UQW1b_Unanc_25_f);
    Img.RAW_UQW1b_PAnc_25    = Img.RAW_UQW1b_PAnc_25    + wr * double(D.RAW_UQW1b_Rel_25_f);
    Img.DSA_PH_PAnc_25       = Img.DSA_PH_PAnc_25       + wd * double(D.DSA_PH_Rel_25_f);
    Img.DSA_IQ70_PAnc_25     = Img.DSA_IQ70_PAnc_25     + wd * double(D.DSA_IQ70_Rel_25_f);
    Img.DSA_W1b_PAnc_25      = Img.DSA_W1b_PAnc_25      + wd * double(D.DSA_W1b_Rel_25_f);
    Img.DSA_UQW1_PAnc_25     = Img.DSA_UQW1_PAnc_25     + wd * double(D.DSA_UQW1_Rel_25_f);
    Img.DSA_UQW1b_PAnc_25    = Img.DSA_UQW1b_PAnc_25    + wd * double(D.DSA_UQW1b_Rel_25_f);
    Img.DSA_UQW1b_PAnc_Amp   = Img.DSA_UQW1b_PAnc_Amp   + wd * double(D.DSA_UQW1b_Rel_Amp_f);
    Img.DSA_UQW1b_PAnc_100   = Img.DSA_UQW1b_PAnc_100   + wd * double(D.DSA_UQW1b_Rel_100_f);
    Img.DSA_UQW1b_EnvAnc_25  = Img.DSA_UQW1b_EnvAnc_25  + we * double(D.DSA_UQW1b_Rel_25_f);
    Img.DSA_UQW1b_EnvAnc_Amp = Img.DSA_UQW1b_EnvAnc_Amp + we * double(D.DSA_UQW1b_Rel_Amp_f);

    spec_res.RAW_R0(idx_f)               = rms(double(D.RAW_R0_f(rm)));
    spec_res.DSA_R0(idx_f)               = rms(double(D.DSA_R0_f(rm)));
    spec_res.DSA_R0_Eq(idx_f)            = rms(wd * double(D.DSA_R0_f(rm)));
    spec_res.RAW_UQW1b_Unanc_100(idx_f)  = rms(double(D.RAW_UQW1b_Unanc_100_f(rm, 3)));
    spec_res.DSA_UQW1b_Unanc_100(idx_f)  = rms(double(D.DSA_UQW1b_Unanc_100_f(rm, 3)));
    spec_res.DSA_UQW1b_Unanc_25(idx_f)   = rms(double(D.DSA_UQW1b_Unanc_25_f(rm, 3)));
    spec_res.RAW_UQW1b_PAnc_25(idx_f)    = rms(wr * double(D.RAW_UQW1b_Rel_25_f(rm, 3)));
    spec_res.DSA_PH_PAnc_25(idx_f)       = rms(wd * double(D.DSA_PH_Rel_25_f(rm, 3)));
    spec_res.DSA_IQ70_PAnc_25(idx_f)     = rms(wd * double(D.DSA_IQ70_Rel_25_f(rm, 3)));
    spec_res.DSA_W1b_PAnc_25(idx_f)      = rms(wd * double(D.DSA_W1b_Rel_25_f(rm, 3)));
    spec_res.DSA_UQW1_PAnc_25(idx_f)     = rms(wd * double(D.DSA_UQW1_Rel_25_f(rm, 3)));
    spec_res.DSA_UQW1b_PAnc_25(idx_f)    = rms(wd * double(D.DSA_UQW1b_Rel_25_f(rm, 3)));
    spec_res.DSA_UQW1b_PAnc_Amp(idx_f)   = rms(wd * double(D.DSA_UQW1b_Rel_Amp_f(rm, 3)));
    spec_res.DSA_UQW1b_PAnc_100(idx_f)   = rms(wd * double(D.DSA_UQW1b_Rel_100_f(rm, 3)));
    spec_res.DSA_UQW1b_EnvAnc_25(idx_f)  = rms(we * double(D.DSA_UQW1b_Rel_25_f(rm, 3)));
    spec_res.DSA_UQW1b_EnvAnc_Amp(idx_f) = rms(we * double(D.DSA_UQW1b_Rel_Amp_f(rm, 3)));
end

r0_repro = NaN; uqw1b_repro = NaN;
if ~smoke
    r0_repro    = norm(Img.RAW_R0 - R0_stage4h) / max(norm(R0_stage4h), 1e-30);
    uqw1b_repro = norm(Img.RAW_UQW1b_Unanc_100(:,3) - UQW1b_stage4h_m1) / max(norm(UQW1b_stage4h_m1), 1e-30);
    fprintf('Exact Stage 4H repro: RAW_R0 rel diff = %.2e | RAW_UQW1b_100(m=1) rel diff = %.2e\n', r0_repro, uqw1b_repro);
end

% 5. Exact Stage 4H Evaluation (tie_fix 8-40 Hz, tie_r0, Zero-Lag, Best-Lag, F-14, Smile_Frac, Res_Slope)
wz_res = (z_in >= cfg.res_z_m(1) & z_in <= cfg.res_z_m(2));
wx_res = (x_in >= cfg.res_x_m(1) & x_in <= cfg.res_x_m(2));
wz_mid = (z_in >= cfg.ovb_z_m(1) & z_in <= cfg.ovb_z_m(2));
vm_in  = P.vmig(P.in{1}, P.in{2});
vr     = median(vm_in(wz_res, wx_res), 'all');
kfix   = 2 * [8.0, 40.0] / vr * 1000.0;

A_r0_prep = core.prep_image(Img.RAW_R0, P.nzi, P.nxi, z_in, P.dx, cfg.bg_sigma_m, cfg.z_mute_m);
[kr0, ~]  = core.res_band(A_r0_prep(wz_res, wx_res), P.dx, 512);

w5  = core.load_well(fullfile(well_dir, 'well_F-5.bin'),  z_in, x_in);
w14 = core.load_well(fullfile(well_dir, 'well_F-14.bin'), z_in, x_in);
rb5_fix  = core.bp_z(w5.rw,  w5.dzw,  kfix(1), kfix(2));
rb5_r0   = core.bp_z(w5.rw,  w5.dzw,  kr0(1),  kr0(2));
rb14_fix = core.bp_z(w14.rw, w14.dzw, kfix(1), kfix(2));

fm = (f_run >= 6.0 & f_run <= 38.0);
res = struct('dataset', 'volve_4_45hz', 'smoke', smoke, 'n_freq', nf_run, 'f', f_run, 'm', ms, ...
             'v_r', vr, 'kfix', kfix, 'kr0', kr0, ...
             'repro_R0_relL2', r0_repro, 'repro_UQW1b_relL2', uqw1b_repro, ...
             'tstar_anc_dsa_ms', tstar_anc_dsa*1e3, 'tstar_anc_flat_ms', tstar_anc_flat*1e3, ...
             'beta_anc_raw', beta_anc_raw, 'beta_anc_dsa', beta_anc_dsa, ...
             'S_anc_raw', S_anc_raw, 'S_anc_dsa', S_anc_dsa, ...
             'snr_med_b1', 1./median(double(INV1(rm, jf_list)), 1), ...
             'snr_med_b300_raw', 1./median(double(INV3(rm, jf_list)), 1), ...
             'snr_med_b300_dsa', snr_med_b300_dsa, ...
             'spec_res', spec_res, 'res_slope', struct(), 'methods', struct());

for ik = 1:numel(mkeys)
    k = mkeys{ik}; sv = spec_res.(k);
    pf = polyfit(f_run(fm), log(max(sv(fm), 1e-30)), 1);
    res.res_slope.(k) = pf(1);
end

fprintf('\n=== TIER-3 VOLVE 4-45 Hz BROADBAND CONFLUENCE (kfix=[%.2f-%.2f] cyc/km) ===\n', kfix(1), kfix(2));
for ik = 1:numel(mkeys)
    k = mkeys{ik}; M = Img.(k); nc = size(M, 2);
    f5_fix = zeros(1,nc); f5_r0 = zeros(1,nc); f5_l0 = zeros(1,nc); f5_lag = zeros(1,nc);
    f14_fix = zeros(1,nc); cen = zeros(1,nc); bw10 = zeros(nc,2); smf = zeros(1,nc); smf25 = zeros(1,nc);
    for c = 1:nc
        A = core.prep_image(M(:,c), P.nzi, P.nxi, z_in, P.dx, cfg.bg_sigma_m, cfg.z_mute_m);
        [f5_fix(c), f5_lag(c), f5_l0(c)] = core.well_tie_full(A, w5,  rb5_fix,  P.dx, kfix(1), kfix(2), 20.0, 30.0);
        [f5_r0(c),  ~,         ~       ] = core.well_tie_full(A, w5,  rb5_r0,   P.dx, kr0(1),  kr0(2),  20.0, 30.0);
        [f14_fix(c), ~,        ~       ] = core.well_tie_full(A, w14, rb14_fix, P.dx, kfix(1), kfix(2), 20.0, 30.0);
        [bw10(c,:), cen(c)] = core.res_band(A(wz_res, wx_res), P.dx, 512);
        smf(c)   = core.steep_smile_frac(A(wz_mid, wx_res), P.dx, 25.0, 3.0, 40.0);
        smf25(c) = core.steep_smile_frac(A(wz_mid, wx_res), P.dx, 25.0, 3.0, 25.0);
    end
    res.methods.(k) = struct('F5_fix', f5_fix, 'F5_fix_mean', mean(f5_fix), ...
        'F5_r0', f5_r0, 'F5_r0_mean', mean(f5_r0), ...
        'F5_ZeroLag', f5_l0, 'F5_BestLag_m', f5_lag, ...
        'F14_fix', f14_fix, 'F14_fix_mean', mean(f14_fix), ...
        'Res_Cen', cen, 'Res_BW10', bw10, 'Smile_Frac', smf, 'Smile_Frac_25', smf25);
    if nc == 1
        fprintf('  %-22s | fix=%.4f (L0=%.3f,%+3dm) | r0=%.4f | F14=%.3f | Sm=%.3f | sl=%+.4f\n', ...
            k, f5_fix(1), f5_l0(1), round(f5_lag(1)), f5_r0(1), f14_fix(1), smf(1), res.res_slope.(k));
    else
        fprintf('  %-22s | fix=%.4f (%.3f/%.3f/%.3f, L0=%.3f,%+3dm) | r0=%.4f | Sm=%.3f | sl=%+.4f\n', ...
            k, mean(f5_fix), f5_fix(1), f5_fix(3), f5_fix(5), f5_l0(3), round(f5_lag(3)), mean(f5_r0), smf(3), res.res_slope.(k));
    end
end
res.elapsed_s = toc(t_all);

tag = ''; if smoke, tag = '_smoke'; end
fn_json = fullfile(res_dir, ['volve45_unified' tag '.json']);
fid = fopen(fn_json, 'w'); fwrite(fid, jsonencode(res, 'PrettyPrint', true), 'char'); fclose(fid);
fprintf('\nSaved %s (total %.1fs)\n', fn_json, res.elapsed_s);

if ~smoke
    Img_s = struct();
    for ik = 1:numel(mkeys), Img_s.(mkeys{ik}) = single(Img.(mkeys{ik})); end
    save(fullfile(scratch_dir, 'volve45_final_imgs.mat'), 'Img_s', 'res', 'z_in', 'x_in', 'phi_eff', '-v7.3');
    fprintf('Saved final 4-45 Hz 2D images snapshot to %s\n', scratch_dir);
end
end

function P = load_fld_ro(cdir)
M = jsondecode(fileread(fullfile(cdir, 'meta.json')));
P.dx = M.dx; P.npml = M.npml; P.nzi = M.nzi; P.nxi = M.nxi; P.nz = M.nz; P.nx = M.nx;
P.x = M.xmin + ((1:P.nx) - P.npml - 1) * P.dx;
P.z = ((1:P.nz) - P.npml - 1) * P.dx;
P.in = {P.npml + (1:P.nzi), P.npml + (1:P.nxi)};
P.f = M.f(:)'; P.nshot = M.nshot; P.ntrace = M.ntrace;
P.shot_of = double(rdc_bin(cdir, 'shot_of.bin', 'int32'));
P.rec_ix  = double(rdc_bin(cdir, 'rec_ix.bin',  'int32'));
P.rec_iz  = double(rdc_bin(cdir, 'rec_iz.bin',  'int32'));
P.src_ix  = double(rdc_bin(cdir, 'src_ix.bin',  'int32'));
P.src_iz  = double(rdc_bin(cdir, 'src_iz.bin',  'int32'));
d = rdc_bin(cdir, 'D.bin', 'single'); d = reshape(d, 2, numel(P.f), P.ntrace);
P.D = single(squeeze(complex(d(1,:,:), d(2,:,:))).');
v = rdc_bin(cdir, 'vp.bin', 'single'); P.vp = double(reshape(v, P.nx, P.nz).');
s = 50 / P.dx; h = exp(-0.5*((-ceil(3*s):ceil(3*s))/s).^2); h = h/sum(h); m = (numel(h)-1)/2;
xp = [repmat(P.vp(1,:), m, 1); P.vp; repmat(P.vp(end,:), m, 1)];
xp = [repmat(xp(:,1), 1, m), xp, repmat(xp(:,end), 1, m)];
P.vmig = conv2(h, h, xp, 'valid');
P.vmig(P.vp < 1500) = P.vp(P.vp < 1500);
f = P.f; fm = max(f); w = ones(size(f));
a = f < 6; w(a) = 0.5*(1 - cos(pi*(f(a)-4)/2));
a = f > fm-2; w(a) = 0.5*(1 + cos(pi*(f(a)-fm+2)/2));
P.W = w;
P.src_x = P.x(P.src_ix);
end

function v = rdc_bin(cdir, n, t)
fid = fopen(fullfile(cdir, n));
v = fread(fid, Inf, ['*' t]);
fclose(fid);
end
