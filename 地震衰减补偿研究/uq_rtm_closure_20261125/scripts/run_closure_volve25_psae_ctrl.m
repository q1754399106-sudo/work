function run_closure_volve25_psae_ctrl(smoke)
% WP3: Volve 4-25 Hz PSAE Attribution Control Experiment
% Tests 4 shallow/spectral routes under unified same-domain DSA B=300:
%   Route A: (2*pi*f)^2 + Unanchored (W=1) + absolute t*_DSA
%   Route B: Remove (2*pi*f)^2 (zero PSAE!), both absolute t*_DSA & relative dt*_rel
%   Route C: (2*pi*f)^2 + 1-param PSAE + relative dt*_rel
%   Route D: (2*pi*f)^2 + 6-Hz smooth shallow spectral deconvolution + relative dt*_rel
%
% Strictly read-only on D:\ChatGPT\研究\地震衰减补偿研究\uq_rtm_20261117\cache\fld.
% Temporary Cs streamed to D:\tmp_uqrtm_closure_20261125\cs_fld25 (auto-cleaned).

if nargin < 1, smoke = false; end
maxNumCompThreads(10);
t_all = tic;

here = fileparts(mfilename('fullpath'));
root = fileparts(here);
res_dir = fullfile(root, 'results');
if ~isfolder(res_dir), mkdir(res_dir); end

wt_root = fileparts(root);
addpath(fullfile(wt_root, 'uq_rtm_20261117', 'scripts'));
core = uq_rtm_unified_core();
cfg  = core.default_cfg('field_obc');

fld_cache   = 'D:\ChatGPT\研究\地震衰减补偿研究\uq_rtm_20261117\cache\fld';
scratch_dir = 'D:\tmp_uqrtm_closure_20261125';
if smoke
    cs_dir = fullfile(scratch_dir, 'cs_fld25_smoke');
else
    cs_dir = fullfile(scratch_dir, 'cs_fld25');
end
if ~isfolder(cs_dir), mkdir(cs_dir); end

P = load_fld_ro(fld_cache);
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

% Blind dip from Stage-4 cached R0
S1 = load(fullfile(fld_cache, 'state1.mat'), 'I');
R0_cached = double(sum(S1.I, 2)); clear S1
R0_2d = core.prep_image(R0_cached, P.nzi, P.nxi, z_in, P.dx, cfg.bg_sigma_m, cfg.z_mute_m);
[phi_eff, c_st] = core.blind_dip(R0_2d, P.dx, cfg.dip_sigma_m, cfg.dip_cap_deg);

% Build DSA mask & TS_dsa
rx = P.x(P.rec_ix);
xr_lo = arrayfun(@(s) min(rx(P.shot_of == s)), 1:ns);
xr_hi = arrayfun(@(s) max(rx(P.shot_of == s)), 1:ns);
TF = load(fullfile(fld_cache, 'tfields.mat'), 'F');
[w_dsa, TS_dsa] = core.build_dsa_full(Xin(:), Zin(:), P.z(P.src_iz), 90.0, P.src_x, ...
    xr_lo, xr_hi, phi_eff(:), TF.F, Zi(:), Xi(:), P.nz, P.nx, cfg);
clear TF

tstar_anc_dsa = median(double(TS_dsa(anc_mask(:), :)), 'all');
TS_dsa_rel = max(TS_dsa - single(tstar_anc_dsa), single(0));
fprintf('=== WP3 Volve 4-25 Hz PSAE Attribution Control (smoke=%d) ===\n', smoke);
fprintf('Blind dip med=%.2f deg (coh=%.3f) | Shallow t*_anc_dsa=%.2f ms\n', ...
    median(abs(rad2deg(phi_eff(wres_mask)))), median(c_st(wres_mask)), tstar_anc_dsa*1e3);

[~, sgn_b300] = core.make_block_signs(P.shot_of, P.rec_ix, ns, 300, cfg.flip_seed);
sgn_b300 = sgn_b300(:).';

if smoke
    jf_list = [9, 43, 73];
else
    jf_list = 1:nf_all;
end
nf_run = numel(jf_list); f_run = P.f(jf_list);
f_mid = 14.5;

% Pass 1: Single Helmholtz pass -> stream Cs + accumulate RAW, DSA, and NoW2 images
I_raw_f  = zeros(nin, nf_run, 'single');
I_dsa_f  = zeros(nin, nf_run, 'single'); I_dsa_fb = I_dsa_f;
I_nw2_f  = zeros(nin, nf_run, 'single'); I_nw2_fb = I_nw2_f;

p1_meta = fullfile(cs_dir, 'pass1_imgs.mat');
done_p1 = false(1, nf_run);
if isfile(p1_meta)
    S = load(p1_meta);
    if numel(S.done_p1) == nf_run
        I_raw_f = S.I_raw_f; I_dsa_f = S.I_dsa_f; I_dsa_fb = S.I_dsa_fb;
        I_nw2_f = S.I_nw2_f; I_nw2_fb = S.I_nw2_fb; done_p1 = S.done_p1;
    end
end

for idx_f = 1:nf_run
    jf = jf_list(idx_f); f = P.f(jf);
    cfile = fullfile(cs_dir, sprintf('cs_%03d.bin', jf));
    if done_p1(idx_f) && isfile(cfile)
        if mod(idx_f, 15) == 0 || idx_f == nf_run
            fprintf('  [Pass1 Cache] f=%5.2f Hz (%2d/%d)\n', f, idx_f, nf_run);
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
    Cs = Us .* single(V(lin, :));
    clear L U Pp Qq R b V Us

    Cr = real(Cs);
    I_raw_f(:, idx_f) = sum(Cr, 2);
    Cr_d = w_dsa .* Cr;
    I_dsa_f(:, idx_f)  = sum(Cr_d, 2);
    I_dsa_fb(:, idx_f) = sum(Cr_d .* sgn_b300, 2);
    nw2_sc = single((f_mid / f)^2);
    I_nw2_f(:, idx_f)  = nw2_sc * I_dsa_f(:, idx_f);
    I_nw2_fb(:, idx_f) = nw2_sc * I_dsa_fb(:, idx_f);
    clear Cr Cr_d

    fid = fopen(cfile, 'w');
    fwrite(fid, [real(Cs(:)), imag(Cs(:))].', 'single');
    fclose(fid);
    clear Cs

    done_p1(idx_f) = true;
    if mod(idx_f, 5) == 0 || idx_f == nf_run
        save(p1_meta, 'I_raw_f', 'I_dsa_f', 'I_dsa_fb', 'I_nw2_f', 'I_nw2_fb', 'done_p1', '-v7.3');
    end
    fprintf('  Pass1 f=%5.2f Hz (%2d/%d): %.1fs\n', f, idx_f, nf_run, toc(t0));
end

% Exact B=300 inverse SNR in BOTH (2*pi*f)^2 DSA domain and NoW2 DSA domain
if smoke
    INV3_DSA = ones(nin, nf_run, 'single') * single(0.2);
    INV3_NW2 = INV3_DSA;
    r0_rel_err = NaN;
else
    INV3_DSA = core.band_inv_snr(I_dsa_f, I_dsa_fb, f_run, P.nzi, P.nxi, P.dx, cfg.snr_band_hz, cfg.snr_box_m, cfg.snr_cap);
    INV3_NW2 = core.band_inv_snr(I_nw2_f, I_nw2_fb, f_run, P.nzi, P.nxi, P.dx, cfg.snr_band_hz, cfg.snr_box_m, cfg.snr_cap);
    r0_rel_err = norm(double(sum(I_raw_f, 2)) - R0_cached) / max(norm(R0_cached), 1e-30);
    fprintf('Exact RAW_R0 repro vs state1.mat: rel L2 diff = %.2e\n', r0_rel_err);
end

% Blind shallow anchor PSAE & 6-Hz smooth spectral decon [1000-1600m]
S_anc_dsa = arrayfun(@(q) rms(double(I_dsa_f(anc_mask(:), q))), 1:nf_run);
S_anc_nw2 = arrayfun(@(q) rms(double(I_nw2_f(anc_mask(:), q))), 1:nf_run);
[W_eq_dsa, beta_anc_dsa, ~] = core.fit_psae(f_run, S_anc_dsa, [6.0, 22.0], f_mid);
[~, beta_anc_nw2, ~]        = core.fit_psae(f_run, S_anc_nw2, [6.0, 22.0], f_mid);
W_decon_dsa = core.fit_env_anc(f_run, S_anc_dsa, P.W(jf_list), 6.0, f_mid);
W_nw2 = (f_mid ./ f_run).^2;

fprintf('Shallow [1000-1600m] slope [6-22 Hz]: with (2*pi*f)^2 beta=%+.5f /Hz | NoW2 beta=%+.5f /Hz\n', ...
    beta_anc_dsa, beta_anc_nw2);

% Pass 2: Evaluate all 14 control & route methods
ms = [0.5, 0.7, 1.0, 1.4, 2.0]; nm = numel(ms);
s70 = single(gain_s(70));
xmax_full = 1.05 * double(max(TS_dsa(:))) / min(ms);
xg_full   = single(linspace(0, xmax_full, 1201)); dxg_full = single(xg_full(2));
xmax_rel  = 1.05 * double(max(TS_dsa_rel(:))) / min(ms) + 1e-6;
xg_rel    = single(linspace(0, xmax_rel, 1201));  dxg_rel  = single(xg_rel(2));
ntab_m2   = single(1201 - 2);

mkeys = {'RAW_R0', 'DSA_R0', 'DSA_R0_Eq', 'DSA_R0_NoW2', 'DSA_R0_Decon', ...
         'DSA_UQW1b_Unanc_25', ...
         'DSA_UQW1b_NoW2_Abs_25', 'DSA_UQW1b_NoW2_Rel_25', ...
         'DSA_PH_NoW2_Rel_25', 'DSA_W1b_NoW2_Rel_25', ...
         'DSA_PH_PAnc_25', 'DSA_IQ70_PAnc_25', 'DSA_W1b_PAnc_25', ...
         'DSA_UQW1b_PAnc_25', 'DSA_UQW1b_Decon_25'};
r0_only = {'RAW_R0', 'DSA_R0', 'DSA_R0_Eq', 'DSA_R0_NoW2', 'DSA_R0_Decon'};

Img = struct(); spec_res = struct();
for ik = 1:numel(mkeys)
    nc = 1 + (nm-1)*~ismember(mkeys{ik}, r0_only);
    Img.(mkeys{ik}) = zeros(nin, nc);
    spec_res.(mkeys{ik}) = zeros(1, nf_run);
end

rm = wres_mask(:);
for idx_f = 1:nf_run
    jf = jf_list(idx_f); f = P.f(jf);
    wd = single(W_eq_dsa(idx_f));
    wn = single(W_nw2(idx_f));
    wdec = single(W_decon_dsa(idx_f));

    fid = fopen(fullfile(cs_dir, sprintf('cs_%03d.bin', jf)));
    a = fread(fid, [2, nin*ns], '*single'); fclose(fid);
    Cs = reshape(complex(a(1,:), a(2,:)), nin, ns); clear a

    iv3d = INV3_DSA(:, idx_f);
    iv3n = INV3_NW2(:, idx_f);

    T1_f25  = uq_tab(f, xg_full, cfg.sig, 25.0);
    P1r_f25 = single(real(T1_f25.P)); P1i_f25 = single(imag(T1_f25.P));
    D1r_f25 = single(real(T1_f25.A - T1_f25.P)); D1i_f25 = single(imag(T1_f25.A - T1_f25.P));
    E1_f25  = single(1 + double(s70)*real(T1_f25.C));

    T0_r25  = uq_tab(f, xg_rel, 0, 25.0);
    P0r_r25 = single(real(T0_r25.P)); P0i_r25 = single(imag(T0_r25.P));
    D0r_r25 = single(real(T0_r25.A - T0_r25.P)); D0i_r25 = single(imag(T0_r25.A - T0_r25.P));
    E0_r25  = single(1 + double(s70)*real(T0_r25.C));
    IQ_r25  = (T0_r25.A + double(s70)*T0_r25.P.*T0_r25.C) ./ (1 + double(s70)*T0_r25.C);
    IQr_r25 = single(real(IQ_r25)); IQi_r25 = single(imag(IQ_r25));

    T1_r25  = uq_tab(f, xg_rel, cfg.sig, 25.0);
    P1r_r25 = single(real(T1_r25.P)); P1i_r25 = single(imag(T1_r25.P));
    D1r_r25 = single(real(T1_r25.A - T1_r25.P)); D1i_r25 = single(imag(T1_r25.A - T1_r25.P));
    E1_r25  = single(1 + double(s70)*real(T1_r25.C));

    acc = struct();
    for ik = 1:numel(mkeys)
        acc.(mkeys{ik}) = zeros(nin, size(Img.(mkeys{ik}), 2), 'single');
    end
    acc.RAW_R0       = I_raw_f(:, idx_f);
    acc.DSA_R0       = I_dsa_f(:, idx_f);
    acc.DSA_R0_Eq    = wd   * I_dsa_f(:, idx_f);
    acc.DSA_R0_NoW2  = wn   * I_dsa_f(:, idx_f);
    acc.DSA_R0_Decon = wdec * I_dsa_f(:, idx_f);

    for c = 1:ceil(ns/31)
        sh = (c-1)*31+1 : min(c*31, ns);
        w_sh = w_dsa(:, sh);
        Cr_d = w_sh .* real(Cs(:, sh));
        Ci_d = w_sh .* imag(Cs(:, sh));
        ts_d  = TS_dsa(:, sh);
        ts_dr = TS_dsa_rel(:, sh);

        for k = 1:nm
            invm = single(1.0 / ms(k));

            % Route A & B1 (absolute TS_dsa)
            [idx, ww] = core.tab_idx(ts_d * invm, dxg_full, ntab_m2); ow = 1 - ww;
            pc = (P1r_f25(idx).*ow + P1r_f25(idx+1).*ww).*Cr_d - (P1i_f25(idx).*ow + P1i_f25(idx+1).*ww).*Ci_d;
            dc = (D1r_f25(idx).*ow + D1r_f25(idx+1).*ww).*Cr_d - (D1i_f25(idx).*ow + D1i_f25(idx+1).*ww).*Ci_d;
            E  = E1_f25(idx).*ow + E1_f25(idx+1).*ww;
            acc.DSA_UQW1b_Unanc_25(:,k)    = acc.DSA_UQW1b_Unanc_25(:,k)    + sum(pc + dc ./ (E + iv3d), 2);
            acc.DSA_UQW1b_NoW2_Abs_25(:,k) = acc.DSA_UQW1b_NoW2_Abs_25(:,k) + wn * sum(pc + dc ./ (E + iv3n), 2);

            % Route B2, C, D (relative TS_dsa_rel)
            [idx, ww] = core.tab_idx(ts_dr * invm, dxg_rel, ntab_m2); ow = 1 - ww;
            p0c = (P0r_r25(idx).*ow + P0r_r25(idx+1).*ww).*Cr_d - (P0i_r25(idx).*ow + P0i_r25(idx+1).*ww).*Ci_d;
            d0c = (D0r_r25(idx).*ow + D0r_r25(idx+1).*ww).*Cr_d - (D0i_r25(idx).*ow + D0i_r25(idx+1).*ww).*Ci_d;
            E0  = E0_r25(idx).*ow + E0_r25(idx+1).*ww;
            iqc = (IQr_r25(idx).*ow + IQr_r25(idx+1).*ww).*Cr_d - (IQi_r25(idx).*ow + IQi_r25(idx+1).*ww).*Ci_d;

            p1c = (P1r_r25(idx).*ow + P1r_r25(idx+1).*ww).*Cr_d - (P1i_r25(idx).*ow + P1i_r25(idx+1).*ww).*Ci_d;
            d1c = (D1r_r25(idx).*ow + D1r_r25(idx+1).*ww).*Cr_d - (D1i_r25(idx).*ow + D1i_r25(idx+1).*ww).*Ci_d;
            E1  = E1_r25(idx).*ow + E1_r25(idx+1).*ww;

            uqb_d = sum(p1c + d1c ./ (E1 + iv3d), 2);
            uqb_n = sum(p1c + d1c ./ (E1 + iv3n), 2);
            w1b_d = sum(p0c + d0c ./ (E0 + iv3d), 2);
            w1b_n = sum(p0c + d0c ./ (E0 + iv3n), 2);
            ph_s  = sum(p0c, 2);

            acc.DSA_UQW1b_NoW2_Rel_25(:,k) = acc.DSA_UQW1b_NoW2_Rel_25(:,k) + wn * uqb_n;
            acc.DSA_PH_NoW2_Rel_25(:,k)    = acc.DSA_PH_NoW2_Rel_25(:,k)    + wn * ph_s;
            acc.DSA_W1b_NoW2_Rel_25(:,k)   = acc.DSA_W1b_NoW2_Rel_25(:,k)   + wn * w1b_n;

            acc.DSA_PH_PAnc_25(:,k)        = acc.DSA_PH_PAnc_25(:,k)        + wd * ph_s;
            acc.DSA_IQ70_PAnc_25(:,k)      = acc.DSA_IQ70_PAnc_25(:,k)      + wd * sum(iqc, 2);
            acc.DSA_W1b_PAnc_25(:,k)       = acc.DSA_W1b_PAnc_25(:,k)       + wd * w1b_d;
            acc.DSA_UQW1b_PAnc_25(:,k)     = acc.DSA_UQW1b_PAnc_25(:,k)     + wd * uqb_d;
            acc.DSA_UQW1b_Decon_25(:,k)    = acc.DSA_UQW1b_Decon_25(:,k)    + wdec * uqb_d;
        end
    end
    clear Cs

    for ik = 1:numel(mkeys)
        k = mkeys{ik};
        Img.(k) = Img.(k) + double(acc.(k));
        c_m1 = min(3, size(acc.(k), 2));
        spec_res.(k)(idx_f) = rms(double(acc.(k)(rm, c_m1)));
    end
    if mod(idx_f, 15) == 0 || idx_f == nf_run
        fprintf('  Pass2 f=%5.2f Hz (%2d/%d)\n', f, idx_f, nf_run);
    end
end

% Blind evaluation on F-5 & F-14 (with full L0 & BestLag for BOTH wells!)
w5  = core.load_well(fullfile(fld_cache, 'well_F-5.bin'),  z_in, x_in);
w14 = core.load_well(fullfile(fld_cache, 'well_F-14.bin'), z_in, x_in);
lo_r0 = 3.515625; hi_r0 = 12.890625; hi_br = 16.0;
rb5_r0  = core.bp_z(w5.rw,  w5.dzw,  lo_r0, hi_r0);
rb5_br  = core.bp_z(w5.rw,  w5.dzw,  lo_r0, hi_br);
rb14_r0 = core.bp_z(w14.rw, w14.dzw, lo_r0, hi_r0);

wz_res = (z_in >= cfg.res_z_m(1) & z_in <= cfg.res_z_m(2));
wx_res = (x_in >= cfg.res_x_m(1) & x_in <= cfg.res_x_m(2));
wz_mid = (z_in >= cfg.ovb_z_m(1) & z_in <= cfg.ovb_z_m(2));

fm = (f_run >= 6.0 & f_run <= 22.0);
res = struct('dataset', 'volve_4_25hz_psae_ctrl', 'smoke', smoke, 'n_freq', nf_run, 'f', f_run, 'm', ms, ...
             'repro_R0_relL2', r0_rel_err, ...
             'tstar_anc_dsa_ms', tstar_anc_dsa*1e3, ...
             'beta_anc_dsa', beta_anc_dsa, 'beta_anc_nw2', beta_anc_nw2, ...
             'S_anc_dsa', S_anc_dsa, 'S_anc_nw2', S_anc_nw2, ...
             'snr_med_b300_dsa', 1./median(double(INV3_DSA(rm,:)), 1), ...
             'snr_med_b300_nw2', 1./median(double(INV3_NW2(rm,:)), 1), ...
             'spec_res', spec_res, 'res_slope', struct(), 'methods', struct());

for ik = 1:numel(mkeys)
    k = mkeys{ik}; sv = spec_res.(k);
    pf = polyfit(f_run(fm), log(max(sv(fm), 1e-30)), 1);
    res.res_slope.(k) = pf(1);
end

fprintf('\n=== WP3 VOLVE 4-25 Hz PSAE vs NoW2 vs Decon RESULTS ===\n');
for ik = 1:numel(mkeys)
    k = mkeys{ik}; M = Img.(k); nc = size(M, 2);
    f5_r0 = zeros(1,nc); f5_br = zeros(1,nc); f5_l0 = zeros(1,nc); f5_lag = zeros(1,nc);
    f14_r0 = zeros(1,nc); f14_l0 = zeros(1,nc); f14_lag = zeros(1,nc);
    cen = zeros(1,nc); bw10 = zeros(nc,2); smf = zeros(1,nc);
    for c = 1:nc
        A = core.prep_image(M(:,c), P.nzi, P.nxi, z_in, P.dx, cfg.bg_sigma_m, cfg.z_mute_m);
        [f5_r0(c),  f5_lag(c),  f5_l0(c) ] = core.well_tie_full(A, w5,  rb5_r0,  P.dx, lo_r0, hi_r0, 20.0, 30.0);
        [f5_br(c),  ~,          ~        ] = core.well_tie_full(A, w5,  rb5_br,  P.dx, lo_r0, hi_br, 20.0, 30.0);
        [f14_r0(c), f14_lag(c), f14_l0(c)] = core.well_tie_full(A, w14, rb14_r0, P.dx, lo_r0, hi_r0, 20.0, 30.0);
        [bw10(c,:), cen(c)] = core.res_band(A(wz_res, wx_res), P.dx, 256);
        smf(c) = core.steep_smile_frac(A(wz_mid, wx_res), P.dx, 25.0, 3.0, 25.0);
    end
    res.methods.(k) = struct('F5_R0band', f5_r0, 'F5_R0band_mean', mean(f5_r0), ...
        'F5_Broad', f5_br, 'F5_Broad_mean', mean(f5_br), ...
        'F5_ZeroLag', f5_l0, 'F5_BestLag_m', f5_lag, ...
        'F14_R0band', f14_r0, 'F14_R0band_mean', mean(f14_r0), ...
        'F14_ZeroLag', f14_l0, 'F14_BestLag_m', f14_lag, ...
        'Res_Cen', cen, 'Res_BW10', bw10, 'Smile_Frac', smf);
    fprintf('  %-24s | R0=%.4f | Br=%.4f | L0=%+.3f(%+3.0fm) | F14=%+.3f | Cen=%.2f | Sm=%.3f | sl=%+.4f\n', ...
        k, mean(f5_r0), mean(f5_br), f5_l0(min(3,nc)), f5_lag(min(3,nc)), mean(f14_r0), cen(min(3,nc)), smf(min(3,nc)), res.res_slope.(k));
end
res.elapsed_s = toc(t_all);

tag = ''; if smoke, tag = '_smoke'; end
fn_json = fullfile(res_dir, ['volve25_psae_ctrl' tag '.json']);
fid = fopen(fn_json, 'w'); fwrite(fid, jsonencode(res, 'PrettyPrint', true), 'char'); fclose(fid);
fprintf('\nSaved %s (total %.1fs)\n', fn_json, res.elapsed_s);

if ~smoke
    Img_s = struct();
    for ik = 1:numel(mkeys), Img_s.(mkeys{ik}) = single(Img.(mkeys{ik})); end
    save(fullfile(scratch_dir, 'volve25_psae_ctrl_imgs.mat'), 'Img_s', 'res', 'z_in', 'x_in', '-v7.3');
    try rmdir(cs_dir, 's'); fprintf('Cleaned temp Cs dir %s\n', cs_dir); catch, end
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
