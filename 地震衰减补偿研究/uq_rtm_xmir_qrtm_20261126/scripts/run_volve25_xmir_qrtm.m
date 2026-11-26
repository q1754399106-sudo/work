function run_volve25_xmir_qrtm(smoke)
% WP4: Volve 4-25 Hz Real OBC Full-Band (85 freqs, df=0.25 Hz, ZERO Aliasing) Showdown:
%   1) Legacy Straight-Ray xmir_str vs Upgraded Bent-Ray xmir_bent (all acoustic UQ/Det methods)
%   2) Head-to-Head vs True Viscoacoustic Wave-Equation Q-RTM (A_comp, 85 freqs x 3 m in {0.5, 1.0, 2.0}):
%      - WE_QRTM_Unreg_Raw     : Textbook unregularized raw true Q-RTM (no DSA, no PSAE)
%      - WE_QRTM_Unreg_PAnc    : Unregularized true Q-RTM + blind shallow [1000-1600m] PSAE
%      - DSA_WE_QRTM_Unreg_PAnc: True Q-RTM + Bent-DSA + PSAE
%      - WE_QRTM_ZhuLP_PAnc    : Zhu et al. (2014) freq-tapered [16-23 Hz] true Q-RTM + PSAE
%      - WE_QRTM_Reg_PAnc      : Wang/Sun-Zhu Tikhonov gain-stabilized true Q-RTM + PSAE
%      - DSA_WE_QRTM_Reg_PAnc  : Gain-stabilized true Q-RTM + Bent-DSA + PSAE
%      - DSA_WE_QRTM_W1b_PAnc  : Noise-adaptive B=300 stabilized true Q-RTM + Bent-DSA + PSAE
%
% Unified m-grid for ALL 20 methods: m in {0.5, 1.0, 2.0} (Q_sed = 60, 120, 240; Q_water = Inf)
% Strictly read-only on D:\ChatGPT\研究\地震衰减补偿研究\uq_rtm_20261117\cache\fld.
% Resumable per frequency in D:\tmp_uqrtm_xmir_qrtm_20261126\volve25 (auto-cleaned on completion).

if nargin < 1, smoke = false; end
maxNumCompThreads(12);
t_all = tic;

here = fileparts(mfilename('fullpath'));
root = fileparts(here);
res_dir = fullfile(root, 'results');
if ~isfolder(res_dir), mkdir(res_dir); end
fn_out = fullfile(res_dir, sprintf('volve25_xmir_qrtm%s.json', ternary(smoke,'_smoke','')));
if isfile(fn_out) && ~smoke
    fprintf('Skipping WP4 Volve 4-25 Hz (already exists: %s)\n', fn_out);
    return;
end

addpath('D:\ChatGPT\研究\地震衰减补偿研究\uq_rtm_20261117\scripts');
core = uq_rtm_xmir_qrtm_core();
cfg  = core.default_cfg('field_obc');

fld_cache = 'D:\ChatGPT\研究\地震衰减补偿研究\uq_rtm_20261117\cache\fld';
scratch_dir = fullfile('D:\tmp_uqrtm_xmir_qrtm_20261126', ternary(smoke,'volve25_smoke','volve25'));
ac_dir = fullfile(scratch_dir, 'ac_cs');
we_dir = fullfile(scratch_dir, 'we_comp');
if ~isfolder(ac_dir), mkdir(ac_dir); end
if ~isfolder(we_dir), mkdir(we_dir); end

P = fld_setup();
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
am = anc_mask(:); rm = wres_mask(:);

% Blind dip from read-only Stage-4 cached R0
S1 = load(fullfile(fld_cache, 'state1.mat'), 'I');
R0_cached = double(sum(S1.I, 2)); clear S1
R0_2d = core.prep_image(R0_cached, P.nzi, P.nxi, z_in, P.dx, cfg.bg_sigma_m, cfg.z_mute_m);
[phi_eff, c_st] = core.blind_dip(R0_2d, P.dx, cfg.dip_sigma_m, cfg.dip_cap_deg);

% Build (or load cached) legacy xmir_str & upgraded xmir_bent geometry
geom_file = fullfile(scratch_dir, 'volve25_geom_str_bent.mat');
if isfile(geom_file)
    fprintf('Loading cached Str/Bent DSA geometry from %s\n', geom_file);
    G = load(geom_file);
    w_str = G.w_str; TS_str = G.TS_str; w_bent = G.w_bent; TS_bent = G.TS_bent; xmir_stats = G.xmir_stats;
    clear G
else
    fprintf('Computing 124-shot Volve 10m eikonal F_T + Str & Bent DSA geometry...\n');
    s_f = 1.0 ./ P.vmig;
    F_T = zeros(P.nz, P.nx, ns, 'single');
    for k = 1:ns
        F_T(:,:,k) = single(eikonal_fs(s_f, P.dx, P.src_iz(k), P.src_ix(k)));
    end
    TF = load(fullfile(fld_cache, 'tfields.mat'), 'F');
    rx = P.x(P.rec_ix);
    xr_lo = arrayfun(@(s) min(rx(P.shot_of == s)), 1:ns);
    xr_hi = arrayfun(@(s) max(rx(P.shot_of == s)), 1:ns);

    [w_str, TS_str, xm_str] = core.build_dsa_str( ...
        Xin(:), Zin(:), P.z(P.src_iz), 90.0, P.src_x, xr_lo, xr_hi, phi_eff(:), ...
        TF.F, Zi(:), Xi(:), P.nz, P.nx, cfg);
    [w_bent, TS_bent, db] = core.build_dsa_bent( ...
        P.vmig, P, P.z(P.src_iz), 90.0, P.src_x, xr_lo, xr_hi, phi_eff, ...
        F_T, TF.F, cfg);
    clear F_T TF

    res_mat = repmat(rm, 1, ns);
    act_e = res_mat & ((w_str > 0.1) | (w_bent > 0.1)) & (~db.turned);
    xmir_stats = struct( ...
        'res_med_dxmir_m', median(abs(double(db.xmir(act_e)) - double(xm_str(act_e)))), ...
        'res_p90_dxmir_m', prctile(abs(double(db.xmir(act_e)) - double(xm_str(act_e))), 90), ...
        'res_gate_flip_pct', mean((w_str(res_mat) > 0.5) ~= (w_bent(res_mat) > 0.5)) * 100, ...
        'res_p90_dts_ms', prctile(abs(double(TS_bent(act_e)) - double(TS_str(act_e))), 90) * 1e3, ...
        'res_med_fold_str', median(sum(double(w_str(rm, :)), 2)), ...
        'res_med_fold_bent', median(sum(double(w_bent(rm, :)), 2)));
    clear xm_str db
    save(geom_file, 'w_str', 'TS_str', 'w_bent', 'TS_bent', 'xmir_stats', '-v7.3');
end

tstar_anc_str  = median(double(TS_str(am, :)), 'all');
tstar_anc_bent = median(double(TS_bent(am, :)), 'all');
TS_str_rel  = max(TS_str  - single(tstar_anc_str),  single(0));
TS_bent_rel = max(TS_bent - single(tstar_anc_bent), single(0));

fprintf('=== WP4 Volve 4-25 Hz Showdown (smoke=%d) ===\n', smoke);
fprintf('  xmir_bent vs xmir_str (Res): med|dx|=%.1fm (p90=%.1fm) | flip=%.1f%% | fold Str=%.1f Bent=%.1f | t*_anc Str=%.2fms Bent=%.2fms\n', ...
    xmir_stats.res_med_dxmir_m, xmir_stats.res_p90_dxmir_m, xmir_stats.res_gate_flip_pct, ...
    xmir_stats.res_med_fold_str, xmir_stats.res_med_fold_bent, tstar_anc_str*1e3, tstar_anc_bent*1e3);

[~, sgn_b300] = core.make_block_signs(P.shot_of, P.rec_ix, ns, 300, cfg.flip_seed);
sgn_b300 = sgn_b300(:).';

if smoke
    jf_list = [9, 43, 73];
else
    jf_list = 1:nf_all;
end
nf_run = numel(jf_list); f_run = P.f(jf_list);
f_mid = 14.5;
ms = [0.5, 1.0, 2.0]; nm = numel(ms); m1_idx = 2;

%% ========================================================================
%  PASS 1: Acoustic Wavefield (85 freqs, df=0.25 Hz) -> Unaliased INV3_Str & INV3_Bent + PSAE
%  ========================================================================
t_ac_start = tic;
I_raw_f = zeros(nin, nf_run, 'single');
I_str_f = I_raw_f; I_str_fb = I_raw_f;
I_bnt_f = I_raw_f; I_bnt_fb = I_raw_f;
p1_meta = fullfile(ac_dir, 'pass1_ac_meta.mat');
done_p1 = false(1, nf_run);
if isfile(p1_meta)
    S = load(p1_meta);
    if numel(S.done_p1) == nf_run
        I_raw_f = S.I_raw_f; I_str_f = S.I_str_f; I_str_fb = S.I_str_fb;
        I_bnt_f = S.I_bnt_f; I_bnt_fb = S.I_bnt_fb; done_p1 = S.done_p1;
    end
end

for idx_f = 1:nf_run
    jf = jf_list(idx_f); f = P.f(jf);
    cfile = fullfile(ac_dir, sprintf('cs_%03d.bin', jf));
    if done_p1(idx_f) && isfile(cfile)
        continue;
    end
    t0 = tic;
    [L,U,Pp,Qq,R] = lu(core.helm_op_visco(P.vmig, P.vmig, f, P));
    bs = zeros(N, ns); bs(sub2ind([N ns], isrc, (1:ns)')) = -P.W(jf) / P.dx^2;
    Us = core.solve_helm_shot(L,U,Pp,Qq,R, bs);
    Us = conj(single(Us(lin, :))) * single((2*pi*f)^2);
    br = zeros(N, ns); br(col) = P.D(:, jf);
    V  = core.solve_helm_adj(L,U,Pp,Qq,R, br, P.dx);
    Cs = Us .* single(V(lin, :));
    clear L U Pp Qq R bs br V Us

    Cr = real(Cs);
    I_raw_f(:, idx_f) = sum(Cr, 2);
    Cr_s = w_str  .* Cr; I_str_f(:, idx_f) = sum(Cr_s, 2); I_str_fb(:, idx_f) = sum(Cr_s .* sgn_b300, 2);
    Cr_b = w_bent .* Cr; I_bnt_f(:, idx_f) = sum(Cr_b, 2); I_bnt_fb(:, idx_f) = sum(Cr_b .* sgn_b300, 2);
    clear Cr Cr_s Cr_b

    fid = fopen(cfile, 'w'); fwrite(fid, [real(Cs(:)), imag(Cs(:))].', 'single'); fclose(fid);
    clear Cs
    done_p1(idx_f) = true;
    if mod(idx_f, 5) == 0 || idx_f == nf_run
        save(p1_meta, 'I_raw_f', 'I_str_f', 'I_str_fb', 'I_bnt_f', 'I_bnt_fb', 'done_p1', '-v7.3');
    end
    if mod(idx_f, 10) == 0 || idx_f <= 3 || idx_f == nf_run
        fprintf('  [Acoustic Pass1] f=%5.2f Hz (%2d/%d): %.1fs\n', f, idx_f, nf_run, toc(t0));
    end
end
time_ac_pass1 = toc(t_ac_start);

if smoke
    INV3_Str  = ones(nin, nf_run, 'single') * single(0.2);
    INV3_Bent = INV3_Str;
    r0_rel_err = NaN;
else
    INV3_Str  = core.band_inv_snr(I_str_f, I_str_fb, f_run, P.nzi, P.nxi, P.dx, cfg.snr_band_hz, cfg.snr_box_m, cfg.snr_cap);
    INV3_Bent = core.band_inv_snr(I_bnt_f, I_bnt_fb, f_run, P.nzi, P.nxi, P.dx, cfg.snr_band_hz, cfg.snr_box_m, cfg.snr_cap);
    r0_rel_err = norm(double(sum(I_raw_f, 2)) - R0_cached) / max(norm(R0_cached), 1e-30);
    fprintf('  Exact RAW_R0 repro vs state1.mat: rel L2 diff = %.2e\n', r0_rel_err);
end

S_anc_str  = arrayfun(@(q) rms(double(I_str_f(am, q))), 1:nf_run);
S_anc_bent = arrayfun(@(q) rms(double(I_bnt_f(am, q))), 1:nf_run);
[W_eq_str,  beta_str,  ~] = core.fit_psae(f_run, S_anc_str,  [6.0, 22.0], f_mid);
[W_eq_bent, beta_bent, ~] = core.fit_psae(f_run, S_anc_bent, [6.0, 22.0], f_mid);
fprintf('  Acoustic PSAE [6-22 Hz]: Str beta=%+.5f /Hz | Bent beta=%+.5f /Hz\n', beta_str, beta_bent);

%% ========================================================================
%  PASS 2: Evaluate All 13 Acoustic Methods (Str vs Bent) + Cache C_ph_bent for WE-QRTM Blending
%  ========================================================================
mkeys = { ...
    'RAW_R0', 'DSA_Str_R0', 'DSA_Str_R0_Eq', 'DSA_Bent_R0', 'DSA_Bent_R0_Eq', ...
    'DSA_Str_PH_PAnc_25', 'DSA_Str_IQ70_PAnc_25', 'DSA_Str_W1b_PAnc_25', 'DSA_Str_UQW1b_PAnc_25', ...
    'DSA_Bent_PH_PAnc_25', 'DSA_Bent_IQ70_PAnc_25', 'DSA_Bent_W1b_PAnc_25', 'DSA_Bent_UQW1b_PAnc_25', ...
    'WE_QRTM_Unreg_Raw', 'WE_QRTM_Unreg_PAnc', 'DSA_WE_QRTM_Unreg_PAnc', ...
    'WE_QRTM_ZhuLP_PAnc', 'WE_QRTM_Reg_PAnc', 'DSA_WE_QRTM_Reg_PAnc', 'DSA_WE_QRTM_W1b_PAnc'};
r0_only = {'RAW_R0', 'DSA_Str_R0', 'DSA_Str_R0_Eq', 'DSA_Bent_R0', 'DSA_Bent_R0_Eq'};

Img = struct(); spec_res = struct();
for ik = 1:numel(mkeys)
    nc = 1 + (nm - 1) * ~ismember(mkeys{ik}, r0_only);
    Img.(mkeys{ik}) = zeros(nin, nc);
    spec_res.(mkeys{ik}) = zeros(1, nf_run);
end

s70 = single(gain_s(70));
xmax_rel = 1.05 * double(max(max(TS_str_rel(:)), max(TS_bent_rel(:)))) / min(ms) + 1e-6;
xg_rel = single(linspace(0, xmax_rel, 1201)); dxg_rel = single(xg_rel(2)); ntab_m2 = single(1201 - 2);

% Keep per-freq acoustic phase-only images (only nin x nm x nf_run? Wait: for WE blending we read cs_%03d.bin directly inside Pass 3!)
t_ac_p2 = tic;
for idx_f = 1:nf_run
    jf = jf_list(idx_f); f = P.f(jf);
    ws = single(W_eq_str(idx_f)); wb = single(W_eq_bent(idx_f));
    fid = fopen(fullfile(ac_dir, sprintf('cs_%03d.bin', jf)));
    a = fread(fid, [2, nin*ns], '*single'); fclose(fid);
    Cs = reshape(complex(a(1,:), a(2,:)), nin, ns); clear a

    iv_s = INV3_Str(:, idx_f); iv_b = INV3_Bent(:, idx_f);
    T0 = uq_tab(f, xg_rel, 0, 25.0); T1 = uq_tab(f, xg_rel, cfg.sig, 25.0);
    P0r = single(real(T0.P)); P0i = single(imag(T0.P));
    D0r = single(real(T0.A - T0.P)); D0i = single(imag(T0.A - T0.P)); E0_t = single(1 + double(s70)*real(T0.C));
    IQ_t = (T0.A + double(s70)*T0.P.*T0.C) ./ (1 + double(s70)*T0.C); IQr = single(real(IQ_t)); IQi = single(imag(IQ_t));
    P1r = single(real(T1.P)); P1i = single(imag(T1.P));
    D1r = single(real(T1.A - T1.P)); D1i = single(imag(T1.A - T1.P)); E1_t = single(1 + double(s70)*real(T1.C));

    acc_r0   = I_raw_f(:, idx_f);
    acc_sr0  = I_str_f(:, idx_f); acc_sre = ws * acc_sr0;
    acc_br0  = I_bnt_f(:, idx_f); acc_bre = wb * acc_br0;
    Img.RAW_R0         = Img.RAW_R0         + double(acc_r0);  spec_res.RAW_R0(idx_f)         = rms(double(acc_r0(rm)));
    Img.DSA_Str_R0     = Img.DSA_Str_R0     + double(acc_sr0); spec_res.DSA_Str_R0(idx_f)     = rms(double(acc_sr0(rm)));
    Img.DSA_Str_R0_Eq  = Img.DSA_Str_R0_Eq  + double(acc_sre); spec_res.DSA_Str_R0_Eq(idx_f)  = rms(double(acc_sre(rm)));
    Img.DSA_Bent_R0    = Img.DSA_Bent_R0    + double(acc_br0); spec_res.DSA_Bent_R0(idx_f)    = rms(double(acc_br0(rm)));
    Img.DSA_Bent_R0_Eq = Img.DSA_Bent_R0_Eq + double(acc_bre); spec_res.DSA_Bent_R0_Eq(idx_f) = rms(double(acc_bre(rm)));

    a_sph = zeros(nin,nm,'single'); a_siq = a_sph; a_sw1 = a_sph; a_suq = a_sph;
    a_bph = zeros(nin,nm,'single'); a_biq = a_bph; a_bw1 = a_bph; a_buq = a_bph;
    for c = 1:ceil(ns/31)
        sh = (c-1)*31+1 : min(c*31, ns);
        Cr = real(Cs(:, sh)); Ci = imag(Cs(:, sh));
        Cr_s = w_str(:, sh)  .* Cr; Ci_s = w_str(:, sh)  .* Ci; ts_sr = TS_str_rel(:, sh);
        Cr_b = w_bent(:, sh) .* Cr; Ci_b = w_bent(:, sh) .* Ci; ts_br = TS_bent_rel(:, sh);
        for k = 1:nm
            invm = single(1.0 / ms(k));
            % Legacy Str
            [idx, ww] = core.tab_idx(ts_sr * invm, dxg_rel, ntab_m2); ow = 1 - ww;
            p0 = (P0r(idx).*ow + P0r(idx+1).*ww).*Cr_s - (P0i(idx).*ow + P0i(idx+1).*ww).*Ci_s;
            d0 = (D0r(idx).*ow + D0r(idx+1).*ww).*Cr_s - (D0i(idx).*ow + D0i(idx+1).*ww).*Ci_s;
            E0 = E0_t(idx).*ow + E0_t(idx+1).*ww;
            iq = (IQr(idx).*ow + IQr(idx+1).*ww).*Cr_s - (IQi(idx).*ow + IQi(idx+1).*ww).*Ci_s;
            p1 = (P1r(idx).*ow + P1r(idx+1).*ww).*Cr_s - (P1i(idx).*ow + P1i(idx+1).*ww).*Ci_s;
            d1 = (D1r(idx).*ow + D1r(idx+1).*ww).*Cr_s - (D1i(idx).*ow + D1i(idx+1).*ww).*Ci_s;
            E1 = E1_t(idx).*ow + E1_t(idx+1).*ww;
            a_sph(:,k) = a_sph(:,k) + ws * sum(p0, 2);
            a_siq(:,k) = a_siq(:,k) + ws * sum(iq, 2);
            a_sw1(:,k) = a_sw1(:,k) + ws * sum(p0 + d0 ./ (E0 + iv_s), 2);
            a_suq(:,k) = a_suq(:,k) + ws * sum(p1 + d1 ./ (E1 + iv_s), 2);

            % Upgraded Bent
            [idx, ww] = core.tab_idx(ts_br * invm, dxg_rel, ntab_m2); ow = 1 - ww;
            p0 = (P0r(idx).*ow + P0r(idx+1).*ww).*Cr_b - (P0i(idx).*ow + P0i(idx+1).*ww).*Ci_b;
            d0 = (D0r(idx).*ow + D0r(idx+1).*ww).*Cr_b - (D0i(idx).*ow + D0i(idx+1).*ww).*Ci_b;
            E0 = E0_t(idx).*ow + E0_t(idx+1).*ww;
            iq = (IQr(idx).*ow + IQr(idx+1).*ww).*Cr_b - (IQi(idx).*ow + IQi(idx+1).*ww).*Ci_b;
            p1 = (P1r(idx).*ow + P1r(idx+1).*ww).*Cr_b - (P1i(idx).*ow + P1i(idx+1).*ww).*Ci_b;
            d1 = (D1r(idx).*ow + D1r(idx+1).*ww).*Cr_b - (D1i(idx).*ow + D1i(idx+1).*ww).*Ci_b;
            E1 = E1_t(idx).*ow + E1_t(idx+1).*ww;
            a_bph(:,k) = a_bph(:,k) + wb * sum(p0, 2);
            a_biq(:,k) = a_biq(:,k) + wb * sum(iq, 2);
            a_bw1(:,k) = a_bw1(:,k) + wb * sum(p0 + d0 ./ (E0 + iv_b), 2);
            a_buq(:,k) = a_buq(:,k) + wb * sum(p1 + d1 ./ (E1 + iv_b), 2);
        end
    end
    Img.DSA_Str_PH_PAnc_25     = Img.DSA_Str_PH_PAnc_25     + double(a_sph); spec_res.DSA_Str_PH_PAnc_25(idx_f)     = rms(double(a_sph(rm, m1_idx)));
    Img.DSA_Str_IQ70_PAnc_25   = Img.DSA_Str_IQ70_PAnc_25   + double(a_siq); spec_res.DSA_Str_IQ70_PAnc_25(idx_f)   = rms(double(a_siq(rm, m1_idx)));
    Img.DSA_Str_W1b_PAnc_25    = Img.DSA_Str_W1b_PAnc_25    + double(a_sw1); spec_res.DSA_Str_W1b_PAnc_25(idx_f)    = rms(double(a_sw1(rm, m1_idx)));
    Img.DSA_Str_UQW1b_PAnc_25  = Img.DSA_Str_UQW1b_PAnc_25  + double(a_suq); spec_res.DSA_Str_UQW1b_PAnc_25(idx_f)  = rms(double(a_suq(rm, m1_idx)));
    Img.DSA_Bent_PH_PAnc_25    = Img.DSA_Bent_PH_PAnc_25    + double(a_bph); spec_res.DSA_Bent_PH_PAnc_25(idx_f)    = rms(double(a_bph(rm, m1_idx)));
    Img.DSA_Bent_IQ70_PAnc_25  = Img.DSA_Bent_IQ70_PAnc_25  + double(a_biq); spec_res.DSA_Bent_IQ70_PAnc_25(idx_f)  = rms(double(a_biq(rm, m1_idx)));
    Img.DSA_Bent_W1b_PAnc_25   = Img.DSA_Bent_W1b_PAnc_25   + double(a_bw1); spec_res.DSA_Bent_W1b_PAnc_25(idx_f)   = rms(double(a_bw1(rm, m1_idx)));
    Img.DSA_Bent_UQW1b_PAnc_25 = Img.DSA_Bent_UQW1b_PAnc_25 + double(a_buq); spec_res.DSA_Bent_UQW1b_PAnc_25(idx_f) = rms(double(a_buq(rm, m1_idx)));
end
time_ac_p2 = toc(t_ac_p2);
fprintf('  Acoustic 13-method Pass2 finished in %.1fs\n', time_ac_p2);

%% ========================================================================
%  PASS 3: True Viscoacoustic Wave-Equation Q-RTM on Volve 4-25 Hz (m in {0.5, 1.0, 2.0})
%  ========================================================================
% Zhu et al. (2014) raised-cosine taper on Volve 4-25 Hz: 1 for f<=16 Hz, cosine [16..23], 0 for f>=23 Hz
w_lp_v = core.cos_taper(23.0 - f_run, 0.0, 7.0);
t_we_start = tic;
we_lu_total = 0;

for im = 1:nm
    mv = ms(im); Q_est = mv * P.Qm;
    fprintf('\n  --- True Viscoacoustic Wave-Equation Q-RTM Pass (m=%.2f, Q_sed=%.0f) ---\n', mv, 120*mv);
    m_dir = fullfile(we_dir, sprintf('m_%03d', round(mv*100)));
    if ~isfolder(m_dir), mkdir(m_dir); end

    % Pass 3a: Compute & checkpoint per-freq 2D images (raw, rel, dsa, reg_raw, reg_dsa, w1b_dsa, zhu_raw)
    I_we_unr_raw = zeros(nin, nf_run, 'single');
    I_we_unr_rel = zeros(nin, nf_run, 'single');
    I_we_unr_dsa = zeros(nin, nf_run, 'single');
    I_we_zhu_raw = zeros(nin, nf_run, 'single');
    I_we_reg_raw = zeros(nin, nf_run, 'single');
    I_we_reg_dsa = zeros(nin, nf_run, 'single');
    I_we_w1b_dsa = zeros(nin, nf_run, 'single');

    for idx_f = 1:nf_run
        jf = jf_list(idx_f); f = P.f(jf);
        ckpt = fullfile(m_dir, sprintf('we_%03d.mat', jf));
        if isfile(ckpt)
            D = load(ckpt);
            I_we_unr_raw(:,idx_f) = D.unr_raw; I_we_unr_rel(:,idx_f) = D.unr_rel;
            I_we_unr_dsa(:,idx_f) = D.unr_dsa; I_we_zhu_raw(:,idx_f) = D.zhu_raw;
            I_we_reg_raw(:,idx_f) = D.reg_raw; I_we_reg_dsa(:,idx_f) = D.reg_dsa;
            I_we_w1b_dsa(:,idx_f) = D.w1b_dsa;
            continue;
        end
        t0 = tic;
        % Genuine True Viscoacoustic Wave-Equation Compensated Propagator A_comp (interior-only Q, real C-PML!)
        vc_c = core.make_visco_vel(P.vmig, Q_est, f, 25.0, P.npml, 'comp', 1.0);
        [L,U,Pp,Qq,R] = lu(core.helm_op_visco(vc_c, P.vmig, f, P));
        we_lu_total = we_lu_total + 1;
        bs = zeros(N, ns); bs(sub2ind([N ns], isrc, (1:ns)')) = -P.W(jf) / P.dx^2;
        Us = core.solve_helm_shot(L,U,Pp,Qq,R, bs);
        Us = conj(single(Us(lin, :))) * single((2*pi*f)^2);
        br = zeros(N, ns); br(col) = P.D(:, jf);
        V  = core.solve_helm_adj(L,U,Pp,Qq,R, br, P.dx);
        Cc = Us .* single(V(lin, :));
        clear L U Pp Qq R bs br V Us

        % Load acoustic Cs for exact dispersion-phase fallback C_ph (same-scale shallow-relative normalization!)
        fid = fopen(fullfile(ac_dir, sprintf('cs_%03d.bin', jf)));
        a = fread(fid, [2, nin*ns], '*single'); fclose(fid);
        Cs_ac = reshape(complex(a(1,:), a(2,:)), nin, ns); clear a

        iv_b = INV3_Bent(:, idx_f);
        % Remove scalar shallow-anchor exponential growth exp(+pi*f*t*_anc/m) so Cc_rel matches C_ph scale in [1000-1600m]!
        g_anc = single(exp(pi * f * tstar_anc_bent / mv));
        inv_g = single(1.0) / g_anc;
        k_disp = single(2.0 * f * log(max(f, 1e-6) / 25.0));
        wl = single(w_lp_v(idx_f));

        unr_raw = zeros(nin, 1, 'single'); unr_rel = unr_raw; unr_dsa = unr_raw;
        zhu_raw = unr_raw; reg_raw = unr_raw; reg_dsa = unr_raw; w1b_dsa = unr_raw;
        for c = 1:ceil(ns/31)
            sh = (c-1)*31+1 : min(c*31, ns);
            cc_sh = real(Cc(:, sh));
            cc_rel = cc_sh * inv_g;
            wb_sh = w_bent(:, sh);
            ts_r  = TS_bent_rel(:, sh) / single(mv);
            ts_a  = TS_bent(:, sh)     / single(mv);
            % Exact phase-dispersion-only cross-correlation matching A_comp's absolute dispersion phase!
            ang = -k_disp * ts_a;
            cph_sh = cos(ang) .* real(Cs_ac(:, sh)) - sin(ang) .* imag(Cs_ac(:, sh));

            % Wave-equation Tikhonov stabilization weights (exact Wang s70 & noise-adaptive B=300 W1b)
            b_rel = exp(-single(pi * f) * ts_r);
            b2 = b_rel .* b_rel;
            ws_reg = b2 ./ (b2 + s70);
            ws_w1b = b2 ./ (b2 + s70 + b2 .* iv_b);

            cr_reg = ws_reg .* cc_rel + (single(1.0) - ws_reg) .* cph_sh;
            cr_w1b = ws_w1b .* cc_rel + (single(1.0) - ws_w1b) .* cph_sh;
            cr_zhu = wl * cc_rel + (single(1.0) - wl) * cph_sh;

            unr_raw = unr_raw + sum(cc_sh, 2);
            unr_rel = unr_rel + sum(cc_rel, 2);
            unr_dsa = unr_dsa + sum(wb_sh .* cc_rel, 2);
            zhu_raw = zhu_raw + sum(cr_zhu, 2);
            reg_raw = reg_raw + sum(cr_reg, 2);
            reg_dsa = reg_dsa + sum(wb_sh .* cr_reg, 2);
            w1b_dsa = w1b_dsa + sum(wb_sh .* cr_w1b, 2);
        end
        clear Cc Cs_ac

        save(ckpt, 'unr_raw', 'unr_rel', 'unr_dsa', 'zhu_raw', 'reg_raw', 'reg_dsa', 'w1b_dsa', '-v6');
        I_we_unr_raw(:,idx_f) = unr_raw; I_we_unr_rel(:,idx_f) = unr_rel;
        I_we_unr_dsa(:,idx_f) = unr_dsa; I_we_zhu_raw(:,idx_f) = zhu_raw;
        I_we_reg_raw(:,idx_f) = reg_raw; I_we_reg_dsa(:,idx_f) = reg_dsa;
        I_we_w1b_dsa(:,idx_f) = w1b_dsa;
        if mod(idx_f, 10) == 0 || idx_f <= 3 || idx_f == nf_run
            fprintf('  [WE-QRTM m=%.2f] f=%5.2f Hz (%2d/%d): %.1fs\n', mv, f, idx_f, nf_run, toc(t0));
        end
    end

    % Fit blind shallow [1000-1600m] PSAE for each WE-QRTM method, matching shallow anchor RMS at f_mid
    % to DSA_Bent_R0_Eq so G_HF measures pure relative reservoir-to-shallow amplification (zero scalar offset)!
    [~, im_mid] = min(abs(f_run - f_mid));
    ref_anc_mid = S_anc_bent(im_mid) * W_eq_bent(im_mid);
    fit_we_psae_norm = @(M) norm_psae(core, f_run, M, am, [6.0, 22.0], f_mid, im_mid, ref_anc_mid);
    Weq_ur = fit_we_psae_norm(I_we_unr_rel);
    Weq_ud = fit_we_psae_norm(I_we_unr_dsa);
    Weq_zh = fit_we_psae_norm(I_we_zhu_raw);
    Weq_rr = fit_we_psae_norm(I_we_reg_raw);
    Weq_rd = fit_we_psae_norm(I_we_reg_dsa);
    Weq_wd = fit_we_psae_norm(I_we_w1b_dsa);

    for idx_f = 1:nf_run
        % Unanchored raw textbook Q-RTM
        v_raw = double(I_we_unr_raw(:, idx_f));
        Img.WE_QRTM_Unreg_Raw(:, im) = Img.WE_QRTM_Unreg_Raw(:, im) + v_raw;
        % Shallow-anchored PSAE True Q-RTM variants
        v_ur = double(I_we_unr_rel(:, idx_f)) * Weq_ur(idx_f);
        v_ud = double(I_we_unr_dsa(:, idx_f)) * Weq_ud(idx_f);
        v_zh = double(I_we_zhu_raw(:, idx_f)) * Weq_zh(idx_f);
        v_rr = double(I_we_reg_raw(:, idx_f)) * Weq_rr(idx_f);
        v_rd = double(I_we_reg_dsa(:, idx_f)) * Weq_rd(idx_f);
        v_wd = double(I_we_w1b_dsa(:, idx_f)) * Weq_wd(idx_f);

        Img.WE_QRTM_Unreg_PAnc(:, im)     = Img.WE_QRTM_Unreg_PAnc(:, im)     + v_ur;
        Img.DSA_WE_QRTM_Unreg_PAnc(:, im) = Img.DSA_WE_QRTM_Unreg_PAnc(:, im) + v_ud;
        Img.WE_QRTM_ZhuLP_PAnc(:, im)     = Img.WE_QRTM_ZhuLP_PAnc(:, im)     + v_zh;
        Img.WE_QRTM_Reg_PAnc(:, im)       = Img.WE_QRTM_Reg_PAnc(:, im)       + v_rr;
        Img.DSA_WE_QRTM_Reg_PAnc(:, im)   = Img.DSA_WE_QRTM_Reg_PAnc(:, im)   + v_rd;
        Img.DSA_WE_QRTM_W1b_PAnc(:, im)   = Img.DSA_WE_QRTM_W1b_PAnc(:, im)   + v_wd;

        if im == m1_idx
            spec_res.WE_QRTM_Unreg_Raw(idx_f)      = rms(v_raw(rm));
            spec_res.WE_QRTM_Unreg_PAnc(idx_f)     = rms(v_ur(rm));
            spec_res.DSA_WE_QRTM_Unreg_PAnc(idx_f) = rms(v_ud(rm));
            spec_res.WE_QRTM_ZhuLP_PAnc(idx_f)     = rms(v_zh(rm));
            spec_res.WE_QRTM_Reg_PAnc(idx_f)       = rms(v_rr(rm));
            spec_res.DSA_WE_QRTM_Reg_PAnc(idx_f)   = rms(v_rd(rm));
            spec_res.DSA_WE_QRTM_W1b_PAnc(idx_f)   = rms(v_wd(rm));
        end
    end
end
time_we_total = toc(t_we_start);

%% ========================================================================
%  PASS 4: Full Zero-Omission Well-Tie & Spectral Evaluation (F-5, F-14, L0, lag, Res_Cen, Smile, Slope, G_HF)
%  ========================================================================
w5  = core.load_well(fullfile(fld_cache, 'well_F-5.bin'),  z_in, x_in);
w14 = core.load_well(fullfile(fld_cache, 'well_F-14.bin'), z_in, x_in);
R0_img = core.prep_image(Img.RAW_R0, P.nzi, P.nxi, z_in, P.dx, cfg.bg_sigma_m, cfg.z_mute_m);
iz_res = find(z_in >= cfg.res_z_m(1) & z_in <= cfg.res_z_m(2));
ix_res = find(x_in >= cfg.res_x_m(1) & x_in <= cfg.res_x_m(2));
iz_ovb = find(z_in >= cfg.ovb_z_m(1) & z_in <= cfg.ovb_z_m(2));
[bw_r0, ~] = core.res_band(R0_img(iz_res, ix_res), P.dx, 256);
rb5_r0  = core.bp_z(w5.rw,  w5.dzw,  bw_r0(1), bw_r0(2));
rb5_br  = core.bp_z(w5.rw,  w5.dzw,  bw_r0(1), 16.0);
rb14_r0 = core.bp_z(w14.rw, w14.dzw, bw_r0(1), bw_r0(2));

out = struct('dataset', 'volve_4_25hz_xmir_qrtm', 'smoke', smoke, 'f', f_run, 'm', ms, ...
             'r0_repro_rel_l2', r0_rel_err, 'bw_r0_cyckm', bw_r0, ...
             'xmir_stats', xmir_stats, ...
             'psae', struct('beta_str', beta_str, 'beta_bent', beta_bent, ...
                            'tstar_anc_str_ms', tstar_anc_str*1e3, 'tstar_anc_bent_ms', tstar_anc_bent*1e3), ...
             'compute_cost', struct('acoustic_LU_total', nf_run, 'we_qrtm_LU_total', we_lu_total, ...
                                    'time_acoustic_p1_s', time_ac_pass1, 'time_acoustic_p2_s', time_ac_p2, ...
                                    'time_we_qrtm_total_s', time_we_total), ...
             'methods', struct(), 'spec_res', spec_res, 'res_slope', struct(), 'G_HF_18_24', struct(), 'elapsed_s', 0);

fm_sl = (f_run >= 6.0 & f_run <= 23.0);
fb_hf = (f_run >= 18.0 & f_run <= 24.0);
sb_eq = spec_res.DSA_Bent_R0_Eq;

fprintf('\n=== Volve 4-25 Hz Full 20-Method Results (bw_r0=[%.2f, %.2f] cyc/km) ===\n', bw_r0(1), bw_r0(2));
saved_2d = struct('z_in', z_in, 'x_in', x_in);
for ik = 1:numel(mkeys)
    key = mkeys{ik}; M = Img.(key); nc = size(M, 2); c_m1 = min(m1_idx, nc);
    f5_r0 = zeros(1, nc); f5_br = f5_r0; f5_l0 = f5_r0; f5_lag = f5_r0;
    f14_r0 = f5_r0; f14_l0 = f5_r0; f14_lag = f5_r0; cen = f5_r0; smf = f5_r0;
    for c = 1:nc
        A = core.prep_image(M(:, c), P.nzi, P.nxi, z_in, P.dx, cfg.bg_sigma_m, cfg.z_mute_m);
        [f5_r0(c), f5_lag(c), f5_l0(c)]   = core.well_tie_full(A, w5,  rb5_r0,  P.dx, bw_r0(1), bw_r0(2), 20.0, 30.0);
        [f5_br(c), ~, ~]                  = core.well_tie_full(A, w5,  rb5_br,  P.dx, bw_r0(1), 16.0,     20.0, 30.0);
        [f14_r0(c), f14_lag(c), f14_l0(c)] = core.well_tie_full(A, w14, rb14_r0, P.dx, bw_r0(1), bw_r0(2), 20.0, 30.0);
        [~, cen(c)] = core.res_band(A(iz_res, ix_res), P.dx, 256);
        smf(c) = core.steep_smile_frac(A(iz_ovb, ix_res), P.dx, 25.0, 3.0, 25.0);
        if c == c_m1 && ismember(key, {'RAW_R0','DSA_Str_R0_Eq','DSA_Bent_R0_Eq','DSA_Str_UQW1b_PAnc_25','DSA_Bent_UQW1b_PAnc_25','WE_QRTM_Unreg_PAnc','DSA_WE_QRTM_Reg_PAnc','DSA_WE_QRTM_W1b_PAnc'})
            saved_2d.([key '_m1']) = single(A);
        end
        if c == 1 && nc > 1 && ismember(key, {'DSA_Bent_UQW1b_PAnc_25','WE_QRTM_Unreg_PAnc','DSA_WE_QRTM_W1b_PAnc'})
            saved_2d.([key '_m05']) = single(A);
        end
    end
    p_sl = polyfit(f_run(fm_sl), log(max(spec_res.(key)(fm_sl), 1e-30)), 1);
    out.res_slope.(key) = p_sl(1);
    ghf = median(spec_res.(key)(fb_hf) ./ max(sb_eq(fb_hf), 1e-30));
    out.G_HF_18_24.(key) = ghf;

    out.methods.(key) = struct('F5_R0band', f5_r0, 'F5_R0band_mean', mean(f5_r0), ...
        'F5_Broad', f5_br, 'F5_Broad_mean', mean(f5_br), ...
        'F5_ZeroLag', f5_l0, 'F5_BestLag_m', f5_lag, ...
        'F14_R0band', f14_r0, 'F14_R0band_mean', mean(f14_r0), ...
        'F14_ZeroLag', f14_l0, 'F14_BestLag_m', f14_lag, ...
        'Res_Cen', cen, 'Smile_Frac', smf);
    fprintf('  %-24s | R0=%.4f (m0.5=%.3f m1=%.3f) | Br=%.4f | L0=%+.3f(%+.0fm) | F14=%+.3f | Cen=%.2f | Sm=%.1f%% | sl=%+.4f | G_HF=%.2fx\n', ...
        key, mean(f5_r0), f5_r0(1), f5_r0(c_m1), mean(f5_br), f5_l0(c_m1), f5_lag(c_m1), ...
        mean(f14_r0), cen(c_m1), smf(c_m1)*100, p_sl(1), ghf);
end

out.elapsed_s = toc(t_all);
fid = fopen(fn_out, 'w'); fwrite(fid, jsonencode(out, 'PrettyPrint', true), 'char'); fclose(fid);
fprintf('\nSaved %s (total %.1fs)\n', fn_out, out.elapsed_s);

if ~smoke
    save('D:\tmp_uqrtm_xmir_qrtm_20261126\volve25_xmir_qrtm_imgs.mat', '-struct', 'saved_2d', '-v7.3');
    try rmdir(ac_dir, 's'); catch, end
    try rmdir(we_dir, 's'); catch, end
    fprintf('Saved final Volve 2D images and cleaned temp slice dirs.\n');
end
end

function W_norm = norm_psae(core, f_run, M, am, f_band, f_mid, im_mid, ref_anc_mid)
S_anc = arrayfun(@(q) rms(double(M(am, q))), 1:numel(f_run));
[W_eq, ~, ~] = core.fit_psae(f_run, S_anc, f_band, f_mid);
scale = ref_anc_mid / max(S_anc(im_mid) * W_eq(im_mid), 1e-30);
W_norm = W_eq * scale;
end

function s = ternary(cond, a, b)
if cond, s = a; else, s = b; end
end
