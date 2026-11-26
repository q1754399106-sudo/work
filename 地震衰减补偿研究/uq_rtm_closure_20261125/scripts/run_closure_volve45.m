function run_closure_volve45(smoke)
% WP2 (P0): Volve 4-45 Hz TRUE Unaliased Same-Domain INV3_DSA Confluence
% Fixes E-1 (wiring bug in prior run_unified_volve45.m) using a lean two-pass
% architecture (< 2.5 GB peak disk in D:\tmp_uqrtm_closure_20261125\, resumable per freq):
%   Pass 1 (165 freqs, df=0.25 Hz, T_max=4.0 s): computes unaliased INV3_DSA,
%          INV3_RAW_FULL, INV3_DSA_NOW2, plus blind shallow [1000-1600m] PSAE & Decon.
%          Persists inv_fullband_45.mat (NEVER discarded!).
%   Pass 2 (165 freqs): injects TRUE unaliased INV3_DSA into all 18 controlled methods,
%          including NoW2 (zero PSAE, no (2*pi*f)^2), SameWeq (exact D2 test),
%          and S4Hiv control (reproducing prior dormant 0.5597).
%
% Strictly read-only on Stage 4H cache & prior fld45_dsa_geom.mat!

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

fld45_cache = 'D:\ChatGPT\研究\.claude\worktrees\study-research-docs-plan-c4d4c6\地震衰减补偿研究\uq_rtm_20261117\cache\fld45';
well_dir    = 'D:\ChatGPT\研究\.claude\worktrees\study-research-docs-plan-c4d4c6\地震衰减补偿研究\uq_rtm_20261117\cache\fld';
old_geom    = 'D:\tmp_uqrtm_unified_20261124\fld45_dsa_geom.mat';

scratch_dir = 'D:\tmp_uqrtm_closure_20261125';
if smoke
    p1_dir = fullfile(scratch_dir, 'p1_fld45_smoke');
    p2_dir = fullfile(scratch_dir, 'p2_fld45_smoke');
else
    p1_dir = fullfile(scratch_dir, 'p1_fld45');
    p2_dir = fullfile(scratch_dir, 'p2_fld45');
end
if ~isfolder(p1_dir), mkdir(p1_dir); end
if ~isfolder(p2_dir), mkdir(p2_dir); end

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

% Load Stage 4H reference R0 & UQW1b for exact repro check
S_im = load(fullfile(fld45_cache, 'images.mat'), 'R0', 'UQW1b');
R0_stage4h = double(S_im.R0);
UQW1b_stage4h_m1 = double(S_im.UQW1b(:, 3));
clear S_im

% Reuse verified read-only DSA geometry
G = load(old_geom);
w_dsa = G.w_dsa; TS_dsa = G.TS_dsa; phi_eff = G.phi_eff; c_st = G.c_st;
tstar_anc_dsa = G.tstar_anc_dsa; tstar_anc_flat = G.tstar_anc_flat;
clear G

TS_flat_S = load(fullfile(fld45_cache, 'ts.mat'), 'TS');
TS_flat = TS_flat_S.TS; clear TS_flat_S

IV_S = load(fullfile(fld45_cache, 'inv.mat'), 'iv', 'iv3');
INV1_S4H = IV_S.iv; INV3_S4H = IV_S.iv3; clear IV_S

[~, sgn_b300] = core.make_block_signs(P.shot_of, P.rec_ix, ns, 300, cfg.flip_seed);
sgn_b300 = sgn_b300(:).';

if smoke
    jf_list = [17, 85, 145];
else
    jf_list = 1:nf_all;
end
nf_run = numel(jf_list); f_run = P.f(jf_list);
f_mid = 25.0;

fprintf('=== WP2 Volve 4-45 Hz TRUE Unaliased INV3_DSA Confluence (smoke=%d, nf=%d) ===\n', smoke, nf_run);

%% ========================================================================
%  PASS 1: Full 165-freq df=0.25 Hz Unaliased B=300 & Shallow Anchoring
%  ========================================================================
inv_file = fullfile(scratch_dir, sprintf('inv_fullband_45%s.mat', ternary(smoke,'_smoke','')));
if ~isfile(inv_file)
    fprintf('\n--- Pass 1: Computing full-band unaliased 2D frequency slices ---\n');
    for idx_f = 1:nf_run
        jf = jf_list(idx_f); f = P.f(jf);
        p1f = fullfile(p1_dir, sprintf('p1_%03d.mat', jf));
        if isfile(p1f)
            if mod(idx_f, 20) == 0 || idx_f == nf_run
                fprintf('  [Pass1 Cache] f=%5.2f Hz (%3d/%d)\n', f, idx_f, nf_run);
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

        I_raw_f = zeros(nin, 1, 'single'); I_raw_flip300_f = I_raw_f;
        I_dsa_f = I_raw_f;                 I_dsa_flip300_f = I_raw_f;
        for c = 1:ceil(ns/31)
            sh = (c-1)*31+1 : min(c*31, ns);
            Cr = real(Us(:, sh) .* V_in(:, sh));
            Cr_d = w_dsa(:, sh) .* Cr;
            I_raw_f         = I_raw_f         + sum(Cr, 2);
            I_raw_flip300_f = I_raw_flip300_f + sum(Cr   .* sgn_b300(sh), 2);
            I_dsa_f         = I_dsa_f         + sum(Cr_d, 2);
            I_dsa_flip300_f = I_dsa_flip300_f + sum(Cr_d .* sgn_b300(sh), 2);
        end
        clear Us V_in Cr Cr_d
        nw2_sc = single((f_mid / f)^2);
        I_dsa_now2_f         = nw2_sc * I_dsa_f;
        I_dsa_now2_flip300_f = nw2_sc * I_dsa_flip300_f;

        save(p1f, 'f', 'I_raw_f', 'I_raw_flip300_f', 'I_dsa_f', 'I_dsa_flip300_f', ...
            'I_dsa_now2_f', 'I_dsa_now2_flip300_f', '-v6');
        if mod(idx_f, 10) == 0 || idx_f == nf_run || idx_f <= 3
            fprintf('  Pass1 f=%5.2f Hz (%3d/%d): %.1fs\n', f, idx_f, nf_run, toc(t0));
        end
    end

    % Assemble full-band 2D matrices & compute unaliased INV3_DSA, INV3_RAW_FULL, INV3_DSA_NOW2
    I_raw_all  = zeros(nin, nf_run, 'single'); I_rawf_all  = I_raw_all;
    I_dsa_all  = I_raw_all;                    I_dsaf_all  = I_raw_all;
    I_nw2_all  = I_raw_all;                    I_nw2f_all  = I_raw_all;
    S_anc_raw  = zeros(1, nf_run);             S_anc_dsa   = zeros(1, nf_run);
    am = anc_mask(:);
    for idx_f = 1:nf_run
        D = load(fullfile(p1_dir, sprintf('p1_%03d.mat', jf_list(idx_f))));
        I_raw_all(:, idx_f)  = D.I_raw_f;        I_rawf_all(:, idx_f) = D.I_raw_flip300_f;
        I_dsa_all(:, idx_f)  = D.I_dsa_f;        I_dsaf_all(:, idx_f) = D.I_dsa_flip300_f;
        I_nw2_all(:, idx_f)  = D.I_dsa_now2_f;   I_nw2f_all(:, idx_f) = D.I_dsa_now2_flip300_f;
        S_anc_raw(idx_f) = rms(double(D.I_raw_f(am)));
        S_anc_dsa(idx_f) = rms(double(D.I_dsa_f(am)));
    end

    if smoke
        INV3_DSA      = INV3_S4H(:, jf_list);
        INV3_RAW_FULL = INV3_S4H(:, jf_list);
        INV3_DSA_NOW2 = INV3_S4H(:, jf_list);
    else
        INV3_DSA      = core.band_inv_snr(I_dsa_all, I_dsaf_all, f_run, P.nzi, P.nxi, P.dx, cfg.snr_band_hz, cfg.snr_box_m, cfg.snr_cap);
        INV3_RAW_FULL = core.band_inv_snr(I_raw_all, I_rawf_all, f_run, P.nzi, P.nxi, P.dx, cfg.snr_band_hz, cfg.snr_box_m, cfg.snr_cap);
        INV3_DSA_NOW2 = core.band_inv_snr(I_nw2_all, I_nw2f_all, f_run, P.nzi, P.nxi, P.dx, cfg.snr_band_hz, cfg.snr_box_m, cfg.snr_cap);
    end
    clear I_raw_all I_rawf_all I_dsa_all I_dsaf_all I_nw2_all I_nw2f_all

    [W_eq_raw, beta_anc_raw, ~] = core.fit_psae(f_run, S_anc_raw, [6.0, 38.0], f_mid);
    [W_eq_dsa, beta_anc_dsa, ~] = core.fit_psae(f_run, S_anc_dsa, [6.0, 38.0], f_mid);
    W_decon_dsa = core.fit_env_anc(f_run, S_anc_dsa, P.W(jf_list), 6.0, f_mid);

    save(inv_file, 'INV3_DSA', 'INV3_RAW_FULL', 'INV3_DSA_NOW2', ...
        'W_eq_raw', 'W_eq_dsa', 'beta_anc_raw', 'beta_anc_dsa', 'W_decon_dsa', ...
        'S_anc_raw', 'S_anc_dsa', '-v7.3');
    fprintf('Saved unaliased Pass-1 priors to %s\n', inv_file);
else
    fprintf('Loading existing Pass-1 priors from %s\n', inv_file);
    IVP = load(inv_file);
    INV3_DSA = IVP.INV3_DSA; INV3_RAW_FULL = IVP.INV3_RAW_FULL; INV3_DSA_NOW2 = IVP.INV3_DSA_NOW2;
    W_eq_raw = IVP.W_eq_raw; W_eq_dsa = IVP.W_eq_dsa;
    beta_anc_raw = IVP.beta_anc_raw; beta_anc_dsa = IVP.beta_anc_dsa;
    W_decon_dsa = IVP.W_decon_dsa; S_anc_raw = IVP.S_anc_raw; S_anc_dsa = IVP.S_anc_dsa;
    clear IVP
end

% Mandatory Pass 1 Self-Check Printout: 6-44 Hz SNR curves
rm = wres_mask(:);
snr_dsa_true = 1 ./ median(double(INV3_DSA(rm, :)), 1);
snr_raw_true = 1 ./ median(double(INV3_RAW_FULL(rm, :)), 1);
snr_nw2_true = 1 ./ median(double(INV3_DSA_NOW2(rm, :)), 1);
snr_s4h_raw  = 1 ./ median(double(INV3_S4H(rm, jf_list)), 1);
snr_b1_s4h   = 1 ./ median(double(INV1_S4H(rm, jf_list)), 1);

fprintf('\n=== Pass-1 Mandatory SNR Self-Check (Reservoir Median SNR in dB) ===\n');
for fchk = [10, 20, 25, 30, 35, 40, 44]
    [~, ic] = min(abs(f_run - fchk));
    fprintf('  f=%5.2f Hz | INV3_DSA=%+6.2f dB (%.2f) | INV3_RAW_FULL=%+6.2f dB | INV3_DSA_NOW2=%+6.2f dB | S4H_aliased=%+6.2f dB\n', ...
        f_run(ic), 10*log10(snr_dsa_true(ic)), snr_dsa_true(ic), ...
        10*log10(snr_raw_true(ic)), 10*log10(snr_nw2_true(ic)), 10*log10(snr_s4h_raw(ic)));
end
fprintf('PSAE [6-38 Hz]: RAW beta=%+.5f /Hz | DSA beta=%+.5f /Hz\n', beta_anc_raw, beta_anc_dsa);

%% ========================================================================
%  PASS 2: Inject TRUE Unaliased INV3_DSA into All 18 Controlled Methods
%  ========================================================================
fprintf('\n--- Pass 2: Full 18-method imaging with TRUE unaliased INV3_DSA ---\n');
ms = [0.6, 0.85, 1.0, 1.2, 1.4]; nm = numel(ms);
s70 = single(gain_s(70));
tanc_d_s = single(tstar_anc_dsa);
tanc_f_s = single(tstar_anc_flat);

xmax_full = 1.05 * double(max(max(TS_flat(:)), max(TS_dsa(:)))) / min(ms);
xg_full   = single(linspace(0, xmax_full, 1201)); dxg_full = single(xg_full(2));
xmax_rel  = 1.05 * max(double(max(TS_flat(:)) - tanc_f_s), double(max(TS_dsa(:)) - tanc_d_s)) / min(ms) + 1e-6;
xg_rel    = single(linspace(0, xmax_rel, 1201));  dxg_rel  = single(xg_rel(2));
ntab_m2   = single(1201 - 2);
W_nw2     = (f_mid ./ f_run).^2;

for idx_f = 1:nf_run
    jf = jf_list(idx_f); f = P.f(jf);
    p2f = fullfile(p2_dir, sprintf('p2_%03d.mat', jf));
    if isfile(p2f)
        if mod(idx_f, 20) == 0 || idx_f == nf_run
            fprintf('  [Pass2 Cache] f=%5.2f Hz (%3d/%d)\n', f, idx_f, nf_run);
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

    % Exact wiring of all noise priors (explicitly tracked for Gate 6!)
    iv1_s4h  = INV1_S4H(:, jf);
    iv3_s4h  = INV3_S4H(:, jf);
    iv3_dsa  = INV3_DSA(:, idx_f);       % TRUE unaliased same-domain DSA B=300!
    iv3_raw  = INV3_RAW_FULL(:, idx_f);  % TRUE unaliased RAW B=300!
    iv3_nw2  = INV3_DSA_NOW2(:, idx_f);  % TRUE unaliased NoW2 DSA B=300!

    T1_f100 = uq_tab(f, xg_full, cfg.sig, 100.0);
    P1r_f100 = single(real(T1_f100.P)); P1i_f100 = single(imag(T1_f100.P));
    D1r_f100 = single(real(T1_f100.A - T1_f100.P)); D1i_f100 = single(imag(T1_f100.A - T1_f100.P));
    E1_f100  = single(1 + double(s70)*real(T1_f100.C));

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

    T1_ramp = uq_tab(f, xg_rel, cfg.sig, f);
    D1_ramp = single(real(T1_ramp.A) - 1);
    E1_ramp = single(1 + double(s70)*real(T1_ramp.C));

    T1_r100 = uq_tab(f, xg_rel, cfg.sig, 100.0);
    P1r_r100 = single(real(T1_r100.P)); P1i_r100 = single(imag(T1_r100.P));
    D1r_r100 = single(real(T1_r100.A - T1_r100.P)); D1i_r100 = single(imag(T1_r100.A - T1_r100.P));
    E1_r100  = single(1 + double(s70)*real(T1_r100.C));

    RAW_R0_f                    = zeros(nin, 1,  'single');
    DSA_R0_f                    = zeros(nin, 1,  'single');
    RAW_UQW1b_S4H_Unanc_100_f   = zeros(nin, nm, 'single');
    DSA_UQW1b_S4Hiv_Rel_25_f    = zeros(nin, nm, 'single');
    DSA_UQW1b_TrueDSAiv_Rel_25_f= zeros(nin, nm, 'single');
    DSA_UQW1b_TrueRAWiv_Rel_25_f= zeros(nin, nm, 'single');
    DSA_UQW1_B1_Rel_25_f        = zeros(nin, nm, 'single');
    DSA_PH_Rel_25_f             = zeros(nin, nm, 'single');
    DSA_IQ70_Rel_25_f           = zeros(nin, nm, 'single');
    DSA_W1b_TrueDSAiv_Rel_25_f  = zeros(nin, nm, 'single');
    DSA_UQW1b_TrueDSAiv_Rel_Amp_f = zeros(nin, nm, 'single');
    DSA_UQW1b_TrueDSAiv_Rel_100_f = zeros(nin, nm, 'single');
    DSA_UQW1b_TrueNW2iv_Rel_25_f  = zeros(nin, nm, 'single');
    RAW_UQW1b_TrueRAWiv_Rel_25_f  = zeros(nin, nm, 'single');

    for c = 1:ceil(ns/31)
        sh = (c-1)*31+1 : min(c*31, ns);
        Cs   = Us(:, sh) .* V_in(:, sh);
        Cr   = real(Cs); Ci = imag(Cs); clear Cs
        w_sh = w_dsa(:, sh);
        Cr_d = w_sh .* Cr; Ci_d = w_sh .* Ci;

        RAW_R0_f = RAW_R0_f + sum(Cr, 2);
        DSA_R0_f = DSA_R0_f + sum(Cr_d, 2);

        ts_f  = TS_flat(:, sh);
        ts_fr = max(ts_f - tanc_f_s, single(0));
        ts_dr = max(TS_dsa(:, sh) - tanc_d_s, single(0));

        for k = 1:nm
            invm = single(1.0 / ms(k));

            % Stage 4H exact repro (TS_flat, fref=100, iv3_s4h)
            [idx, ww] = core.tab_idx(ts_f * invm, dxg_full, ntab_m2); ow = 1 - ww;
            pc = (P1r_f100(idx).*ow + P1r_f100(idx+1).*ww).*Cr - (P1i_f100(idx).*ow + P1i_f100(idx+1).*ww).*Ci;
            dc = (D1r_f100(idx).*ow + D1r_f100(idx+1).*ww).*Cr - (D1i_f100(idx).*ow + D1i_f100(idx+1).*ww).*Ci;
            E  = E1_f100(idx).*ow + E1_f100(idx+1).*ww;
            RAW_UQW1b_S4H_Unanc_100_f(:,k) = RAW_UQW1b_S4H_Unanc_100_f(:,k) + sum(pc + dc ./ (E + iv3_s4h), 2);

            % RAW with true unaliased INV3_RAW_FULL (TS_flat_rel, fref=25)
            [idx, ww] = core.tab_idx(ts_fr * invm, dxg_rel, ntab_m2); ow = 1 - ww;
            pc = (P1r_r25(idx).*ow + P1r_r25(idx+1).*ww).*Cr - (P1i_r25(idx).*ow + P1i_r25(idx+1).*ww).*Ci;
            dc = (D1r_r25(idx).*ow + D1r_r25(idx+1).*ww).*Cr - (D1i_r25(idx).*ow + D1i_r25(idx+1).*ww).*Ci;
            E  = E1_r25(idx).*ow + E1_r25(idx+1).*ww;
            RAW_UQW1b_TrueRAWiv_Rel_25_f(:,k) = RAW_UQW1b_TrueRAWiv_Rel_25_f(:,k) + sum(pc + dc ./ (E + iv3_raw), 2);

            % All DSA relative-dt* variants (shared idx lookup!)
            [idx, ww] = core.tab_idx(ts_dr * invm, dxg_rel, ntab_m2); ow = 1 - ww;
            p0c = (P0r_r25(idx).*ow + P0r_r25(idx+1).*ww).*Cr_d - (P0i_r25(idx).*ow + P0i_r25(idx+1).*ww).*Ci_d;
            d0c = (D0r_r25(idx).*ow + D0r_r25(idx+1).*ww).*Cr_d - (D0i_r25(idx).*ow + D0i_r25(idx+1).*ww).*Ci_d;
            E0  = E0_r25(idx).*ow + E0_r25(idx+1).*ww;
            iqc = (IQr_r25(idx).*ow + IQr_r25(idx+1).*ww).*Cr_d - (IQi_r25(idx).*ow + IQi_r25(idx+1).*ww).*Ci_d;

            p1c = (P1r_r25(idx).*ow + P1r_r25(idx+1).*ww).*Cr_d - (P1i_r25(idx).*ow + P1i_r25(idx+1).*ww).*Ci_d;
            d1c = (D1r_r25(idx).*ow + D1r_r25(idx+1).*ww).*Cr_d - (D1i_r25(idx).*ow + D1i_r25(idx+1).*ww).*Ci_d;
            E1  = E1_r25(idx).*ow + E1_r25(idx+1).*ww;

            DSA_PH_Rel_25_f(:,k)             = DSA_PH_Rel_25_f(:,k)             + sum(p0c, 2);
            DSA_IQ70_Rel_25_f(:,k)           = DSA_IQ70_Rel_25_f(:,k)           + sum(iqc, 2);
            DSA_W1b_TrueDSAiv_Rel_25_f(:,k)  = DSA_W1b_TrueDSAiv_Rel_25_f(:,k)  + sum(p0c + d0c ./ (E0 + iv3_dsa), 2);
            DSA_UQW1_B1_Rel_25_f(:,k)        = DSA_UQW1_B1_Rel_25_f(:,k)        + sum(p1c + d1c ./ (E1 + iv1_s4h), 2);
            DSA_UQW1b_S4Hiv_Rel_25_f(:,k)    = DSA_UQW1b_S4Hiv_Rel_25_f(:,k)    + sum(p1c + d1c ./ (E1 + iv3_s4h), 2);
            DSA_UQW1b_TrueDSAiv_Rel_25_f(:,k)= DSA_UQW1b_TrueDSAiv_Rel_25_f(:,k)+ sum(p1c + d1c ./ (E1 + iv3_dsa), 2);
            DSA_UQW1b_TrueRAWiv_Rel_25_f(:,k)= DSA_UQW1b_TrueRAWiv_Rel_25_f(:,k)+ sum(p1c + d1c ./ (E1 + iv3_raw), 2);
            DSA_UQW1b_TrueNW2iv_Rel_25_f(:,k)= DSA_UQW1b_TrueNW2iv_Rel_25_f(:,k)+ sum(p1c + d1c ./ (E1 + iv3_nw2), 2);

            Da = D1_ramp(idx).*ow + D1_ramp(idx+1).*ww; Ea = E1_ramp(idx).*ow + E1_ramp(idx+1).*ww;
            DSA_UQW1b_TrueDSAiv_Rel_Amp_f(:,k) = DSA_UQW1b_TrueDSAiv_Rel_Amp_f(:,k) + sum((1 + Da ./ (Ea + iv3_dsa)) .* Cr_d, 2);

            phc = (P1r_r100(idx).*ow + P1r_r100(idx+1).*ww).*Cr_d - (P1i_r100(idx).*ow + P1i_r100(idx+1).*ww).*Ci_d;
            dhc = (D1r_r100(idx).*ow + D1r_r100(idx+1).*ww).*Cr_d - (D1i_r100(idx).*ow + D1i_r100(idx+1).*ww).*Ci_d;
            Eh  = E1_r100(idx).*ow + E1_r100(idx+1).*ww;
            DSA_UQW1b_TrueDSAiv_Rel_100_f(:,k) = DSA_UQW1b_TrueDSAiv_Rel_100_f(:,k) + sum(phc + dhc ./ (Eh + iv3_dsa), 2);
        end
    end
    clear Us V_in Cr Ci Cr_d Ci_d

    save(p2f, 'f', 'RAW_R0_f', 'DSA_R0_f', 'RAW_UQW1b_S4H_Unanc_100_f', ...
        'DSA_UQW1b_S4Hiv_Rel_25_f', 'DSA_UQW1b_TrueDSAiv_Rel_25_f', 'DSA_UQW1b_TrueRAWiv_Rel_25_f', ...
        'DSA_UQW1_B1_Rel_25_f', 'DSA_PH_Rel_25_f', 'DSA_IQ70_Rel_25_f', ...
        'DSA_W1b_TrueDSAiv_Rel_25_f', 'DSA_UQW1b_TrueDSAiv_Rel_Amp_f', 'DSA_UQW1b_TrueDSAiv_Rel_100_f', ...
        'DSA_UQW1b_TrueNW2iv_Rel_25_f', 'RAW_UQW1b_TrueRAWiv_Rel_25_f', '-v6');
    if mod(idx_f, 10) == 0 || idx_f == nf_run || idx_f <= 3
        fprintf('  Pass2 f=%5.2f Hz (%3d/%d): total=%.1fs (helm=%.1fs, op+io=%.1fs)\n', ...
            f, idx_f, nf_run, toc(t0), t_helm, toc(t0) - t_helm);
    end
end
clear w_dsa TS_flat TS_dsa

%% ========================================================================
%  STACK ALL 18 METHODS & EVALUATE (Full F-5, F-14, L0, Lag, G_HF, Slope)
%  ========================================================================
mkeys = {'RAW_R0', 'DSA_R0', 'DSA_R0_Eq', 'DSA_R0_NoW2', ...
         'RAW_UQW1b_S4H_Unanc_100', 'DSA_UQW1b_S4Hiv_PAnc_25', ...
         'DSA_UQW1b_TrueDSAiv_PAnc_25', 'DSA_UQW1b_TrueRAWiv_PAnc_25', 'DSA_UQW1_B1_PAnc_25', ...
         'DSA_PH_PAnc_25', 'DSA_IQ70_PAnc_25', 'DSA_W1b_TrueDSAiv_PAnc_25', ...
         'DSA_UQW1b_TrueDSAiv_PAnc_Amp', 'DSA_UQW1b_TrueDSAiv_PAnc_100', ...
         'DSA_UQW1b_TrueDSAiv_NoW2_25', 'DSA_UQW1b_TrueDSAiv_Decon_25', ...
         'RAW_UQW1b_TrueRAWiv_PAnc_25', 'RAW_UQW1b_SameWeq_PAnc_25'};
r0_only = {'RAW_R0', 'DSA_R0', 'DSA_R0_Eq', 'DSA_R0_NoW2'};

Img = struct(); spec_res = struct();
for ik = 1:numel(mkeys)
    nc = 1 + (nm-1)*~ismember(mkeys{ik}, r0_only);
    Img.(mkeys{ik}) = zeros(nin, nc);
    spec_res.(mkeys{ik}) = zeros(1, nf_run);
end

for idx_f = 1:nf_run
    D = load(fullfile(p2_dir, sprintf('p2_%03d.mat', jf_list(idx_f))));
    wr = W_eq_raw(idx_f); wd = W_eq_dsa(idx_f); wn = W_nw2(idx_f); wdec = W_decon_dsa(idx_f);

    Img.RAW_R0                       = Img.RAW_R0                       + double(D.RAW_R0_f);
    Img.DSA_R0                       = Img.DSA_R0                       + double(D.DSA_R0_f);
    Img.DSA_R0_Eq                    = Img.DSA_R0_Eq                    + wd   * double(D.DSA_R0_f);
    Img.DSA_R0_NoW2                  = Img.DSA_R0_NoW2                  + wn   * double(D.DSA_R0_f);
    Img.RAW_UQW1b_S4H_Unanc_100      = Img.RAW_UQW1b_S4H_Unanc_100      + double(D.RAW_UQW1b_S4H_Unanc_100_f);
    Img.DSA_UQW1b_S4Hiv_PAnc_25      = Img.DSA_UQW1b_S4Hiv_PAnc_25      + wd   * double(D.DSA_UQW1b_S4Hiv_Rel_25_f);
    Img.DSA_UQW1b_TrueDSAiv_PAnc_25  = Img.DSA_UQW1b_TrueDSAiv_PAnc_25  + wd   * double(D.DSA_UQW1b_TrueDSAiv_Rel_25_f);
    Img.DSA_UQW1b_TrueRAWiv_PAnc_25  = Img.DSA_UQW1b_TrueRAWiv_PAnc_25  + wd   * double(D.DSA_UQW1b_TrueRAWiv_Rel_25_f);
    Img.DSA_UQW1_B1_PAnc_25          = Img.DSA_UQW1_B1_PAnc_25          + wd   * double(D.DSA_UQW1_B1_Rel_25_f);
    Img.DSA_PH_PAnc_25               = Img.DSA_PH_PAnc_25               + wd   * double(D.DSA_PH_Rel_25_f);
    Img.DSA_IQ70_PAnc_25             = Img.DSA_IQ70_PAnc_25             + wd   * double(D.DSA_IQ70_Rel_25_f);
    Img.DSA_W1b_TrueDSAiv_PAnc_25    = Img.DSA_W1b_TrueDSAiv_PAnc_25    + wd   * double(D.DSA_W1b_TrueDSAiv_Rel_25_f);
    Img.DSA_UQW1b_TrueDSAiv_PAnc_Amp = Img.DSA_UQW1b_TrueDSAiv_PAnc_Amp + wd   * double(D.DSA_UQW1b_TrueDSAiv_Rel_Amp_f);
    Img.DSA_UQW1b_TrueDSAiv_PAnc_100 = Img.DSA_UQW1b_TrueDSAiv_PAnc_100 + wd   * double(D.DSA_UQW1b_TrueDSAiv_Rel_100_f);
    Img.DSA_UQW1b_TrueDSAiv_NoW2_25  = Img.DSA_UQW1b_TrueDSAiv_NoW2_25  + wn   * double(D.DSA_UQW1b_TrueNW2iv_Rel_25_f);
    Img.DSA_UQW1b_TrueDSAiv_Decon_25 = Img.DSA_UQW1b_TrueDSAiv_Decon_25 + wdec * double(D.DSA_UQW1b_TrueDSAiv_Rel_25_f);
    Img.RAW_UQW1b_TrueRAWiv_PAnc_25  = Img.RAW_UQW1b_TrueRAWiv_PAnc_25  + wr   * double(D.RAW_UQW1b_TrueRAWiv_Rel_25_f);
    Img.RAW_UQW1b_SameWeq_PAnc_25    = Img.RAW_UQW1b_SameWeq_PAnc_25    + wd   * double(D.RAW_UQW1b_TrueRAWiv_Rel_25_f);

    spec_res.RAW_R0(idx_f)                       = rms(double(D.RAW_R0_f(rm)));
    spec_res.DSA_R0(idx_f)                       = rms(double(D.DSA_R0_f(rm)));
    spec_res.DSA_R0_Eq(idx_f)                    = rms(wd   * double(D.DSA_R0_f(rm)));
    spec_res.DSA_R0_NoW2(idx_f)                  = rms(wn   * double(D.DSA_R0_f(rm)));
    spec_res.RAW_UQW1b_S4H_Unanc_100(idx_f)      = rms(double(D.RAW_UQW1b_S4H_Unanc_100_f(rm, 3)));
    spec_res.DSA_UQW1b_S4Hiv_PAnc_25(idx_f)      = rms(wd   * double(D.DSA_UQW1b_S4Hiv_Rel_25_f(rm, 3)));
    spec_res.DSA_UQW1b_TrueDSAiv_PAnc_25(idx_f)  = rms(wd   * double(D.DSA_UQW1b_TrueDSAiv_Rel_25_f(rm, 3)));
    spec_res.DSA_UQW1b_TrueRAWiv_PAnc_25(idx_f)  = rms(wd   * double(D.DSA_UQW1b_TrueRAWiv_Rel_25_f(rm, 3)));
    spec_res.DSA_UQW1_B1_PAnc_25(idx_f)          = rms(wd   * double(D.DSA_UQW1_B1_Rel_25_f(rm, 3)));
    spec_res.DSA_PH_PAnc_25(idx_f)               = rms(wd   * double(D.DSA_PH_Rel_25_f(rm, 3)));
    spec_res.DSA_IQ70_PAnc_25(idx_f)             = rms(wd   * double(D.DSA_IQ70_Rel_25_f(rm, 3)));
    spec_res.DSA_W1b_TrueDSAiv_PAnc_25(idx_f)    = rms(wd   * double(D.DSA_W1b_TrueDSAiv_Rel_25_f(rm, 3)));
    spec_res.DSA_UQW1b_TrueDSAiv_PAnc_Amp(idx_f) = rms(wd   * double(D.DSA_UQW1b_TrueDSAiv_Rel_Amp_f(rm, 3)));
    spec_res.DSA_UQW1b_TrueDSAiv_PAnc_100(idx_f) = rms(wd   * double(D.DSA_UQW1b_TrueDSAiv_Rel_100_f(rm, 3)));
    spec_res.DSA_UQW1b_TrueDSAiv_NoW2_25(idx_f)  = rms(wn   * double(D.DSA_UQW1b_TrueNW2iv_Rel_25_f(rm, 3)));
    spec_res.DSA_UQW1b_TrueDSAiv_Decon_25(idx_f) = rms(wdec * double(D.DSA_UQW1b_TrueDSAiv_Rel_25_f(rm, 3)));
    spec_res.RAW_UQW1b_TrueRAWiv_PAnc_25(idx_f)  = rms(wr   * double(D.RAW_UQW1b_TrueRAWiv_Rel_25_f(rm, 3)));
    spec_res.RAW_UQW1b_SameWeq_PAnc_25(idx_f)    = rms(wd   * double(D.RAW_UQW1b_TrueRAWiv_Rel_25_f(rm, 3)));
end

r0_repro = NaN; uqw1b_repro = NaN;
if ~smoke
    r0_repro    = norm(Img.RAW_R0 - R0_stage4h) / max(norm(R0_stage4h), 1e-30);
    uqw1b_repro = norm(Img.RAW_UQW1b_S4H_Unanc_100(:,3) - UQW1b_stage4h_m1) / max(norm(UQW1b_stage4h_m1), 1e-30);
    fprintf('Exact Stage 4H repro: RAW_R0=%.2e | RAW_UQW1b_S4H_Unanc_100(m=1)=%.2e\n', r0_repro, uqw1b_repro);
end

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
f_hf = (f_run >= 28.0 & f_run <= 40.0);

res = struct('dataset', 'volve_4_45hz_closure', 'smoke', smoke, 'n_freq', nf_run, 'f', f_run, 'm', ms, ...
             'v_r', vr, 'kfix', kfix, 'kr0', kr0, ...
             'repro_R0_relL2', r0_repro, 'repro_UQW1b_relL2', uqw1b_repro, ...
             'tstar_anc_dsa_ms', tstar_anc_dsa*1e3, 'tstar_anc_flat_ms', tstar_anc_flat*1e3, ...
             'beta_anc_raw', beta_anc_raw, 'beta_anc_dsa', beta_anc_dsa, ...
             'S_anc_raw', S_anc_raw, 'S_anc_dsa', S_anc_dsa, ...
             'snr_med_b1_s4h', snr_b1_s4h, 'snr_med_b300_s4h', snr_s4h_raw, ...
             'snr_med_b300_raw_full', snr_raw_true, 'snr_med_b300_dsa_full', snr_dsa_true, ...
             'snr_med_b300_nw2_full', snr_nw2_true, ...
             'spec_res', spec_res, 'res_slope', struct(), 'G_HF_28_40', struct(), 'methods', struct());

for ik = 1:numel(mkeys)
    k = mkeys{ik}; sv = spec_res.(k);
    pf = polyfit(f_run(fm), log(max(sv(fm), 1e-30)), 1);
    res.res_slope.(k) = pf(1);
    res.G_HF_28_40.(k) = median(sv(f_hf) ./ max(spec_res.DSA_R0_Eq(f_hf), 1e-30));
end

fprintf('\n=== WP2 VOLVE 4-45 Hz TRUE CONFLUENCE RESULTS (18 Methods) ===\n');
for ik = 1:numel(mkeys)
    k = mkeys{ik}; M = Img.(k); nc = size(M, 2);
    f5_fix = zeros(1,nc); f5_r0 = zeros(1,nc); f5_l0 = zeros(1,nc); f5_lag = zeros(1,nc);
    f14_fix = zeros(1,nc); f14_l0 = zeros(1,nc); f14_lag = zeros(1,nc);
    cen = zeros(1,nc); bw10 = zeros(nc,2); smf = zeros(1,nc); smf25 = zeros(1,nc);
    for c = 1:nc
        A = core.prep_image(M(:,c), P.nzi, P.nxi, z_in, P.dx, cfg.bg_sigma_m, cfg.z_mute_m);
        [f5_fix(c),  f5_lag(c),  f5_l0(c) ] = core.well_tie_full(A, w5,  rb5_fix,  P.dx, kfix(1), kfix(2), 20.0, 30.0);
        [f5_r0(c),   ~,          ~        ] = core.well_tie_full(A, w5,  rb5_r0,   P.dx, kr0(1),  kr0(2),  20.0, 30.0);
        [f14_fix(c), f14_lag(c), f14_l0(c)] = core.well_tie_full(A, w14, rb14_fix, P.dx, kfix(1), kfix(2), 20.0, 30.0);
        [bw10(c,:), cen(c)] = core.res_band(A(wz_res, wx_res), P.dx, 512);
        smf(c)   = core.steep_smile_frac(A(wz_mid, wx_res), P.dx, 25.0, 3.0, 40.0);
        smf25(c) = core.steep_smile_frac(A(wz_mid, wx_res), P.dx, 25.0, 3.0, 25.0);
    end
    res.methods.(k) = struct('F5_fix', f5_fix, 'F5_fix_mean', mean(f5_fix), ...
        'F5_r0', f5_r0, 'F5_r0_mean', mean(f5_r0), ...
        'F5_ZeroLag', f5_l0, 'F5_BestLag_m', f5_lag, ...
        'F14_fix', f14_fix, 'F14_fix_mean', mean(f14_fix), ...
        'F14_ZeroLag', f14_l0, 'F14_BestLag_m', f14_lag, ...
        'Res_Cen', cen, 'Res_BW10', bw10, 'Smile_Frac', smf, 'Smile_Frac_25', smf25);
    c3 = min(3, nc);
    fprintf('  %-29s | fix=%.4f (L0=%+.3f,%+3.0fm) | r0=%.4f | F14=%+.3f | Sm=%.3f | sl=%+.4f | G_HF=%.2fx\n', ...
        k, mean(f5_fix), f5_l0(c3), f5_lag(c3), mean(f5_r0), mean(f14_fix), smf(c3), res.res_slope.(k), res.G_HF_28_40.(k));
end
res.elapsed_s = toc(t_all);

tag = ''; if smoke, tag = '_smoke'; end
fn_json = fullfile(res_dir, ['volve45_closure' tag '.json']);
fid = fopen(fn_json, 'w'); fwrite(fid, jsonencode(res, 'PrettyPrint', true), 'char'); fclose(fid);
fprintf('\nSaved %s (total %.1fs)\n', fn_json, res.elapsed_s);

if ~smoke
    Img_s = struct();
    for ik = 1:numel(mkeys), Img_s.(mkeys{ik}) = single(Img.(mkeys{ik})); end
    save(fullfile(scratch_dir, 'volve45_closure_imgs.mat'), 'Img_s', 'res', 'z_in', 'x_in', 'phi_eff', '-v7.3');
    fprintf('Saved final 4-45 Hz closure 2D images to %s\n', scratch_dir);
    try rmdir(p1_dir, 's'); rmdir(p2_dir, 's'); fprintf('Cleaned temp Pass1/Pass2 slice dirs.\n'); catch, end
end
end

function s = ternary(cond, a, b)
if cond, s = a; else, s = b; end
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
