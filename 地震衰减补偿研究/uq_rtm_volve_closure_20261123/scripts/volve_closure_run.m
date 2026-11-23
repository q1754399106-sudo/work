function volve_closure_run(smoke)
% Pre-registered Real Volve North Sea OBC Closure Experiment (uq_rtm_volve_closure_20261123/PLAN_ZH.md)
% Strictly read-only on uq_rtm_20261117/cache/fld (never touches fld45).
% 100% blind to well logs (F-5 / F-14 are loaded ONLY at the very end for scoring!).
% Single 85-frequency Helmholtz pass (~35-45 min on 8 threads, < 2.2 GB RAM), resumable per frequency.
%
% Usage:
%   volve_closure_run(true)   % 3-frequency syntax & timing smoke test
%   volve_closure_run(false)  % full 85-frequency pre-registered production run

if nargin < 1, smoke = false; end
maxNumCompThreads(8);
t_all = tic;

here = fileparts(mfilename('fullpath'));
root = fileparts(here);
res_dir = fullfile(root, 'results');
if ~isfolder(res_dir), mkdir(res_dir); end

% Keep all heavy intermediate slices strictly outside git worktree!
scratch_base = 'C:\Users\qq175\.gemini\antigravity\brain\c0dc2750-de05-424f-adb4-e7ef332a1950\scratch\volve_closure';
if smoke
    slice_dir = fullfile(scratch_base, 'slices_smoke_v2');
else
    slice_dir = fullfile(scratch_base, 'slices_prod');
end
if ~isfolder(slice_dir), mkdir(slice_dir); end

addpath('D:\ChatGPT\研究\地震衰减补偿研究\uq_rtm_20261117\scripts');
P = fld_setup();
ns = P.nshot; N = P.nz * P.nx;
[Zi, Xi] = ndgrid(P.in{1}, P.in{2});
lin = sub2ind([P.nz P.nx], Zi(:), Xi(:));
nin = numel(lin);
irec = sub2ind([P.nz P.nx], P.rec_iz, P.rec_ix);
isrc = sub2ind([P.nz P.nx], P.src_iz, P.src_ix);
col  = sub2ind([N ns], irec, P.shot_of);
z_in = P.z(P.in{1}); x_in = P.x(P.in{2});
[Zin, Xin] = ndgrid(z_in, x_in);

% 1. Load Stage-4 band-summed inv.mat (exact apples-to-apples noise prior) & cached R0 for blind dip
IV_S = load(fullfile(P.cache, 'inv.mat'), 'iv');
INV = IV_S.iv; clear IV_S

S1 = load(fullfile(P.cache, 'state1.mat'), 'I');
R0_cached_vec = double(sum(S1.I, 2));
clear S1
R0_2d = prep_volve(R0_cached_vec, P.nzi, P.nxi, z_in);
[phi_eff, c_st] = blind_structure_dip(R0_2d, P.dx, 6.0);

wres_mask = (Zin >= 2600 & Zin <= 3200) & (Xin >= -2500 & Xin <= 1500);
anc_mask  = (Zin >= 1000 & Zin <= 1600) & (Xin >= -2500 & Xin <= 1500);

fprintf('=== Volve Closure (smoke=%d, Stage-4 band-summed inv.mat) ===\n', smoke);
fprintf('Blind dip reservoir: med=%.2f deg, p90=%.2f deg, med_coh=%.3f\n', ...
    median(abs(rad2deg(phi_eff(wres_mask)))), prctile(abs(rad2deg(phi_eff(wres_mask))), 90), median(c_st(wres_mask)));

% 2. Build pre-stack DSA mask (theta_max=50 deg, Nyquist shot-aliasing limit) & dip-aware TS_dsa from cached tfields.mat
rx = P.x(P.rec_ix);
xr_lo = arrayfun(@(s) min(rx(P.shot_of == s)), 1:ns);
xr_hi = arrayfun(@(s) max(rx(P.shot_of == s)), 1:ns);

TF = load(fullfile(P.cache, 'tfields.mat'), 'F');
TS_flat_S = load(fullfile(P.cache, 'ts.mat'), 'TS');
TS_flat = TS_flat_S.TS; clear TS_flat_S

[xs_sort, ord] = sort(P.src_x(:));
w_dsa  = zeros(nin, ns, 'single');
TS_dsa = zeros(nin, ns, 'single');
xv = Xin(:); zv = Zin(:); phiv = phi_eff(:);
zi_v = Zi(:); xi_v = Xi(:);
dz_r = max(zv - 90.0, 10.0);

for k = 1:ns
    dz_s = max(zv - P.z(P.src_iz(k)), 10.0);
    th_s = atan2(xv - P.src_x(k), dz_s);
    th_r = max(min(th_s + 2*phiv, deg2rad(78)), -deg2rad(78));
    xmir = xv + dz_r .* tan(th_r);
    d_in = min(xmir - xr_lo(k), xr_hi(k) - xmir);
    t_arr = cos_taper(d_in, -150.0, 250.0);
    t_gs  = cos_taper(50.0 - abs(rad2deg(th_s)), 0.0, 10.0);
    t_gr  = cos_taper(50.0 - abs(rad2deg(th_r)), 0.0, 10.0);
    w_dsa(:, k) = single(max(t_arr .* t_gs .* t_gr, 1e-4));

    q = interp1(xs_sort, 1:ns, min(max(xmir, xs_sort(1)), xs_sort(end)));
    k0 = min(floor(q), ns - 1); ww = single(q - k0);
    a = sub2ind([P.nz P.nx ns], zi_v, xi_v, ord(k0));
    b = sub2ind([P.nz P.nx ns], zi_v, xi_v, ord(k0 + 1));
    TS_dsa(:, k) = TF.F(sub2ind([P.nz P.nx ns], zi_v, xi_v, k*ones(nin, 1))) + (1 - ww).*TF.F(a) + ww.*TF.F(b);
end
clear TF

tstar_anc_dsa  = median(double(TS_dsa(anc_mask(:), :)), 'all');
tstar_anc_flat = median(double(TS_flat(anc_mask(:), :)), 'all');
TS_dsa_rel  = max(TS_dsa  - single(tstar_anc_dsa),  single(0));
TS_flat_rel = max(TS_flat - single(tstar_anc_flat), single(0));
fprintf('Shallow anchor [1000-1600m] med t*: DSA=%.2f ms, flat=%.2f ms | Res fold med=%.1f / %d\n', ...
    tstar_anc_dsa*1e3, tstar_anc_flat*1e3, median(sum(double(w_dsa(wres_mask(:), :)), 2)), ns);

% 3. Single-pass frequency loop (resumable per frequency, fast -v6 slice cache!)
if smoke
    jf_list = [9, 43, 73];  % 6.0, 14.5, 22.0 Hz
else
    jf_list = 1:numel(P.f);
end
nf_run = numel(jf_list);
ms = [0.5, 0.7, 1.0, 1.4, 2.0]; nm = numel(ms);
s70 = single(gain_s(70));
f_mid = 14.5;

xmax_full = 1.05 * double(max(max(TS_flat(:)), max(TS_dsa(:)))) / min(ms);
xg_full   = single(linspace(0, xmax_full, 1201)); dxg_full = single(xg_full(2));
xmax_rel  = 1.05 * double(max(max(TS_flat_rel(:)), max(TS_dsa_rel(:)))) / min(ms) + 1e-6;
xg_rel    = single(linspace(0, xmax_rel, 1201));  dxg_rel  = single(xg_rel(2));
ntab_m2   = single(1201 - 2);

for idx_f = 1:nf_run
    jf = jf_list(idx_f);
    f = P.f(jf);
    sfile = fullfile(slice_dir, sprintf('slice_%03d.mat', jf));
    if isfile(sfile)
        if mod(idx_f, 15) == 0 || idx_f == nf_run
            fprintf('  [Cache hit] f=%5.2f Hz (%d/%d)\n', f, idx_f, nf_run);
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

    iv_f = INV(:, jf);

    % Operator lookup tables (fref=f amp-only for primary; fref=14.5 for Disp15 check)
    T1_full = uq_tab(f, xg_full, P.sig, f);
    A1_f = single(real(T1_full.A)); C1_f = single(real(T1_full.C));

    T0_rel  = uq_tab(f, xg_rel, 0,     f);
    T1_rel  = uq_tab(f, xg_rel, P.sig, f);
    A0_r = single(real(T0_rel.A)); C0_r = single(real(T0_rel.C));
    A1_r = single(real(T1_rel.A)); C1_r = single(real(T1_rel.C));

    T1_d15  = uq_tab(f, xg_rel, P.sig, f_mid);
    A1_d = single(T1_d15.A); C1_d = single(T1_d15.C); P1_d = single(T1_d15.P);

    RAW_R0_f              = zeros(nin, 1,  'single');
    DSA_R0_f              = zeros(nin, 1,  'single');
    RAW_UQW1_Unanc_f      = zeros(nin, nm, 'single');
    DSA_UQW1_Unanc_f      = zeros(nin, nm, 'single');
    RAW_UQW1_Rel_f        = zeros(nin, nm, 'single');
    DSA_IQ70_Rel_f        = zeros(nin, nm, 'single');
    DSA_W1_Rel_f          = zeros(nin, nm, 'single');
    DSA_UQW1_Rel_f        = zeros(nin, nm, 'single');
    DSA_UQW1_Rel_Disp15_f = zeros(nin, nm, 'single');

    % Single shot-chunk pass (note: for fref=f, P0=P1=1 real, so H*real(Cs) == real(H*Cs) -> 2x faster!)
    for c = 1:ceil(ns/31)
        sh = (c-1)*31+1 : min(c*31, ns);
        Cs   = Us(:, sh) .* V_in(:, sh);
        Cr   = real(Cs);
        w_sh = w_dsa(:, sh);
        Cr_d = w_sh .* Cr;
        Cs_d = w_sh .* Cs;

        RAW_R0_f = RAW_R0_f + sum(Cr, 2);
        DSA_R0_f = DSA_R0_f + sum(Cr_d, 2);

        ts_f_sh  = TS_flat(:, sh);
        ts_d_sh  = TS_dsa(:, sh);
        ts_fr_sh = TS_flat_rel(:, sh);
        ts_dr_sh = TS_dsa_rel(:, sh);

        for k = 1:nm
            inv_m = single(1.0 / ms(k));

            % 1) RAW_UQW1_Unanc (real H, P=1)
            [idx, ww] = tab_idx(ts_f_sh * inv_m, dxg_full, ntab_m2);
            A = A1_f(idx).*(1-ww) + A1_f(idx+1).*ww;
            C = C1_f(idx).*(1-ww) + C1_f(idx+1).*ww;
            H = (A + s70*C + iv_f) ./ (1 + s70*C + iv_f);
            RAW_UQW1_Unanc_f(:, k) = RAW_UQW1_Unanc_f(:, k) + sum(H .* Cr, 2);

            % 2) DSA_UQW1_Unanc (real H, P=1)
            [idx, ww] = tab_idx(ts_d_sh * inv_m, dxg_full, ntab_m2);
            A = A1_f(idx).*(1-ww) + A1_f(idx+1).*ww;
            C = C1_f(idx).*(1-ww) + C1_f(idx+1).*ww;
            H = (A + s70*C + iv_f) ./ (1 + s70*C + iv_f);
            DSA_UQW1_Unanc_f(:, k) = DSA_UQW1_Unanc_f(:, k) + sum(H .* Cr_d, 2);

            % 3) RAW_UQW1_Rel (real H, P=1)
            [idx, ww] = tab_idx(ts_fr_sh * inv_m, dxg_rel, ntab_m2);
            A = A1_r(idx).*(1-ww) + A1_r(idx+1).*ww;
            C = C1_r(idx).*(1-ww) + C1_r(idx+1).*ww;
            H = (A + s70*C + iv_f) ./ (1 + s70*C + iv_f);
            RAW_UQW1_Rel_f(:, k) = RAW_UQW1_Rel_f(:, k) + sum(H .* Cr, 2);

            % 4-7) Shared index lookup for TS_dsa_rel: IQ70, W1, UQW1, and Disp15!
            [idx, ww] = tab_idx(ts_dr_sh * inv_m, dxg_rel, ntab_m2);
            ow = 1 - ww;
            A0 = A0_r(idx).*ow + A0_r(idx+1).*ww;
            C0 = C0_r(idx).*ow + C0_r(idx+1).*ww;
            A1 = A1_r(idx).*ow + A1_r(idx+1).*ww;
            C1 = C1_r(idx).*ow + C1_r(idx+1).*ww;

            DSA_IQ70_Rel_f(:, k) = DSA_IQ70_Rel_f(:, k) + sum(((A0 + s70*C0) ./ (1 + s70*C0)) .* Cr_d, 2);
            DSA_W1_Rel_f(:, k)   = DSA_W1_Rel_f(:, k)   + sum(((A0 + s70*C0 + iv_f) ./ (1 + s70*C0 + iv_f)) .* Cr_d, 2);
            DSA_UQW1_Rel_f(:, k) = DSA_UQW1_Rel_f(:, k) + sum(((A1 + s70*C1 + iv_f) ./ (1 + s70*C1 + iv_f)) .* Cr_d, 2);

            Ad = A1_d(idx).*ow + A1_d(idx+1).*ww;
            Cd = C1_d(idx).*ow + C1_d(idx+1).*ww;
            Pd = P1_d(idx).*ow + P1_d(idx+1).*ww;
            H_d15 = (Ad + s70*Pd.*Cd + Pd.*iv_f) ./ (1 + s70*Cd + iv_f);
            DSA_UQW1_Rel_Disp15_f(:, k) = DSA_UQW1_Rel_Disp15_f(:, k) + real(sum(H_d15 .* Cs_d, 2));
        end
    end
    clear Us V_in

    save(sfile, 'f', 'RAW_R0_f', 'DSA_R0_f', 'RAW_UQW1_Unanc_f', 'DSA_UQW1_Unanc_f', ...
        'RAW_UQW1_Rel_f', 'DSA_IQ70_Rel_f', 'DSA_W1_Rel_f', 'DSA_UQW1_Rel_f', ...
        'DSA_UQW1_Rel_Disp15_f', '-v6');
    fprintf('  f=%5.2f Hz (%2d/%d): total=%.1fs (helm=%.1fs, op+io=%.1fs)\n', ...
        f, idx_f, nf_run, toc(t0), t_helm, toc(t0) - t_helm);
end

clear w_dsa TS_flat TS_dsa TS_flat_rel TS_dsa_rel INV

% 4. Blind Parametric Smooth Shallow-Slope Anchoring (PSAE) over [6, 22] Hz in [1000, 1600] m
f_run = P.f(jf_list);
S_anc_raw = zeros(1, nf_run); S_anc_dsa = zeros(1, nf_run);
for idx_f = 1:nf_run
    jf = jf_list(idx_f);
    D = load(fullfile(slice_dir, sprintf('slice_%03d.mat', jf)), 'RAW_R0_f', 'DSA_R0_f');
    S_anc_raw(idx_f) = rms(double(D.RAW_R0_f(anc_mask(:))));
    S_anc_dsa(idx_f) = rms(double(D.DSA_R0_f(anc_mask(:))));
end

fm = (f_run >= 6.0 & f_run <= 22.0);
p_raw = polyfit(f_run(fm), log(max(S_anc_raw(fm), 1e-30)), 1); beta_anc_raw = p_raw(1);
p_dsa = polyfit(f_run(fm), log(max(S_anc_dsa(fm), 1e-30)), 1); beta_anc_dsa = p_dsa(1);
W_eq_raw = exp(-beta_anc_raw * (f_run - f_mid));
W_eq_dsa = exp(-beta_anc_dsa * (f_run - f_mid));
fprintf('\n=== Blind PSAE Shallow [1000-1600m] Slopes ([6-22 Hz]) ===\n');
fprintf('  RAW beta_anc = %+.5f /Hz (eff_t* = %.2f ms)\n', beta_anc_raw, -beta_anc_raw/pi*1e3);
fprintf('  DSA beta_anc = %+.5f /Hz (eff_t* = %.2f ms)\n', beta_anc_dsa, -beta_anc_dsa/pi*1e3);

% 5. Stack final 2D images across frequency & record reservoir spectra
Img = struct();
Img.RAW_R0               = zeros(nin, 1);
Img.DSA_R0               = zeros(nin, 1);
Img.DSA_R0_Eq            = zeros(nin, 1);
Img.RAW_UQW1_Unanc       = zeros(nin, nm);
Img.DSA_UQW1_Unanc       = zeros(nin, nm);
Img.RAW_UQW1_PAnc        = zeros(nin, nm);
Img.DSA_IQ70_PAnc        = zeros(nin, nm);
Img.DSA_W1_PAnc          = zeros(nin, nm);
Img.DSA_UQW1_PAnc        = zeros(nin, nm);
Img.DSA_UQW1_PAnc_Disp15 = zeros(nin, nm);

spec_res = struct();
spec_keys = {'RAW_R0', 'DSA_R0', 'DSA_R0_Eq', 'RAW_UQW1_Unanc_m1', 'DSA_UQW1_Unanc_m1', ...
             'RAW_UQW1_PAnc_m1', 'DSA_IQ70_PAnc_m1', 'DSA_W1_PAnc_m1', 'DSA_UQW1_PAnc_m1', 'DSA_UQW1_PAnc_Disp15_m1'};
for sk = spec_keys, spec_res.(sk{1}) = zeros(1, nf_run); end
spec_anc_eq = zeros(1, nf_run);

for idx_f = 1:nf_run
    jf = jf_list(idx_f);
    D = load(fullfile(slice_dir, sprintf('slice_%03d.mat', jf)));
    wr = W_eq_raw(idx_f); wd = W_eq_dsa(idx_f);

    Img.RAW_R0               = Img.RAW_R0               + double(D.RAW_R0_f);
    Img.DSA_R0               = Img.DSA_R0               + double(D.DSA_R0_f);
    Img.DSA_R0_Eq            = Img.DSA_R0_Eq            + wd * double(D.DSA_R0_f);
    Img.RAW_UQW1_Unanc       = Img.RAW_UQW1_Unanc       + double(D.RAW_UQW1_Unanc_f);
    Img.DSA_UQW1_Unanc       = Img.DSA_UQW1_Unanc       + double(D.DSA_UQW1_Unanc_f);
    Img.RAW_UQW1_PAnc        = Img.RAW_UQW1_PAnc        + wr * double(D.RAW_UQW1_Rel_f);
    Img.DSA_IQ70_PAnc        = Img.DSA_IQ70_PAnc        + wd * double(D.DSA_IQ70_Rel_f);
    Img.DSA_W1_PAnc          = Img.DSA_W1_PAnc          + wd * double(D.DSA_W1_Rel_f);
    Img.DSA_UQW1_PAnc        = Img.DSA_UQW1_PAnc        + wd * double(D.DSA_UQW1_Rel_f);
    Img.DSA_UQW1_PAnc_Disp15 = Img.DSA_UQW1_PAnc_Disp15 + wd * double(D.DSA_UQW1_Rel_Disp15_f);

    rm = wres_mask(:);
    spec_anc_eq(idx_f)                = rms(wd * double(D.DSA_R0_f(anc_mask(:))));
    spec_res.RAW_R0(idx_f)            = rms(double(D.RAW_R0_f(rm)));
    spec_res.DSA_R0(idx_f)            = rms(double(D.DSA_R0_f(rm)));
    spec_res.DSA_R0_Eq(idx_f)         = rms(wd * double(D.DSA_R0_f(rm)));
    spec_res.RAW_UQW1_Unanc_m1(idx_f) = rms(double(D.RAW_UQW1_Unanc_f(rm, 3)));
    spec_res.DSA_UQW1_Unanc_m1(idx_f) = rms(double(D.DSA_UQW1_Unanc_f(rm, 3)));
    spec_res.RAW_UQW1_PAnc_m1(idx_f)  = rms(wr * double(D.RAW_UQW1_Rel_f(rm, 3)));
    spec_res.DSA_IQ70_PAnc_m1(idx_f)  = rms(wd * double(D.DSA_IQ70_Rel_f(rm, 3)));
    spec_res.DSA_W1_PAnc_m1(idx_f)    = rms(wd * double(D.DSA_W1_Rel_f(rm, 3)));
    spec_res.DSA_UQW1_PAnc_m1(idx_f)  = rms(wd * double(D.DSA_UQW1_Rel_f(rm, 3)));
    spec_res.DSA_UQW1_PAnc_Disp15_m1(idx_f) = rms(wd * double(D.DSA_UQW1_Rel_Disp15_f(rm, 3)));
end

if ~smoke
    r0_rel_err = norm(Img.RAW_R0 - R0_cached_vec) / max(norm(R0_cached_vec), 1e-30);
    fprintf('Exact RAW_R0 reproducibility check vs cache/fld/state1.mat: rel L2 diff = %.2e\n', r0_rel_err);
end

% 6. Blind evaluation on F-5 (exact fld_eval.m R0-band & broadband), F-14, spectral slope, smile energy, bandwidth
w5  = load_well(fullfile(P.cache, 'well_F-5.bin'),  z_in, x_in);
w14 = load_well(fullfile(P.cache, 'well_F-14.bin'), z_in, x_in);
lo_r0 = 3.515625; hi_r0 = 12.890625; hi_br = 16.0;
rb5_r0  = bp_z(w5.rw,  P.dx, lo_r0, hi_r0);
rb5_br  = bp_z(w5.rw,  P.dx, lo_r0, hi_br);
rb14_r0 = bp_z(w14.rw, P.dx, lo_r0, hi_r0);

wz_res = (z_in >= 2600 & z_in <= 3200);
wx_res = (x_in >= -2500 & x_in <= 1500);
wz_mid = (z_in >= 1600 & z_in <= 2400);

res = struct();
res.smoke = smoke; res.n_freq = nf_run; res.f = f_run; res.m = ms;
res.tstar_anc_dsa_ms = tstar_anc_dsa * 1e3;
res.tstar_anc_flat_ms = tstar_anc_flat * 1e3;
res.beta_anc_raw = beta_anc_raw;
res.beta_anc_dsa = beta_anc_dsa;
res.S_anc_raw = S_anc_raw;
res.S_anc_dsa = S_anc_dsa;
res.spec_anc_eq = spec_anc_eq;
res.spec_res = spec_res;

res.res_slope = struct();
for sk = spec_keys
    sv = spec_res.(sk{1});
    pf = polyfit(f_run(fm), log(max(sv(fm), 1e-30)), 1);
    res.res_slope.(sk{1}) = pf(1);
end
pf_aeq = polyfit(f_run(fm), log(max(spec_anc_eq(fm), 1e-30)), 1);
res.anc_eq_slope = pf_aeq(1);

mkeys = fieldnames(Img)';
fprintf('\n=== VOLVE CLOSURE RESULTS (F-5 R0-band [3.52-12.89] & Broad [3.52-16.0]) ===\n');
for ik = 1:numel(mkeys)
    k = mkeys{ik};
    M = Img.(k); nc = size(M, 2);
    f5_r0 = zeros(1, nc); f5_br = zeros(1, nc); f14_r0 = zeros(1, nc);
    cen = zeros(1, nc); bw10 = zeros(nc, 2); sm_fr = zeros(1, nc);
    for c = 1:nc
        A = prep_volve(M(:, c), P.nzi, P.nxi, z_in);
        f5_r0(c)  = best_cc_well(A, w5,  rb5_r0,  P.dx, lo_r0, hi_r0);
        f5_br(c)  = best_cc_well(A, w5,  rb5_br,  P.dx, lo_r0, hi_br);
        f14_r0(c) = best_cc_well(A, w14, rb14_r0, P.dx, lo_r0, hi_r0);
        [bw10(c,:), cen(c)] = res_band(A(wz_res, wx_res), P.dx);
        sm_fr(c)  = steep_smile_frac(A(wz_mid, wx_res), P.dx, 25.0);
    end
    res.methods.(k) = struct( ...
        'F5_R0band', f5_r0, 'F5_R0band_mean', mean(f5_r0), ...
        'F5_Broad',  f5_br, 'F5_Broad_mean',  mean(f5_br), ...
        'F14_R0band', f14_r0, 'F14_R0band_mean', mean(f14_r0), ...
        'Res_Cen', cen, 'Res_BW10', bw10, 'Smile_Frac', sm_fr);
    if nc == 1
        fprintf('  %-21s | F5_R0=%.4f | F5_Br=%.4f | F14=%.3f | Cen=%.2f | Smile=%.3f\n', ...
            k, f5_r0(1), f5_br(1), f14_r0(1), cen(1), sm_fr(1));
    else
        fprintf('  %-21s | F5_R0=%.4f (m0.5=%.4f,m1=%.4f,m2=%.4f) | F5_Br=%.4f (m1=%.4f) | Cen=%.2f | Smile=%.3f\n', ...
            k, mean(f5_r0), f5_r0(1), f5_r0(3), f5_r0(5), mean(f5_br), f5_br(3), cen(3), sm_fr(3));
    end
end
res.elapsed_s = toc(t_all);

tag = ''; if smoke, tag = '_smoke'; end
fn_json = fullfile(res_dir, ['volve_closure_results' tag '.json']);
fid = fopen(fn_json, 'w'); fwrite(fid, jsonencode(res, 'PrettyPrint', true), 'char'); fclose(fid);
fprintf('\nSaved %s (total %.1fs)\n', fn_json, res.elapsed_s);

if ~smoke
    Img_s = struct();
    for ik = 1:numel(mkeys), Img_s.(mkeys{ik}) = single(Img.(mkeys{ik})); end
    save(fullfile(scratch_base, 'volve_closure_final_imgs.mat'), 'Img_s', 'res', 'z_in', 'x_in', 'phi_eff', '-v7.3');
    fprintf('Saved final 2D images snapshot to scratch.\n');
end
end

%% ========================================================================
%% HELPERS
%% ========================================================================
function [idx, ww] = tab_idx(ts_m, dxg, ntab_m2)
u = ts_m / dxg;
i0 = min(max(floor(u), 0), ntab_m2);
ww = u - i0;
idx = int32(i0) + 1;
end

function [phi_eff, c_st] = blind_structure_dip(R0_2d, dx, sig)
[gx, gz] = gradient(R0_2d, dx);
Jxx = smooth2d(gx.^2, sig); Jzz = smooth2d(gz.^2, sig); Jxz = smooth2d(gx.*gz, sig);
phi_raw = max(min(-0.5 * atan2(2*Jxz, Jzz - Jxx), deg2rad(30)), -deg2rad(30));
c_st = sqrt((Jzz - Jxx).^2 + 4*Jxz.^2) ./ max(Jzz + Jxx, 1e-30);
phi_eff = c_st .* phi_raw;
end

function y = prep_volve(v, nzi, nxi, z_in)
y = reshape(double(v), nzi, nxi);
ker = exp(-0.5 * ((-12:12)/4).^2); ker = ker / sum(ker);
y = y - conv2(ker', 1, [repmat(y(1,:), 12, 1); y; repmat(y(end,:), 12, 1)], 'valid');
y(z_in < 150, :) = 0;
end

function w = load_well(fn, z, x)
fid = fopen(fn); a = fread(fid, Inf, 'single'); fclose(fid);
a = reshape(a, 4, [])'; w.zw = a(:,1); w.sw = a(:,2); w.rw = double(a(:,4));
[~, w.iz] = min(abs(z(:) - w.zw'), [], 1);
w.ix = arrayfun(@(xv) find(abs(x - xv) == min(abs(x - xv)), 1), w.sw);
end

function best = best_cc_well(A, w, rb, dx, lo, hi)
nw = numel(w.zw); tr = zeros(nw, 1);
for i = 1:nw
    tr(i) = mean(A(w.iz(i), max(w.ix(i)-2, 1):min(w.ix(i)+2, size(A, 2))));
end
tr_b = bp_z(tr, dx, lo, hi);
best = -Inf;
for L = -3:3
    a = circshift(tr_b, L); a = a - mean(a); b = rb - mean(rb);
    c = (a' * b) / max(sqrt((a' * a) * (b' * b)), 1e-30);
    if c > best, best = c; end
end
end

function [bw10, cen] = res_band(A_res, dx)
n = 256; kz = (0:n/2) / (n * dx) * 1000;
s = mean(abs(fft(A_res .* hann(size(A_res, 1)), n, 1)), 2);
s = s(1:n/2+1)';
[mx, i] = max(s); t = mx * 10^(-10/20);
l = find(s(1:i) < t, 1, 'last'); if isempty(l), l = 1; end
h = i - 1 + find(s(i:end) < t, 1); if isempty(h), h = numel(s); end
bw10 = [kz(l), kz(h)];
cen = sum(kz .* (s.^2)) / max(sum(s.^2), 1e-30);
end

function fr = steep_smile_frac(A_mid, dx, deg_cut)
[nz, nx] = size(A_mid);
win = hann(nz) * hann(nx)';
F2 = abs(fftshift(fft2(A_mid .* win))).^2;
kz = ((-floor(nz/2)):(ceil(nz/2)-1))' / (nz * dx);
kx = ((-floor(nx/2)):(ceil(nx/2)-1))  / (nx * dx);
[KZ, KX] = ndgrid(kz, kx);
k_mag = sqrt(KX.^2 + KZ.^2);
in_band = (k_mag >= 3.0/1000) & (k_mag <= 25.0/1000);
steep = in_band & (abs(KX) > tand(deg_cut) * max(abs(KZ), 1e-9));
fr = sum(F2(steep), 'all') / max(sum(F2(in_band), 'all'), 1e-30);
end

function y = bp_z(r, dx, lo, hi)
n = 2^nextpow2(4 * numel(r)); R = fft(r, n);
k = abs([0:n/2, -n/2+1:-1]') / (n * dx) * 1000;
y = real(ifft(R .* (k >= lo & k <= hi))); y = y(1:numel(r));
end

function t = cos_taper(x, x0, x1)
u = min(max((x - x0) / max(x1 - x0, 1e-12), 0), 1);
t = 0.5 - 0.5 * cos(pi * u);
end

function B = smooth2d(A, sig)
r = max(1, ceil(3*sig)); g = exp(-0.5*((-r:r)/sig).^2); g = g/sum(g);
Ap = [repmat(A(1,:), r, 1); A; repmat(A(end,:), r, 1)];
B1 = conv2(g', 1, Ap, 'valid');
B1p = [repmat(B1(:,1), 1, r), B1, repmat(B1(:,end), 1, r)];
B = conv2(1, g, B1p, 'valid');
end
