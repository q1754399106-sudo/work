function wsuq_rtm2d_bench(smoke)
% T2.1 Part B: True 2D Full-Wave-Equation Viscoacoustic Helmholtz Multi-Shot RTM Benchmark
% comparing R0, IQ70_old, UQW1_old(25), VA_old, WS_UQ1(25), and WS_VA_UQ1 across m=[0.5,0.7,1,1.4,2]
% and two noise levels (clean 15 dB, noisy 3 dB).
%
% Zero inverse crime:
%   - Forward data: exact 2D viscoacoustic Helmholtz with complex velocity v_c(z,x,f) (4-layer + Q=20 gas chimney)
%   - Reference: exact 2D elastic (Q=inf) Helmholtz data migrated with smooth v_mig
%   - Migration: smooth Gaussian background v_mig + FMM eikonal ray-traced t*_s + blind shot-flip iv estimate
% Peak RAM < 0.6 GB, 2 threads.

if nargin < 1, smoke = false; end
maxNumCompThreads(2);
rng(20261124, 'twister');

here = fileparts(mfilename('fullpath'));
root = fileparts(here);
out_dir = fullfile(root, 'results');
if ~isfolder(out_dir), mkdir(out_dir); end
addpath('D:\ChatGPT\研究\地震衰减补偿研究\uq_rtm_20261117\scripts');
addpath(here);

% ---- 1. Compact 2D Model Geometry (136 x 196 incl. 18 PML => < 0.4 GB RAM, ~6 min full 45 freqs) ----
P.dx = 10.0;
P.npml = 18;
P.nzi = 100; P.nxi = 160;     % Interior: 1000 m x 1600 m
P.nz = P.nzi + 2 * P.npml;
P.nx = P.nxi + 2 * P.npml;
P.z = ((1:P.nz)' - P.npml - 1) * P.dx;
P.x = ((1:P.nx)  - P.npml - 1) * P.dx;
P.in = {P.npml + (1:P.nzi), P.npml + (1:P.nxi)};
P.izs = P.npml + 3;                            % z = 20 m
P.ixs = P.npml + 1 + (8:6:152);                % 25 shots, 60 m spacing (x = 80..1520 m)
P.ixr = P.npml + 1 + (0:2:P.nxi-1);            % 80 receivers, 20 m spacing
ns = numel(P.ixs);
nr = numel(P.ixr);
N = P.nz * P.nx;
[Zi, Xi] = ndgrid(P.in{1}, P.in{2});
lin = sub2ind([P.nz P.nx], Zi(:), Xi(:));
nin = numel(lin);
isrc = sub2ind([P.nz P.nx], P.izs * ones(1, ns), P.ixs);
irec = sub2ind([P.nz P.nx], P.izs * ones(1, nr), P.ixr);

[Zfull, Xfull] = ndgrid(P.z, P.x);
z_in = P.z(P.in{1}); x_in = P.x(P.in{2});

% True velocity v_true (at f_ref = 25 Hz): 5 layers with gentle anticline under gas chimney
z_anc = 240.0;                                 % Shallow anchor interface (240 m)
z_r1  = 640.0 - 25.0 * exp(-((Xfull - 800.0) / 300.0).^2);
z_r2  = 770.0 - 15.0 * exp(-((Xfull - 800.0) / 300.0).^2);
z_r3  = 890.0 * ones(size(Xfull));

v_true = 1850.0 * ones(P.nz, P.nx);
v_true(Zfull >= z_anc) = 2150.0;
v_true(Zfull >= z_r1)  = 2550.0;
v_true(Zfull >= z_r2)  = 2900.0;
v_true(Zfull >= z_r3)  = 3300.0;

% Smooth migration velocity (Gaussian sigma = 25 m = 2.5 grid cells -> ZERO sharp reflectors in v_mig!)
v_mig = smooth2d_gauss(v_true, 2.5);

% True Q model: background Q = 90 below z=80 m + realistic North Sea shallow gas cloud Q = 12 at x=800, z in [240, 620]
invQ_true = (1.0 / 90.0) * ones(P.nz, P.nx);
invQ_true(Zfull < 80.0) = 1.0 / 1000.0;
chim_env = exp(-(((Xfull - 800.0) / 210.0).^2 + ((Zfull - 430.0) / 160.0).^2).^1.5);
invQ_true = invQ_true + (1.0 / 12.0 - 1.0 / 90.0) * chim_env;
Q_true = 1.0 ./ invQ_true;

% Anchored Q model (1/Q = 0 above z_anc = 240 m)
Q_anc = Q_true;
Q_anc(Zfull < z_anc) = Inf;

% ---- 2. FMM Eikonal Two-Way t* Fields (exact same Fermat central-diff recipe as tstar_fields.m) ----
t0_all = tic;
s_mig = 1.0 ./ v_mig;
eps_q = 0.5;
F_abs.T = zeros(P.nz, P.nx, ns, 'single'); F_abs.ts = F_abs.T;
F_anc.T = F_abs.T; F_anc.ts = F_abs.T;
for k = 1:ns
    F_abs.T(:, :, k)  = eikonal_fs(s_mig, P.dx, P.izs, P.ixs(k));
    F_abs.ts(:, :, k) = (eikonal_fs(s_mig .* (1 + eps_q ./ Q_true), P.dx, P.izs, P.ixs(k)) - ...
                         eikonal_fs(s_mig .* (1 - eps_q ./ Q_true), P.dx, P.izs, P.ixs(k))) / (2 * eps_q);
    F_anc.T(:, :, k)  = F_abs.T(:, :, k);
    F_anc.ts(:, :, k) = (eikonal_fs(s_mig .* (1 + eps_q ./ Q_anc),  P.dx, P.izs, P.ixs(k)) - ...
                         eikonal_fs(s_mig .* (1 - eps_q ./ Q_anc),  P.dx, P.izs, P.ixs(k))) / (2 * eps_q);
end
TS_abs = zeros(nin, ns);
TS_anc = zeros(nin, ns);
for k = 1:ns
    [~, tt] = two_way(F_abs, k, P); TS_abs(:, k) = double(tt(lin));
    [~, ta] = two_way(F_anc, k, P); TS_anc(:, k) = max(double(ta(lin)), 0);
end
clear F_abs F_anc

% Quick audit of inter-shot t* contrast under the chimney!
[Zin_g, Xin_g] = ndgrid(z_in, x_in);
m_chim_pix = (Zin_g >= 600 & Zin_g <= 940) & (Xin_g >= 580 & Xin_g <= 1020);
ts_std_chim = median(std(TS_abs(m_chim_pix(:), :), 0, 2)) * 1e3;
ts_rng_chim = median(max(TS_abs(m_chim_pix(:), :), [], 2) - min(TS_abs(m_chim_pix(:), :), [], 2)) * 1e3;
fprintf('2D Gas-Chimney Eikonal ready (%.1fs): W_chim t* std=%.2f ms, range=%.2f ms\n', toc(t0_all), ts_std_chim, ts_rng_chim);

% ---- 3. Frequency Loop: True 2D Forward + Direct-Wave Removal + Multi-Shot RTM ----
if smoke
    freqs = [10.0, 22.0, 34.0];
else
    freqs = 4.0:0.75:37.0;  % 45 frequencies, df = 0.75 Hz -> wraparound 1.33s > T_max 0.95s (ZERO aliasing!)
end
nf = numel(freqs);
fref = 25.0;
fv_band = [15.0, 35.0];
fpeak = 22.0;
sig = 0.447;
s70 = gain_s(70);
ms = [0.5, 0.7, 1.0, 1.4, 2.0];
nm = numel(ms);

% Random shot-flip signs for blind image-domain noise estimation (fixed across freqs like est_noise2)
flip_sgn = 2 * (rand(1, ns) > 0.5) - 1;

% Allocate image accumulators for ref (elastic Q=inf), clean (15 dB), noisy (3 dB), and exact signal/noise splits!
meths = {'R0', 'IQ70_old', 'UQW1_old', 'VA_old', 'WS_UQ0', 'WS_UQ1', 'WS_VA_UQ1'};
I_ref = zeros(nin, 1);
for rg = {'clean', 'noisy'}
    rgn = rg{1};
    for k = 1:numel(meths)
        nc = 1 + (nm - 1) * ~strcmp(meths{k}, 'R0');
        IS.(rgn).(meths{k}) = zeros(nin, nc);
        IN.(rgn).(meths{k}) = zeros(nin, nc);
    end
end

% Precompute reference signal RMS scale at f=22 Hz to calibrate 15 dB & 3 dB receiver noise cleanly
[~, ~, D_cal] = forward_refl(v_true, v_mig, Q_true, 22.0, fref, fpeak, P, isrc, irec, ns, N);
d_rms_peak = sqrt(mean(abs(D_cal(:)).^2));
% Receiver noise amplitudes for 15 dB (clean) and 3 dB (noisy) relative to peak-band reflected data
sig_r_clean = d_rms_peak * 10^(-15.0 / 20.0);
sig_r_noisy = d_rms_peak * 10^(-3.0  / 20.0);
fprintf('Calibrated peak reflection RMS=%.3e (clean noise=%.3e, noisy=%.3e)\n', d_rms_peak, sig_r_clean, sig_r_noisy);

for jf = 1:nf
    t1 = tic;
    f = freqs(jf);
    w = 2 * pi * f;
    W_src = (f / fpeak)^2 * exp(-(f / fpeak)^2);

    % 3a. Smooth migration operator LU (used for direct-wave subtraction AND all RTM backpropagations!)
    A_mig = helm_op(v_mig, f, P);
    [Lm, Um, Pm, Qm, Rm] = lu(A_mig); clear A_mig
    b_src = zeros(N, ns);
    b_src(sub2ind([N ns], isrc, 1:ns)) = -W_src / P.dx^2;
    U_bg = Qm * (Um \ (Lm \ (Pm * (Rm \ b_src))));
    D_bg = U_bg(irec, :);
    Us_conj = conj(U_bg(lin, :)) * (w^2);
    clear U_bg

    % 3b. True elastic (Q=inf) reflection data -> Elastic Reference Image I_ref
    A_el = helm_op(v_true, f, P);
    [Le, Ue, Pe, Qe, Re] = lu(A_el); clear A_el
    U_el = Qe * (Ue \ (Le \ (Pe * (Re \ b_src))));
    D_el = U_el(irec, :) - D_bg;
    clear U_el Le Ue Pe Qe Re

    % 3c. True viscoacoustic (Azimi complex velocity) reflection data
    vc = v_true ./ (1.0 - (1.0 ./ (pi * Q_true)) * log(max(f, 1e-6) / fref) - 1i ./ (2.0 * Q_true));
    A_vis = helm_op(vc, f, P);
    [Lv, Uv, Pv, Qv, Rv] = lu(A_vis); clear A_vis vc
    U_vis = Qv * (Uv \ (Lv \ (Pv * (Rv \ b_src))));
    D_vis = U_vis(irec, :) - D_bg;
    clear U_vis Lv Uv Pv Qv Rv b_src

    % 3d. Backpropagate elastic ref, visco signal, and unit white receiver noise using ALREADY-FACTORED Lm,Um!
    N_unit = (1.0 / sqrt(2.0)) * (randn(nr, ns) + 1i * randn(nr, ns));

    backprop = @(Dmat) extract_backprop(Dmat, Lm, Um, Pm, Qm, Rm, irec, lin, N, ns, P.dx);
    V_el  = backprop(D_el);
    V_vis = backprop(D_vis);
    V_nu  = backprop(N_unit);
    clear Lm Um Pm Qm Rm

    Cs_ref = Us_conj .* V_el;
    I_ref = I_ref + real(sum(Cs_ref, 2));
    clear Cs_ref V_el

    Cs_sig = Us_conj .* V_vis;
    Cs_nu  = Us_conj .* V_nu;
    % Exact source-illumination weights a_s = |U_s|^2 / sum |U_k|^2 for exact MRC in 2D RTM!
    shot_ill = abs(Us_conj).^2;
    clear Us_conj V_vis V_nu

    % 3e. Apply all operators to BOTH clean (15 dB) and noisy (3 dB) regimes
    for rg = {'clean', 'noisy'}
        rgn = rg{1};
        s_nr = sig_r_clean;
        if strcmp(rgn, 'noisy'), s_nr = sig_r_noisy; end
        Cs_noi = s_nr * Cs_nu;
        Cs_tot = Cs_sig + Cs_noi;

        % Blind image-domain 1/SNR estimate via random-polarity shot-flip (same principle as est_noise2!)
        stk_tot = abs(sum(Cs_tot, 2)).^2;
        stk_flp = abs(sum(Cs_tot .* flip_sgn, 2)).^2;
        p_tot = smooth2d_gauss(reshape(stk_tot, P.nzi, P.nxi), 5.0);
        p_flp = smooth2d_gauss(reshape(stk_flp, P.nzi, P.nxi), 5.0);
        iv_est = min(max(p_flp(:) ./ max(p_tot(:) - p_flp(:), 0.05 * p_flp(:)), 1e-3), 1e3);

        % R0
        IS.(rgn).R0 = IS.(rgn).R0 + real(sum(Cs_sig, 2));
        IN.(rgn).R0 = IN.(rgn).R0 + real(sum(Cs_noi, 2));

        for im = 1:nm
            m = ms(im);
            TS_m = TS_abs / m;
            TS_am = TS_anc / m;

            [~, uq_s, ws0_s, ws1_s, ~, uq_n, ws0_n, ws1_n] = ws_uq_apply(Cs_sig, TS_m, f, sig, fref, iv_est, s70, [], shot_ill, Cs_noi);

            [~, iq_s, ~, ~, ~, iq_n, ~, ~] = ws_uq_apply(Cs_sig, TS_m, f, 0.0, fref, zeros(nin, 1), s70, [], shot_ill, Cs_noi);

            [~, va_s, ~, wva_s, ~, va_n, ~, wva_n] = ws_uq_apply(Cs_sig, TS_am, f, sig, fref, iv_est, s70, fv_band, shot_ill, Cs_noi);

            IS.(rgn).IQ70_old(:, im)  = IS.(rgn).IQ70_old(:, im)  + iq_s;  IN.(rgn).IQ70_old(:, im)  = IN.(rgn).IQ70_old(:, im)  + iq_n;
            IS.(rgn).UQW1_old(:, im)  = IS.(rgn).UQW1_old(:, im)  + uq_s;  IN.(rgn).UQW1_old(:, im)  = IN.(rgn).UQW1_old(:, im)  + uq_n;
            IS.(rgn).VA_old(:, im)    = IS.(rgn).VA_old(:, im)    + va_s;  IN.(rgn).VA_old(:, im)    = IN.(rgn).VA_old(:, im)    + va_n;
            IS.(rgn).WS_UQ0(:, im)    = IS.(rgn).WS_UQ0(:, im)    + ws0_s; IN.(rgn).WS_UQ0(:, im)    = IN.(rgn).WS_UQ0(:, im)    + ws0_n;
            IS.(rgn).WS_UQ1(:, im)    = IS.(rgn).WS_UQ1(:, im)    + ws1_s; IN.(rgn).WS_UQ1(:, im)    = IN.(rgn).WS_UQ1(:, im)    + ws1_n;
            IS.(rgn).WS_VA_UQ1(:, im) = IS.(rgn).WS_VA_UQ1(:, im) + wva_s; IN.(rgn).WS_VA_UQ1(:, im) = IN.(rgn).WS_VA_UQ1(:, im) + wva_n;
        end
    end

    if mod(jf, 10) == 0 || jf == 1 || jf == nf
        fprintf('rtm2d f=%4.1f Hz (%2d/%d): %.1fs (elapsed %.0fs)\n', f, jf, nf, toc(t1), toc(t0_all));
    end
end

% ---- 4. Non-Saturating 4-Metric Evaluation (E1, E1_HF, E_amp_chim, SNR_img_dB) ----
ref_2d = prep_img(I_ref, P.nzi, P.nxi, z_in);
ref_hf = highpass_kz(ref_2d, P.dx, 10.0);   % kz >= 10 cycles/km (~26+ Hz!)

W_chim  = {[580 1020], [580 950]};
W_flank = {[180 500],  [580 950]};
W_deep  = {[180 1420], [580 950]};

rms_win = @(A, W) sqrt(mean(A(z_in>=W{2}(1)&z_in<=W{2}(2), x_in>=W{1}(1)&x_in<=W{1}(2)).^2, 'all'));
ref_ratio_cf = log10(rms_win(ref_2d, W_chim) / max(rms_win(ref_2d, W_flank), 1e-30));

res = struct('smoke', smoke, 'n_freq', nf, 'm', ms, 'elapsed_s', toc(t0_all), ...
             'ts_std_chim_ms', ts_std_chim, 'ts_rng_chim_ms', ts_rng_chim);

for rg = {'clean', 'noisy'}
    rgn = rg{1};
    fprintf('\n=== 2D Full-Wave Gas-Chimney RTM [%s] (W_chim) ===\n', rgn);
    for wn = {'W_chim', 'W_deep'}
        wname = wn{1};
        W = W_chim; if strcmp(wname, 'W_deep'), W = W_deep; end
        iz_w = find(z_in >= W{2}(1) & z_in <= W{2}(2));
        ix_w = find(x_in >= W{1}(1) & x_in <= W{1}(2));
        for k = 1:numel(meths)
            mk = meths{k};
            Sm = IS.(rgn).(mk); Nm = IN.(rgn).(mk);
            nc = size(Sm, 2);
            e1 = zeros(1, nc); e1_hf = zeros(1, nc); eamp = zeros(1, nc); snr_db = zeros(1, nc);
            for c = 1:nc
                S2 = prep_img(Sm(:, c), P.nzi, P.nxi, z_in);
                N2 = prep_img(Nm(:, c), P.nzi, P.nxi, z_in);
                T2 = S2 + N2;
                Th = highpass_kz(T2, P.dx, 10.0);
                e1(c) = score_e1(T2, ref_2d, iz_w, ix_w);
                e1_hf(c) = score_e1(Th, ref_hf, iz_w, ix_w);
                snr_db(c) = 10 * log10(sum(S2(iz_w, ix_w).^2, 'all') / max(sum(N2(iz_w, ix_w).^2, 'all'), 1e-30));
                rat_cf = log10(rms_win(T2, W_chim) / max(rms_win(T2, W_flank), 1e-30));
                eamp(c) = abs(rat_cf - ref_ratio_cf);
            end
            res.(rgn).(wname).(mk) = struct( ...
                'E1', e1, 'E1_mean', mean(e1), ...
                'E1_HF', e1_hf, 'E1_HF_mean', mean(e1_hf), ...
                'E_amp', eamp, 'E_amp_mean', mean(eamp), ...
                'SNR_img_dB', snr_db, 'SNR_img_dB_mean', mean(snr_db));
            if strcmp(wname, 'W_chim')
                fprintf('  %-10s | E1=%.4f (m0.5=%.4f, m1=%.4f) | E1_HF=%.4f (m0.5=%.4f) | E_amp=%.3f | SNR=%.2f dB (m0.5=%.2f)\n', ...
                        mk, mean(e1), e1(1), e1(min(3, nc)), mean(e1_hf), e1_hf(1), mean(eamp), mean(snr_db), snr_db(1));
            end
        end
    end
end

tag = ''; if smoke, tag = '_smoke'; end
fn = fullfile(out_dir, ['wsuq_rtm2d_bench' tag '.json']);
fid = fopen(fn, 'w'); fwrite(fid, jsonencode(res, 'PrettyPrint', true), 'char'); fclose(fid);
fprintf('\nSaved %s (total %.1fs)\n', fn, toc(t0_all));
end

function [D_el, D_bg, D_vis] = forward_refl(v_true, v_mig, Q_true, f, fref, fpeak, P, isrc, irec, ns, N)
W_src = (f / fpeak)^2 * exp(-(f / fpeak)^2);
b_src = zeros(N, ns); b_src(sub2ind([N ns], isrc, 1:ns)) = -W_src / P.dx^2;
Am = helm_op(v_mig, f, P); [L, U, Pp, Qq, R] = lu(Am); clear Am
Ub = Qq * (U \ (L \ (Pp * (R \ b_src)))); D_bg = Ub(irec, :); clear Ub L U Pp Qq R
vc = v_true ./ (1.0 - (1.0 ./ (pi * Q_true)) * log(max(f, 1e-6) / fref) - 1i ./ (2.0 * Q_true));
Av = helm_op(vc, f, P); [L, U, Pp, Qq, R] = lu(Av); clear Av
Uv = Qq * (U \ (L \ (Pp * (R \ b_src)))); D_vis = Uv(irec, :) - D_bg; clear Uv L U Pp Qq R
D_el = D_vis;
end

function V = extract_backprop(Dmat, L, U, Pp, Qq, R, irec, lin, N, ns, dx)
b = zeros(N, ns);
b(irec, :) = Dmat;
Vf = -(R \ (Pp' * (L' \ (U' \ (Qq' * b))))) / dx^2;
V = Vf(lin, :);
end

function B = smooth2d_gauss(A, sig)
r = max(1, ceil(3 * sig));
g = exp(-0.5 * ((-r:r) / sig).^2); g = g / sum(g);
Ap = [repmat(A(1, :), r, 1); A; repmat(A(end, :), r, 1)];
B1 = conv2(g', 1, Ap, 'valid');
B1p = [repmat(B1(:, 1), 1, r), B1, repmat(B1(:, end), 1, r)];
B = conv2(1, g, B1p, 'valid');
end

function y = prep_img(v, nzi, nxi, z_in)
y = reshape(double(v), nzi, nxi);
y = y - smooth2d_gauss(y, 5.0);
y(z_in < 160, :) = 0;
end

function yh = highpass_kz(y, dx, kz_cut)
n = 2^nextpow2(2 * size(y, 1));
Y = fft(y, n, 1);
kz = abs([0:n/2, -n/2+1:-1]') / (n * dx) * 1000.0;
g = min(max((kz - (kz_cut - 2)) / 4.0, 0), 1);
g = 0.5 - 0.5 * cos(pi * g);
yh = real(ifft(Y .* g, n, 1));
yh = yh(1:size(y, 1), :);
end

function s = score_e1(a, r, iz, ix)
rr = r(iz, ix);
c0 = cc(a(iz, ix), rr);
cl = c0;
nz = size(a, 1);
for L = [-2 -1 1 2]
    iz2 = min(max(iz + L, 1), nz);
    cl = max(cl, cc(a(iz2, ix), rr));
end
s = 0.5 * (c0 + cl);
end

function c = cc(a, b)
a = a(:) - mean(a(:)); b = b(:) - mean(b(:));
c = (a' * b) / max(sqrt((a' * a) * (b' * b)), 1e-30);
end
