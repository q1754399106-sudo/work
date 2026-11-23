function dsa_confirm_bench(models_req, smoke)
% Pre-registered Confirmation Benchmark for DSA-UQW1 (uq_rtm_dsa_confirm_20261123/PLAN_ZH.md)
% Part 1: 4 brand-new unseen 2D full-wave confirmation models (H1-H4, seeds 20261201-20261204)
% Part 2: Volve-exact flat-W (2*pi*f)^2 real-data well-tie replica H5_VolveAnc (seed 20261205)
% Peak RAM < 0.6 GB, 4 threads, resumable per model (results/<name>.json).

if nargin < 1 || isempty(models_req), models_req = {'H1', 'H2', 'H3', 'H4', 'H5_VolveAnc'}; end
if ischar(models_req) || isstring(models_req), models_req = {char(models_req)}; end
if nargin < 2, smoke = false; end

maxNumCompThreads(4);
here = fileparts(mfilename('fullpath'));
root = fileparts(here);
out_dir = fullfile(root, 'results');
if ~isfolder(out_dir), mkdir(out_dir); end
addpath('D:\ChatGPT\研究\地震衰减补偿研究\uq_rtm_20261117\scripts');

for imod = 1:numel(models_req)
    mname = models_req{imod};
    tag = ''; if smoke, tag = '_smoke'; end
    fn = fullfile(out_dir, [mname tag '.json']);
    if isfile(fn) && ~smoke
        fprintf('Skipping %s (already exists: %s)\n', mname, fn);
        continue;
    end
    if strcmp(mname, 'H5_VolveAnc')
        res = run_h5_volve_anc(smoke);
    else
        res = run_confirm_model(mname, smoke);
    end
    fid = fopen(fn, 'w'); fwrite(fid, jsonencode(res, 'PrettyPrint', true), 'char'); fclose(fid);
    fprintf('Saved %s (%.1fs)\n\n', fn, res.elapsed_s);
end
end

function res = run_confirm_model(mname, smoke)
t0_all = tic;
P.dx = 10.0; P.npml = 18; P.nzi = 100; P.nxi = 160;
P.nz = P.nzi + 2*P.npml; P.nx = P.nxi + 2*P.npml;
P.z = ((1:P.nz)' - P.npml - 1) * P.dx;
P.x = ((1:P.nx)  - P.npml - 1) * P.dx;
P.in = {P.npml + (1:P.nzi), P.npml + (1:P.nxi)};
P.izs = P.npml + 3;
P.ixs = P.npml + 1 + (8:6:152);      % 25 shots, 60 m spacing (x = 80..1520 m)
P.ixr = P.npml + 1 + (0:2:P.nxi-1);  % 80 receivers, 20 m spacing (x = 0..1580 m)
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
corr_noise_sig = 0.0;

switch mname
    case 'H1'
        % Brand-new Syncline + Twin off-center low-Q lenses (Q=14 & Q=16)
        seed = 20261201;
        z1 = 230.0 + 15.0 * exp(-((Xfull - 780.0)/350.0).^2);
        z2 = 620.0 + 40.0 * exp(-((Xfull - 820.0)/320.0).^2);  % Syncline!
        z3 = 755.0 + 30.0 * exp(-((Xfull - 820.0)/340.0).^2);
        z4 = 885.0 + 15.0 * exp(-((Xfull - 820.0)/380.0).^2);
        v_true = 1840.0 * ones(P.nz, P.nx);
        v_true(Zfull >= z1) = 2140.0;
        v_true(Zfull >= z2) = 2540.0;
        v_true(Zfull >= z3) = 2900.0;
        v_true(Zfull >= z4) = 3300.0;
        invQ = (1/90)*ones(P.nz, P.nx); invQ(Zfull < 80) = 1/1000;
        L1 = exp(-(((Xfull - 560)/180).^2 + ((Zfull - 410)/105).^2).^1.5);
        L2 = exp(-(((Xfull - 1060)/170).^2 + ((Zfull - 460)/110).^2).^1.5);
        invQ = invQ + (1/14 - 1/90)*L1 + (1/16 - 1/90)*L2;
        Q_true = 1./invQ;
        W_targ = {[420 720], [560 950]}; W_flank = {[740 900], [560 950]};

    case 'H2'
        % Reverse fault + Steeper opposing dips (+22 deg / -16 deg) + shallow low-Q wedge
        seed = 20261202;
        Zg = Zfull;
        % Reverse fault: hanging wall shifted UP by 30 m
        hw = (Xfull < 840.0 - (Zfull - 200.0) / tand(64.0)) & (Zfull > 200.0);
        Zg(hw) = Zg(hw) + 30.0;
        z1 = 235.0 - 0.06 * (Xfull - 800.0);
        z2 = 660.0 + 0.40 * (Xfull - 800.0);   % +21.8 deg steep dip!
        z3 = 790.0 - 0.29 * (Xfull - 800.0);   % -16.2 deg opposite steep dip!
        z4 = 900.0 - 0.12 * (Xfull - 800.0);
        v_true = 1860.0 * ones(P.nz, P.nx);
        v_true(Zg >= z1) = 2160.0;
        v_true(Zg >= z2) = 2620.0;
        v_true(Zg >= z3) = 2980.0;
        v_true(Zg >= z4) = 3350.0;
        invQ = (1/85)*ones(P.nz, P.nx); invQ(Zfull < 80) = 1/1000;
        wedge = exp(-(((Xfull - 980)/260).^2 + ((Zfull - 360)/95).^2).^1.5);
        invQ = invQ + (1/16 - 1/85)*wedge;
        Q_true = 1./invQ;
        W_targ = {[780 1200], [560 950]}; W_flank = {[220 560], [560 950]};

    case 'H3'
        % Volve-like One-Sided Trailing Marine Streamer (reverse sail: [xs-860, xs+240]) + Chimney + 12-deg dip
        seed = 20261203;
        for s = 1:ns
            rec_mask(:, s) = (xr >= xs(s) - 860.0) & (xr <= xs(s) + 240.0);
        end
        z1 = 240.0 - 0.05 * (Xfull - 800.0);
        z2 = 640.0 - 0.21 * (Xfull - 800.0);   % -11.9 deg dipping reservoir!
        z3 = 770.0 - 0.18 * (Xfull - 800.0);   % -10.2 deg dipping base!
        z4 = 890.0 + 0.04 * (Xfull - 800.0);
        v_true = 1830.0 * ones(P.nz, P.nx);
        v_true(Zfull >= z1) = 2130.0;
        v_true(Zfull >= z2) = 2550.0;
        v_true(Zfull >= z3) = 2900.0;
        v_true(Zfull >= z4) = 3260.0;
        invQ = (1/90)*ones(P.nz, P.nx); invQ(Zfull < 80) = 1/1000;
        chim = exp(-(((Xfull - 760)/195).^2 + ((Zfull - 420)/150).^2).^1.5);
        invQ = invQ + (1/15 - 1/90)*chim;
        Q_true = 1./invQ;
        W_targ = {[560 960], [560 950]}; W_flank = {[1040 1380], [560 950]};

    case 'H4'
        % Severe laterally correlated noise (sigma=3 tr = 60 m) + -3% velocity error + lateral Q gradient + folds
        seed = 20261204;
        vscale_mig = 0.97;      % Opposite sign (-3%) vs G4 (+3%)!
        corr_noise_sig = 3.0;   % 60 m lateral correlation length
        z1 = 235.0 + 12.0 * sin(2*pi*Xfull/900.0);
        z2 = 635.0 + 22.0 * sin(2*pi*(Xfull - 200)/850.0);
        z3 = 765.0 + 18.0 * cos(2*pi*(Xfull - 150)/950.0);
        z4 = 885.0 - 12.0 * sin(2*pi*Xfull/1100.0);
        v_true = 1850.0 * ones(P.nz, P.nx);
        v_true(Zfull >= z1) = 2150.0;
        v_true(Zfull >= z2) = 2560.0;
        v_true(Zfull >= z3) = 2920.0;
        v_true(Zfull >= z4) = 3300.0;
        % Strong lateral Q gradient: 110 on left -> 28 on right
        Q_true = 110.0 - 82.0 * min(max(Xfull / 1600.0, 0), 1);
        Q_true(Zfull < 80) = 1000.0;
        W_targ = {[980 1400], [560 950]}; W_flank = {[200 600], [560 950]};

    otherwise
        error('Unknown model %s', mname);
end

rng(seed, 'twister');
v_mig_exact = smooth2d(v_true, 2.5);
v_mig = vscale_mig * v_mig_exact;

% Eikonal one-way t* fields F.ts (kept for BOTH flat two_way and dip-aware two_way!)
s_mig = 1.0 ./ v_mig; eps_q = 0.5;
F.T = zeros(P.nz, P.nx, ns, 'single'); F.ts = F.T;
for k = 1:ns
    F.T(:,:,k) = eikonal_fs(s_mig, P.dx, P.izs, P.ixs(k));
    F.ts(:,:,k) = (eikonal_fs(s_mig.*(1+eps_q./Q_true), P.dx, P.izs, P.ixs(k)) - ...
                   eikonal_fs(s_mig.*(1-eps_q./Q_true), P.dx, P.izs, P.ixs(k))) / (2*eps_q);
end
TS_flat = zeros(nin, ns);
Fts_in  = zeros(nin, ns);
for k = 1:ns
    [~, tt] = two_way(F, k, P); TS_flat(:, k) = double(tt(lin));
    tk = F.ts(:,:,k); Fts_in(:, k) = double(tk(lin));
end
clear F

% Live receiver bounds per shot
xr_lo = zeros(1, ns); xr_hi = zeros(1, ns);
for k = 1:ns
    live_xr = xr(rec_mask(:, k));
    xr_lo(k) = min(live_xr); xr_hi(k) = max(live_xr);
end

% Flat-mirror SA mask (ablation baseline phi = 0)
w_sa_flat = build_dsa_mask(Xin(:), Zin(:), xs, xr_lo, xr_hi, zeros(nin, 1));
w_deep_mask = (Zin(:) >= 560 & Zin(:) <= 950) & (Xin(:) >= 180 & Xin(:) <= 1420);
out_pairs = 0; tot_pairs = 0;
for k = 1:ns
    xm = 2*Xin(w_deep_mask) - xs(k);
    out_pairs = out_pairs + sum(xm < xr_lo(k) | xm > xr_hi(k));
    tot_pairs = tot_pairs + numel(xm);
end
out_frac = out_pairs / max(tot_pairs, 1);
fprintf('=== Confirm Model %s (seed=%d, out-of-array=%.1f%%, vscale=%.2f, corr_sig=%.1f tr) ===\n', ...
    mname, seed, 100*out_frac, vscale_mig, corr_noise_sig);

if smoke, freqs = [10.0, 22.0, 34.0]; else, freqs = 4.0:1.0:37.0; end
nf = numel(freqs);
fref = 25.0; fpeak = 22.0; sig = 0.447; s70 = gain_s(70);
ms = [0.5, 0.7, 1.0, 1.4, 2.0]; nm = numel(ms);
flip_sgn = 2 * (rand(1, ns) > 0.5) - 1;

% Calibrate 15 dB and 3 dB noise on live reflections at f=22 Hz
W22 = (22/fpeak)^2 * exp(-(22/fpeak)^2);
b22 = zeros(N, ns); b22(sub2ind([N ns], isrc, 1:ns)) = -W22 / P.dx^2;
[Lm,Um,Pm,Qm,Rm] = lu(helm_op(v_mig_exact, 22, P)); Ub = Qm*(Um\(Lm\(Pm*(Rm\b22)))); Dbg22 = Ub(irec,:); clear Ub Lm Um Pm Qm Rm
vc22 = v_true ./ (1 - (1./(pi*Q_true))*log(22/fref) - 1i./(2*Q_true));
[Lv,Uv,Pv,Qv,Rv] = lu(helm_op(vc22, 22, P)); Uv22 = Qv*(Uv\(Lv\(Pv*(Rv\b22)))); Dcal = (Uv22(irec,:) - Dbg22) .* rec_mask; clear Uv22 Lv Uv Pv Qv Rv
d_rms = sqrt(mean(abs(Dcal(rec_mask)).^2));
sig_clean = d_rms * 10^(-15/20);
sig_noisy = d_rms * 10^(-3/20);

if corr_noise_sig > 0
    rc = ceil(3 * corr_noise_sig);
    g_nc = exp(-0.5 * ((-rc:rc)' / corr_noise_sig).^2);
    g_nc = g_nc / sqrt(sum(g_nc.^2));
end

% Pass 1: Compute & store single-freq Cs_sig, Cs_nu (only ~108 MB in single!) + I_ref + raw R0 for blind dip tensor!
Cs_sig_all = zeros(nin, ns, nf, 'single');
Cs_nu_all  = zeros(nin, ns, nf, 'single');
I_ref = zeros(nin, 1);
I_R0_sig = zeros(nin, 1);
I_R0_nu  = zeros(nin, 1);

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
    if mod(jf, 10) == 0 || jf == nf
        fprintf('  %s Pass-1 f=%4.1f Hz (%2d/%d), elapsed %.1fs\n', mname, f, jf, nf, toc(t0_all));
    end
end

% Pass 2: Evaluate RAW, SA_flat, and DSA (coherence-weighted dip estimated BLINDLY from each noise regime's R0!)
modes = {'RAW', 'SA_flat', 'DSA'};
ops   = {'R0', 'IQ70', 'W1', 'UQW1'};
for rg = {'clean', 'noisy'}
    for imd = 1:3
        for io = 1:4
            key = [modes{imd} '_' ops{io}];
            nc = 1 + (nm-1)*~strcmp(ops{io}, 'R0');
            IS.(rg{1}).(key) = zeros(nin, nc);
            IN.(rg{1}).(key) = zeros(nin, nc);
        end
    end
end

dip_stats = struct();
for rg = {'clean', 'noisy'}
    rgn = rg{1};
    s_nr = sig_clean; if strcmp(rgn, 'noisy'), s_nr = sig_noisy; end
    R0_tot_2d = prep_img(I_R0_sig + s_nr * I_R0_nu, P.nzi, P.nxi, z_in);
    [phi_eff, c_st] = blind_structure_dip(R0_tot_2d, P.dx);
    dip_stats.(rgn) = struct('med_abs_phi_deg', median(abs(rad2deg(phi_eff(w_deep_mask)))), ...
                             'p90_abs_phi_deg', prctile(abs(rad2deg(phi_eff(w_deep_mask))), 90), ...
                             'med_coh', median(c_st(w_deep_mask)));

    [w_dsa, TS_dsa] = build_dsa_mask_and_tstar(Xin(:), Zin(:), xs, xr_lo, xr_hi, phi_eff(:), Fts_in);
    W_list  = {ones(nin, ns), w_sa_flat, w_dsa};
    TS_list = {TS_flat,       TS_flat,   TS_dsa};
    xmax = 1.05 * max(max(TS_flat(:)), max(TS_dsa(:))) / min(ms);
    xg = linspace(0, xmax, 1201); dxg = xg(2);

    for jf = 1:nf
        f = freqs(jf);
        Cs_s = double(Cs_sig_all(:,:,jf));
        Cs_n = s_nr * double(Cs_nu_all(:,:,jf));
        Cs_t = Cs_s + Cs_n;
        T0 = uq_tab(f, xg, 0, fref);
        T1 = uq_tab(f, xg, sig, fref);

        for imd = 1:3
            Wm = W_list{imd}; TSm = TS_list{imd}; mstr = modes{imd};
            IS.(rgn).([mstr '_R0']) = IS.(rgn).([mstr '_R0']) + real(sum(Wm .* Cs_s, 2));
            IN.(rgn).([mstr '_R0']) = IN.(rgn).([mstr '_R0']) + real(sum(Wm .* Cs_n, 2));

            % Mode-consistent blind shot-flip iv
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
                IS.(rgn).([mstr '_IQ70'])(:,im) = IS.(rgn).([mstr '_IQ70'])(:,im) + real(sum((Wm.*H_iq).*Cs_s, 2));
                IN.(rgn).([mstr '_IQ70'])(:,im) = IN.(rgn).([mstr '_IQ70'])(:,im) + real(sum((Wm.*H_iq).*Cs_n, 2));
                IS.(rgn).([mstr '_W1'])(:,im)   = IS.(rgn).([mstr '_W1'])(:,im)   + real(sum((Wm.*H_w1).*Cs_s, 2));
                IN.(rgn).([mstr '_W1'])(:,im)   = IN.(rgn).([mstr '_W1'])(:,im)   + real(sum((Wm.*H_w1).*Cs_n, 2));
                IS.(rgn).([mstr '_UQW1'])(:,im) = IS.(rgn).([mstr '_UQW1'])(:,im) + real(sum((Wm.*H_uq).*Cs_s, 2));
                IN.(rgn).([mstr '_UQW1'])(:,im) = IN.(rgn).([mstr '_UQW1'])(:,im) + real(sum((Wm.*H_uq).*Cs_n, 2));
            end
        end
    end
end

% Evaluate all 12 combinations on W_deep and W_targ
ref_2d = prep_img(I_ref, P.nzi, P.nxi, z_in);
ref_hf = highpass_kz(ref_2d, P.dx, 10.0);
W_deep = {[180 1420], [560 950]};
rms_w = @(A, W) sqrt(mean(A(z_in>=W{2}(1)&z_in<=W{2}(2), x_in>=W{1}(1)&x_in<=W{1}(2)).^2, 'all'));
ref_rcf = log10(rms_w(ref_2d, W_targ) / max(rms_w(ref_2d, W_flank), 1e-30));

res = struct('model', mname, 'seed', seed, 'smoke', smoke, 'n_freq', nf, 'm', ms, ...
             'out_of_array_frac', out_frac, 'vscale_mig', vscale_mig, ...
             'corr_noise_sig_tr', corr_noise_sig, 'dip_stats', dip_stats, 'elapsed_s', toc(t0_all));

for rg = {'clean', 'noisy'}
    rgn = rg{1};
    for wn = {'W_deep', 'W_targ'}
        wname = wn{1}; W = W_deep; if strcmp(wname, 'W_targ'), W = W_targ; end
        iz_w = find(z_in >= W{2}(1) & z_in <= W{2}(2));
        ix_w = find(x_in >= W{1}(1) & x_in <= W{1}(2));
        if strcmp(wname, 'W_deep')
            fprintf('\n  [%s %s W_deep] (blind dip med=%.1f deg, p90=%.1f deg)\n', ...
                mname, rgn, dip_stats.(rgn).med_abs_phi_deg, dip_stats.(rgn).p90_abs_phi_deg);
        end
        for imd = 1:3
            for io = 1:4
                key = [modes{imd} '_' ops{io}];
                Sm = IS.(rgn).(key); Nm = IN.(rgn).(key); nc = size(Sm, 2);
                e1 = zeros(1,nc); e1h = zeros(1,nc); eamp = zeros(1,nc); snr = zeros(1,nc);
                for c = 1:nc
                    S2 = prep_img(Sm(:,c), P.nzi, P.nxi, z_in);
                    N2 = prep_img(Nm(:,c), P.nzi, P.nxi, z_in);
                    T2 = S2 + N2; Th = highpass_kz(T2, P.dx, 10.0);
                    e1(c)  = score_e1(T2, ref_2d, iz_w, ix_w);
                    e1h(c) = score_e1(Th, ref_hf, iz_w, ix_w);
                    snr(c) = 10*log10(sum(S2(iz_w,ix_w).^2,'all') / max(sum(N2(iz_w,ix_w).^2,'all'), 1e-30));
                    eamp(c)= abs(log10(rms_w(T2, W_targ)/max(rms_w(T2, W_flank),1e-30)) - ref_rcf);
                end
                res.(rgn).(wname).(key) = struct( ...
                    'E1', e1, 'E1_mean', mean(e1), ...
                    'E1_HF', e1h, 'E1_HF_mean', mean(e1h), ...
                    'E_amp', eamp, 'E_amp_mean', mean(eamp), ...
                    'SNR_img_dB', snr, 'SNR_img_dB_mean', mean(snr));
                if strcmp(wname, 'W_deep')
                    fprintf('    %-12s | E1=%.4f (m0.5=%.4f,m1=%.4f,m2=%.4f) | HF=%.4f | Eamp=%.3f | SNR=%.2f dB\n', ...
                        key, mean(e1), e1(1), e1(min(3,nc)), e1(end), mean(e1h), mean(eamp), mean(snr));
                end
            end
        end
    end
end
end

function res = run_h5_volve_anc(smoke)
% Part 2: Volve-exact real-data gap & Shallow-Anchored DSA-UQW1-Anc verification!
% Forward data uses realistic Ricker + Q=60 overburden + Q=35 reservoir + one-sided trailing marine streamer,
% while Migration does NOT know the source wavelet -> uses flat boxcar W_mig = 1 & (2*pi*f)^2 exactly like fld_run.m!
% Evaluates 5-well flat-bandpass vertical well-tie (exact fld_eval.m metric!) against true vertical reflectivity rw(z).
t0_all = tic;
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

% One-sided trailing marine streamer (asymmetric like Volve!)
rec_mask = false(nr, ns);
xr_lo = zeros(1, ns); xr_hi = zeros(1, ns);
for s = 1:ns
    rec_mask(:, s) = (xr >= xs(s) - 240.0) & (xr <= xs(s) + 860.0);
    xr_lo(s) = min(xr(rec_mask(:,s))); xr_hi(s) = max(xr(rec_mask(:,s)));
end

% Realistic finely layered earth (like Volve well logs: shallow reflectors + reservoir stack in [560, 920] m)
v_bg = 1800.0 + 1.45 * max(Zfull, 0);
% Add 25 random thin-layer velocity impedance steps with gentle 6-deg dip!
z_dip = Zfull + 0.10 * (Xfull - 800.0);
dv_layer = zeros(size(Zfull));
z_edges = [180, 230, 290, 350, 420, 490, 560, 595, 635, 670, 710, 750, 790, 830, 865, 905, 945];
for il = 1:numel(z_edges)-1
    step = 0.05 * randn();
    step = max(min(step, 0.09), -0.09);
    dv_layer(z_dip >= z_edges(il) & z_dip < z_edges(il+1)) = step;
end
v_true = v_bg .* (1 + dv_layer);
v_mig  = smooth2d(v_true, 3.0);

% True Q model: Q = 55 overburden (below z=80m) + Q = 35 reservoir below z=520m (strong Volve-like attenuation!)
Q_true = 55.0 * ones(P.nz, P.nx);
Q_true(Zfull < 80) = 1000.0;
Q_true(Zfull >= 520) = 35.0;

% Compute eikonal F.ts and flat TS
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

if smoke, freqs = [10.0, 22.0, 34.0]; else, freqs = 4.0:1.0:37.0; end
nf = numel(freqs);
fref = 25.0; fpeak_true = 20.0; sig = 0.447; s70 = gain_s(70);
ms = [0.5, 0.7, 1.0, 1.4, 2.0]; nm = numel(ms);
flip_sgn = 2 * (rand(1, ns) > 0.5) - 1;

% Volve-like flat boxcar migration source wavelet W_mig(f) (cosine tails 4-6 & 34-37 Hz)
W_mig = ones(1, nf);
a1 = freqs < 6; W_mig(a1) = 0.5 * (1 - cos(pi*(freqs(a1)-4)/2));
a2 = freqs > 34; W_mig(a2) = 0.5 * (1 + cos(pi*(freqs(a2)-34)/3));

% Calibrate 6 dB field-realistic noise at f=20 Hz
W20 = (20/fpeak_true)^2 * exp(-(20/fpeak_true)^2);
b20 = zeros(N, ns); b20(sub2ind([N ns], isrc, 1:ns)) = -W20 / P.dx^2;
[Lm,Um,Pm,Qm,Rm] = lu(helm_op(v_mig, 20, P)); Ub = Qm*(Um\(Lm\(Pm*(Rm\b20)))); Dbg20 = Ub(irec,:); clear Ub Lm Um Pm Qm Rm
vc20 = v_true ./ (1 - (1./(pi*Q_true))*log(20/fref) - 1i./(2*Q_true));
[Lv,Uv,Pv,Qv,Rv] = lu(helm_op(vc20, 20, P)); Uv20 = Qv*(Uv\(Lv\(Pv*(Rv\b20)))); Dcal = (Uv20(irec,:) - Dbg20) .* rec_mask; clear Uv20 Lv Uv Pv Qv Rv
sig_field = sqrt(mean(abs(Dcal(rec_mask)).^2)) * 10^(-6/20);

% Forward data uses TRUE Ricker W_true(f); Migration source Us uses FLAT W_mig(f) * (2*pi*f)^2 (EXACTLY like fld_run.m!)
Cs_tot_all = zeros(nin, ns, nf, 'single');
I_RAW_R0 = zeros(nin, 1);
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

% Blind dip from R0 -> DSA mask & dip-aware t*
R0_2d = prep_img(I_RAW_R0, P.nzi, P.nxi, z_in);
[phi_eff, ~] = blind_structure_dip(R0_2d, P.dx);
[w_dsa, TS_dsa] = build_dsa_mask_and_tstar(Xin(:), Zin(:), xs, xr_lo, xr_hi, phi_eff(:), Fts_in);

% Measure shallow-anchor RMS spectrum S_anc(f) of DSA_R0 in [190, 320] m and anchor t*_anc!
anc_mask = (Zin(:) >= 190 & Zin(:) <= 320) & (Xin(:) >= 240 & Xin(:) <= 1360);
S_anc = zeros(1, nf);
for jf = 1:nf
    img_f = real(sum(w_dsa .* double(Cs_tot_all(:,:,jf)), 2));
    S_anc(jf) = rms(img_f(anc_mask));
end
% Smooth shallow reference spectrum over 3 bins & form target shaping factor W_shape(f) = W_mig(f) / S_anc_smooth(f)
S_anc_s = movmean(S_anc, 3);
W_shape = W_mig ./ max(S_anc_s, 0.10 * max(S_anc_s));
W_shape = W_shape / median(W_shape(freqs >= 10 & freqs <= 26));
% Anchor t*_anc = median of TS_dsa in shallow anchor window
tstar_anc = median(TS_dsa(anc_mask, :), 'all');
TS_rel = max(TS_dsa - tstar_anc, 0.0);

% Accumulate RAW_R0, DSA_R0, RAW_UQW1 (unanchored surface t*!), DSA_UQW1 (unanchored), DSA_UQW1_Anc (shallow-anchored + spectrum balanced!)
I_DSA_R0     = zeros(nin, 1);
I_RAW_UQW1   = zeros(nin, nm);
I_DSA_UQW1   = zeros(nin, nm);
I_DSA_Anc    = zeros(nin, nm);
I_DSA_Anc_W1 = zeros(nin, nm);

xmax = 1.05 * max(max(TS_flat(:)), max(TS_dsa(:))) / min(ms);
xg = linspace(0, xmax, 1201); dxg = xg(2);

for jf = 1:nf
    f = freqs(jf);
    Cs_t = double(Cs_tot_all(:,:,jf));
    I_DSA_R0 = I_DSA_R0 + real(sum(w_dsa .* Cs_t, 2));

    % Blind shot-flip iv for RAW and DSA
    pt_r = smooth2d(reshape(abs(sum(Cs_t, 2)).^2, P.nzi, P.nxi), 5.0);
    pf_r = smooth2d(reshape(abs(sum(Cs_t .* flip_sgn, 2)).^2, P.nzi, P.nxi), 5.0);
    iv_r = min(max(pf_r(:) ./ max(pt_r(:) - pf_r(:), 0.05*pf_r(:)), 1e-3), 1e3);

    pt_d = smooth2d(reshape(abs(sum(w_dsa .* Cs_t, 2)).^2, P.nzi, P.nxi), 5.0);
    pf_d = smooth2d(reshape(abs(sum(w_dsa .* Cs_t .* flip_sgn, 2)).^2, P.nzi, P.nxi), 5.0);
    iv_d = min(max(pf_d(:) ./ max(pt_d(:) - pf_d(:), 0.05*pf_d(:)), 1e-3), 1e3);

    T0 = uq_tab(f, xg, 0, fref);
    T1 = uq_tab(f, xg, sig, fref);

    for im = 1:nm
        % RAW_UQW1 (unanchored surface t*, exactly like Volve fld_run.m!)
        u = (TS_flat / ms(im)) / dxg; i0 = min(floor(u), numel(xg)-2); ww = u - i0; i0 = i0 + 1;
        it = @(tb) tb(i0).*(1-ww) + tb(i0+1).*ww;
        A1 = it(T1.A); C1 = it(T1.C); P1 = it(T1.P);
        Huq = (A1 + s70*P1.*C1 + P1.*iv_r) ./ (1 + s70*C1 + iv_r);
        I_RAW_UQW1(:, im) = I_RAW_UQW1(:, im) + real(sum(Huq .* Cs_t, 2));

        % DSA_UQW1 (unanchored surface t*)
        u = (TS_dsa / ms(im)) / dxg; i0 = min(floor(u), numel(xg)-2); ww = u - i0; i0 = i0 + 1;
        it = @(tb) tb(i0).*(1-ww) + tb(i0+1).*ww;
        A1 = it(T1.A); C1 = it(T1.C); P1 = it(T1.P);
        Huq = (A1 + s70*P1.*C1 + P1.*iv_d) ./ (1 + s70*C1 + iv_d);
        I_DSA_UQW1(:, im) = I_DSA_UQW1(:, im) + real(sum((w_dsa .* Huq) .* Cs_t, 2));

        % DSA_UQW1_Anc & DSA_W1_Anc (shallow-anchored relative dt* + shallow-anchor spectrum normalization!)
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

% Exact fld_eval.m 5-well flat-bandpass vertical well-tie in reservoir [540, 940] m!
well_xs = [400, 600, 800, 1000, 1200];
iz_res = find(z_in >= 540 & z_in <= 940);
v_in = v_true(P.in{1}, P.in{2});
% True vertical reflectivity rw = d(ln v)/dz
rw_map = zeros(size(v_in));
rw_map(2:end, :) = diff(log(v_in), 1, 1);

% Fixed flat well-tie wavenumber band [4.0, 18.0] cyc/km (like Stage 4H fixed-band well tie!)
lo_k = 4.0; hi_k = 18.0;
eval_wells = @(Ivec) score_wells(prep_img(Ivec, P.nzi, P.nxi, z_in), rw_map, x_in, iz_res, well_xs, P.dx, lo_k, hi_k);

tie_raw_r0 = eval_wells(I_RAW_R0);
tie_dsa_r0 = eval_wells(I_DSA_R0);
tie_raw_uq = arrayfun(@(im) eval_wells(I_RAW_UQW1(:,im)), 1:nm);
tie_dsa_uq = arrayfun(@(im) eval_wells(I_DSA_UQW1(:,im)), 1:nm);
tie_anc_uq = arrayfun(@(im) eval_wells(I_DSA_Anc(:,im)),  1:nm);
tie_anc_w1 = arrayfun(@(im) eval_wells(I_DSA_Anc_W1(:,im)), 1:nm);

fprintf('=== Part 2 H5_VolveAnc (5-well flat-band [4-18 cyc/km] well-tie, Volve flat-W replica) ===\n');
fprintf('  RAW_R0        : %.4f\n', tie_raw_r0);
fprintf('  DSA_R0        : %.4f\n', tie_dsa_r0);
fprintf('  RAW_UQW1      : mean=%.4f (m0.5=%.4f, m1.0=%.4f, m2.0=%.4f)\n', mean(tie_raw_uq), tie_raw_uq(1), tie_raw_uq(3), tie_raw_uq(5));
fprintf('  DSA_UQW1      : mean=%.4f (m0.5=%.4f, m1.0=%.4f, m2.0=%.4f)\n', mean(tie_dsa_uq), tie_dsa_uq(1), tie_dsa_uq(3), tie_dsa_uq(5));
fprintf('  DSA_W1_Anc    : mean=%.4f (m0.5=%.4f, m1.0=%.4f, m2.0=%.4f)\n', mean(tie_anc_w1), tie_anc_w1(1), tie_anc_w1(3), tie_anc_w1(5));
fprintf('  DSA_UQW1_Anc  : mean=%.4f (m0.5=%.4f, m1.0=%.4f, m2.0=%.4f)\n', mean(tie_anc_uq), tie_anc_uq(1), tie_anc_uq(3), tie_anc_uq(5));

res = struct('model', 'H5_VolveAnc', 'seed', seed, 'smoke', smoke, 'm', ms, ...
             'tstar_anc_ms', tstar_anc*1e3, 'band_cyckm', [lo_k hi_k], ...
             'RAW_R0', tie_raw_r0, 'DSA_R0', tie_dsa_r0, ...
             'RAW_UQW1', struct('tie', tie_raw_uq, 'mean', mean(tie_raw_uq)), ...
             'DSA_UQW1', struct('tie', tie_dsa_uq, 'mean', mean(tie_dsa_uq)), ...
             'DSA_W1_Anc', struct('tie', tie_anc_w1, 'mean', mean(tie_anc_w1)), ...
             'DSA_UQW1_Anc', struct('tie', tie_anc_uq, 'mean', mean(tie_anc_uq)), ...
             'elapsed_s', toc(t0_all));
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
    for L = -3:3
        best = max(best, cc(circshift(tr_b, L), rw_b));
    end
    cs(iw) = best;
end
c_mean = mean(cs);
end

function [phi_eff, c_st] = blind_structure_dip(R0_2d, dx)
[gx, gz] = gradient(R0_2d, dx);
Jxx = smooth2d(gx.^2, 4.0);
Jzz = smooth2d(gz.^2, 4.0);
Jxz = smooth2d(gx.*gz, 4.0);
phi_raw = -0.5 * atan2(2*Jxz, Jzz - Jxx);
phi_raw = max(min(phi_raw, deg2rad(30)), -deg2rad(30));
c_st = sqrt((Jzz - Jxx).^2 + 4*Jxz.^2) ./ max(Jzz + Jxx, 1e-30);
phi_eff = c_st .* phi_raw;
end

function w_sa = build_dsa_mask(xv, zv, xs, xr_lo, xr_hi, phi_v)
nin = numel(xv); ns = numel(xs);
w_sa = zeros(nin, ns);
dz = max(zv - 20.0, 10.0);
for k = 1:ns
    th_s = atan2(xv - xs(k), dz);
    th_r = max(min(th_s + 2*phi_v, deg2rad(78)), -deg2rad(78));
    xmir = xv + dz .* tan(th_r);
    d_in = min(xmir - xr_lo(k), xr_hi(k) - xmir);
    t_arr = cos_taper(d_in, -80.0, 120.0);
    t_gs  = cos_taper(80.0 - abs(rad2deg(th_s)), 0.0, 10.0);
    t_gr  = cos_taper(80.0 - abs(rad2deg(th_r)), 0.0, 10.0);
    w_sa(:, k) = max(t_arr .* t_gs .* t_gr, 1e-4);
end
end

function [w_dsa, TS_dsa] = build_dsa_mask_and_tstar(xv, zv, xs, xr_lo, xr_hi, phi_v, Fts_in)
nin = numel(xv); ns = numel(xs);
w_dsa  = zeros(nin, ns);
TS_dsa = zeros(nin, ns);
dz = max(zv - 20.0, 10.0);
row_idx = (1:nin)';
for k = 1:ns
    th_s = atan2(xv - xs(k), dz);
    th_r = max(min(th_s + 2*phi_v, deg2rad(78)), -deg2rad(78));
    xmir = xv + dz .* tan(th_r);
    d_in = min(xmir - xr_lo(k), xr_hi(k) - xmir);
    t_arr = cos_taper(d_in, -80.0, 120.0);
    t_gs  = cos_taper(80.0 - abs(rad2deg(th_s)), 0.0, 10.0);
    t_gr  = cos_taper(80.0 - abs(rad2deg(th_r)), 0.0, 10.0);
    w_dsa(:, k) = max(t_arr .* t_gs .* t_gr, 1e-4);

    q = (xmir - xs(1)) / (xs(2) - xs(1)) + 1;
    q = min(max(q, 1), ns); k0 = floor(q); k0(k0 == ns) = ns - 1; ww = q - k0;
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
