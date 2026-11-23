function spg_rtm_bench(models_req, smoke)
% Pre-registered SPG-UQ-RTM 2D Full-Wave Multi-Geometry Suite (G1-G4)
% Compares RAW, SA (live-array specular mask), and SPG (SA x Poynting structural-dip taper)
% fairly across R0, IQ70, W1, and UQW1 on 4 brand-new unseen models:
%   G1: Gas Chimney (Q=12) + Anticline (seed 20261131)
%   G2: Normal Fault + Steep 18-deg Dipping Layers + Deep Lens (seed 20261132)
%   G3: Volve-like One-Sided Trailing Marine Streamer + Lateral Q Gradient (seed 20261133)
%   G4: Laterally Correlated Noise (sigma=60m) + Gas Chimney (Q=15) + +3% v_mig error (seed 20261134)
% Peak RAM < 0.5 GB, 2 threads, resumable per model (results/<name>.json).

if nargin < 1 || isempty(models_req), models_req = {'G1', 'G2', 'G3', 'G4'}; end
if ischar(models_req) || isstring(models_req), models_req = {char(models_req)}; end
if nargin < 2, smoke = false; end

maxNumCompThreads(2);
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
    res = run_model(mname, smoke);
    fid = fopen(fn, 'w'); fwrite(fid, jsonencode(res, 'PrettyPrint', true), 'char'); fclose(fid);
    fprintf('Saved %s (%.1fs)\n\n', fn, res.elapsed_s);
end
end

function res = run_model(mname, smoke)
t0_all = tic;
P.dx = 10.0; P.npml = 18; P.nzi = 100; P.nxi = 160;
P.nz = P.nzi + 2*P.npml; P.nx = P.nxi + 2*P.npml;
P.z = ((1:P.nz)' - P.npml - 1) * P.dx;
P.x = ((1:P.nx)  - P.npml - 1) * P.dx;
P.in = {P.npml + (1:P.nzi), P.npml + (1:P.nxi)};
P.izs = P.npml + 3;
P.ixs = P.npml + 1 + (8:6:152);      % 25 shots, 60 m spacing (x = 80..1520 m)
P.ixr = P.npml + 1 + (0:2:P.nxi-1);  % 80 receiver grid slots, 20 m spacing (x = 0..1580 m)
ns = numel(P.ixs); nr = numel(P.ixr); N = P.nz * P.nx;
[Zi, Xi] = ndgrid(P.in{1}, P.in{2}); lin = sub2ind([P.nz P.nx], Zi(:), Xi(:)); nin = numel(lin);
isrc = sub2ind([P.nz P.nx], P.izs*ones(1,ns), P.ixs);
irec = sub2ind([P.nz P.nx], P.izs*ones(1,nr), P.ixr);
[Zfull, Xfull] = ndgrid(P.z, P.x);
z_in = P.z(P.in{1}); x_in = P.x(P.in{2});
[Zin, Xin] = ndgrid(z_in, x_in);

xs = P.x(P.ixs);
xr = P.x(P.ixr);

% Default receiver mask (all 80 receivers live for every shot; overridden in G3!)
rec_mask = true(nr, ns);
vscale_mig = 1.0;
corr_noise_sig = 0.0; % in receiver traces (dx_r = 20 m)
has_chimney = false;

switch mname
    case 'G1'
        seed = 20261131;
        has_chimney = true;
        z1 = 230.0 + 0*Xfull;
        z2 = 635.0 - 30.0 * exp(-((Xfull - 820.0)/280.0).^2);
        z3 = 765.0 - 20.0 * exp(-((Xfull - 820.0)/300.0).^2);
        z4 = 885.0 + 10.0 * sin(2*pi*(Xfull - 200)/1400);
        v_true = 1850.0 * ones(P.nz, P.nx);
        v_true(Zfull >= z1) = 2120.0;
        v_true(Zfull >= z2) = 2520.0;
        v_true(Zfull >= z3) = 2880.0;
        v_true(Zfull >= z4) = 3280.0;
        invQ = (1/95)*ones(P.nz, P.nx); invQ(Zfull < 80) = 1/1000;
        chim = exp(-(((Xfull - 820)/200).^2 + ((Zfull - 425)/155).^2).^1.5);
        invQ = invQ + (1/12 - 1/95)*chim;
        Q_true = 1./invQ;
        W_chim = {[600 1040], [560 950]}; W_flank = {[180 500], [560 950]};

    case 'G2'
        % Normal fault + steep 18-deg dipping layers (tan 18 = 0.325) + deep low-Q lens
        seed = 20261132;
        Zg = Zfull;
        hw = (Xfull > 760.0 + (Zfull - 200.0) / tand(68.0)) & (Zfull > 200.0);
        Zg(hw) = Zg(hw) - 35.0; % 35 m normal fault throw
        z1 = 240.0 + 0.08 * (Xfull - 800.0);
        z2 = 650.0 + 0.32 * (Xfull - 800.0);  % ~18 deg steep dip!
        z3 = 790.0 - 0.25 * (Xfull - 800.0);  % ~-14 deg opposite steep dip!
        z4 = 900.0 + 0.10 * (Xfull - 800.0);
        v_true = 1850.0 * ones(P.nz, P.nx);
        v_true(Zg >= z1) = 2150.0;
        v_true(Zg >= z2) = 2600.0;
        v_true(Zg >= z3) = 2950.0;
        v_true(Zg >= z4) = 3320.0;
        invQ = (1/85)*ones(P.nz, P.nx); invQ(Zfull < 80) = 1/1000;
        lens = exp(-(((Xfull - 560)/210).^2 + ((Zfull - 480)/90).^2).^2);
        invQ = invQ + (1/18 - 1/85)*lens;
        Q_true = 1./invQ;
        W_chim = {[350 770], [560 950]}; W_flank = {[950 1350], [560 950]};

    case 'G3'
        % Volve-like One-Sided Trailing Marine Streamer (960 m cable on right/left depending on shot,
        % here vessel sails left->right, streamer trails to the LEFT [xs-960, xs-40] or right where needed ->
        % realistic fixed trailing streamer [xs - 960, xs - 40] + mirror at left edge? No, exact marine sail:
        % trailing cable [xs - 980, xs - 40] for right-going line, active receivers inside [0, 1580] ->
        % wait: if shots span [80..1520] and streamer is trailing left [xs-980, xs-40], left shots have few recs!
        % Just like Volve (src_x=[-3600..2550], rec_x=[-3430..2550], offset max 5000m!), give each shot an
        % asymmetric live spread [xs - 300, xs + 900] intersected with [0, 1580] -> ~40% out-of-array mirror shots!
        seed = 20261133;
        for s = 1:ns
            rec_mask(:, s) = (xr >= xs(s) - 240.0) & (xr <= xs(s) + 860.0);
        end
        z1 = 240.0 + 0.04 * (Xfull - 800.0);
        z2 = 620.0 - 20.0 * cos(pi * (Xfull - 800.0) / 700.0);
        z3 = 760.0 + 0.06 * (Xfull - 800.0);
        z4 = 885.0 - 0.03 * (Xfull - 800.0);
        v_true = 1820.0 * ones(P.nz, P.nx);
        v_true(Zfull >= z1) = 2120.0;
        v_true(Zfull >= z2) = 2500.0;
        v_true(Zfull >= z3) = 2860.0;
        v_true(Zfull >= z4) = 3240.0;
        % Strong lateral Q gradient: Q = 35 on left -> 115 on right
        Q_grad = 35.0 + 80.0 * min(max(Xfull / 1600.0, 0), 1);
        Q_true = Q_grad; Q_true(Zfull < 80) = 1000.0;
        W_chim = {[180 620], [560 950]}; W_flank = {[980 1420], [560 950]};

    case 'G4'
        % Laterally correlated noise (sigma_tr = 3 traces = 60 m) + Q=15 chimney + +3% velocity error!
        seed = 20261134;
        has_chimney = true;
        vscale_mig = 1.03;
        corr_noise_sig = 3.0;  % 3 traces * 20 m = 60 m lateral correlation length!
        z1 = 235.0 + 0*Xfull;
        z2 = 630.0 - 20.0 * exp(-((Xfull - 780.0)/290.0).^2);
        z3 = 770.0 - 15.0 * exp(-((Xfull - 780.0)/310.0).^2);
        z4 = 890.0 + 0*Xfull;
        v_true = 1860.0 * ones(P.nz, P.nx);
        v_true(Zfull >= z1) = 2160.0;
        v_true(Zfull >= z2) = 2560.0;
        v_true(Zfull >= z3) = 2920.0;
        v_true(Zfull >= z4) = 3300.0;
        invQ = (1/90)*ones(P.nz, P.nx); invQ(Zfull < 80) = 1/1000;
        chim = exp(-(((Xfull - 780)/210).^2 + ((Zfull - 430)/160).^2).^1.5);
        invQ = invQ + (1/15 - 1/90)*chim;
        Q_true = 1./invQ;
        W_chim = {[560 1000], [560 950]}; W_flank = {[180 480], [560 950]};

    otherwise
        error('Unknown model %s', mname);
end

rng(seed, 'twister');
v_mig_exact = smooth2d(v_true, 2.5);
v_mig = vscale_mig * v_mig_exact;

% Eikonal two-way t* (using v_mig & standard two_way.m!) + source ray directions (px_s, pz_s)
s_mig = 1.0 ./ v_mig; eps_q = 0.5;
F.T = zeros(P.nz, P.nx, ns, 'single'); F.ts = F.T;
for k = 1:ns
    F.T(:,:,k) = eikonal_fs(s_mig, P.dx, P.izs, P.ixs(k));
    F.ts(:,:,k) = (eikonal_fs(s_mig.*(1+eps_q./Q_true), P.dx, P.izs, P.ixs(k)) - ...
                   eikonal_fs(s_mig.*(1-eps_q./Q_true), P.dx, P.izs, P.ixs(k))) / (2*eps_q);
end
TS = zeros(nin, ns);
px_s = zeros(nin, ns); pz_s = zeros(nin, ns);
for k = 1:ns
    [~, tt] = two_way(F, k, P); TS(:, k) = double(tt(lin));
    Tk = double(F.T(P.in{1}, P.in{2}, k));
    [gx, gz] = gradient(Tk, P.dx);
    gn = max(sqrt(gx.^2 + gz.^2), 1e-12);
    px_s(:, k) = reshape(gx ./ gn, nin, 1);
    pz_s(:, k) = reshape(gz ./ gn, nin, 1);
end
clear F

% Build Live-Array Specular Mask w_sa (nin x ns):
% Uses EACH SHOT'S ACTUAL LIVE RECEIVER INTERVAL [xr_lo(s), xr_hi(s)]!
% NO premature 55-deg opening-angle cutoff (only grazing-wave mute at 70-80 deg)!
w_sa = zeros(nin, ns);
out_pairs = 0; tot_pairs = 0;
w_deep_mask = (Zin(:) >= 560 & Zin(:) <= 950) & (Xin(:) >= 180 & Xin(:) <= 1420);
for k = 1:ns
    live_xr = xr(rec_mask(:, k));
    xr_lo = min(live_xr); xr_hi = max(live_xr);
    xmir = 2.0 * Xin(:) - xs(k);
    out_pairs = out_pairs + sum((xmir(w_deep_mask) < xr_lo) | (xmir(w_deep_mask) > xr_hi));
    tot_pairs = tot_pairs + sum(w_deep_mask);
    d_in = min(xmir - xr_lo, xr_hi - xmir);
    t_arr = cos_taper(d_in, -80.0, 120.0);
    th_s = acosd(min(max(pz_s(:, k), -1), 1));
    t_graz = cos_taper(80.0 - th_s, 0.0, 10.0);
    w_sa(:, k) = max(t_arr .* t_graz, 1e-4);
end
out_frac = out_pairs / max(tot_pairs, 1);
fprintf('=== Model %s (seed=%d, out-of-array=%.1f%%, vscale=%.2f, corr_sig=%.1f tr) ===\n', ...
    mname, seed, 100*out_frac, vscale_mig, corr_noise_sig);

if smoke
    freqs = [10.0, 22.0, 34.0];
else
    freqs = 4.0:1.0:37.0;  % 34 freqs, df = 1.0 Hz -> wraparound 1.0s > T_max
end
nf = numel(freqs);
fref = 25.0; fpeak = 22.0; sig = 0.447; s70 = gain_s(70);
ms = [0.5, 0.7, 1.0, 1.4, 2.0]; nm = numel(ms);
flip_sgn = 2 * (rand(1, ns) > 0.5) - 1;

% Calibrate 15 dB and 3 dB noise on LIVE masked reflections at f = 22 Hz
W22 = (22/fpeak)^2 * exp(-(22/fpeak)^2);
b22 = zeros(N, ns); b22(sub2ind([N ns], isrc, 1:ns)) = -W22 / P.dx^2;
[Lm,Um,Pm,Qm,Rm] = lu(helm_op(v_mig_exact, 22, P)); Ub = Qm*(Um\(Lm\(Pm*(Rm\b22)))); Dbg22 = Ub(irec,:); clear Ub Lm Um Pm Qm Rm
vc22 = v_true ./ (1 - (1./(pi*Q_true))*log(22/fref) - 1i./(2*Q_true));
[Lv,Uv,Pv,Qv,Rv] = lu(helm_op(vc22, 22, P)); Uv22 = Qv*(Uv\(Lv\(Pv*(Rv\b22)))); Dcal = (Uv22(irec,:) - Dbg22) .* rec_mask; clear Uv22 Lv Uv Pv Qv Rv
d_rms = sqrt(mean(abs(Dcal(rec_mask)).^2));
sig_clean = d_rms * 10^(-15/20);
sig_noisy = d_rms * 10^(-3/20);

modes = {'RAW', 'SA', 'SPG'};
ops   = {'R0', 'IQ70', 'W1', 'UQW1'};
I_ref = zeros(nin, 1);
for rg = {'clean', 'noisy'}
    for imd = 1:numel(modes)
        for io = 1:numel(ops)
            key = [modes{imd} '_' ops{io}];
            nc = 1 + (nm-1)*~strcmp(ops{io}, 'R0');
            IS.(rg{1}).(key) = zeros(nin, nc);
            IN.(rg{1}).(key) = zeros(nin, nc);
        end
    end
end

xmax = 1.05 * max(TS(:)) / min(ms);
xg = linspace(0, xmax, 1201); dxg = xg(2);

% Build 1D receiver correlation kernel if corr_noise_sig > 0
if corr_noise_sig > 0
    rc = ceil(3 * corr_noise_sig);
    g_nc = exp(-0.5 * ((-rc:rc)' / corr_noise_sig).^2);
    g_nc = g_nc / sqrt(sum(g_nc.^2)); % unit variance preservation!
end

for jf = 1:nf
    f = freqs(jf); w = 2*pi*f;
    W_src = (f/fpeak)^2 * exp(-(f/fpeak)^2);
    b_src = zeros(N, ns); b_src(sub2ind([N ns], isrc, 1:ns)) = -W_src / P.dx^2;

    % Direct wave background uses exact smooth velocity v_mig_exact (so true data D_el, D_vis are physical!)
    [Lb,Ub_,Pb,Qb,Rb] = lu(helm_op(v_mig_exact, f, P));
    U_bg_ex = Qb*(Ub_\(Lb\(Pb*(Rb\b_src)))); D_bg = U_bg_ex(irec,:);
    if vscale_mig == 1.0
        Lm = Lb; Um = Ub_; Pm = Pb; Qm = Qb; Rm = Rb;
        Us_c = conj(U_bg_ex(lin,:)) * w^2;
        clear U_bg_ex
    else
        clear U_bg_ex Lb Ub_ Pb Qb Rb
        [Lm,Um,Pm,Qm,Rm] = lu(helm_op(v_mig, f, P));
        U_mig = Qm*(Um\(Lm\(Pm*(Rm\b_src))));
        Us_c = conj(U_mig(lin,:)) * w^2;
        clear U_mig
    end

    % Elastic reference data (Q = Inf)
    [Le,Ue,Pe,Qe,Re] = lu(helm_op(v_true, f, P));
    U_el = Qe*(Ue\(Le\(Pe*(Re\b_src)))); D_el = (U_el(irec,:) - D_bg) .* rec_mask; clear U_el Le Ue Pe Qe Re

    % Viscoacoustic true data
    vc = v_true ./ (1 - (1./(pi*Q_true))*log(max(f,1e-6)/fref) - 1i./(2*Q_true));
    [Lv,Uv,Pv,Qv,Rv] = lu(helm_op(vc, f, P));
    U_vi = Qv*(Uv\(Lv\(Pv*(Rv\b_src)))); D_vi = (U_vi(irec,:) - D_bg) .* rec_mask; clear U_vi Lv Uv Pv Qv Rv b_src

    % Receiver noise (white or laterally correlated, masked to live receivers!)
    N_unit = (1/sqrt(2)) * (randn(nr, ns) + 1i*randn(nr, ns));
    if corr_noise_sig > 0
        N_unit = conv2(g_nc, 1, N_unit, 'same');
    end
    N_unit = N_unit .* rec_mask;

    bp = @(D) extract_bp(D, Lm, Um, Pm, Qm, Rm, irec, lin, N, ns, P.dx);
    V_el = bp(D_el); V_vi = bp(D_vi); V_nu = bp(N_unit);
    clear Lm Um Pm Qm Rm

    % Elastic reference uses exact migration velocity? Wait: I_ref is always migrated with same v_mig_exact?
    % Keep consistent with uq_rtm_20261117: ref is elastic data migrated with v_mig (exact when vscale=1, or true ref)!
    % Better: when vscale_mig ~= 1, migrate elastic ref with v_mig_exact? No, exact smooth v_mig_exact would shift
    % depth by 3% (30 m > ±2 samples E1 lag!), whereas uq_rtm_20261117 migrates ref with the SAME vmig so E1 measures
    % Q-compensation quality rather than kinematic shift! (migrate2.m line 53 migrates ref with vmig!)
    I_ref = I_ref + real(sum(Us_c .* V_el, 2)); clear V_el

    Cs_sig = Us_c .* V_vi; Cs_nu = Us_c .* V_nu;
    T0 = uq_tab(f, xg, 0, fref);
    T1 = uq_tab(f, xg, sig, fref);

    for rg = {'clean', 'noisy'}
        rgn = rg{1};
        s_nr = sig_clean; if strcmp(rgn, 'noisy'), s_nr = sig_noisy; end
        Cs_noi = s_nr * Cs_nu;
        Cs_tot = Cs_sig + Cs_noi;
        V_tot  = V_vi + s_nr * V_nu;

        % Strictly blind Poynting structural-dip taper (30-50 deg taper -> preserves <=30 deg steep dips & 60 deg wide angles!)
        w_spg = compute_poynting_spg(V_tot, px_s, pz_s, P.nzi, P.nxi, P.dx, w, w_sa);
        W_list = {ones(nin, ns), w_sa, w_spg};

        % Mode-consistent blind shot-flip (B=all receivers/shot -> unbiased even for l30 laterally correlated noise!)
        iv_list = cell(1, 3);
        for imd = 1:3
            Wm = W_list{imd};
            stk_tot = abs(sum(Wm .* Cs_tot, 2)).^2;
            stk_flp = abs(sum(Wm .* Cs_tot .* flip_sgn, 2)).^2;
            p_tot = smooth2d(reshape(stk_tot, P.nzi, P.nxi), 5.0);
            p_flp = smooth2d(reshape(stk_flp, P.nzi, P.nxi), 5.0);
            iv_list{imd} = min(max(p_flp(:) ./ max(p_tot(:) - p_flp(:), 0.05*p_flp(:)), 1e-3), 1e3);
            key0 = [modes{imd} '_R0'];
            IS.(rgn).(key0) = IS.(rgn).(key0) + real(sum(Wm .* Cs_sig, 2));
            IN.(rgn).(key0) = IN.(rgn).(key0) + real(sum(Wm .* Cs_noi, 2));
        end

        for im = 1:nm
            u = (TS / ms(im)) / dxg; i0 = min(floor(u), numel(xg)-2); ww = u - i0; i0 = i0 + 1;
            it = @(tb) tb(i0).*(1-ww) + tb(i0+1).*ww;
            A0 = it(T0.A); C0 = it(T0.C); P0 = it(T0.P);
            A1 = it(T1.A); C1 = it(T1.C); P1 = it(T1.P);
            H_iq = (A0 + s70*P0.*C0) ./ (1 + s70*C0);

            for imd = 1:3
                Wm = W_list{imd}; iv = iv_list{imd};
                H_w1 = (A0 + s70*P0.*C0 + P0.*iv) ./ (1 + s70*C0 + iv);
                H_uq = (A1 + s70*P1.*C1 + P1.*iv) ./ (1 + s70*C1 + iv);
                mstr = modes{imd};
                IS.(rgn).([mstr '_IQ70'])(:,im) = IS.(rgn).([mstr '_IQ70'])(:,im) + real(sum((Wm.*H_iq).*Cs_sig, 2));
                IN.(rgn).([mstr '_IQ70'])(:,im) = IN.(rgn).([mstr '_IQ70'])(:,im) + real(sum((Wm.*H_iq).*Cs_noi, 2));
                IS.(rgn).([mstr '_W1'])(:,im)   = IS.(rgn).([mstr '_W1'])(:,im)   + real(sum((Wm.*H_w1).*Cs_sig, 2));
                IN.(rgn).([mstr '_W1'])(:,im)   = IN.(rgn).([mstr '_W1'])(:,im)   + real(sum((Wm.*H_w1).*Cs_noi, 2));
                IS.(rgn).([mstr '_UQW1'])(:,im) = IS.(rgn).([mstr '_UQW1'])(:,im) + real(sum((Wm.*H_uq).*Cs_sig, 2));
                IN.(rgn).([mstr '_UQW1'])(:,im) = IN.(rgn).([mstr '_UQW1'])(:,im) + real(sum((Wm.*H_uq).*Cs_noi, 2));
            end
        end
    end
    if mod(jf, 10) == 0 || jf == nf
        fprintf('  %s f=%4.1f Hz (%2d/%d), elapsed %.1fs\n', mname, f, jf, nf, toc(t0_all));
    end
end

% Evaluate all 12 methods on W_deep and W_chim
ref_2d = prep_img(I_ref, P.nzi, P.nxi, z_in);
ref_hf = highpass_kz(ref_2d, P.dx, 10.0);
W_deep = {[180 1420], [560 950]};
rms_w = @(A, W) sqrt(mean(A(z_in>=W{2}(1)&z_in<=W{2}(2), x_in>=W{1}(1)&x_in<=W{1}(2)).^2, 'all'));
ref_rcf = log10(rms_w(ref_2d, W_chim) / max(rms_w(ref_2d, W_flank), 1e-30));

res = struct('model', mname, 'seed', seed, 'smoke', smoke, 'n_freq', nf, 'm', ms, ...
             'out_of_array_frac', out_frac, 'vscale_mig', vscale_mig, ...
             'corr_noise_sig_tr', corr_noise_sig, 'has_chimney', has_chimney, ...
             'elapsed_s', toc(t0_all));

for rg = {'clean', 'noisy'}
    rgn = rg{1};
    for wn = {'W_deep', 'W_chim'}
        wname = wn{1}; W = W_deep; if strcmp(wname, 'W_chim'), W = W_chim; end
        iz_w = find(z_in >= W{2}(1) & z_in <= W{2}(2));
        ix_w = find(x_in >= W{1}(1) & x_in <= W{1}(2));
        if strcmp(wname, 'W_deep')
            fprintf('\n  [%s %s W_deep]\n', mname, rgn);
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
                    eamp(c)= abs(log10(rms_w(T2, W_chim)/max(rms_w(T2, W_flank),1e-30)) - ref_rcf);
                end
                res.(rgn).(wname).(key) = struct( ...
                    'E1', e1, 'E1_mean', mean(e1), ...
                    'E1_HF', e1h, 'E1_HF_mean', mean(e1h), ...
                    'E_amp', eamp, 'E_amp_mean', mean(eamp), ...
                    'SNR_img_dB', snr, 'SNR_img_dB_mean', mean(snr));
                if strcmp(wname, 'W_deep')
                    fprintf('    %-9s | E1=%.4f (m0.5=%.4f,m1=%.4f,m2=%.4f) | HF=%.4f | Eamp=%.3f | SNR=%.2f dB\n', ...
                        key, mean(e1), e1(1), e1(min(3,nc)), e1(end), mean(e1h), mean(eamp), mean(snr));
                end
            end
        end
    end
end
end

function w_spg = compute_poynting_spg(V_tot, px_s, pz_s, nzi, nxi, dx, w, w_sa)
% Structural dip taper ONLY (30-50 deg smooth cosine taper):
% Preserves steep reflectors up to 30 deg AND wide-angle flat undershooting rays (dip=0 deg),
% while suppressing steep >50 deg migration smiles and coherent noise arcs.
[nin, ns] = size(V_tot);
w_spg = zeros(nin, ns);
for k = 1:ns
    Vk = reshape(V_tot(:, k), nzi, nxi);
    [Vx, Vz] = gradient(Vk, dx);
    num_x = smooth2d(real(conj(Vk) .* Vx / (1i * w)), 1.5);
    num_z = smooth2d(real(conj(Vk) .* Vz / (1i * w)), 1.5);
    pnorm = max(sqrt(num_x.^2 + num_z.^2), 1e-30);
    px_r = num_x(:) ./ pnorm(:); pz_r = num_z(:) ./ pnorm(:);
    sgn = sign(pz_r); sgn(sgn == 0) = 1;
    px_r = px_r .* sgn; pz_r = pz_r .* sgn;
    dip_deg = atan2d(abs(px_s(:, k) + px_r), max(pz_s(:, k) + pz_r, 1e-6));
    t_dip = cos_taper(50.0 - dip_deg, 0.0, 20.0);  % 1 for dip<=30 deg, 0 for dip>=50 deg
    w_spg(:, k) = max(w_sa(:, k) .* t_dip, 1e-4);
end
end

function t = cos_taper(x, x0, x1)
u = min(max((x - x0) / max(x1 - x0, 1e-12), 0), 1);
t = 0.5 - 0.5 * cos(pi * u);
end

function V = extract_bp(Dmat, L, U, Pp, Qq, R, irec, lin, N, ns, dx)
b = zeros(N, ns); b(irec, :) = Dmat;
Vf = -(R \ (Pp' * (L' \ (U' \ (Qq' * b))))) / dx^2;
V = Vf(lin, :);
end

function B = smooth2d(A, sig)
r = max(1, ceil(3*sig)); g = exp(-0.5*((-r:r)/sig).^2); g = g/sum(g);
Ap = [repmat(A(1,:), r, 1); A; repmat(A(end,:), r, 1)];
B1 = conv2(g', 1, Ap, 'valid');
B1p = [repmat(B1(:,1), 1, r), B1, repmat(B1(:,end), 1, r)];
B = conv2(1, g, B1p, 'valid');
end

function y = prep_img(v, nzi, nxi, z_in)
y = reshape(double(v), nzi, nxi);
y = y - smooth2d(y, 5.0);
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
rr = r(iz, ix); c0 = cc(a(iz, ix), rr); cl = c0; nz = size(a, 1);
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
